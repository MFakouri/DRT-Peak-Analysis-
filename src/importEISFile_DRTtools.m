function [data_from_file, cancelled] = importEISFile_DRTtools(fullFileName)
% IMPORT EIS FILE FOR DRTTOOLS
%
% Output:
%   column 1 = Frequency
%   column 2 = Zreal
%   column 3 = Zimag
%
% If requested by the user, columns 2 and 3 are multiplied
% by the active area in cm^2 before the data are passed to DRTtools.

    cancelled = false;

    [~,~,ext] = fileparts(fullFileName);

    % ---------------------------------------------------------
    % 1. Read / standardize input
    % ---------------------------------------------------------

    switch lower(ext)

        case '.mat'

            storedStructure = load(fullFileName);

            data_from_file = [ ...
                storedStructure.freq, ...
                storedStructure.Z_prime, ...
                storedStructure.Z_double_prime];

        case {'.txt','.csv','.z','.dat'}

            try
                % First try the new intelligent column detector
                [data_from_file, ~, ~] = readEISFile(fullFileName);

            catch conversionError

                % -------------------------------------------------
                % Backward compatibility with normal DRTtools files
                % -------------------------------------------------

                switch lower(ext)

                    case '.txt'

                        try
                            fid = fopen(fullFileName, 'r');

                            if fid == -1
                                error('Could not open the selected TXT file.');
                            end

                            originalText = fread(fid, '*char')';
                            fclose(fid);

                            correctedText = strrep(originalText, ', ', '.');

                            fid = fopen(fullFileName, 'w');
                            fprintf(fid, '%s', correctedText);
                            fclose(fid);

                            try
                                data_from_file = dlmread(fullFileName);
                            catch ME
                                % Restore original file before reporting error
                                fid = fopen(fullFileName, 'w');
                                fprintf(fid, '%s', originalText);
                                fclose(fid);
                                rethrow(ME);
                            end

                            % Restore original file
                            fid = fopen(fullFileName, 'w');
                            fprintf(fid, '%s', originalText);
                            fclose(fid);

                        catch
                            rethrow(conversionError);
                        end

                    case '.csv'

                        try
                            data_from_file = csvread(fullFileName);
                        catch
                            rethrow(conversionError);
                        end

                    otherwise

                        rethrow(conversionError);
                end
            end

        otherwise

            error('Unsupported file type: %s', ext);
    end


    % ---------------------------------------------------------
    % 2. Ask about active area
    % ---------------------------------------------------------

    [applyArea, activeArea, cancelled] = askActiveArea();

    if cancelled
        data_from_file = [];
        return;
    end


    % ---------------------------------------------------------
    % 3. Apply active area
    % ---------------------------------------------------------

    if applyArea

        if size(data_from_file,2) < 3
            error('The imported data must contain Frequency, Zreal, and Zimag.');
        end

        data_from_file(:,2:3) = ...
            data_from_file(:,2:3) .* activeArea;
    end

end

function [applyArea, activeArea, cancelled] = askActiveArea()

    applyArea = false;
    activeArea = 1;
    cancelled = true;

    d = dialog( ...
        'Name', 'Active Area', ...
        'Position', [500 400 470 220], ...
        'WindowStyle', 'modal', ...
        'CloseRequestFcn', @cancelDialog);


    uicontrol(d, ...
        'Style', 'text', ...
        'String', 'Active Area Correction', ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', ...
        'Position', [25 170 250 25]);


    uicontrol(d, ...
        'Style', 'text', ...
        'String', ...
        'Check the box only if active area has NOT already been applied to the impedance.', ...
        'HorizontalAlignment', 'left', ...
        'Position', [25 135 420 35]);


    areaCheck = uicontrol(d, ...
        'Style', 'checkbox', ...
        'String', 'Active area has NOT already been applied', ...
        'Value', 0, ...
        'Position', [25 100 300 25], ...
        'Callback', @toggleArea);


    areaLabel = uicontrol(d, ...
        'Style', 'text', ...
        'String', 'Active area (cm^2):', ...
        'HorizontalAlignment', 'left', ...
        'Enable', 'off', ...
        'Position', [45 65 135 25]);


    areaEdit = uicontrol(d, ...
        'Style', 'edit', ...
        'String', '', ...
        'BackgroundColor', 'white', ...
        'Enable', 'off', ...
        'Position', [185 66 100 25]);


    uicontrol(d, ...
        'Style', 'pushbutton', ...
        'String', 'OK', ...
        'Position', [275 20 75 30], ...
        'Callback', @acceptDialog);


    uicontrol(d, ...
        'Style', 'pushbutton', ...
        'String', 'Cancel', ...
        'Position', [365 20 75 30], ...
        'Callback', @cancelDialog);


    uiwait(d);

    if ishghandle(d)
        delete(d);
    end


    function toggleArea(~,~)

        if get(areaCheck,'Value') == 1

            set(areaLabel,'Enable','on');
            set(areaEdit,'Enable','on');

        else

            set(areaLabel,'Enable','off');
            set(areaEdit,'Enable','off');

        end
    end


    function acceptDialog(~,~)

        if get(areaCheck,'Value') == 1

            value = str2double( ...
                strtrim(get(areaEdit,'String')));

            if isnan(value) || ...
                    ~isfinite(value) || ...
                    value <= 0

                errordlg( ...
                    'Please enter a valid positive Active area (cm^2).', ...
                    'Active Area Error', ...
                    'modal');

                return;
            end

            applyArea = true;
            activeArea = value;

        else

            applyArea = false;
            activeArea = 1;

        end

        cancelled = false;
        uiresume(d);
    end


    function cancelDialog(~,~)

        cancelled = true;
        uiresume(d);
    end

end

function [data, detectedHeaders, sourceWasNegativeImag] = readEISFile(fileName)
% Locate the actual header row, identify the three required EIS columns,
% and read only numeric rows beneath that header.

    rawText = fileread(fileName);

    % Remove UTF-8/Unicode BOM if MATLAB decoded it into the text.
    rawText = strrep(rawText, char(65279), '');

    lines = regexp(rawText, '\r\n|\n|\r', 'split');

    headerLine = 0;
    freqIdx = [];
    realIdx = [];
    imagIdx = [];
    sourceWasNegativeImag = false;
    delimiterType = '';
    detectedHeaders = {};

    for i = 1:numel(lines)
        [tokens, thisDelimiterType] = splitCandidateHeader(lines{i});

        if numel(tokens) < 3
            continue;
        end

        tempFreq = [];
        tempReal = [];
        tempImag = [];
        tempNegativeImag = false;

        for j = 1:numel(tokens)
            [kind, isNegativeImag] = classifyHeader(tokens{j});

            switch kind
                case 'freq'
                    if isempty(tempFreq)
                        tempFreq = j;
                    end

                case 'real'
                    if isempty(tempReal)
                        tempReal = j;
                    end

                case 'imag'
                    if isempty(tempImag)
                        tempImag = j;
                        tempNegativeImag = isNegativeImag;
                    end
            end
        end

        if ~isempty(tempFreq) && ~isempty(tempReal) && ~isempty(tempImag)
            headerLine = i;
            freqIdx = tempFreq;
            realIdx = tempReal;
            imagIdx = tempImag;
            sourceWasNegativeImag = tempNegativeImag;
            delimiterType = thisDelimiterType;

            detectedHeaders = { ...
                strtrim(tokens{freqIdx}), ...
                strtrim(tokens{realIdx}), ...
                strtrim(tokens{imagIdx})};

            break;
        end
    end

    if headerLine == 0
        error(['Could not find the required EIS columns. Supported labels include ', ...
               'Freq/Frequency, Zreal/Z''/impedance'', and ', ...
               'Zimag/Z''''/impedance'''' (including -Z'''').']);
    end

    maxIdx = max([freqIdx, realIdx, imagIdx]);
    data = zeros(0, 3);

    for i = headerLine + 1:numel(lines)
        line = strtrim(lines{i});

        if isempty(line)
            continue;
        end

        tokens = splitDataLine(lines{i}, delimiterType);

        if numel(tokens) < maxIdx
            continue;
        end

        f = parseNumber(tokens{freqIdx});
        zr = parseNumber(tokens{realIdx});
        zi = parseNumber(tokens{imagIdx});

        if ~isnan(f) && ~isnan(zr) && ~isnan(zi)
            if sourceWasNegativeImag
                zi = -zi;
            end

            data(end + 1, :) = [f, zr, zi]; %#ok<AGROW>
        end
    end

    if isempty(data)
        error('The required columns were found, but no numeric EIS rows could be read.');
    end
end


function [tokens, delimiterType] = splitCandidateHeader(line)
% Detect a likely column delimiter without assuming a specific instrument.

    if contains(line, sprintf('\t'))
        delimiterType = 'tab';
        tokens = regexp(line, '\t', 'split');

    elseif countCharacter(line, ',') >= 2
        delimiterType = 'comma';
        tokens = strsplit(line, ',');

    elseif countCharacter(line, ';') >= 2
        delimiterType = 'semicolon';
        tokens = strsplit(line, ';');

    else
        % Useful for files whose columns are separated by multiple spaces.
        delimiterType = 'multispace';
        tokens = regexp(strtrim(line), '\s{2,}', 'split');
    end
end


function tokens = splitDataLine(line, delimiterType)

    switch delimiterType
        case 'tab'
            tokens = regexp(line, '\t', 'split');

        case 'comma'
            tokens = strsplit(line, ',');

        case 'semicolon'
            tokens = strsplit(line, ';');

        otherwise
            tokens = regexp(strtrim(line), '\s+', 'split');
    end
end


function [kind, isNegativeImag] = classifyHeader(token)
% Convert many header styles to a common form and classify them.

    kind = '';
    isNegativeImag = false;

    token = strtrim(token);
    token = strrep(token, char(65279), '');
    token = lower(token);

    % Normalize common prime characters.
    token = strrep(token, char(8243), 'primeprime'); % Unicode double prime
    token = strrep(token, char(8242), 'prime');      % Unicode prime
    token = strrep(token, char(8217), 'prime');      % right single quote
    token = strrep(token, char(39), 'prime');        % apostrophe
    token = strrep(token, '"', 'primeprime');

    % Remove units or instrument suffixes such as (Hz), (Ohm), (a), (b).
    token = regexprep(token, '\([^)]*\)', '');
    token = regexprep(token, '\[[^\]]*\]', '');

    token = strtrim(token);

    if ~isempty(token) && token(1) == '-'
        isNegativeImag = true;
        token = token(2:end);
    elseif ~isempty(token) && token(1) == '+'
        token = token(2:end);
    end

    % Remove separators but preserve the normalized words.
    token = regexprep(token, '[\s_\-]', '');

    % Frequency
    if strcmp(token, 'freq') || strcmp(token, 'frequency') || ...
            startsWith(token, 'freq')
        kind = 'freq';
        return;
    end

    % Real impedance
    realAliases = { ...
        'zreal', ...
        'zre', ...
        'zprime', ...
        'impedanceprime', ...
        'realz', ...
        'rez'};

    if any(strcmp(token, realAliases))
        kind = 'real';
        return;
    end

    % Imaginary impedance
    imagAliases = { ...
        'zimag', ...
        'zimaginary', ...
        'zim', ...
        'zprimeprime', ...
        'impedanceprimeprime', ...
        'imagz', ...
        'imz'};

    if any(strcmp(token, imagAliases))
        kind = 'imag';
    end
end


function value = parseNumber(token)

    token = strtrim(token);
    token = strrep(token, '"', '');
    token = strrep(token, char(160), ''); % non-breaking space

    value = str2double(token);
end


function n = countCharacter(text, ch)
% Compatibility helper so the code does not depend on MATLAB count().

    n = sum(text == ch);
end
