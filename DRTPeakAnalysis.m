function fig2csv_gui
% FIG -> CSV extractor with a tiny GUI.
%
% IMPORTANT: This tool is designed to work with outputs from DRTtools.
% Link: https://github.com/ciuccislab/DRTtools
%
% STEPS TO USE:
% 1. Perform deconvolution on γ(lnτ)-τ curve in the 'Peak Analysis' section of DRTtools.
% 2. Save the resulting figure using 'Export Results' in DRTtools.
% 3. Use this module to 'Load a single .fig file' OR a 'folder of .fig files'.
%
%   Features:
% - Extract all Line/Scatter series (XData/YData)
% - Preview curves (only the last loaded file's data is previewed)
% - Save CSV: Creates individual CSV files for the full data (all curves).
% - Save Parameters CSV: Creates CSV files for peak parameters (R, Tau, Capacitance, Frequency).
%
% Author: Masood Fakouri Hasanabadi
% Note: Any use of this code requires proper citation of the related publication.
    % --- App state container
    % app.fileData: cell array of structs, one per loaded file.
    % Each struct contains:
    %   .filename: Base name of the .fig file
    %   .series:   Cell array of extracted series: struct('name',string,'x',vec,'y',vec)
    %   .summary:  Peak parameters table (excluding Total)
    app.fileData = {};      % Stores data and summary for ALL loaded files
    app.data = [];          % Used temporarily for preview/single file load
    app.seriesCount = 0;    % Total series in the file being previewed/last loaded
    app.summary = table();  % Summary table for the file being previewed/last loaded
    
    % --- UI
    app.ui = uifigure('Name','FIG → CSV Exporter','Position',[100 100 980 560]);
    uilabel(app.ui,'Text','1) Load Figures/Folder  2) Save Full Data CSVs  3) Save Parameters CSVs',...
        'FontWeight','bold','Position',[20 525 700 22]);
    
    % NEW Button: Load Folder
    app.btnLoadFolder = uibutton(app.ui,'Text','Load Folder','Position',[150 490 120 30],...
        'ButtonPushedFcn',@onLoadFolder);
        
    app.btnLoad = uibutton(app.ui,'Text','Load .fig','Position',[20 490 120 30],...
        'ButtonPushedFcn',@onLoad);
    
    % Renamed/Modified Save Buttons
    app.btnSave = uibutton(app.ui,'Text','Save Full Data CSVs','Position',[280 490 150 30],...
        'ButtonPushedFcn',@onSave,'Enable','off');
    
    app.btnSaveParams = uibutton(app.ui,'Text','Save Params CSVs (P)','Position',[440 490 170 30],...
        'ButtonPushedFcn',@onSaveParams,'Enable','off');
    
    app.btnClear = uibutton(app.ui,'Text','Clear','Position',[620 490 90 30],...
        'ButtonPushedFcn',@onClear);
    
    % Axes (right side)
    app.ax = uiaxes(app.ui,'Position',[360 60 600 450]);
    title(app.ax,'Preview of extracted series');
    xlabel(app.ax,'\tau / s','Interpreter','tex','FontSize',14);
    ylabel(app.ax,'\gamma(\ln\tau) / \Omega','Interpreter','tex','FontSize',14);
    set(app.ax,'XScale','log'); % log scale on X axis
    grid(app.ax,'on'); hold(app.ax,'on');
    
    % Parameters table (left side)
    uilabel(app.ui,'Text','Peak Parameters (excluding Total):','Position',[20 460 320 20],...
        'FontWeight','bold');
    app.tbl = uitable(app.ui,'Position',[20 60 320 390],...
        'ColumnEditable',false(1,6));
    % --- Core Function to Load and Process a Single Figure ---
    function fileStruct = loadAndProcessFigure(fullPath, fn)
        fileStruct = struct('filename', fn, 'series', [], 'summary', table(), 'R_correction_ratio', 1.0); % Added R_correction_ratio field
        
        try
            % Try to open invisibly.
            try
                f = openfig(fullPath,'invisible');
            catch
                try
                    f = hgload(fullPath, struct('Visible','off'));
                catch ME
                    uialert(app.ui, sprintf('Could not open figure "%s":\n%s', fn, ME.message), 'Open Error');
                    return;
                end
            end
            
            % Find series: line & scatter
            lines = findall(f,'Type','line');
            scat  = findall(f,'Type','scatter');
            series = {};
            
            % Extract data
            for k = numel(lines):-1:1
                L = lines(k); x = get(L,'XData'); y = get(L,'YData');
                nm = getOrDefaultName(get(L,'DisplayName'), k, 'Line');
                series{end+1} = struct('name',nm,'x',x(:),'y',y(:)); %#ok<AGROW>
            end
            for k = numel(scat):-1:1
                S = scat(k); x = get(S,'XData'); y = get(S,'YData');
                nm = getOrDefaultName(get(S,'DisplayName'), k, 'Scatter');
                series{end+1} = struct('name',nm,'x',x(:),'y',y(:)); %#ok<AGROW>
            end
            
            try, close(f); end
            
            if isempty(series)
                warning('FIG2CSV:NoData', 'No Line/Scatter series found in figure: %s', fn);
                return;
            end
            
            % Compute summary
            [summary, R_correction_ratio] = computePeakSummary(series, true); % Get R_correction_ratio
            
            % Populate struct
            fileStruct.series = series;
            fileStruct.summary = summary;
            fileStruct.R_correction_ratio = R_correction_ratio; % Store the ratio
            
        catch ME
            uialert(app.ui, sprintf('Error processing %s: %s', fn, ME.message), 'Processing Error');
            fileStruct = struct('filename', fn, 'series', [], 'summary', table(), 'R_correction_ratio', 1.0); % Return empty on failure
        end
    end
    function updateUI(fileStruct)
        % Updates the UI elements based on the data from ONE fileStruct
        
        series = fileStruct.series;
        summary = fileStruct.summary;
        R_corr_ratio = fileStruct.R_correction_ratio; % Get the ratio
        
        % Update state for preview/table display
        app.data = series;
        app.seriesCount = numel(series);
        
        % --- START MODIFICATION for summary table ---
        if height(summary) > 0
            % Calculate corrected R_Ohm values for display
            summary.corrected_R_Ohm = summary.R_Ohm * R_corr_ratio;
            % Select and reorder columns for display (keeping original R_Ohm and adding corrected)
            summaryDisplay = summary(:, {'ID','Series_Name','R_Ohm','corrected_R_Ohm','Tau_s','Capacitance_F','Frequency_Hz'});
            app.summary = summaryDisplay; % Use the enhanced table for UI
        else
            app.summary = table();
        end
        % --- END MODIFICATION ---
        
        % Preview
        cla(app.ax); hold(app.ax,'on');
        if app.seriesCount >= 1
            plot(app.ax, app.data{1}.x, app.data{1}.y, 'r-', 'LineWidth', 2.2, 'DisplayName', sprintf('%s (Total)', app.data{1}.name));
        end
        for i = 2:app.seriesCount
            plot(app.ax, app.data{i}.x, app.data{i}.y, 'LineWidth', 1.2, 'DisplayName', app.data{i}.name);
        end
        legend(app.ax,'show','Interpreter','none','Location','best');
        title(app.ax, sprintf('Extracted %d series from: %s', app.seriesCount, fileStruct.filename));
        drawnow;
        
        % Populate table
        if height(app.summary) > 0
            % 1. Sort and re-index (This ensures the UI table is sorted)
            % Sort by the original Frequency_Hz, which is present in the full summary table (if available), 
            % but since app.summary is the modified table, we use its current sort.
            % We will rely on the sort done inside loadAndProcessFigure for the full summary. 
            % For display, let's sort based on Frequency_Hz which is retained.
            app.summary = sortrows(app.summary, 'Frequency_Hz', 'descend');
            app.summary.ID = (1:height(app.summary)).'; 
            
            % 2. Prepare data for display
            % Updated to include corrected_R_Ohm column
            displayTableData = cell(height(app.summary), width(app.summary));
            displayTableData(:, 1) = num2cell(app.summary.ID);
            displayTableData(:, 2) = cellstr(app.summary.Series_Name);
            % Updated numeric columns to include 'corrected_R_Ohm'
            numericCols = {'R_Ohm', 'corrected_R_Ohm', 'Tau_s', 'Capacitance_F', 'Frequency_Hz'};
            
            for colIdx = 1:length(numericCols)
                varName = numericCols{colIdx};
                formattedStrings = arrayfun(@(x) sprintf('%.4e', x), app.summary.(varName), 'UniformOutput', false);
                displayTableData(:, colIdx + 2) = formattedStrings;
            end
            
            app.tbl.Data = displayTableData;
            % Update column names for display
            colNames = {'ID','Series_Name','R_Ohm','corrected R_Ohm','Tau_s','Capacitance_F','Frequency_Hz'};
            app.tbl.ColumnName = colNames;
            app.tbl.ColumnWidth = {'auto','auto','auto','auto','auto','auto','auto'}; % Added one 'auto'
            app.btnSave.Enable = 'on';
            app.btnSaveParams.Enable = 'on';
        else
            app.tbl.Data = {};
            app.tbl.ColumnName = {};
            app.btnSave.Enable = 'off';
            app.btnSaveParams.Enable = 'off';
        end
    end
    
    % --- Callbacks
    
    function onLoad(~,~)
        % Load single file: clears app.fileData first
        onClear(); 
        [fn,fp] = uigetfile({'*.fig','MATLAB Figure (*.fig)'},'Select a figure');
        if isequal(fn,0); return; end
        full = fullfile(fp,fn);
        
        fileStruct = loadAndProcessFigure(full, fn);
        
        if ~isempty(fileStruct.series)
            app.fileData = {fileStruct};
            updateUI(fileStruct);
        end
    end
    function onLoadFolder(~,~)
        % Load folder of files: clears app.fileData first
        onClear(); 
        folderPath = uigetdir(pwd,'Select a folder containing .fig files');
        if isequal(folderPath,0); return; end
        
        figFiles = dir(fullfile(folderPath, '*.fig'));
        if isempty(figFiles)
            uialert(app.ui, 'No .fig files found in the selected folder.', 'No Files');
            return;
        end
        
        loadedCount = 0;
        newFileData = {};
        for i = 1:numel(figFiles)
            fn = figFiles(i).name;
            full = fullfile(folderPath, fn);
            
            fileStruct = loadAndProcessFigure(full, fn);
            
            if ~isempty(fileStruct.series)
                newFileData{end+1} = fileStruct; %#ok<AGROW>
                loadedCount = loadedCount + 1;
            end
        end
        
        if loadedCount > 0
            app.fileData = newFileData;
            % Update UI with the last loaded file's data
            updateUI(newFileData{end});
        else
             uialert(app.ui, 'Failed to load any series from the .fig files.', 'Load Error');
        end
    end
    
    function onSave(~,~)
        % Creates individual CSV files for the full data (all curves).
        if isempty(app.fileData)
            uialert(app.ui,'No data loaded. Click "Load .fig" or "Load Folder" first.','Nothing to save');
            return;
        end
        
        savePath = uigetdir(pwd,'Select a folder to save individual Full Data CSV files');
        if isequal(savePath,0); return; end
        successCount = 0;
        for i = 1:numel(app.fileData)
            fileStruct = app.fileData{i};
            
            % Determine the output file name (replace .fig with .csv)
            [~, name, ~] = fileparts(fileStruct.filename);
            saveName = fullfile(savePath, [name '.csv']);
            
            % Build a padded matrix and headers
            % PASS THE CORRECTION RATIO TO THE FUNCTION
            [M, headers] = buildPaddedMatrix(fileStruct.series, fileStruct.R_correction_ratio);
            
            try
                T = array2table(M, 'VariableNames', matlab.lang.makeValidName(headers));
                writetable(T, saveName);
                successCount = successCount + 1;
            catch ME
                warning('FIG2CSV:SaveError', 'Failed to save %s: %s', saveName, ME.message);
            end
        end
        
        if successCount == 0
            uialert(app.ui, 'Failed to save any Full Data CSV files.', 'Save Error','Icon','error');
        end
    end
    
    function onSaveParams(~,~)
        % MODIFIED: Creates individual CSV files for the parameters with 'P' suffix, sorted by Frequency_Hz (high to low).
        if isempty(app.fileData)
            uialert(app.ui,'No data loaded. Click "Load .fig" or "Load Folder" first.','Nothing to save');
            return;
        end
        
        savePath = uigetdir(pwd,'Select a folder to save individual Parameters CSV files');
        if isequal(savePath,0); return; end
        successCount = 0;
        for i = 1:numel(app.fileData)
            fileStruct = app.fileData{i};
            
            % --- START MODIFICATION to get summary for saving ---
            % Recompute the summary table with the corrected R_Ohm values for saving
            [fullSummary, R_corr_ratio] = computePeakSummary(fileStruct.series, true);
            
            if height(fullSummary) > 0
                % 1. Get the current summary table
                summaryTable = fullSummary;
                
                % 2. Calculate the corrected R_Ohm and add to the table
                summaryTable.corrected_R_Ohm = summaryTable.R_Ohm * R_corr_ratio;
                
                % 3. Sort the table by Frequency_Hz (high to low)
                summaryTable = sortrows(summaryTable, 'Frequency_Hz', 'descend');
                
                % 4. Re-index the 'ID' column after sorting (for clarity in the file)
                summaryTable.ID = (1:height(summaryTable)).';
                
                % 5. Reorder columns to place corrected R_Ohm next to R_Ohm
                summaryTable = summaryTable(:, {'ID','Series_Name','R_Ohm','corrected_R_Ohm','Tau_s','Capacitance_F','Frequency_Hz'});
                
                % --- START NEW MODIFICATION: Add Total Row ---
                % Calculate Sums
                sum_R_Ohm = sum(summaryTable.R_Ohm, 'omitnan');
                sum_corrected_R_Ohm = sum(summaryTable.corrected_R_Ohm, 'omitnan');
                
                % Create the Total row
                T_total = table(height(summaryTable) + 1, ...
                                "Total Sum", ...
                                sum_R_Ohm, ...
                                sum_corrected_R_Ohm, ...
                                NaN, NaN, NaN, ...
                                'VariableNames', {'ID','Series_Name','R_Ohm','corrected_R_Ohm','Tau_s','Capacitance_F','Frequency_Hz'});

                % Append the Total row
                summaryTable = [summaryTable; T_total];
                % --- END NEW MODIFICATION ---
                
                % Determine the output file name (e.g. 'figure_nameP.csv')
                [~, name, ~] = fileparts(fileStruct.filename);
                saveName = fullfile(savePath, [name 'P.csv']);
                
                try
                    % Write the sorted table (now including the total row)
                    writetable(summaryTable, saveName);
                    successCount = successCount + 1;
                catch ME
                    warning('FIG2CSV:ParamSaveError', 'Failed to save parameters for %s: %s', fileStruct.filename, ME.message);
                end
            end
        end
        
        if successCount == 0
            uialert(app.ui, 'Failed to save any Parameters CSV files.', 'Save Error','Icon','error');
        end
    end
    
    function onClear(~,~)
        app.fileData = {}; % Clear all loaded files
        app.data = [];
        app.seriesCount = 0;
        app.summary = table();
        cla(app.ax); legend(app.ax,'off');
        title(app.ax,'Preview of extracted series');
        app.tbl.Data = {};
        app.tbl.ColumnName = {};
        app.btnSave.Enable = 'off';
        app.btnSaveParams.Enable = 'off';
    end
end
% ---------- helpers (MODIFIED) ----------
function nm = getOrDefaultName(dispName, k, prefix)
    if isstring(dispName) || (ischar(dispName) && ~isempty(dispName))
        s = string(dispName);
        if strlength(s) > 0
            nm = char(s);
            return;
        end
    end
    nm = sprintf('%s_%d', prefix, k);
end

% The buildPaddedMatrix function remains logically unchanged as it handles X/Y data, not the R parameters.
function [M, headers] = buildPaddedMatrix(seriesCell, R_correction_ratio)
    if nargin < 2, R_correction_ratio = 1.0; end % Unused in this function
    
    n = numel(seriesCell);
    lens = cellfun(@(s) numel(s.x), seriesCell);
    L = max(lens);
    % Preallocate with NaNs
    M = nan(L, 2*n);
    headers = cell(1, 2*n);
    for i = 1:n
        xi = seriesCell{i}.x;
        yi = seriesCell{i}.y;
        % Ensure column vectors & finite values aligned
        mask = isfinite(xi) & isfinite(yi);
        xi = xi(mask); yi = yi(mask);
        m = numel(xi);
        M(1:m, 2*i-1) = xi;
        M(1:m, 2*i)   = yi;
        base = seriesCell{i}.name;
        headers{2*i-1} = sprintf('%s_X', base);
        headers{2*i}   = sprintf('%s_Y', base);
    end
end

% MODIFIED: Function signature and logic updated to return correction ratio
function [T, R_correction_ratio] = computePeakSummary(seriesCell, ignoreFirst)
% Computes, for each series (assumed tau=X, gamma(ln_tau)=Y):
%   R = INT [gamma(ln_tau) d(ln_tau)]
%   tau_peak = x at max(gamma)
%   f = 1/(2*pi*tau_peak)
%   C = tau_peak / R
% If ignoreFirst==true, skip seriesCell{1} (treat as Total).

    if nargin < 2, ignoreFirst = false; end
    startIdx = 1 + double(ignoreFirst);
    nTot = numel(seriesCell);
    
    % --- START MODIFICATION for R correction (steps 1, 2, 3) ---
    R_total = 0;
    R_cell_1 = 0;
    
    % First, pre-calculate all R values (Total and individual)
    all_R = nan(nTot, 1);
    for i = 1:nTot
        xi = seriesCell{i}.x;
        yi = seriesCell{i}.y;
        mask = isfinite(xi) & isfinite(yi) & (xi > 0);
        xi = xi(mask);
        yi = yi(mask);
        if numel(xi) >= 2
            [xs, idx] = sort(xi);
            ys = yi(idx);
            ln_xs = log(xs); 
            all_R(i) = trapz(ln_xs, ys); % R = INT [gamma(ln_tau) d(ln_tau)]
        end
    end
    
    % 1. Calculate the sum of all R as R_total (excluding seriesCell{1})
    R_total = sum(all_R(startIdx:end), 'omitnan');

    % 2. calculate the area under the curve seriesCell{1} as R_cell.
    if nTot >= 1
        R_cell_1 = all_R(1); 
    end
    
    % 3. Divide R-cell by R_total to obtain correction ratio.
    R_correction_ratio = 1.0;
    if R_total ~= 0 && isfinite(R_cell_1) && isfinite(R_total)
        R_correction_ratio = R_cell_1 / R_total;
    end
    % --- END MODIFICATION ---
    
    
    n = max(0, nTot - (startIdx-1));
    ID = (1:n).';
    Series_Name = strings(n,1);
    R_Ohm = nan(n,1);
    Tau_s = nan(n,1);
    Cap_F = nan(n,1);
    Freq_Hz = nan(n,1);
    row = 0;
    for i = startIdx:nTot
        row = row + 1;
        R = all_R(i); % Use the pre-calculated R value
        
        xi = seriesCell{i}.x;
        yi = seriesCell{i}.y;
        mask = isfinite(xi) & isfinite(yi) & (xi > 0);
        xi = xi(mask);
        yi = yi(mask);
        
        Series_Name(row) = string(seriesCell{i}.name);
        
        if numel(xi) >= 2
            [xs, idx] = sort(xi); % Need to sort for peak finding
            ys = yi(idx);
            
            % R_Ohm calculation is done outside the loop; just assign it.
            R_Ohm(row) = R; 
            
            % Peak tip location (tau)
            [~, imax] = max(ys);
            tau_pk = xs(imax);
            Tau_s(row) = tau_pk;
            
            % Derived
            if isfinite(tau_pk) && tau_pk > 0
                Freq_Hz(row) = 1/(2*pi*tau_pk);
            end
            if isfinite(R) && R ~= 0
                Cap_F(row) = tau_pk / R;
            end
        end
    end
    
    % T still contains the UNCORRECTED R_Ohm values for now, but 
    % R_correction_ratio is returned to be used elsewhere.
    T = table(ID, Series_Name, R_Ohm, Tau_s, Cap_F, Freq_Hz, ...
              'VariableNames', {'ID','Series_Name','R_Ohm','Tau_s','Capacitance_F','Frequency_Hz'});
end