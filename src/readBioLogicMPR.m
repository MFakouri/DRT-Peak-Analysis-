function data = readBioLogicMPR(fullFileName)
%READBIOLOGICMPR Read BioLogic EC-Lab .mpr EIS data without Python.
%
% Output columns:
%   1 = Frequency / Hz
%   2 = Zreal / Ohm
%   3 = Zimag / Ohm
%
% This reader targets EC-Lab impedance data modules and extracts the native
% EIS columns (frequency, real impedance and -imaginary impedance). If the
% real/imaginary columns are absent but magnitude and phase are present,
% they are used as a fallback.

    fid = fopen(fullFileName, 'r', 'ieee-le');
    if fid == -1
        error('Could not open the selected BioLogic MPR file.');
    end
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    magic = fread(fid, 52, '*uint8')';
    if numel(magic) ~= 52 || ...
            ~isequal(magic(1:22), uint8('BIO-LOGIC MODULAR FILE')) || ...
            magic(23) ~= 26
        error('The selected file is not a valid BioLogic MPR file.');
    end

    dataModuleFound = false;
    dataModuleBytes = [];
    dataModuleVersion = [];

    while ~feof(fid)
        moduleMagic = fread(fid, 6, '*char')';
        if isempty(moduleMagic)
            break;
        end
        if numel(moduleMagic) ~= 6 || ~strcmp(moduleMagic, 'MODULE')
            error('Invalid BioLogic MPR module structure.');
        end

        header = fread(fid, 51, '*uint8')';
        if numel(header) ~= 51
            error('Unexpected end of BioLogic MPR module header.');
        end

        shortName = char(header(1:10));
        shortName(shortName == 0) = ' ';
        shortName = strtrim(shortName);

        % EC-Lab >= 11.50 uses the extended 59-byte module header. In the
        % short header view, bytes 36:39 contain 0xFFFFFFFF as a marker.
        if all(header(36:39) == 255)
            extra = fread(fid, 8, '*uint8')';
            if numel(extra) ~= 8
                error('Unexpected end of BioLogic MPR extended module header.');
            end
            header = [header, extra]; %#ok<AGROW>
            moduleLength = readU32LE(header(40:43));
            moduleVersion = readU32LE(header(44:47));
        else
            moduleLength = readU32LE(header(36:39));
            moduleVersion = readU32LE(header(40:43));
        end

        if moduleLength < 0
            error('Invalid BioLogic MPR module length.');
        end

        if strcmp(shortName, 'VMP data')
            dataModuleBytes = fread(fid, moduleLength, '*uint8')';
            if numel(dataModuleBytes) ~= moduleLength
                error('Unexpected end of BioLogic MPR data module.');
            end
            dataModuleVersion = moduleVersion;
            dataModuleFound = true;
            break;
        else
            status = fseek(fid, moduleLength, 'cof');
            if status ~= 0
                error('Could not seek to the next BioLogic MPR module.');
            end
        end
    end

    if ~dataModuleFound
        error('BioLogic MPR file does not contain a VMP data module.');
    end

    data = parseMPRDataModule(dataModuleBytes, dataModuleVersion);
end


function data = parseMPRDataModule(bytes, version)

    if numel(bytes) < 6
        error('BioLogic MPR data module is incomplete.');
    end

    numberOfPoints = readU32LE(bytes(1:4));
    numberOfColumns = double(bytes(5));

    if numberOfPoints <= 0 || numberOfColumns <= 0
        error('BioLogic MPR data module contains no EIS data.');
    end

    switch version
        case {2, 3}
            requiredHeaderBytes = 5 + 2 * numberOfColumns;
            if numel(bytes) < requiredHeaderBytes
                error('BioLogic MPR column header is incomplete.');
            end

            columnIds = zeros(1, numberOfColumns);
            p = 6;
            for k = 1:numberOfColumns
                columnIds(k) = readU16LE(bytes(p:p+1));
                p = p + 2;
            end

            if version == 3
                mainDataOffset = 406;
            else
                mainDataOffset = 405;
            end

        case 0
            % Older data modules use one-byte column identifiers. Some
            % EC-Lab versions interleave each identifier with a zero byte.
            if bytes(6) ~= 0
                lastHeaderByte = 5 + numberOfColumns;
                if numel(bytes) < lastHeaderByte
                    error('BioLogic MPR column header is incomplete.');
                end
                columnIds = double(bytes(6:lastHeaderByte));
                mainDataOffset = 100;
            else
                lastHeaderByte = 5 + 2 * numberOfColumns;
                if numel(bytes) < lastHeaderByte
                    error('BioLogic MPR column header is incomplete.');
                end
                rawIds = bytes(6:lastHeaderByte);
                columnIds = double(rawIds(2:2:end));
                mainDataOffset = 1007;
            end

        otherwise
            error('Unsupported BioLogic MPR data-module version: %d', version);
    end

    if numel(bytes) <= mainDataOffset
        error('BioLogic MPR data module does not contain a data array.');
    end

    payloadBytes = numel(bytes) - mainDataOffset;
    recordSize = payloadBytes / numberOfPoints;
    if abs(recordSize - round(recordSize)) > 1e-9
        error('BioLogic MPR data records have an unexpected size.');
    end
    recordSize = round(recordSize);

    freqIndex = find(columnIds == 32, 1, 'first');
    realIndex = find(columnIds == 37, 1, 'first');
    negImagIndex = find(columnIds == 38, 1, 'first');

    useMagnitudePhase = false;
    if isempty(realIndex) || isempty(negImagIndex)
        magnitudeIndex = find(columnIds == 36, 1, 'first');
        phaseIndex = find(columnIds == 35, 1, 'first');
        if ~isempty(magnitudeIndex) && ~isempty(phaseIndex)
            useMagnitudePhase = true;
        else
            error(['BioLogic MPR file does not contain the required EIS ', ...
                   'impedance columns (Re(Z), -Im(Z), or |Z| and Phase(Z)).']);
        end
    else
        magnitudeIndex = [];
        phaseIndex = [];
    end

    if isempty(freqIndex)
        error('BioLogic MPR file does not contain the frequency EIS column.');
    end

    neededIndices = freqIndex;
    if useMagnitudePhase
        neededIndices = [neededIndices, magnitudeIndex, phaseIndex];
    else
        neededIndices = [neededIndices, realIndex, negImagIndex];
    end

    columnOffsets = computeColumnOffsets(columnIds, max(neededIndices));

    freqOffset = columnOffsets(freqIndex);
    if useMagnitudePhase
        magnitudeOffset = columnOffsets(magnitudeIndex);
        phaseOffset = columnOffsets(phaseIndex);
    else
        realOffset = columnOffsets(realIndex);
        negImagOffset = columnOffsets(negImagIndex);
    end

    mainStart = mainDataOffset + 1;
    frequency = zeros(numberOfPoints, 1);
    zReal = zeros(numberOfPoints, 1);
    zImag = zeros(numberOfPoints, 1);

    for row = 1:numberOfPoints
        recordStart = mainStart + (row - 1) * recordSize;

        frequency(row) = readF32LE(bytes(recordStart + freqOffset : ...
                                               recordStart + freqOffset + 3));

        if useMagnitudePhase
            magnitude = readF32LE(bytes(recordStart + magnitudeOffset : ...
                                             recordStart + magnitudeOffset + 3));
            phaseDeg = readF32LE(bytes(recordStart + phaseOffset : ...
                                            recordStart + phaseOffset + 3));
            zReal(row) = magnitude * cosd(phaseDeg);
            zImag(row) = magnitude * sind(phaseDeg);
        else
            zReal(row) = readF32LE(bytes(recordStart + realOffset : ...
                                             recordStart + realOffset + 3));
            storedNegativeImag = readF32LE(bytes(recordStart + negImagOffset : ...
                                                      recordStart + negImagOffset + 3));
            % EC-Lab stores the Nyquist ordinate -Im(Z). DRTtools expects
            % the actual imaginary component Zimag.
            zImag(row) = -storedNegativeImag;
        end
    end

    data = [frequency, zReal, zImag];
end


function offsets = computeColumnOffsets(columnIds, lastNeededIndex)
% Compute byte offsets within one EC-Lab data record up to the last EIS
% column needed by this importer.

    offsets = nan(size(columnIds));
    currentOffset = 0;
    flagByteAlreadyAdded = false;

    for k = 1:lastNeededIndex
        offsets(k) = currentOffset;
        [width, isPackedFlag, known] = mprColumnWidth(columnIds(k));

        if ~known
            error(['Unsupported BioLogic MPR column ID %d occurs before ', ...
                   'the required EIS columns.'], columnIds(k));
        end

        if isPackedFlag
            if ~flagByteAlreadyAdded
                currentOffset = currentOffset + 1;
                flagByteAlreadyAdded = true;
            end
        else
            currentOffset = currentOffset + width;
        end
    end
end


function [width, isPackedFlag, known] = mprColumnWidth(columnId)
% Byte widths for EC-Lab columns commonly encountered before/within EIS
% data. Flag columns share one packed byte.

    isPackedFlag = any(columnId == [1 2 3 21 31 65]);
    if isPackedFlag
        width = 0;
        known = true;
        return;
    end

    % 8-byte floating-point columns.
    eightByteIds = [4 7 11 13 23 24 74 123 124 125 126 182 211 ...
                    438 467 498 499 500 501 502];

    % 2-byte unsigned integer columns.
    twoByteIds = [39 131];

    % 1-byte scalar columns.
    oneByteIds = 509;

    % Known 4-byte columns used by EC-Lab, including the impedance fields.
    fourByteIds = [5 6 8 9 16 17 19 20 26 27 32 33 34 35 36 37 38 ...
                   69 70 75 76 77 78 96 98 99 100 101 163 168 169 ...
                   172 173 174 178 179 212 213 217 218 220 221 223 224 ...
                   230 231 232 233 234 235 236 237 238 239 240 241 242 ...
                   271 272 301 302 331 332 361 362 391 392 422 423 424 ...
                   425 426 430 431 432 433 434 435 441 462 468 469 471 ...
                   473 474 476 477 479 480 486 487 488 489 490 491 492 ...
                   493 494 495 496 497 505];

    if any(columnId == eightByteIds)
        width = 8;
        known = true;
    elseif any(columnId == fourByteIds)
        width = 4;
        known = true;
    elseif any(columnId == twoByteIds)
        width = 2;
        known = true;
    elseif any(columnId == oneByteIds)
        width = 1;
        known = true;
    else
        width = 0;
        known = false;
    end
end


function value = readU16LE(bytes)
    value = double(bytes(1)) + 256 * double(bytes(2));
end


function value = readU32LE(bytes)
    value = double(bytes(1)) + ...
            256 * double(bytes(2)) + ...
            65536 * double(bytes(3)) + ...
            16777216 * double(bytes(4));
end


function value = readF32LE(bytes)
% Convert four little-endian bytes to IEEE-754 single, independent of the
% host machine byte order used by MATLAB.

    bytes = uint8(bytes(:)');
    if ~isLittleEndianHost()
        bytes = fliplr(bytes);
    end
    value = double(typecast(bytes, 'single'));
end


function tf = isLittleEndianHost()
    persistent cachedValue
    if isempty(cachedValue)
        testValue = typecast(uint16(1), 'uint8');
        cachedValue = (testValue(1) == 1);
    end
    tf = cachedValue;
end
