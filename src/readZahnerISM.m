function data = readZahnerISM(fullFileName)
%READZAHNERISM Read a Zahner Thales .ism impedance file without Python.
%
% Output columns:
%   1 = Frequency / Hz
%   2 = Zreal / Ohm
%   3 = Zimag / Ohm
%
% The ISM data arrays are stored as big-endian values. Zahner files may
% contain an overlapping part of the frequency sweep; this reader keeps the
% largest non-overlapping range, matching the normal Thales/Zahner import
% behaviour.

    fid = fopen(fullFileName, 'r', 'ieee-be');
    if fid == -1
        error('Could not open the selected Zahner ISM file.');
    end
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    % File version is stored as a signed 6-byte integer. It is not required
    % for extracting the impedance arrays, but reading it advances the file
    % pointer to the sample count.
    readInt48BE(fid);

    numberOfSamples = readInt48BE(fid) + 1;
    if numberOfSamples <= 0 || numberOfSamples > 1e7
        error('Invalid Zahner ISM sample count.');
    end

    frequency = fread(fid, numberOfSamples, 'double=>double', 0, 'ieee-be');
    magnitude = fread(fid, numberOfSamples, 'double=>double', 0, 'ieee-be');
    phase = fread(fid, numberOfSamples, 'double=>double', 0, 'ieee-be');

    if numel(frequency) ~= numberOfSamples || ...
            numel(magnitude) ~= numberOfSamples || ...
            numel(phase) ~= numberOfSamples
        error('Unexpected end of Zahner ISM file while reading EIS data.');
    end

    % ISM phase is stored in radians.
    zReal = magnitude .* cos(phase);
    zImag = magnitude .* sin(phase);

    % Some ISM sweeps include an overlapping frequency section. Keep the
    % continuous range between the global frequency extrema and orient it
    % from high to low frequency.
    [~, minIndex] = min(frequency);
    [~, maxIndex] = max(frequency);

    firstIndex = min(minIndex, maxIndex);
    lastIndex = max(minIndex, maxIndex);

    frequency = frequency(firstIndex:lastIndex);
    zReal = zReal(firstIndex:lastIndex);
    zImag = zImag(firstIndex:lastIndex);

    if minIndex > maxIndex
        frequency = flipud(frequency);
        zReal = flipud(zReal);
        zImag = flipud(zImag);
    end

    data = [frequency(:), zReal(:), zImag(:)];
end


function value = readInt48BE(fid)
% Read a signed 48-bit big-endian integer.

    bytes = fread(fid, 6, '*uint8');
    if numel(bytes) ~= 6
        error('Unexpected end of Zahner ISM file.');
    end

    unsignedValue = 0;
    for k = 1:6
        unsignedValue = unsignedValue * 256 + double(bytes(k));
    end

    if bitand(bytes(1), uint8(128)) ~= 0
        value = unsignedValue - 2^48;
    else
        value = unsignedValue;
    end
end
