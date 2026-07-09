function createLcrBedFile(fastaPath, outputDir, finalBedPath)
% createLcrBedFile
%
% Converts LCRFinder outputs into a BED-like file:
%
% ProteinID    Start    End
%
% Coordinates remain 1-based.
%
% Usage:
% createLcrBedFile( ...
%     'mouse_proteome.fasta', ...
%     'LCRFinder_Output', ...
%     'lcrfinder.bed');

    fprintf('========================================\n');
    fprintf('LCRFinder BED File Generator\n');
    fprintf('========================================\n');

    %% ----------------------------------------------------
    % Create output directory if needed
    %% ----------------------------------------------------
    outFolder = fileparts(finalBedPath);

    if ~isempty(outFolder) && ~exist(outFolder,'dir')
        mkdir(outFolder);
    end

    %% ----------------------------------------------------
    % Read FASTA headers
    %% ----------------------------------------------------
    fprintf('Reading FASTA file...\n');

    if ~exist(fastaPath,'file')
        error('FASTA file not found:\n%s', fastaPath);
    end

    [HEADER, ~] = fastaread(fastaPath);

    fprintf('Loaded %d FASTA entries\n', numel(HEADER));

    %% ----------------------------------------------------
    % Locate LCRFinder output files
    %% ----------------------------------------------------
    searchPattern = fullfile(outputDir,'*_');

    files = dir(searchPattern);

    files = files(~[files.isdir]);

    if isempty(files)
        error('No LCRFinder output files found in:\n%s', outputDir);
    end

    fprintf('Found %d LCRFinder files\n', length(files));

    %% ----------------------------------------------------
    % Open BED output
    %% ----------------------------------------------------
    fid_out = fopen(finalBedPath,'w');

    if fid_out == -1
        error('Cannot create output file:\n%s', finalBedPath);
    end

    %% ----------------------------------------------------
    % Process files
    %% ----------------------------------------------------
    totalRegions = 0;

    for f = 1:length(files)

        fname = files(f).name;

        % Example:
        % 14_  -> index 14

        idx = str2double(regexprep(fname,'_',''));

        if isnan(idx)
            warning('Skipping file: %s', fname);
            continue;
        end

        if idx < 1 || idx > numel(HEADER)
            warning('Index %d outside FASTA range. Skipping.', idx);
            continue;
        end

        %% Get protein ID from FASTA header

        if isstruct(HEADER)
            headerLine = HEADER(idx).Header;
        else
            headerLine = HEADER{idx};
        end

        % Extract UniProt accession
        %
        % >sp|Q9XYZ1|PROT_MOUSE ...
        %
        % becomes:
        % Q9XYZ1

        parts = split(headerLine,'|');

        if numel(parts) >= 3

            protein_id = strtrim(parts{2});

        else

            protein_id = strtok(headerLine,' ');

        end

        %% Open LCRFinder file

        fullFile = fullfile(outputDir,fname);

        fid = fopen(fullFile,'r');

        if fid == -1
            warning('Cannot open %s', fullFile);
            continue;
        end

        %% Read coordinates

        while ~feof(fid)

            line = fgetl(fid);

            if ~ischar(line)
                continue;
            end

            % Skip FASTA headers

            if startsWith(line,'>')
                continue;
            end

            % Example:
            %
            % (142)...(165)

            tokens = regexp( ...
                line, ...
                '\((\d+)\).*?\((\d+)\)', ...
                'tokens');

            if isempty(tokens)
                continue;
            end

            start_pos = str2double(tokens{1}{1});
            end_pos   = str2double(tokens{1}{2});

            fprintf( ...
                fid_out, ...
                '%s\t%d\t%d\n', ...
                protein_id, ...
                start_pos, ...
                end_pos);

            totalRegions = totalRegions + 1;

        end

        fclose(fid);

    end

    %% ----------------------------------------------------
    % Finish
    %% ----------------------------------------------------
    fclose(fid_out);

    fprintf('\n');
    fprintf('Finished successfully.\n');
    fprintf('Total LCRs written: %d\n', totalRegions);
    fprintf('Output BED file:\n%s\n', finalBedPath);
    fprintf('========================================\n');

end