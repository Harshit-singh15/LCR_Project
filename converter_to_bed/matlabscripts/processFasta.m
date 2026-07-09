function processFasta(fastaPath, targetFolder, outputName)
    % 1. Read the specific FASTA file
    [HEADER, SEQ] = fastaread(fastaPath);

    % 2. Path to where run_LCRFinder.m lives (stays constant)
    toolsFolder = '\\wsl.localhost\ubuntu\home\harshit_15\Lcr_Project\tools\lcr_2finder';
    addpath(toolsFolder);

    % 3. Create the unique target folder if it doesn't exist
    if ~exist(targetFolder, 'dir')
        mkdir(targetFolder);
    end

    % 4. Save current directory, then switch to the unique target folder
    originalFolder = pwd; 
    cd(targetFolder);

    % 5. Run the tool
    RES = run_LCRFinder(SEQ, 'AA', 'Medium');

    % 6. Switch back to your original working directory
    cd(originalFolder);

    % 7. Save inside the unique target folder with the unique filename
    save(fullfile(targetFolder, outputName), 'RES'); 
end