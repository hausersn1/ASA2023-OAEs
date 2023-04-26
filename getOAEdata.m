subjects = {'SH_R'};


mainDir = pwd;
%% Get result file for each one.
for s = 1:length(subjects)
    subj = string(subjects(s));
    
    %% Run DPOAE analysis
    folder = ['sweptDPOAE-main/Results/' subj];
    cd(folder)
    
    file = sprintf('DPOAEswept_%s_*.mat', subj);
    windowdur = 0.04;
    offsetwin = 0.0;
    npoints = 1024;
    
    DP = DPanalysis(file, windowdur, offsetwin, npoints);
    
    cd(mainDir)
    
    %% Run SFOAE analysis
    folder = ['sweptDPOAE-main/Results/' subj];
    cd(folder)
    
    file = sprintf('DPOAEswept_%s_%s_*.mat', subj, ear);
    windowdur = 0.04;
    offsetwin = 0.0;
    npoints = 1024;
    SF = DPanalysis(input, windowdur, offsetwin, npoints);
    
    %% Get Calibration data
    
    res.SF = SF;
    res.DP = DP;
    res.calib = calib;
end
end