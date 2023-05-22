function [res] = DPanalysis(subj, ear)

% Analyze swept tone DPOAE data using least-squares fit of chirp model

%%%%%%%%% Set these parameters %%%%%%%%%%%%%%%%%%

windowdur = 0.25;
offsetwin = 0.0; % not finding additional delay
npoints = 512;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load in the data

maindir = pwd;

% Select correct subject file for chins or humans
if subj(1) == 'Q'
    datadir = [maindir '/sweptDPOAE/Results/'];
else
    datadir = [maindir '/sweptDPOAE-main/Results/'];
end
subjfile = sprintf('DPOAEswept_%s_%s*.mat', subj, ear);
subjdir = [datadir char(subj)];
cd(subjdir)

checkDIR=dir(subjfile);
if isempty(checkDIR)
    fprintf('No file for subject %s %s ear\n', subj, ear)
    cd(maindir)
    return
elseif size(checkDIR,1) > 1
    fprintf('Too many files for subject %s %s ear', subj, ear)
    cd(maindir)
    return
    % error('No such OAEs for subject %s',subj);
end

load(checkDIR.name);
cd(maindir)

%% Set variables from the stim
phi1_inst = 2 * pi * stim.phi1_inst;
phi2_inst = 2 * pi * stim.phi2_inst;
phi_dp_inst = (2.*stim.phi1_inst - stim.phi2_inst) * 2 * pi;
rdp = 2 / stim.ratio - 1;    % f_dp = f2 * rdp

trials = size(stim.resp,1); 

t = stim.t;
if stim.speed < 0 % downsweep
    f_start = stim.fmax;
    f_end = stim.fmin;
else
    f_start = stim.fmin;
    f_end = stim.fmax;
end

% set freq we're testing and the timepoints when they happen.
if abs(stim.speed) < 20         % then in octave scaling
    freq_f2 = 2 .^ linspace(log2(f_start), log2(f_end), npoints);
    freq_f1 = freq_f2 ./ stim.ratio;
    freq_dp = 2.*freq_f1 - freq_f2;
    t_freq = log2(freq_f2/f_start)/stim.speed + stim.buffdur;
else                            % otherwise linear scaling
    freq_f2 = linspace(f_start, f_end, npoints);
    freq_f1 = freq_f2 ./ stim.ratio;
    freq_dp = 2.*freq_f1 - freq_f2;
    t_freq = (freq_f2-f_start)/stim.speed + stim.buffdur;
end

%% Artifact Rejection
% high pass filter the response (can also be done on ER10X hardware) 
%filtcutoff = 300;
%b = fir1(1000, filtcutoff*2/stim.Fs, 'high');
%DPOAEtrials= filtfilt(b, 1, stim.resp')';
DPOAEtrials = stim.resp;

% Set empty matricies for next steps
coeffs = zeros(npoints, 6);
a_temp = zeros(trials, npoints);
b_temp = zeros(trials, npoints);

% Least Squares fit of DP Only for AR
for x = 1:trials
    DPOAE = DPOAEtrials(x, :);
    fprintf(1, 'Checking trial %d / %d for artifact\n', x, (trials));
    
    for k = 1:npoints
        win = find( (t > (t_freq(k) - windowdur/2)) & ...
            (t < (t_freq(k) + windowdur/2)));
        taper = hanning(numel(win))';
        
        model_dp = [cos(phi_dp_inst(win)) .* taper;
            -sin(phi_dp_inst(win)) .* taper];
        
        resp = DPOAE(win) .* taper;
        
        coeffs(k, 1:2) = model_dp' \ resp';
    end
    a_temp(x,:) = coeffs(:, 1);
    b_temp(x,:) = coeffs(:, 2);
end

oae = abs(complex(a_temp, b_temp));

median_oae = median(oae);
std_oae = std(oae);
resp_AR = DPOAEtrials;
for j = 1:trials
    for k = 1:npoints
        if oae(j,k) > median_oae(1,k) + 4*std_oae(1,k)
            win = find( (t > (t_freq(k) - windowdur.*.1)) & ...
                (t < (t_freq(k) + windowdur.*.1)));
            resp_AR(j,win) = NaN;
        end
    end
end

%% Calculate Noise Floor

% First way to calculate noise floor, just subtracting alternate trials
numOfTrials = floor(trials/2)*2; % need even number of trials
noise = zeros(numOfTrials/2, size(resp_AR, 2));
for x = 1:2:numOfTrials
    noise(ceil(x/2),:) = (resp_AR(x,:) - resp_AR(x+1,:)) / 2;
end

DPOAE = mean(resp_AR, 1, "omitNaN");
NOISE = mean(noise,1, "omitNaN");


%% LSF analysis

% Set empty matricies for next steps
maxoffset = ceil(stim.Fs * offsetwin);
coeffs = zeros(npoints, 2);
coeffs_n = zeros(npoints, 2);
tau_dp = zeros(npoints, 1); % delay if offset > 0
coeffs_noise = zeros(npoints,8);
%durs = -.5*(2.^(0.003*t_freq)-1)/ (0.003*log(2)) + 0.5; 


% Least Squares fit of Chirp model (stimuli, DP, noise two ways)
for k = 1:npoints
    
    fprintf(1, 'Running window %d / %d\n', k, npoints);
    %windowdur = durs(k); 
    
    win = find( (t > (t_freq(k) - windowdur/2)) & ...
        (t < (t_freq(k) + windowdur/2)));
    taper = hanning(numel(win))';
    
    % set the models
    model_dp = [cos(phi_dp_inst(win)) .* taper;
        -sin(phi_dp_inst(win)) .* taper];
    model_f1 = [cos(phi1_inst(win)) .* taper;
        -sin(phi1_inst(win)) .* taper];
    model_f2 = [cos(phi2_inst(win)) .* taper;
        -sin(phi2_inst(win)) .* taper];
    model_noise = ...
        [cos(0.9*phi_dp_inst(win)) .* taper;
        -sin(0.9*phi_dp_inst(win)) .* taper;
        cos(0.88*phi_dp_inst(win)) .* taper;
        -sin(0.88*phi_dp_inst(win)) .* taper;
        cos(0.86*phi_dp_inst(win)) .* taper;
        -sin(0.86*phi_dp_inst(win)) .* taper;
        cos(0.84*phi_dp_inst(win)) .* taper;
        -sin(0.84*phi_dp_inst(win)) .* taper];
    
    % zero out variables for offset calc
    coeff = zeros(maxoffset, 6);
    coeff_n = zeros(maxoffset, 6);
    resid = zeros(maxoffset, 3);
    
    for offset = 0:maxoffset
        resp = DPOAE(win+offset) .* taper;
        resp_n = NOISE(win+offset) .* taper;
        
        % for model_dp
        coeff(offset + 1, 1:2) = model_dp' \ resp';
        coeff_n(offset + 1, 1:2) = model_dp' \ resp_n';
        resid(offset + 1, 1) = sum( (resp  - coeff(offset + 1, 1:2) * model_dp).^2);
    end
    
    resp = DPOAE(win) .* taper;
    resp_n = NOISE(win) .* taper;
    
    % for model_f1
    coeff(1, 3:4) = model_f1' \ resp';
    coeff_n(1, 3:4) = model_f1' \ resp_n';
    resid(1, 2) = sum( (resp  - coeff(1, 3:4) * model_f1).^2);
    
    % for model_f2
    coeff(1, 5:6) = model_f2' \ resp';
    coeff_n(1, 5:6) = model_f2' \ resp_n';
    resid(1, 3) = sum( (resp  - coeff(1, 5:6) * model_f2).^2);
    
    % for model_noise
    coeffs_noise(k,:) = model_noise' \ resp';
    
    [~, ind] = min(resid(:,1));
    coeffs(k, 1:2) = coeff(ind, 1:2);
    coeffs_n(k, 1:2) = coeff_n(ind, 1:2);
    coeffs(k, 3:6) = coeff(1,3:6);
    
    tau_dp(k) = (ind(1) - 1) * 1/stim.Fs; % delay in sec
end

%% Amplitude and Delay calculations
a_dp = coeffs(:, 1);
b_dp = coeffs(:, 2);
a_f1 = coeffs(:, 3);
b_f1 = coeffs(:, 4);
a_f2 = coeffs(:, 5);
b_f2 = coeffs(:, 6);
a_n = coeffs_n(:, 1);
b_n = coeffs_n(:, 2);

% for noise
noise2 = zeros(npoints,4);
for i = 1:2:8
    noise2(:,ceil(i/2)) = complex(coeffs_noise(:,i), coeffs_noise(:,i+1));
end

phi_dp = tau_dp.*freq_dp'; % cycles (from delay/offset)
phasor_dp = exp(-1j * phi_dp * 2 * pi);

oae_complex = complex(a_dp, b_dp);
noise_complex2 = complex(a_n, b_n);
noise_complex = mean(noise2,2);
res.multiplier = stim.VoltageToPascal.* stim.PascalToLinearSPL;

%% Plot Results Figure
figure;
plot(freq_f2/1000, db(abs(oae_complex).*res.multiplier), 'linew', 1.75);
hold on;
plot(freq_f2/1000, db(abs(noise_complex).*res.multiplier), '--', 'linew', 1.5);
%plot(freq_f2/1000, db(abs(noise_complex2).*res.multiplier));
%plot(freq_f2/1000, db(abs(complex(a_f2,b_f2)).*res.multiplier));
%plot(freq_f1/1000, db(abs(complex(a_f1, b_f1)).*res.multiplier));
title(sprintf('DPOAE Subj: %s, Ear: %s', string(subj), string(ear)))
set(gca, 'XScale', 'log', 'FontSize', 14)
xlim([.5, 16])
ylim([-40, 40])
xticks([.5, 1, 2, 4, 8, 16])
xlabel('F_2 Frequency (kHz)')
legend('DPOAE', 'NF')
drawnow; 

%% Apply EPL

if subj(1) == 'S'
    calib1 = 0; 
    calib2 = 0; 
    % find calib file
    calibdir = [maindir '/EARCAL/'];
    subjdir = [calibdir char(subj)];
    cd(subjdir)
    
    checkDIR=dir(sprintf('Calib_Ph1ER-10X_%s%s*.mat', subj, ear));
    if isempty(checkDIR)
        fprintf('No such Calibs for subject %s\n',subj);
    else
        if size(checkDIR, 1) > 1
            checkDIR = checkDIR(size(checkDIR,1));
        end
        load(checkDIR.name);
        res.calib.Ph1 = calib;
        calib1 = 1;
        clear calib;
    end
    
    checkDIR=dir(sprintf('Calib_Ph2ER-10X_%s%s*.mat', subj, ear));
    if isempty(checkDIR)
        fprintf('No such Calibs for subject %s\n',subj);
    else
        if size(checkDIR, 1) > 1
            checkDIR = checkDIR(size(checkDIR,1));
        end
        load(checkDIR.name);
        res.calib.Ph2 = calib;
        calib2 = 1;
    end
    
    cd(maindir)
    
    % if calibration is here
    if calib1 == 1 
        
        % Get EPL units
        [DP] = calc_EPL(freq_dp, oae_complex.*res.multiplier, res.calib.Ph1);
        res.complex.dp_epl = DP.P_epl;
        res.f_epl = DP.f;
        res.dbEPL_dp = db(abs(DP.P_epl));
        
        [NF] = calc_EPL(freq_dp, noise_complex.*res.multiplier, res.calib.Ph1);
        res.complex.nf_epl = NF.P_epl;
        res.f_epl = NF.f;
        res.dbEPL_nf = db(abs(NF.P_epl));
        
        %         [F1] = calc_FPL(res.f.f1, res.complex_f1, res.calib.Ph1);
        %         res.complex_f1_fpl = F1.P_fpl;
        %         res.f1_fpl = F1.f;
        %         if exist('res.calib.Ph2', 'var')
        %             [F2] = calc_FPL(res.f.f2, res.complex_f2, res.calib.Ph2);
        %         else
        %             [F2] = calc_FPL(res.f.f2, res.complex_f2, res.calib.Ph1);
        %         end
        %         res.complex_f2_fpl = F2.P_fpl;
        %         res.f2_fpl = F2.f;
        
        % plot figure again
        figure;
        plot(freq_f2/1000, db(res.dbEPL_dp));
        hold on;
        plot(freq_f2/1000, db(res.dbEPL_nf));
        plot(freq_f2/1000, db(abs(complex(a_f2,b_f2)).*stim.VoltageToPascal.*stim.PascalToLinearSPL));
        plot(freq_f1/1000, db(abs(complex(a_f1, b_f1)).*stim.VoltageToPascal.*stim.PascalToLinearSPL));
        %title(sprintf('Subj: %s, Ear: %s', string(subj), string(ear)))
        set(gca, 'XScale', 'log', 'FontSize', 14)
        xlim([.5, 16])
        ylim([-30, 80])
        xticks([.5, 1, 2, 4, 8, 16])
        ylabel('Amplitude db EPL')
        xlabel('F_2 Frequency (kHz)')
        legend('DP', 'NF', 'F2', 'F1')
    end
    
end


%% Save result function
res.windowdur = windowdur;
res.offsetwin = offsetwin;
res.npoints = npoints;
res.avgDPOAEresp = DPOAE;   % average mic response
res.avgNOISEresp = NOISE;
res.t_freq = t_freq;
res.f.f2 = freq_f2;         % frequency vectors
res.f.f1 = freq_f1;
res.f.dp = freq_dp;
res.a.dp = a_dp;            % coefficients
res.b.dp = b_dp;
res.a.f1 = a_f1;
res.b.f1 = b_f1;
res.a.f2 = a_f2;
res.b.f2 = b_f2;
res.a.n = a_n; % subtraction method
res.b.n = b_n;
res.tau.dp = tau_dp;
res.stim = stim;
res.subj = subj;
res.ear = ear;
res.complex.oae = oae_complex; 
res.complex.nf = noise_complex; 
res.complex.nf2 = noise_complex2; 

% Save
respdir = [maindir '/Results/DP/'];
fname = strcat(respdir, 'DP_', stim.subj, '_', stim.ear,'.mat');
save(fname,'res');


end

