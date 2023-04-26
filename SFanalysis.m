function [res] = SFanalysis(subj, ear, windowdur, offsetwin, npoints)

maindir = pwd;

if subj(1) == 'Q'
    datadir = [maindir '/swept_wCalib/Results/'];
else
    datadir = [maindir '/sweptSFOAE-main/Results/'];
end

subjdir = [datadir char(subj)];
cd(subjdir)

% Check if subject exists
checkDIR=dir(sprintf('SFOAE_log_%s%s*.mat', subj, ear));
if isempty(checkDIR)
    error('No such OAEs for subject %s',subj);
end

load(checkDIR.name);
cd(maindir)

if nargin == 2
    windowdur = 0.060; % 40ms in paper
    offsetwin = 0.020; % 20ms in paper
    npoints = 1024;
    plot = 0;
end

hpfilter = 1;

% setting stuff from stim
phiProbe_inst = stim.phiProbe_inst;
t = stim.t;

if stim.speed < 0
    f1 = stim.fmax;
    f2 = stim.fmin;
else
    f1 = stim.fmin;
    f2 = stim.fmax;
end

% set SFOAE
SFOAEtrials = stim.ProbeBuffs + stim.SuppBuffs - stim.BothBuffs;

% HP filter
if hpfilter == 1
    filtcutoff = 300;
    b = fir1(999, filtcutoff*2/stim.Fs, 'high');
    SFOAEtrials = filtfilt(b, 1, SFOAEtrials')';
end


% Alternate Artifact rejection and set mean SFOAE and NOISE
% do twice, once for pos and once for neg suppressor
energy = squeeze(sum(SFOAEtrials.^2, 2)); % same cut off for both trial types
good = energy < median(energy) + 1.2*mad(energy);

count = 0;
trial_2x = floor(size(SFOAEtrials,1)/2)*2;
for y = 1:2:trial_2x
    if good(y) == 1 && good(y+1) == 1
        count = count +1;
        pos_SFclean(count,:) = SFOAEtrials(y,:);
        neg_SFclean(count,:) = SFOAEtrials(y+1,:);
    end
end

count_2x = floor(count/2)*2;
SFOAEclean = [pos_SFclean(1:count_2x,:); neg_SFclean(1:count_2x,:) ];
SFOAE = mean(SFOAEclean, 1);

pos_noise = zeros(count_2x, size(pos_SFclean, 2));
neg_noise = zeros(count_2x, size(pos_SFclean, 2));
count = 0;
for x = 1:2:count_2x
    count = count + 1;
    pos_noise(count,:) = (pos_SFclean(x, :) - pos_SFclean(x+1, :)) / 2;
    neg_noise(count,:) = (neg_SFclean(x, :) - neg_SFclean(x+1, :)) / 2;
end
noise = [pos_noise; neg_noise ];
NOISE = mean(noise,1);

% set freq we're testing and the timepoints when they happen.

if stim.speed < 20
    testfreq = 2 .^ linspace(log2(f1), log2(f2), npoints);
    t_freq = log2(testfreq/f1)/stim.speed + stim.buffdur;
else
    testfreq = linspace(f1, f2, npoints);
    t_freq = (testfreq-f1)/stim.speed + stim.buffdur;
end

% Set empty matricies for next steps
maxoffset = ceil(stim.Fs * offsetwin);
coeffs = zeros(npoints, 2);
coeffs_n = zeros(npoints, 2);
tau = zeros(npoints, 1);

% Generate model of chirp and test against response
for k = 1:npoints
    fprintf(1, 'Running window %d / %d\n', k, (npoints));
    
    win = find( (t > (t_freq(k)-windowdur/2)) & ...
        (t < (t_freq(k)+windowdur/2)));
    taper = hanning(numel(win))';
    resp = SFOAE(win) .* taper;
    resp_n = NOISE(win) .* taper;
    
    % zero out variables for offset calc
    coeff = zeros(maxoffset, 2);
    coeff_n = zeros(maxoffset, 2);
    resid = zeros(maxoffset, 1);
    resid_n = zeros(maxoffset, 1);
    
    for offset = 0:maxoffset
        model = [cos(2*pi*phiProbe_inst(win-offset)) .* taper;
            -sin(2*pi*phiProbe_inst(win-offset)) .* taper];
        
        coeff(offset + 1, :) = model' \ resp';
        coeff_n(offset + 1, :) = model' \ resp_n';
        
        resid(offset +1) = sum( (resp - coeff(offset+1, :) * model).^2);
        
    end
    
    [~, ind] = min(resid);
    
    coeffs(k, :) = coeff(ind, :);
    coeffs_n(k, :) = coeff_n(ind, :);
    
    tau(k) = (ind - 1) * (1/stim.Fs); % delay in sec
    
end

a = coeffs(:, 1);
b = coeffs(:, 2);
a_n = coeffs_n(:, 1);
b_n = coeffs_n(:, 2);
t_freq_tau = t_freq' - tau;
f_tau = ((t_freq_tau - stim.buffdur).*stim.speed)+f1;
phi = tau .* f_tau; % cycles
phasor = exp(-1j * phi * 2 *pi);
oae = abs(complex(a, b) .* phasor) .* stim.VoltageToPascal .* stim.PascalToLinearSPL;
nf = abs(complex(a_n, b_n)) .* stim.VoltageToPascal .* stim.PascalToLinearSPL;

theta = unwrap(angle(complex(a,b).*phasor))/(2*pi); % cycles
tau_pg = -diff(theta)./diff((f_tau./1000));
x = (f_tau(2:end) + f_tau(1:end-1))/2;

%% Save result function
res.windowdur = windowdur;
res.offsetwin = offsetwin;
res.npoints = npoints;
res.avgSFOAEresp = SFOAE;   % average mic response
res.avgNOISEresp = NOISE;
res.t_freq = t_freq;
res.f = testfreq;           % frequency vectors
res.a = a;                  % coefficients
res.b = b;
res.a_n = a_n;
res.b_n = b_n;
res.stim = stim;
res.tau = tau;
res.complex_sf = complex(a,b) .* phasor .* stim.VoltageToPascal .* stim.PascalToLinearSPL;
res.complex_nf = complex(a_n, b_n) .* stim.VoltageToPascal .* stim.PascalToLinearSPL;

%% Calib

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
    clear calib;
end

cd(maindir)

%% Get EPL units
[SF] = calc_EPL(res.f, res.complex_sf, res.calib.Ph1);
res.complex_sf_epl = SF.P_epl;
res.f_epl = SF.f;
[NF] = calc_EPL(res.f, res.complex_nf, res.calib.Ph1);
res.complex_nf_epl = NF.P_epl;
res.subj = subj;
res.ear = ear;

%% Save
respdir = [maindir '/Results/'];
fname = strcat(respdir, 'SF_', res.subj, '_', res.ear,'.mat');
save(fname,'res');

%% plot
SFplot(res)

end

