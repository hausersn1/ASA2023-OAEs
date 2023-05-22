t = stim.t;
testfreq = [.75, 1, 1.5, 2, 3, 4, 6, 8, 12].* 1000;
windowdur = 0.5;

if stim.speed < 0
    f1 = stim.fmax;
    f2 = stim.fmin;
else
    f1 = stim.fmin;
    f2 = stim.fmax;
end

if stim.speed < 20
    t_freq = log2(testfreq/f1)/stim.speed + stim.buffdur;
else
    t_freq = (testfreq-f1)/stim.speed + stim.buffdur;
end
figure; 

for k = [6, 12, 23]
    OAEtrials = stim.resp(1:k, :);
    OAE = median(OAEtrials,1);
    coeffs_temp = zeros(length(testfreq), 2);
    coeffs_noise = zeros(length(testfreq), 8);
    for m = 1:length(testfreq)
        win = find( (t > (t_freq(m)-windowdur/2)) & ...
            (t < (t_freq(m)+windowdur/2)));
        taper = hanning(numel(win))';
        
        oae_win = OAE(win) .* taper;
        
        phi_dp_inst = (2.*stim.phi1_inst - stim.phi2_inst) * 2 * pi;
        phiProbe_inst = phi_dp_inst;
        model_dp = [cos(phiProbe_inst(win)) .* taper;
            -sin(phiProbe_inst(win)) .* taper];
        
        model_noise = ...
            [cos(0.9*phiProbe_inst(win)) .* taper;
            -sin(0.9*phiProbe_inst(win)) .* taper;
            cos(0.88*phiProbe_inst(win)) .* taper;
            -sin(0.88*phiProbe_inst(win)) .* taper;
            cos(0.86*phiProbe_inst(win)) .* taper;
            -sin(0.86*phiProbe_inst(win)) .* taper;
            cos(0.84*phiProbe_inst(win)) .* taper;
            -sin(0.84*phiProbe_inst(win)) .* taper];
        
        coeffs_temp(m,:) = model_dp' \ oae_win';
        coeffs_noise(m,:) = model_noise' \ oae_win';
    end
    
    % for noise
    noise2 = zeros(length(testfreq),4);
    for i = 1:2:8
        noise2(:,ceil(i/2)) = abs(complex(coeffs_noise(:,i), coeffs_noise(:,i+1)));
    end
    noise = mean(noise2, 2);
    
    oae = abs(complex(coeffs_temp(:,1), coeffs_temp(:,2)));
    
    SNR_temp = db(oae) - db(noise);
   
    figure_prop_name = {'PaperPositionMode', 'units', 'Position'};
    figure_prop_val = {'auto', 'inches', [1 1 6.5 3]}; % xcor, ycor, xwid, yheight
    set(gcf, figure_prop_name, figure_prop_val);
    subplot(1, 3, ceil(k/11))
    plot(testfreq./1000,db(oae.*10000), 'o', 'linew', 2);
    hold on;
    plot(testfreq./1000,db(noise.*10000), 'x', 'linew', 2);
    if k == 23
    lgd = legend('DPOAE', 'NOISE', 'location', 'Northwest');
    lgd.FontSize = 8; 
    end
    title(sprintf('Trials complete: %d', k))
    xlabel('Frequency (kHz)')
    ylabel('Median Amplitude (dB)')
    set(gca, 'XScale', 'log', 'FontSize', 10)
    xticks([.5, 1, 2, 4, 8, 16])
    xlim([0.5, 16])
    ylim([-20, 40])
    
end
print -djpeg -r600 SNRexample
