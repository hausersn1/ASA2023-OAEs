
maindir = pwd; 
tts = {'Q406_R', 'Q406_L', 'Q407_R', 'Q407_L'};
carbo = {'Q403_R', 'Q403_L', 'Q381_R'};
ONH = {'Q298_L', 'Q298_R', 'Q351_R', 'Q351_L', 'Q368_R', 'Q368_L'};
YNH = {'Q421_L', 'Q421_R', 'Q412_L', 'Q412_R'}; 
groups = {tts, carbo, ONH, YNH};
labels = {'TTS', 'carbo', 'ONH', 'YNH'};

avg_oae = zeros(length(groups), 512); 
avg_nf = zeros(length(groups), 512); 
for x = 1:length(groups)
    figure(x+10);
    figure_prop_name = {'PaperPositionMode', 'units', 'Position'};
    figure_prop_val = {'auto', 'inches', [1 1 5 4]}; % xcor, ycor, xwid, yheight
    set(gcf, figure_prop_name, figure_prop_val);
    title(sprintf('Group: %s', string(labels(x))))
    xlabel('Frequency (kHz)')
    ylabel('SF Amplitude (dB SPL)')
    set(gca, 'FontSize', 14, 'XScale', 'log')
    group = groups{x};
    oae_group = zeros(length(group), 512); 
    nf_group = zeros(length(group), 512);
    hold on; 
    for y = 1:length(group)
        hold on;
        load(sprintf('%s/Results/SF/SF_%s.mat', maindir, string(group(y))))
        
        oae = db(abs(res.complex.oae).*res.multiplier);
        nf = db(abs(res.complex.nf).*res.multiplier);
        f2 = res.f/1000;
        
        plot(f2, oae, 'linew', 1.5)
        plot(f2, nf, '--', 'linew', 1.25)
        
        oae_group(y,:) = oae';
        nf_group(y,:) = nf';
        
        set(gca, 'FontSize', 14, 'XScale', 'log')
        xlim([0.5, 16])
        ylim([-40, 60])
        
    end
    avg_oae(x,:) = mean(oae_group,1);
    avg_nf(x,:) = mean(nf_group,1);
    
end

figure(20);
figure_prop_name = {'PaperPositionMode', 'units', 'Position'};
figure_prop_val = {'auto', 'inches', [1 1 6 5]}; % xcor, ycor, xwid, yheight
set(gcf, figure_prop_name, figure_prop_val);
hold on;
plot(f2, avg_oae(1,:), 'color', [0.6350 0.0780 0.1840], 'linew', 2)
plot(f2, avg_oae(2,:), 'color', [0 0.447 0.741], 'linew', 2)
plot(f2, avg_oae(3,:), 'color', [.466 0.674 0.188], 'linew', 2)
plot(f2, avg_oae(4,:), 'k', 'linew', 2)

plot(f2, avg_nf(1,:), '--','color', [0.6350 0.0780 0.1840], 'linew', 1.5)
plot(f2, avg_nf(2,:), '--','color', [0 0.447 0.741],  'linew', 1.5)
plot(f2, avg_nf(3,:), '--', 'color', [.466 0.674 0.188], 'linew', 1.5)
plot(f2, avg_nf(4,:),'k--', 'linew', 1.5)
xticks([0.5, 1, 2, 4, 8, 16])
xlim([0.5, 16])
ylim([-40, 50])
yticks([-30, -15 0, 15, 30, 45])
xlabel('Frequency (kHz)')
ylabel('Amplitude (dB SPL)')
title('Group Average SFOAE')
set(gca, 'XScale', 'log', 'FontSize', 14)
lgd = legend('TTS', 'Carbo', 'ONH', 'YNH', 'location', 'Northwest'); 
lgd.FontSize = 12; 
print -djpeg -r600 groupSFOAE


