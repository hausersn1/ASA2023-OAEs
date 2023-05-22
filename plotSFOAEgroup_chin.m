
maindir = pwd; 
tts = {'Q402_R', 'Q402_L', 'Q406_R', 'Q406_L', 'Q407_R', 'Q407_L'};
carbo = {'Q405_R', 'Q405_L', 'Q403_R', 'Q403_L'};  %, 'Q381_R'};
ONH = {'Q298_L', 'Q298_R', 'Q351_R', 'Q351_L', 'Q368_R', 'Q368_L'};
YNH = {'Q421_L', 'Q421_R', 'Q412_L', 'Q412_R', 'Q414_L', 'Q414_R', 'Q415_L'}; 
groups = {tts, carbo, ONH, YNH};
labels = {'TTS', 'carbo', 'ONH', 'YNH'};

color1 = [140,81,10; 191,129,45; 223,194,125; % gold to teal
    246,232,195; 245,245,245; 199,234,229; 
    128,205,193; 53,151,143; 1,102,94]./255; 

color2 = [178,24,43; 214,96,77; 244,165,130; % red to blue
    253,219,199; 247,247,247; 209,229,240; 
    146,197,222; 67,147,195; 33,102,172]./255; 

c3 = [158,1,66; 213,62,79; 244,109,67; 
    253,174,97; 254,224,139; 255,255,191; 
    230,245,152; 171,221,164; 102,194,165; 
    50,136,189; 94,79,162]./255; 

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
%% Plot 
figure(20);
figure_prop_name = {'PaperPositionMode', 'units', 'Position'};
figure_prop_val = {'auto', 'inches', [1 1 7.5 6]}; % xcor, ycor, xwid, yheight
set(gcf, figure_prop_name, figure_prop_val);
hold on;
plot(f2, avg_oae(4,:), 'color', c3(9,:), 'linew', 3) % young
plot(f2, avg_oae(1,:), 'color', c3(2,:), 'linew', 3) % tts
plot(f2, avg_oae(2,:), 'color', c3(11,:), 'linew', 3) % carbo
plot(f2, avg_oae(3,:), 'color', c3(5,:), 'linew', 3) % older

plot(f2, avg_nf(4,:), '--', 'color', c3(9,:), 'linew', 2)
plot(f2, avg_nf(1,:), '--', 'color', c3(2,:), 'linew', 2)
plot(f2, avg_nf(2,:), '--', 'color', c3(11,:),  'linew', 2)
plot(f2, avg_nf(3,:), '--', 'color', c3(5,:), 'linew', 2)

xticks([0.5, 1, 2, 4, 8, 16])
xlim([0.5, 16])
ylim([-40, 50])
yticks([-30, -15 0, 15, 30, 45])
xlabel('Frequency (kHz)')
ylabel('Amplitude (dB SPL)')
title('Group Average SFOAE')
lgd = legend('Young, unexposed','TTS', 'Carboplatin', 'Older', 'location', 'Northwest'); 
lgd.FontSize = 16; 
set(gca, 'XScale', 'log', 'FontSize', 18, 'FontName', 'Franklin Gothic')
legend boxoff
print -djpeg -r600 groupSFOAE


