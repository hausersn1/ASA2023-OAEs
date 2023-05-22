%% Binning and further analysis in 1/3 octave bins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takes analyzed data and
% 1. Further bins into 1/2 octaves
% 2.
% 3.
%
% Samantha Hauser, May 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subjects to set in format 'SubjID_ear')
tts = {'Q402_R', 'Q402_L', 'Q406_R', 'Q406_L', 'Q407_R', 'Q407_L'};
carbo = {'Q405_R', 'Q405_L', 'Q403_R', 'Q403_L'};  %, 'Q381_R'};
ONH = {'Q298_L', 'Q298_R', 'Q351_R', 'Q351_L', 'Q368_R', 'Q368_L'};
YNH = {'Q421_L', 'Q421_R', 'Q412_L', 'Q412_R', 'Q414_L', 'Q414_R', 'Q415_L'}; 
%human = {'S356_R', 'S343_R', 'S354_R', 'S357_R', 'SH1_R', 'SBC_R', ...
%    'S356_L', 'S343_L', 'S354_L', 'S357_L', 'SH1_L', 'SBC_L'}; 
% Set groups to analyze
groups = {tts, carbo, ONH, YNH};
labels = {'TTS', 'Carbo', 'ONH', 'YNH'};

% Colors
% colors light to dark
blue=[247,251,255;222,235,247;198,219,239;
    158,202,225; 107,174,214; 66,146,198 ; 
    33,113,181; 8,81,156; 8,48,107]./255;
orange = [255,245,235; 254,230,206; 253,208,162;
    253,174,107; 253,141,60; 241,105,19; 
    217,72,1; 166,54,3; 127,39,4]./255; 
purple = [252,251,253; 239,237,245; 218,218,235; 
    188,189,220; 158,154,200; 128,125,186; 
    106,81,163; 84,39,143; 63,0,125]./255; 
black = [255,255,255; 240,240,240; 217,217,217; 
    189,189,189; 150,150,150; 115,115,115; 
    82,82,82; 37,37,37; 10,10,10]./255;

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

m={'+', 'o', '*', '.', 'x', 'square', 'diamond', 'v', '^', '>', '<', 'pentagram'}; 


% Set frequency bands to analyze
fmin = 0.5;
fmax = 16;
edges = 2 .^ linspace(log2(fmin), log2(fmax), 21);
bandEdges = edges(2:2:end-1);
centerFreqs = edges(3:2:end-2);

% load data
maindir = pwd;
avg_oae = zeros(length(groups), length(centerFreqs));
avg_nf = zeros(length(groups), length(centerFreqs));
for x = 1:length(groups)
    
    group = groups{x};
    oae_group = zeros(length(group), 9);
    nf_group = zeros(length(group), 9);
    
    for y = 1:length(group)
        
        % For DPOAE
        load(sprintf('%s/Results/DP/DP_%s.mat', maindir, string(group(y))))
        
        dpoae_full = db(abs(res.complex.oae).*res.multiplier);
        dpnf_full = db(abs(res.complex.nf).*res.multiplier);
        f2 = res.f.f2/1000;
        
        dpoae = zeros(length(centerFreqs),1);
        dpnf = zeros(length(centerFreqs),1);
        dpoae_w = zeros(length(centerFreqs),1);
        dpnf_w = zeros(length(centerFreqs),1);
        % resample / average to 9 center frequencies
        for z = 1:length(centerFreqs)
            band = find( f2 >= bandEdges(z) & f2 < bandEdges(z+1));
            
            % Do some weighting by SNR
            
            % TO DO: NF from which SNR was calculated included median of 7 points
            % nearest the target frequency. 
            SNR = dpoae_full(band) - dpnf_full(band);
            weight = (10.^(SNR./10)).^2;
            
            dpoae(z, 1) = mean(dpoae_full(band));
            dpnf(z,1) = mean(dpnf_full(band));
            
            dpoae_w(z,1) = sum(weight.*dpoae_full(band))/sum(weight); 
            dpnf_w(z,1) = sum(weight.*dpnf_full(band))/sum(weight); 
            
        end
        
        % Same for SFOAE
        load(sprintf('%s/Results/SF/SF_%s.mat', maindir, string(group(y))))
        
        sfoae_full = db(abs(res.complex.oae).*res.multiplier);
        sfnf_full = db(abs(res.complex.nf).*res.multiplier);
        f = res.f/1000;
        
        sfoae = zeros(length(centerFreqs),1);
        sfnf = zeros(length(centerFreqs),1);
        sfoae_w = zeros(length(centerFreqs),1);
        sfnf_w = zeros(length(centerFreqs),1);
        % resample / average to 9 center frequencies
        for z = 1:length(centerFreqs)
            band = find( f >= bandEdges(z) & f < bandEdges(z+1));
            
            % Do some weighting by SNR
            SNR = sfoae_full(band) - sfnf_full(band);
            weight = (10.^(SNR./10)).^2;
            
            sfoae(z, 1) = mean(sfoae_full(band));
            sfnf(z,1) = mean(sfnf_full(band));
            
            sfoae_w(z,1) = sum(weight.*sfoae_full(band))/sum(weight); 
            sfnf_w(z,1) = sum(weight.*sfnf_full(band))/sum(weight); 
            
        end
        
%         figure;
%         hold on;
%         plot(centerFreqs, dpoae, 'ob', 'linew', 1.5)
%         %plot(centerFreqs, dpnf, '--', 'linew', 1.25)
%         plot(centerFreqs, dpoae_w, '*b', 'linew', 1.5)
%         %plot(centerFreqs, dpnf_w, '--', 'linew', 1.25)
%         plot(centerFreqs, sfoae, 'or', 'linew', 1.5)
%         plot(centerFreqs, sfoae_w, '*r', 'linew', 1.5)
%                 set(gca, 'FontSize', 14, 'XScale', 'log')
%         xlim([0.5, 16])
%         ylim([-40, 60])
        
        
        dpoae_group(y,:) = dpoae';
        dpnf_group(y,:) = dpnf';
        dpoae_w_group(y,:) = dpoae_w';
        dpnf_w_group(y,:) = dpnf_w';
        sfoae_group(y,:) = sfoae';
        sfnf_group(y,:) = sfnf';
        sfoae_w_group(y,:) = sfoae_w';
        sfnf_w_group(y,:) = sfnf_w';

        
    end
    avg_dpoae(x,:) = mean(dpoae_w_group,1);
    avg_dpnf(x,:) = mean(dpnf_w_group,1);
 
    avg_sfoae(x,:) = mean(sfoae_w_group,1);
    avg_sfnf(x,:) = mean(sfnf_w_group,1);
    
    if x == 1
        tts_dpoae = dpoae_w_group; 
        tts_sfoae = sfoae_w_group; 
    elseif x == 2
        carbo_dpoae = dpoae_w_group; 
        carbo_sfoae = sfoae_w_group; 
    elseif x == 3
        onh_dpoae = dpoae_w_group; 
        onh_sfoae = sfoae_w_group;   
    elseif x == 4
        ynh_dpoae = dpoae_w_group; 
        ynh_sfoae = sfoae_w_group; 
    else
        human_dpoae = dpoae_w_group; 
        human_sfoae = sfoae_w_group; 
    end
end


%% LME 

%make the table
points = length(centerFreqs); 

for i = 1:size(carbo_dpoae,1)
    amp_carbo((1:points)+((i-1).*points), 1) = carbo_dpoae(i, :); 
    group_carbo((1:points)+((i-1).*points), 1) = {'carbo'};
    freq_carbo((1:points)+((i-1).*points), 1) = centerFreqs;
    type_carbo((1:points)+((i-1).*points), 1) = {'dp'}; 
    chin_carbo((1:points)+((i-1).*points), 1) = i; 
    amp_carbo2((1:points)+((i-1).*points), 1) = carbo_sfoae(i, :); 
    group_carbo2((1:points)+((i-1).*points), 1) = {'carbo'};
    freq_carbo2((1:points)+((i-1).*points), 1) = centerFreqs;
    type_carbo2((1:points)+((i-1).*points), 1) = {'sf'}; 
    chin_carbo2((1:points)+((i-1).*points), 1) = i; 
end 
j = i; 

for i = 1:size(tts_dpoae,1)
    amp_tts((1:points)+((i-1).*points), 1) = tts_dpoae(i, :); 
    group_tts((1:points)+((i-1).*points), 1) = {'tts'};
    freq_tts((1:points)+((i-1).*points), 1) = centerFreqs;
    type_tts((1:points)+((i-1).*points), 1) = {'dp'}; 
    amp_tts2((1:points)+((i-1).*points), 1) = tts_sfoae(i, :); 
    group_tts2((1:points)+((i-1).*points), 1) = {'tts'};
    freq_tts2((1:points)+((i-1).*points), 1) = centerFreqs;
    type_tts2((1:points)+((i-1).*points), 1) = {'sf'}; 
    chin_tts((1:points)+((i-1).*points), 1) = j+i; 
    chin_tts2((1:points)+((i-1).*points), 1) = j+i; 
end 
k = j+i; 
for i = 1:size(onh_dpoae,1)
    amp_onh((1:points)+((i-1).*points), 1) = onh_dpoae(i, :); 
    group_onh((1:points)+((i-1).*points), 1) = {'onh'};
    freq_onh((1:points)+((i-1).*points), 1) = centerFreqs;
    type_onh((1:points)+((i-1).*points), 1) = {'dp'}; 
    amp_onh2((1:points)+((i-1).*points), 1) = onh_sfoae(i, :); 
    group_onh2((1:points)+((i-1).*points), 1) = {'onh'};
    freq_onh2((1:points)+((i-1).*points), 1) = centerFreqs;
    type_onh2((1:points)+((i-1).*points), 1) = {'sf'}; 
    chin_onh((1:points)+((i-1).*points), 1) = i+k; 
    chin_onh2((1:points)+((i-1).*points), 1) = i+k; 
end 
l = i+k; 

for i = 1:size(ynh_dpoae,1)
    amp_ynh((1:points)+((i-1).*points), 1) = ynh_dpoae(i, :); 
    group_ynh((1:points)+((i-1).*points), 1) = {'ynh'};
    freq_ynh((1:points)+((i-1).*points), 1) = centerFreqs;
    type_ynh((1:points)+((i-1).*points), 1) = {'dp'}; 
    amp_ynh2((1:points)+((i-1).*points), 1) = ynh_sfoae(i, :); 
    group_ynh2((1:points)+((i-1).*points), 1) = {'ynh'};
    freq_ynh2((1:points)+((i-1).*points), 1) = centerFreqs;
    type_ynh2((1:points)+((i-1).*points), 1) = {'sf'}; 
    chin_ynh((1:points)+((i-1).*points), 1) = l+i; 
    chin_ynh2((1:points)+((i-1).*points), 1) = l+i; 
end 

amp = [amp_carbo; amp_carbo2; amp_tts; amp_tts2;amp_onh; amp_onh2;amp_ynh; amp_ynh2]; 
group = [group_carbo; group_carbo2; group_tts; group_tts2;group_onh; group_onh2;group_ynh; group_ynh2]; 
freq = [freq_carbo; freq_carbo2; freq_tts; freq_tts2;freq_onh; freq_onh2;freq_ynh; freq_ynh2]; 
type = [type_carbo; type_carbo2; type_tts; type_tts2;type_onh; type_onh2;type_ynh; type_ynh2]; 
chin = [chin_carbo; chin_carbo2; chin_tts; chin_tts2;chin_onh; chin_onh2;chin_ynh; chin_ynh2];

table = table(chin, type, freq, group, amp,'VariableNames',{'chin', 'OAEtype','Frequency','Group','dBSPL'});

table.OAEtype = categorical(table.OAEtype); 
table.Group = categorical(table.Group); 
table.chin = categorical(table.chin); 

writetable(table, 'data.csv')

%% Plotting Results -- x and o for different half octave bands
figure;
figure_prop_name = {'PaperPositionMode', 'units', 'Position'};
figure_prop_val = {'auto', 'inches', [1 1 5 4]}; % xcor, ycor, xwid, yheight
set(gcf, figure_prop_name, figure_prop_val);
hold on;
plot(centerFreqs, avg_dpoae(1,:), 'x', 'color', [0.6350 0.0780 0.1840], 'linew', 1.5)
plot(centerFreqs, avg_dpoae(2,:), 'x','color', [0 0.447 0.741], 'linew', 1.5)
plot(centerFreqs, avg_dpoae(3,:), 'x','color', [.466 0.674 0.188], 'linew', 1.5)
plot(centerFreqs, avg_dpoae(4,:), 'xk', 'linew', 2)

plot(centerFreqs, avg_sfoae(1,:), 'o','color', [0.6350 0.0780 0.1840], 'linew', 1.5)
plot(centerFreqs, avg_sfoae(2,:), 'o','color', [0 0.447 0.741],  'linew', 1.5)
plot(centerFreqs, avg_sfoae(3,:), 'o', 'color', [.466 0.674 0.188], 'linew', 1.5)
plot(centerFreqs, avg_sfoae(4,:),'ok', 'linew', 1.5)
xticks([0.5, 1, 2, 4, 8, 16])
xlim([0.5, 16])
ylim([-10, 50])
yticks([0, 10, 20, 30, 40, 50])
title('Group Averages')
ylabel('Amplitude (dB SPL)')
xlabel('Frequency (kHz)')
legend('TTS', 'Carbo', 'ONH', 'YNH', 'location', 'Northwest')
set(gca, 'XScale', 'log', 'FontSize', 14)

print -djpeg -r600 groupDPvsSFbyFreq
%% all data plot ---- don't need, by frequency 

figure; 
figure_prop_name = {'PaperPositionMode', 'units', 'Position'};
figure_prop_val = {'auto', 'inches', [1 1 5 4]}; % xcor, ycor, xwid, yheight
set(gcf, figure_prop_name, figure_prop_val);
hold on;
sz = 50; 
plot([-10, 50], [-10, 50], '--', 'color', black(6, :), 'linew', 2)
for i=1:size(tts_sfoae,2)
    scatter(tts_sfoae(:,i), tts_dpoae(:,i), sz, 'Marker', m{2},...
        'MarkerFaceColor', orange(2*ceil(i/3),:), 'MarkerEdgeColor', orange(9,:))
    scatter(avg_sfoae(1,i), avg_dpoae(1,i), sz, 'Marker', m{7},...
        'MarkerFaceColor', orange(6+ceil(i/3),:), 'MarkerEdgeColor', orange(9,:))
end
for i=1:size(carbo_sfoae,2)
    scatter(carbo_sfoae(:,i), carbo_dpoae(:,i), sz, 'Marker', m{2},...
        'MarkerFaceColor', blue(2*ceil(i/3),:), 'MarkerEdgeColor', blue(9,:))
    scatter(avg_sfoae(2,i), avg_dpoae(2,i), sz, 'Marker', m{7},...
        'MarkerFaceColor', blue(6+ceil(i/3),:), 'MarkerEdgeColor', blue(9,:))
end
for i=1:size(onh_sfoae,2)
    scatter(onh_sfoae(:,i), onh_dpoae(:,i), sz, 'Marker', m{2},...
        'MarkerFaceColor', purple(2*ceil(i/3),:), 'MarkerEdgeColor', purple(9,:))
    scatter(avg_sfoae(3,i), avg_dpoae(3,i), sz, 'Marker', m{7},...
        'MarkerFaceColor', purple(6+ceil(i/3),:), 'MarkerEdgeColor', purple(9,:))
end

title('SF vs DP amplitude')
xlabel('SFOAE dB SPL')
ylabel('DPOAE dB SPL')
%legend('TTS', 'Carbo', 'ONH', 'location', 'best')
set(gca, 'FontSize', 14)
xlim([-10, 50])
ylim([-10, 50])

print -djpeg -r600 groupSFvsDP

%% %% 
figure; 
figure_prop_name = {'PaperPositionMode', 'units', 'Position'};
figure_prop_val = {'auto', 'inches', [1 1 6 5]}; % xcor, ycor, xwid, yheight
set(gcf, figure_prop_name, figure_prop_val);
hold on;
sz = 50; 
plot([-10, 50], [-10, 50], '--', 'color', black(6, :), 'linew', 2)
for i=1:size(tts_sfoae,2)
    scatter(tts_sfoae(:,i), tts_dpoae(:,i), sz, 'Marker', m{2},...
        'MarkerFaceColor', color2(2,:), 'MarkerEdgeColor', 'k')
    p1 = scatter(avg_sfoae(1,i), avg_dpoae(1,i), sz, 'Marker', m{7},...
        'MarkerFaceColor', color2(1,:), 'MarkerEdgeColor', 'k');
end
for i=1:size(carbo_sfoae,2)
    scatter(carbo_sfoae(:,i), carbo_dpoae(:,i), sz, 'Marker', m{2},...
        'MarkerFaceColor', color2(8,:), 'MarkerEdgeColor', 'k')
    p2 = scatter(avg_sfoae(2,i), avg_dpoae(2,i), sz, 'Marker', m{7},...
        'MarkerFaceColor', color2(9,:), 'MarkerEdgeColor', 'k');
end
for i=1:size(onh_sfoae,2)
    scatter(onh_sfoae(:,i), onh_dpoae(:,i), sz, 'Marker', m{2},...
        'MarkerFaceColor', color1(3,:), 'MarkerEdgeColor', 'k')
    p3 = scatter(avg_sfoae(3,i), avg_dpoae(3,i), sz, 'Marker', m{7},...
        'MarkerFaceColor', color1(2,:), 'MarkerEdgeColor', 'k');
end
for i=1:size(ynh_sfoae,2)
    scatter(ynh_sfoae(:,i), ynh_dpoae(:,i), sz, 'Marker', m{2},...
        'MarkerFaceColor', color1(7,:), 'MarkerEdgeColor', 'k')
    p4 = scatter(avg_sfoae(4,i), avg_dpoae(4,i), sz, 'Marker', m{7},...
        'MarkerFaceColor', color1(8,:), 'MarkerEdgeColor', 'k'); 
end
title('SF vs DP amplitude')
xlabel('SFOAE (dB SPL)')
ylabel('DPOAE (dB SPL)')
legend([p1 p2 p3 p4], {'TTS', 'Carbo', 'ONH', 'YNH'}, 'location', 'Southeast')
set(gca, 'FontSize', 16, 'FontName', 'Franklin Gothic')
xlim([-10, 50])
ylim([-10, 50])

print -djpeg -r600 allDvR

%% Young vs Old 
figure; 
figure_prop_name = {'PaperPositionMode', 'units', 'Position'};
figure_prop_val = {'auto', 'inches', [1 1 5 4.5]}; % xcor, ycor, xwid, yheight
set(gcf, figure_prop_name, figure_prop_val);
hold on;
sz = 70; 
plot([-10, 50], [-10, 50], '--', 'color', black(6, :), 'linew', 2)
for i=1:size(ynh_sfoae,2)
    if i == 1
    p2 = scatter(ynh_sfoae(:,i), ynh_dpoae(:,i), sz, 'Marker', m{2},...
        'MarkerFaceColor', c3(9,:), 'MarkerEdgeColor', 'k'); 
    else 
        scatter(ynh_sfoae(:,i), ynh_dpoae(:,i), sz, 'Marker', m{2},...
        'MarkerFaceColor', c3(9,:), 'MarkerEdgeColor', 'k'); 
    end
    %p4 = scatter(avg_sfoae(4,i), avg_dpoae(4,i), sz, 'Marker', m{7},...
    %    'MarkerFaceColor', color1(8,:), 'MarkerEdgeColor', 'k'); 
end
for i=1:size(onh_sfoae,2)
    if i ==1
    p1 = scatter(onh_sfoae(:,i), onh_dpoae(:,i), sz, 'Marker', m{2},...
        'MarkerFaceColor', c3(5,:), 'MarkerEdgeColor', 'k'); 
    else
        scatter(onh_sfoae(:,i), onh_dpoae(:,i), sz, 'Marker', m{2},...
        'MarkerFaceColor', c3(5,:), 'MarkerEdgeColor', 'k')
    end
    %p3 = scatter(avg_sfoae(3,i), avg_dpoae(3,i), sz, 'Marker', m{7},...
    %    'MarkerFaceColor', color1(2,:), 'MarkerEdgeColor', 'k');
end

title('Age')
xlabel('SFOAE (dB SPL)')
ylabel('DPOAE (dB SPL)')
legend([p1 p2], {'Older', 'Younger'}, 'location', 'Southeast')
xlim([-10, 50])
ylim([-10, 50])
xticks([-10:10:50])
yticks([0:10:50])
set(gca, 'FontSize', 18, 'FontName', 'Franklin Gothic')

print -djpeg -r600 youngOldDvR

%% Young vs TTS 
figure; 
figure_prop_name = {'PaperPositionMode', 'units', 'Position'};
figure_prop_val = {'auto', 'inches', [1 1 5 4.5]}; % xcor, ycor, xwid, yheight
set(gcf, figure_prop_name, figure_prop_val);
hold on;
sz = 70; 
plot([-10, 50], [-10, 50], '--', 'color', black(6, :), 'linew', 2)
for i=1:size(ynh_sfoae,2)
    if i == 1
    p2 = scatter(ynh_sfoae(:,i), ynh_dpoae(:,i), sz, 'Marker', m{2},...
        'MarkerFaceColor', c3(9,:), 'MarkerEdgeColor', 'k'); 
    else 
        scatter(ynh_sfoae(:,i), ynh_dpoae(:,i), sz, 'Marker', m{2},...
        'MarkerFaceColor', c3(9,:), 'MarkerEdgeColor', 'k'); 
    end
    %p4 = scatter(avg_sfoae(4,i), avg_dpoae(4,i), sz, 'Marker', m{7},...
    %    'MarkerFaceColor', color1(8,:), 'MarkerEdgeColor', 'k'); 
end
for i=1:size(tts_sfoae,2)
    if i ==1
    p1 = scatter(tts_sfoae(:,i), tts_dpoae(:,i), sz, 'Marker', m{2},...
        'MarkerFaceColor', c3(2,:), 'MarkerEdgeColor', 'k'); 
    else
        scatter(tts_sfoae(:,i), tts_dpoae(:,i), sz, 'Marker', m{2},...
        'MarkerFaceColor', c3(2,:), 'MarkerEdgeColor', 'k')
    end
    %p3 = scatter(avg_sfoae(3,i), avg_dpoae(3,i), sz, 'Marker', m{7},...
    %    'MarkerFaceColor', color1(2,:), 'MarkerEdgeColor', 'k');
end

title('Synaptopathy')
xlabel('SFOAE (dB SPL)')
ylabel('DPOAE (dB SPL)')
legend([p1 p2], {'TTS', 'Unexposed'}, 'location', 'Southeast')
xlim([-10, 50])
ylim([-10, 50])
xticks([-10:10:50])
yticks([0:10:50])
set(gca, 'FontSize', 18, 'FontName', 'Franklin Gothic')

print -djpeg -r600 youngTTSDvR

%% Young vs Carbo 
figure; 
figure_prop_name = {'PaperPositionMode', 'units', 'Position'};
figure_prop_val = {'auto', 'inches', [1 1 5 4.5]}; % xcor, ycor, xwid, yheight
set(gcf, figure_prop_name, figure_prop_val);
hold on;
sz = 70; 
plot([-10, 50], [-10, 50], '--', 'color', black(6, :), 'linew', 2)
for i=1:size(ynh_sfoae,2)
    if i == 1
    p2 = scatter(ynh_sfoae(:,i), ynh_dpoae(:,i), sz, 'Marker', m{2},...
        'MarkerFaceColor', c3(9,:), 'MarkerEdgeColor', 'k'); 
    else 
        scatter(ynh_sfoae(:,i), ynh_dpoae(:,i), sz, 'Marker', m{2},...
        'MarkerFaceColor', c3(9,:), 'MarkerEdgeColor', 'k'); 
    end
    %p4 = scatter(avg_sfoae(4,i), avg_dpoae(4,i), sz, 'Marker', m{7},...
    %    'MarkerFaceColor', color1(8,:), 'MarkerEdgeColor', 'k'); 
end
for i=1:size(carbo_sfoae,2)
    if i ==1
    p1 = scatter(carbo_sfoae(:,i), carbo_dpoae(:,i), sz, 'Marker', m{2},...
        'MarkerFaceColor', c3(11,:), 'MarkerEdgeColor', 'k'); 
    else
        scatter(carbo_sfoae(:,i), carbo_dpoae(:,i), sz, 'Marker', m{2},...
        'MarkerFaceColor', c3(11,:), 'MarkerEdgeColor', 'k')
    end
    %p3 = scatter(avg_sfoae(3,i), avg_dpoae(3,i), sz, 'Marker', m{7},...
    %    'MarkerFaceColor', color1(2,:), 'MarkerEdgeColor', 'k');
end

title('Carboplatin')
xlabel('SFOAE (dB SPL)')
ylabel('DPOAE (dB SPL)')
legend([p1 p2], {'Carboplatin', 'Unexposed'}, 'location', 'Southeast')
xlim([-10, 50])
ylim([-10, 50])
xticks([-10:10:50])
yticks([0:10:50])
set(gca, 'FontSize', 18, 'FontName', 'Franklin Gothic')

print -djpeg -r600 youngCarboDvR