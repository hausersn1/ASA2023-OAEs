
% Subjects to set in format 'SubjID_ear')

nh = { 'S343_R',  'S357_R', 'SH1_R', 'SBC_L', 'SBC_R', 'S343_L', 'S357_L', 'SH1_L'};
hl = { 'S354_R', 'S354_L', 'S359_L', 'S359_R', 'S356_R', 'S356_L'};

% Set groups to analyze
groups = {hl, nh};
labels = {'HL', 'NH'};

% Colors
% colors light to dark
c3 = [158,1,66; 213,62,79; 244,109,67;
    253,174,97; 254,224,139; 255,255,191;
    230,245,152; 171,221,164; 102,194,165;
    50,136,189; 94,79,162]./255;
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
        
        if res.stim.speed < 0
            dpoae_full = db(abs(res.complex.oae(end:-1:1)).*res.multiplier);
            dpnf_full = db(abs(res.complex.nf(end:-1:1)).*res.multiplier);
            f2 = res.f.f2(end:-1:1)/1000;
        else
            dpoae_full = db(abs(res.complex.oae).*res.multiplier);
            dpnf_full = db(abs(res.complex.nf).*res.multiplier);
            f2 = res.f.f2/1000;
        end
        
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
        hl_dpoae = dpoae_w_group;
        hl_sfoae = sfoae_w_group;
    elseif x == 2
        nh_dpoae = dpoae_w_group;
        nh_sfoae = sfoae_w_group;
    end
end

%% plot human
figure;
figure_prop_name = {'PaperPositionMode', 'units', 'Position'};
figure_prop_val = {'auto', 'inches', [1 1 6 5]}; % xcor, ycor, xwid, yheight
set(gcf, figure_prop_name, figure_prop_val);
hold on;
sz = 80;
plot([-40, 50], [-40, 50], '--', 'color', black(6, :), 'linew', 2)
for i=1:size(nh_sfoae,1)
    scatter(nh_sfoae(i,:), nh_dpoae(i,:), sz, 'Marker', m{2},...
        'MarkerFaceColor', blue(i+1,:), 'MarkerEdgeColor', 'k')
    if i == size(nh_sfoae,1)-2
        p1 = scatter(nh_sfoae(i,:), nh_dpoae(i,:), sz, 'Marker', m{2},...
            'MarkerFaceColor', blue(i+1,:), 'MarkerEdgeColor', 'k');
    end
end
for i=1:size(hl_sfoae,1)
    scatter(hl_sfoae(i,:), hl_dpoae(i,:), sz, 'Marker', m{2},...
        'MarkerFaceColor', orange(i+2,:), 'MarkerEdgeColor', 'k')
    if i == size(hl_sfoae,1)-2
        p2= scatter(hl_sfoae(i,:), hl_dpoae(i,:), sz, 'Marker', m{2},...
            'MarkerFaceColor', orange(i+2,:), 'MarkerEdgeColor', 'k');
    end
end
title('Human')
xlabel('SFOAE (dB EPL)')
ylabel('DPOAE (dB EPL)')
legend([p1, p2], {'Normal Hearing', 'Hearing Loss'}, 'location', 'Northwest')
set(gca, 'FontSize', 18, 'FontName', 'Franklin Gothic')
xlim([-40, 40])
ylim([-40, 40])
xticks([-40:20:40])
yticks([-40:20:40])
print -djpeg -r600 HumanDvR


%% Plotting Results -- x and o for different half octave bands
figure;
figure_prop_name = {'PaperPositionMode', 'units', 'Position'};
figure_prop_val = {'auto', 'inches', [1 1 6.5 4.5]}; % xcor, ycor, xwid, yheight
set(gcf, figure_prop_name, figure_prop_val);
hold on;

for i=1:size(hl_dpoae,1)
    plot(centerFreqs, hl_dpoae(i,:), 'x-', 'color', orange(2,:), 'linew', 2, 'MarkerSize', 10)
end

for i=1:size(nh_dpoae,1)
    plot(centerFreqs, nh_dpoae(i,:), 'x-', 'color', blue(2,:), 'linew', 2, 'MarkerSize', 10)
end

p1 = plot(centerFreqs, avg_dpoae(2,:), 'x-', 'color', blue(6,:), 'linew', 2, 'MarkerSize', 10)
p2 = plot(centerFreqs, avg_dpoae(1,:), 'x-','color', orange(6,:), 'linew', 2, 'MarkerSize', 10)

xticks([0.5, 1, 2, 4, 8, 16])
xlim([0.5, 16])
ylim([-40, 40])
yticks([-40:20:40])
title('Group Average DPOAE')
ylabel('Amplitude (dB EPL)')
xlabel('Frequency (kHz)')
legend([p1, p2], {'NH-DP', 'HL-DP'},'location', 'Northwest')
set(gca, 'XScale', 'log', 'FontSize', 18)

print -djpeg -r600 humanDPbyFreq


figure;
figure_prop_name = {'PaperPositionMode', 'units', 'Position'};
figure_prop_val = {'auto', 'inches', [1 1 6.5 4.5]}; % xcor, ycor, xwid, yheight
set(gcf, figure_prop_name, figure_prop_val);
hold on;


for i=1:size(hl_sfoae,1)
    plot(centerFreqs, hl_sfoae(i,:), 'o-', 'color', orange(3,:), 'linew', 2, 'MarkerSize', 10)
end

for i=1:size(nh_sfoae,1)
    plot(centerFreqs, nh_sfoae(i,:), 'o-', 'color', blue(3,:), 'linew', 2, 'MarkerSize', 10)
end

p3 = plot(centerFreqs, avg_sfoae(2,:), 'o-','color', blue(7,:), 'linew', 2, 'MarkerSize', 10)
p4 = plot(centerFreqs, avg_sfoae(1,:), 'o-','color', orange(7,:),  'linew', 2, 'MarkerSize', 10)

xticks([0.5, 1, 2, 4, 8, 16])
xlim([0.5, 16])
ylim([-40, 40])
yticks([-40:20:40])
title('Group Average SFOAE')
ylabel('Amplitude (dB EPL)')
xlabel('Frequency (kHz)')
legend([p3, p4], {'NH-SF', 'HL-SF'},'location', 'Southwest')
set(gca, 'XScale', 'log', 'FontSize', 18)

print -djpeg -r600 humanSFbyFreq