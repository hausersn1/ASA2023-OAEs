
maindir = pwd;

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
orange = [255,245,235; 254,230,206; 253,208,162;
    253,174,107; 253,141,60; 241,105,19;
    217,72,1; 166,54,3; 127,39,4]./255;

nh = { 'S343_R',  'S357_R', 'SH1_R', 'SBC_L', 'SBC_R', 'S343_L', 'S357_L', 'SH1_L'};
hl = { 'S354_R', 'S354_L', 'S359_L', 'S359_R', 'S356_R', 'S356_L'};

% Set groups to analyze
groups = {nh, hl};
labels = { 'NH', 'HL'};

avg_oae = zeros(length(groups), 512);
avg_nf = zeros(length(groups), 512);

figure;
figure_prop_name = {'PaperPositionMode', 'units', 'Position'};
figure_prop_val = {'auto', 'inches', [1 1 7 4]}; % xcor, ycor, xwid, yheight
set(gcf, figure_prop_name, figure_prop_val);
hold on;

for x = 1:length(groups)
    group = groups{x};
    oae_group = zeros(length(group), 512);
    nf_group = zeros(length(group), 512);
    
    for y = 1:length(group)
        hold on;
        load(sprintf('%s/Results/DP/DP_%s.mat', maindir, string(group(y))))
        
        if res.stim.speed < 0
            oae = db(abs(res.complex.oae(end:-1:1)).*res.multiplier);
            nf = db(abs(res.complex.nf(end:-1:1)).*res.multiplier);
            f2 = res.f.f2(end:-1:1)/1000;
        else
            oae = db(abs(res.complex.oae).*res.multiplier);
            nf = db(abs(res.complex.nf).*res.multiplier);
            f2 = res.f.f2/1000;
        end
        
        if x == 1 % NH
            oae_nh(y,:) = oae';
            nf_nh(y,:) = nf';
            color = color2(6,:);
        else
            oae_HL(y,:) = oae';
            nf_HL(y,:) = nf';
            color = orange(2,:);
        end
        plot(f2, oae, 'linew', 2, 'color', color)
        plot(f2, nf, '--', 'linew', 2,  'color', color)
        
        set(gca, 'FontSize', 14, 'XScale', 'log')
        xticks([0.5, 1, 2, 4, 8, 16])
        xlim([0.5, 16])
        ylim([-40, 30])
        yticks([-30, -15 0, 15, 30])
        xlabel('Frequency (kHz)')
        ylabel('Amplitude (dB EPL)')
        title('DPOAE')
        set(gca, 'XScale', 'log', 'FontSize', 18, 'FontName', 'Franklin Gothic')
        
        
    end
    
end

avg_oae(1,:) = mean(oae_nh,1);
avg_oae(2,:) = mean(oae_HL, 1);
avg_nf(1,:) = mean(nf_nh,1);
avg_nf(2,:) = mean(nf_HL, 1);

p1 = plot(f2, avg_oae(1,:),  'color', c3(10,:), 'linew', 3);
p2 = plot(f2, avg_oae(2,:), 'color', c3(3,:), 'linew', 3);

plot(f2, avg_nf(1,:), '--','color', c3(10,:), 'linew', 2)
plot(f2, avg_nf(2,:), '--','color', c3(3,:), 'linew', 2)

legend([p1,p2], {'Normal Hearing', 'Hearing Loss'}, 'location', 'Northeast')
print -djpeg -r600 humanDPOAE
%% Same for SF


avg_oae = zeros(length(groups), 512);
avg_nf = zeros(length(groups), 512);

figure;
figure_prop_name = {'PaperPositionMode', 'units', 'Position'};
figure_prop_val = {'auto', 'inches', [1 1 7 4]}; % xcor, ycor, xwid, yheight
set(gcf, figure_prop_name, figure_prop_val);
hold on;

for x = 1:length(groups)
    group = groups{x};
    oae_group = zeros(length(group), 512);
    nf_group = zeros(length(group), 512);
    
    for y = 1:length(group)
        hold on;
        load(sprintf('%s/Results/SF/SF_%s.mat', maindir, string(group(y))))
        
        oae = db(abs(res.complex.oae).*res.multiplier);
        nf = db(abs(res.complex.nf).*res.multiplier);
        f2 = res.f/1000;
        
        if x == 1 % NH
            color = color2(6,:);
        else
            color = orange(2,:);
        end
        plot(f2, oae, 'linew', 2, 'color', color)
        plot(f2, nf, '--', 'linew', 2,  'color', color)
        
        oae_group(y,:) = oae';
        nf_group(y,:) = nf';
        
        set(gca, 'FontSize', 14, 'XScale', 'log')
        xticks([0.5, 1, 2, 4, 8, 16])
        xlim([0.5, 16])
        ylim([-40, 30])
        yticks([-30, -15 0, 15, 30])
        %xlabel('Frequency (kHz)')
        ylabel('Amplitude (dB EPL)')
        title('SFOAE')
        set(gca, 'XScale', 'log', 'FontSize', 18, 'FontName', 'Franklin Gothic')
        
        
    end
    avg_oae(x,:) = mean(oae_group,1);
    std_oae(x,:) = std(oae_group,1);
    avg_nf(x,:) = mean(nf_group,1);
    std_nf(x,:) = std(nf_group,1);
    
end

p1 = plot(f2, avg_oae(1,:), 'color', c3(10,:), 'linew', 3);
p2 = plot(f2, avg_oae(2,:), 'color', c3(3,:), 'linew', 3);

plot(f2, avg_nf(1,:), '--','color', c3(10,:), 'linew', 2)
plot(f2, avg_nf(2,:), '--','color', c3(3,:), 'linew', 2)

print -djpeg -r600 humanSFOAE