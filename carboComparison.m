%% Carbo increases OAEs
maindir = pwd;


purple = [252,251,253; 239,237,245; 218,218,235;
    188,189,220; 158,154,200; 128,125,186;
    106,81,163; 84,39,143; 63,0,125]./255;
color1 = [140,81,10; 191,129,45; 223,194,125; % gold to teal
    246,232,195; 245,245,245; 199,234,229;
    128,205,193; 53,151,143; 1,102,94]./255;

carbo = {'Q405_R', 'Q405_L', 'Q403_R', 'Q403_L'};  %, 'Q381_R'};
YNH = {'Q421_L', 'Q421_R', 'Q412_L', 'Q412_R', 'Q414_L', 'Q414_R', 'Q415_L'};

groups = {YNH, carbo};
labels = {'YNH','carbo'};

avg_oae = zeros(length(groups), 512);
avg_nf = zeros(length(groups), 512);

figure;
figure_prop_name = {'PaperPositionMode', 'units', 'Position'};
figure_prop_val = {'auto', 'inches', [1 1 7.5 4.3]}; % xcor, ycor, xwid, yheight
set(gcf, figure_prop_name, figure_prop_val);
hold on;

for x = 1:length(groups)
    group = groups{x};
    oae_group = zeros(length(group), 512);
    nf_group = zeros(length(group), 512);
    
    for y = 1:length(group)
        hold on;
        load(sprintf('%s/Results/DP/DP_%s.mat', maindir, string(group(y))))
        
        oae = db(abs(res.complex.oae).*res.multiplier);
        nf = db(abs(res.complex.nf).*res.multiplier);
        f2 = res.f.f2/1000;
        
        if x == 1
            color = color1(6,:);
        else
            color = purple(3,:);
        end
        plot(f2, oae, 'linew', 3, 'color', color)
        plot(f2, nf, '--', 'linew', 2, 'color', color)
        
        oae_group(y,:) = oae';
        nf_group(y,:) = nf';
        
        set(gca, 'FontSize', 14, 'XScale', 'log')
        xticks([0.5, 1, 2, 4, 8, 16])
        xlim([0.5, 16])
        ylim([-40, 50])
        yticks([-30, -15 0, 15, 30, 45])
        xlabel('Frequency (kHz)')
        ylabel('Amplitude (dB SPL)')
        title('DPOAE')
        set(gca, 'XScale', 'log', 'FontSize', 18, 'FontName', 'Franklin Gothic')
        
        
    end
    avg_oae(x,:) = mean(oae_group,1);
    std_oae(x,:) = std(oae_group,1);
    avg_nf(x,:) = mean(nf_group,1);
    std_nf(x,:) = std(nf_group,1);
    
end

p1 = plot(f2, avg_oae(1,:), 'color', color1(8,:), 'linew', 3);
p2 = plot(f2, avg_oae(2,:), 'color', purple(8,:), 'linew', 3);

plot(f2, avg_nf(1,:), '--', 'color', color1(8,:), 'linew', 2)
plot(f2, avg_nf(2,:), '--', 'color', purple(8,:), 'linew', 2)
%legend([p1, p2], {'Unexposed', 'Carboplatin'}, 'location', 'Northwest');
print -djpeg -r600 carboDPOAE

%% Carbo increases OAEs
maindir = pwd;

carbo = {'Q405_R', 'Q405_L', 'Q403_R', 'Q403_L'};  %, 'Q381_R'};
YNH = {'Q421_L', 'Q421_R', 'Q412_L', 'Q412_R', 'Q414_L', 'Q414_R', 'Q415_L'};

groups = {YNH, carbo};
labels = {'YNH','carbo'};

avg_oae = zeros(length(groups), 512);
avg_nf = zeros(length(groups), 512);

figure;
figure_prop_name = {'PaperPositionMode', 'units', 'Position'};
figure_prop_val = {'auto', 'inches', [1 1 7.5 4.3]}; % xcor, ycor, xwid, yheight
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
        
        if x == 1
            color = color1(6,:);
        else
            color = purple(3,:);
        end
        plot(f2, oae, 'linew', 3, 'color', color)
        plot(f2, nf, '--', 'linew', 2, 'color', color)
        
        oae_group(y,:) = oae';
        nf_group(y,:) = nf';
        
        set(gca, 'FontSize', 14, 'XScale', 'log')
        xticks([0.5, 1, 2, 4, 8, 16])
        xlim([0.5, 16])
        ylim([-40, 50])
        yticks([-30, -15 0, 15, 30, 45])
        xlabel('Frequency (kHz)')
        ylabel('Amplitude (dB SPL)')
        title('SFOAE')
        set(gca, 'XScale', 'log', 'FontSize', 18, 'FontName', 'Franklin Gothic')
        
        
    end
    avg_oae(x,:) = mean(oae_group,1);
    std_oae(x,:) = std(oae_group,1);
    avg_nf(x,:) = mean(nf_group,1);
    std_nf(x,:) = std(nf_group,1);
    
end

p1 = plot(f2, avg_oae(1,:), 'color', color1(8,:), 'linew', 3);
p2 = plot(f2, avg_oae(2,:), 'color', purple(8,:), 'linew', 3);

plot(f2, avg_nf(1,:), '--', 'color', color1(8,:), 'linew', 2)
plot(f2, avg_nf(2,:), '--', 'color', purple(8,:), 'linew', 2)
legend([p1, p2], {'Unexposed', 'Carboplatin'}, 'location', 'Northwest');
print -djpeg -r600 carboSFOAE