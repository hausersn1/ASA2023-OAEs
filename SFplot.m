function SFplot(res)

%% Plot 
figure; 
hold on;
plot(res.f_epl, db(abs(res.complex_sf_epl)), 'linew', 1.5)
plot(res.f_epl, db(abs(res.complex_nf_epl)), '--', 'linew', 1.5)
xticks([500, 1000, 2000, 4000, 8000, 16000]); 
legend('SF', 'NF', 'location', 'Northwest')
title('SFOAEs', 'FontSize', 14)
subtitle(sprintf('Subject: %s, Ear: %s', res.subj, res.ear))
ylabel('Amplitude (dB FPL/EPL)')
xlabel ('Frequency (Hz)')
set(gca, 'XScale', 'log', 'FontSize', 14); 

end