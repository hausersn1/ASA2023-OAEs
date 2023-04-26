function DPplot(res)

%% Plot 
figure; 
hold on;
plot(res.f2_fpl, db(abs(res.complex_f2_fpl)), 'linew', 2)
plot(res.f1_fpl, db(abs(res.complex_f1_fpl)), 'linew', 2);  
plot(res.f2_fpl, db(abs(res.complex_dp_epl)), 'linew', 1.5)
plot(res.f2_fpl, db(abs(res.complex_nf_epl)), '--', 'linew', 1.5)
xticks([500, 1000, 2000, 4000, 8000, 16000]); 
legend('F2', 'F1', 'DP', 'NF', 'location', 'Northwest')
title('DPOAEs', 'FontSize', 14)
subtitle(sprintf('Subject: %s, Ear: %s', res.subj, res.ear))
ylabel('Amplitude (dB FPL/EPL)')
xlabel ('F_2 Frequency (Hz)')
set(gca, 'XScale', 'log', 'FontSize', 14); 



end
