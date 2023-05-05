%% Converts Raw data to Analyzed format Needed for Further Analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enter the subjects you need to run 
% Will automatically check for both SF and DPOAEs
% 
% Samantha Hauser, May 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Enter subject (human or chin) here: 
%subjects = { 'Q381' 'Q406', 'Q351', 'Q298', 'Q403', 'Q422', 'Q407', 'Q368', 'Q412'}; 

%subjects = { 'S343', 'S354', 'S356', 'S357', 'S359' }; 
subjects = { 'Q421' }; 

ears = { 'R', 'L'}; 
for i = 1:length(subjects)
    for j = 1:length(ears)
        s = char(subjects(i)); 
        e = char(ears(j)); 
        %DPanalysis(s, e); 
        SFanalysis(s,e); 
    end
end



