%% Converts Raw data to Analyzed format Needed for Further Analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enter the subjects you need to run
% Will automatically check for both SF and DPOAEs
%
% Samantha Hauser, May 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Enter subject (human or chin) here:
%subjects = { 'Q381', 'Q406', 'Q351', 'Q298', 'Q421', 'Q414', 'Q415', 'Q402', ...
%    'Q403', 'Q422', 'Q407', 'Q368', 'Q412', 'Q405'};
% subjects = {'S356', 'S343', 'S354', 'S357'};
subjects = {'S359'};

%subjects = { 'S343', 'S354', 'S356', 'S357', 'S359' };
%subjects = {'Q415', 'Q414'};

ears = { 'L', 'R'};
for i = 1:length(subjects)
    for j = 1:length(ears)
        s = char(subjects(i));
        e = char(ears(j));
        DPanalysis(s, e);
        SFanalysis(s,e);
    end
end



