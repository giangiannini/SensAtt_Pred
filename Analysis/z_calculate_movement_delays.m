%% TAKE OUT DIFFERENCE FROM GO and Touch
clc;clear;
addpath(genpath('C:/Users/nnu02/Documents/MATLAB/spm12'));
addpath('C:/Users/nnu02/Documents/MATLAB/fieldtrip-20220827');
addpath('E:/02Data/03Utils/Functions/');
ft_defaults

subjects = ["04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "24" "25" "26" "27" "28" "29"];
folder = 'E:/Gian/GG_SensAtt_Prediction/02Data';

touch_events = [124 125 252]; %ACHTUNG! PReGO EVENTS INSTEAD OF TOUCH

diffs = []; 
for sogg = 1:length(subjects)
    ID = char(subjects(sogg));
    
    %Stuff for SPM
    subj_folder = strcat(folder, '/ID', ID, '/01EEG/');
    cd(subj_folder)
    load("rejected_trials_manual.mat")
    load('rejected_trials_RT.mat')
    
    %copy paste the loc file
    loc_dir = strcat(folder, '/ID', ID, '/00Behavioural/');
    bdf_file = strcat(folder, '/ID', ID, '/01EEG/*ID', ID, '*.bdf');
    bdf_file = dir(bdf_file); 
    bdf_file = strcat(bdf_file.folder, '\', bdf_file.name); 
    
    %create trl for Touch
    trl = [];
    cfg                         = [];
    cfg.dataset                 = bdf_file;
    cfg.trialfun                = 'ft_trialfun_general'; % this is the default
    cfg.trialdef.eventtype      = 'STATUS';
    cfg.trialdef.eventvalue     = touch_events; % the values of the stimulus trigger for the three conditions
    cfg.trialdef.prestim        = 0; % in seconds
    cfg.trialdef.poststim       = 1; % in seconds
    cfg = ft_definetrial(cfg);
    trl_Touch = cfg.trl;
    trl_Touch(unique([rejected_trials rejected_trials_RT]),:) = []; %eye blinks will be looked for only in the "clean" channels. 

    trl_GO = []; 
    samples_Touch = [cfg.event.sample];
    for i = 1:length(trl_Touch)
        trl_GO = [trl_GO samples_Touch(find(samples_Touch == trl_Touch(i,1)) - 1)];
    end

    diff = []; 
    diff = trl_Touch(:,1) - trl_GO';
    save("movement_delays.mat", "diff")
    diffs = [diffs mean(diff)];
end