clc;clear;
addpath(genpath('C:/Users/nnu02/Documents/MATLAB/spm12'));
addpath('C:/Users/nnu02/Documents/MATLAB/fieldtrip-20220827');
addpath('E:/02Data/03Utils/Functions/');
ft_defaults

subjects = ["04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "24" "25" "26" "27" "28" "29"];

folder = 'E:/Gian/GG_SensAtt_Prediction/02Data';

touch_events = [125 124 252]; %markers inside the EEG
touch_labels = ["noTouch" "Touch"];

move_events = [1, 3; 2, 4]; %labels of the conditions as saved in the log files
move_labels = ["Move", "Stay"];

EOG_modality = "SPM_vEOG_hEOG";%correct for both hEOG and vEOG in a reiterative manner

restart_preprocessing = 1; %set to 1 for re-running all the preprocessing even when files are already saved
skip_spm = 0;

for ID = subjects
    ID = char(ID);
    
    verbose = 0; %specify whether you want each preproc step printed in a ps file (yes = 1; no = 0);
    notch = 0;
    
    %Stuff for SPM
    subj_folder = strcat(folder, '/ID', ID, '/01EEG/');
    
    %copy paste the loc file
    loc_dir = strcat(folder, '/ID', ID, '/00Behavioural/');
    bdf_file = strcat(folder, '/ID', ID, '/01EEG/*ID', ID, '*.bdf');
    bdf_file = dir(bdf_file); 
    bdf_file = strcat(bdf_file.folder, '\', bdf_file.name); 
    
    out_folder = strcat(folder, '/ID', ID, '/01EEG/spm/');
    if ~exist(out_folder, 'dir')
        mkdir(out_folder)
    end
    out_file = strcat(folder, '/ID', ID, '/01EEG/spm/ID');
    
    %Stuff for FT
    caplocation = 'E:/02Data/03Utils/biosemi64.lay';
    neighbourslocation = 'E:/02Data/03Utils/biosemi64_neighb.mat';
    
    if EOG_modality == "uncorr"
        prefix = 'uncorr';
    elseif EOG_modality == "SPM_vEOG_hEOG"
        thresh_file = 'EOGthresh_vEOG_hEOG';
        eogchan = ["vEOG", "hEOG"];
        components_file = 'num_components_vEOG_hEOG';
        prefix = ["vEOG_t", "hEOG_t"];
    end

    cd(subj_folder) %jump in the right folder
    
    %% PART 2 Import data from .bdf datastructure into ft .mat format
    if EOG_modality == "SPM_vEOG_hEOG"
        if (double(~exist(strcat(subj_folder, strjoin(eogchan, ''), '_filtered_ID', ID, '.mat'))) + restart_preprocessing) > 0
            if ~exist(strcat(subj_folder, '/spm/hdMID', ID, '.mat')) | skip_spm == 0
                addpath(genpath('C:/Users/nnu02/Documents/MATLAB/spm12'));
                cd('spm')

                load('E:/Gian/GG_SensAtt_Prediction/02Data/ID100/avref_combEOG.mat');
        
                %% FIRST MAKE SURE THAT THERE ARE NO BAD CHANS BEFORE MONTAGE (reref)
                cfg = []; 
                cfg.dataset = bdf_file; 
                origin_file = ft_preprocessing(cfg);
                labels = origin_file.label;
        
                if isfile(strcat(subj_folder, 'bad_chan_manual.mat'))
                    load(strcat(subj_folder, 'bad_chan_manual.mat'));
                else
                    origin_file.trial{1,1}(69,:) = origin_file.trial{1,1}(65,:) - origin_file.trial{1,1}(66,:);
                    origin_file.trial{1,1}(70,:) = origin_file.trial{1,1}(67,:) - origin_file.trial{1,1}(68,:);
    
                    trl = [];
                    cfg                         = [];
                    cfg.dataset                 = bdf_file;
                    cfg.trialfun                = 'ft_trialfun_general'; % this is the default
                    cfg.trialdef.eventtype      = 'STATUS';
                    cfg.trialdef.eventvalue     = touch_events; % the values of the stimulus trigger for the three conditions
                    cfg.trialdef.prestim        = 3; % in seconds
                    cfg.trialdef.poststim       = 3; % in seconds
                    cfg = ft_definetrial(cfg);
            %         cfg = check_extra_trial(cfg);
            %         cfg = check_outliers(cfg);
                    trl = cfg.trl;
                    cfg = [];
                    cfg.trl = trl;
                    origin_epoched = ft_redefinetrial(cfg, origin_file);
             
                    cfg = []; 
                    cfg.preproc.demean = 'yes';
                    cfg.preproc.lpfilter = 'yes'; 
                    cfg.preproc.lpfreq = 45; 
                    cfg.ylim = [-20 20];
                    ft_databrowser(cfg, origin_file); 
        
                    cfg          = [];
                    cfg.method   = 'summary';
                    cfg.layout   = caplocation;  % for plotting individual trials
                    hand_cleaned   = ft_rejectvisual(cfg, origin_epoched);
                
                    bad_chan_manual = setdiff(labels, hand_cleaned.label);
                %     noisy_channels_manual = labels(noisy_channels_manual);
                    save(strcat(subj_folder, 'bad_chan_manual.mat'), "bad_chan_manual");
                end
        
                load(neighbourslocation)
                if length(bad_chan_manual) > 0
                    cfg               = [];
                    cfg.method = 'average';
                    cfg.badchannel    = bad_chan_manual;
                    cfg.neighbours = neighbours;
                %     cfg.neighbourdist = 4;
                %     cfg.elec = elec;
                    origin_file_corr = ft_channelrepair(cfg, origin_file);
                else
                    origin_file_corr = origin_file;
                end
        
                %D = spm_eeg_ft2spm(origin_file_corr, strcat('ID', ID));

                %% RUN MOST OF PREPROCESSING IN FIELDTRIP
                cfg = []; 
                cfg.resamplefs = 512;
                %cfg.detrend = 'no'; 
                origin_file_corr = ft_resampledata(cfg, origin_file_corr); 
                origin_file_corr.hdr.Fs = cfg.resamplefs; 

                cfg = []; 
                cfg.hpfilter = 'yes'; 
                cfg.hpfilttype = 'firws';
                cfg.hpfreq = 0.01; 
                origin_file_corr = ft_preprocessing(cfg, origin_file_corr); %HPF is done in fieldtrip 

                % cfg = [];
                % cfg.reref = 'yes';
                % cfg.refchannel = 'all';
                % cfg.refmethod = 'avg';
                % origin_file_corr = ft_preprocessing(cfg, origin_file_corr);
                
                %% THEN CONVERT OUR CORRECTED FILE INTO SPM FORMAT
                S = [];
                S.dataset = bdf_file;
                S.outfile = strcat(out_file, ID); 
                D = spm_eeg_convert(S);

                %downsample
                S = []; 
                S.D = D; 
                S.fsample_new = 512;
                S.prefix = 'dM';
                D = spm_eeg_downsample(S); 

                D(1:64,:) = origin_file_corr.trial{1,1}(1:64,:); %push the HPF data from fieldtrip into SPM structure and continue preprocessing from there. 

                %montage 
                S = []; 
                S.D = D; 
                S.mode = 'write';
                S.montage = montage;
                S.prefix = 'h';
                D = spm_eeg_montage(S); 

                %insert eeg default locations
                S = [];
                S.D = D;
                S.task = 'defaulteegsens';
                S.save = 1;
                D = spm_eeg_prep(S);
        
            else %if file already exists, then simply load it. 
                addpath(genpath('C:/Users/nnu02/Documents/MATLAB/spm12'));
                cd('spm')
                D = spm_eeg_load(strcat('hdMID', ID, '.mat'));
            end

            %% EOG Stuff: reiteratively correct for vEOG and hEOG in this order. 
            % this part of code also samples a group of frontal electrodes
            % before and after correction to check that it went well. 
            thresh = []; 
            num_components = []; 
            for eyes = 1:2
               if eyes == 1
                   D = spm_eeg_load(strcat('hdMID', ID, '.mat'));
               else
                   D = spm_eeg_load(strcat('vEOG_thdMID', ID, '.mat'));
               end
               if isfile(strcat(folder, '/ID', ID, '/01EEG/',thresh_file,'.mat'))
                    load(strcat(folder, '/ID', ID, '/01EEG/',thresh_file,'.mat'))
        
                    S = [];
                    S.D = D;
                    S.eogchan = {char(eogchan(eyes))};
                    S.stdthresh= thresh(eyes);
                    S.overwrite = 1;
                    S.direction_h_movement = 2;
                    D_ebf = spm_eeg_detect_eyeblinks_GIAN(S);
                else
                    thresh(eyes) = 2;
    
                    S = [];
                    S.D = D;
                    S.eogchan = {char(eogchan(eyes))};
                    S.stdthresh= thresh(eyes);
                    S.overwrite = 1;
                    S.direction_h_movement = 2;
                    D_ebf = spm_eeg_detect_eyeblinks_GIAN(S);
        
                    if eyes == 2
                        save(strcat(folder, '/ID', ID, '/01EEG/',thresh_file,'.mat'), 'thresh');
                    end
                end
                %export img
                exportgraphics(gcf, strcat(folder, '/EEG_group_data/ID', ID, '_', eogchan(eyes), 'epochs.png'))
                close all
                if eyes == 1
                    before_corr = squeeze(D_ebf([1:5 33:38 65 66],:,:)); 
                end
        
                %% Remove eye blinks
                if isfile(strcat(folder, '/ID', ID, '/01EEG/',components_file,'.mat'))
                    load(strcat(folder, '/ID', ID, '/01EEG/',components_file,'.mat'))
                    compute_eye_blink_components(D_ebf, num_components(eyes)); 
                else
                    num_components(eyes) = 2; 
                    compute_eye_blink_components(D_ebf, num_components(eyes)); 
                    if eyes == 2
                        save(strcat(folder, '/ID', ID, '/01EEG/', components_file,'.mat'), 'num_components')
                    end
                end

                exportgraphics(gcf, strcat(folder, '/EEG_group_data/ID', ID, '_', eogchan(eyes), 'components_topo.png'))
                close
                exportgraphics(gcf, strcat(folder, '/EEG_group_data/ID', ID, '_', eogchan(eyes), 'components_ERP.png'))
                close
        
                try
                    % remove any spatial confounds file if present the meeg object
                    S           = [];
                    S.D         = D;
                    S.method    = 'CLEAR';
                    D = spm_eeg_spatial_confounds(S);
                catch 
                end
                % add the spatial confound to the meeg object
                S           = [];
                S.D         = D;
                S.method    = 'SPMEEG';
                S.conffile  = 'ebf_conf.mat';
                D = spm_eeg_spatial_confounds(S);
                
                % correct for the spatial confounds (Berg and Scherg)
                S               = [];
                S.D             = D;
                S.correction    = 'Berg';
                S.prefix        = char(prefix(eyes));
                D = spm_eeg_correct_sensor_data_pia(S);
                
                close all

                if eyes == 1
                    after_vEOG_corr = squeeze(D([1:5 33:38 65 66],:,:));
                elseif eyes == 2
                    after_hEOG_corr = squeeze(D([1:5 33:38 65 66],:,:));
                end
            end

            %% plot before and after imgs
            load('vEOG_spikes.mat')
            vEOG_spikes = spikes; 
            load('hEOG_spikes.mat')
            hEOG_spikes = spikes; 

            %select a random point in the hEOG spikes
            if exist("spike_pos.mat")
                load("spike_pos.mat")
            else
                h_EOG_pos = 5;
                save("spike_pos.mat", "h_EOG_pos")
            end
            %select a matching vEOG spike that is close to that
            [value v_EOG_pos] = min(abs(vEOG_spikes - hEOG_spikes(h_EOG_pos)));

            center = round((vEOG_spikes(v_EOG_pos) + hEOG_spikes(h_EOG_pos))/2);
            window_spikes = [center-5000 : center+5000];

            %plot before correction
            multichanplot(before_corr(:,window_spikes)', 10000, 'ylim', [-60 60])
            exportgraphics(gcf, strcat(folder, '/EEG_group_data/ID', ID, '_before_EOGCorr.png'))
            close
            %plot after vEOG corr
            multichanplot(after_vEOG_corr(:,window_spikes)', 10000, 'ylim', [-60 60])
            exportgraphics(gcf, strcat(folder, '/EEG_group_data/ID', ID, '_after_vEOGCorr.png'))
            close
            %plot after hEOG corr
            multichanplot(after_hEOG_corr(:,window_spikes)', 10000, 'ylim', [-60 60])
            exportgraphics(gcf, strcat(folder, '/EEG_group_data/ID', ID, '_after_hEOGCorr.png'))
            close
           
            %% Convert spm D object into fieldtrip data format
            %Remove Matlab path otherwise it gets in contrast with fieldtrip proper
            %functions
            rmpath(genpath('C:/Users/nnu02/Documents/MATLAB/spm12'));
            
            cfg = [];
            cfg.dataset = char(strcat(strjoin([prefix(2) prefix(1)], ''), 'hdMID', ID, '.mat'));
            raw = ft_preprocessing(cfg);
            labels = raw.label(1:64);
            
            cfg = [];
            cfg.channel = labels;
            raw = ft_selectdata(cfg, raw);

            %delete superfluous spm files
            delete thdMID* thdMID* fvEOG_thdMID* fvEOG_thdMID* efvEOG_thdMID* comb* ecomb* noise_tool* enoise_tool* ID* MID* spike_pos.mat spikemat.mat dMID* hdMID*

            cd .. %go back in the main folder and save data
            save(strcat(strjoin(eogchan, ''), '_filtered_ID', ID, '.mat'), "raw");
        else %if file exists, simply load it
            load(strcat(subj_folder, strjoin(eogchan, ''), '_filtered_ID', ID, '.mat'))
        end
    elseif EOG_modality == "uncorr"
        cd('spm')
        rmpath(genpath('C:/Users/nnu02/Documents/MATLAB/spm12'));
        
        cfg = [];
        cfg.dataset = strcat('hdMID', ID, '.mat');
        raw = ft_preprocessing(cfg);
        labels = raw.label(1:64);
        
        cd .. %go back in the main folder
    end

    % %% HPF 0.05
    % cfg = [];
    % cfg.hpfilter = 'yes';
    % cfg.hpfreq = 0.5;
    % raw = ft_preprocessing(cfg, raw); 
    
    %% INSERT EVENTS (stimulation, not divided per conditions)
    cfg                         = [];
    cfg.dataset                 = bdf_file;
    cfg.trialfun                = 'ft_trialfun_general'; % this is the default
    cfg.trialdef.eventtype      = 'STATUS';
    cfg.trialdef.eventvalue     = touch_events; % the values of the stimulus trigger for the three conditions
    cfg.trialdef.prestim        = 3; % in seconds
    cfg.trialdef.poststim       = 3; % in seconds
    cfg = ft_definetrial(cfg);
    events_list_EEG = cfg.event;
    trl = cfg.trl;

    trl(:,1:3) = round(trl(:,1:3)/4); 

    cfg = [];
    cfg.trl = trl;
    data_to_clean = ft_redefinetrial(cfg, raw);

    %% REMOVE OUTLIERS BEHAVIOURAL TRIALS
    if ~exist(strcat(subj_folder, 'rejected_trials_RT.mat'))
        rejected_trials_RT = catch_RT_outliers(ID, events_list_EEG, data_to_clean);
        save("rejected_trials_RT.mat", "rejected_trials_RT")
    end

    %% ARTEFACTS REMOVAL
    if isfile(strcat(subj_folder, 'rejected_trials_manual.mat'))
        load(strcat(subj_folder, 'rejected_channels.mat'))
        load(strcat(subj_folder, 'rejected_trials_manual.mat'))

        tot_rejected_trials = unique([rejected_trials_RT rejected_trials]);

        cfg = []; 
        cfg.trials = 1:numel(data_to_clean.trialinfo);
        cfg.trials(tot_rejected_trials) = []; 
        cfg.channels = 1:64;
        cfg.channels(rejected_channels) = []; 
        data = ft_selectdata(cfg, data_to_clean); 
    else
        % old_reject = load(strcat(subj_folder, 'rejected_trials.mat'));
        % old_reject = old_reject.rejected_trials; 
        rejected_channels = [];
        rejected_trials = [];

        adjusted_trl = trl;
        adjusted_trl(:,1) = adjusted_trl(:,1) + 1536 - round(0.1*512);
        adjusted_trl(:,2) = adjusted_trl(:,2) - 1536 + round(0.5*512);
        adjusted_trl(:,3) = round(0.1*(-512)) * ones(length(trl),1);

        %just have a look at the data
        cfg = [];
        cfg.trl = adjusted_trl; 
        cfg.preproc.lpfilter = 'yes'; 
        cfg.preproc.lpfreq = 45; 
        artf = ft_databrowser(cfg, data_to_clean);

        %after you had a look at the data, repeat the procedure that you did
        %before. 
        cfg          = [];
        cfg.method   = 'summary';
        cfg.layout   = caplocation;  % for plotting individual trials
        data_clean   = ft_rejectvisual(cfg, data_to_clean);

        %same procedure as before, simply update the previously created vectors
        rejchan = setdiff(data_to_clean.label, data_clean.label);
        rejtrl = setdiff(data_to_clean.sampleinfo(:,1), data_clean.sampleinfo(:,1));
        if isempty(rejchan) == 0
            for j = 1:size(rejchan,1)
                rejected_channels = [rejected_channels string(rejchan{j,1})];
            end
        end
        if isempty(rejtrl) == 0
            for j=1:size(rejtrl,1)
                rejected_trials = [rejected_trials find(data_to_clean.sampleinfo(:,1) == rejtrl(j,1))];
            end
        end

        %perform proper trial rejection here
        load(strcat(subj_folder, 'rejected_trials_RT.mat'))
        tot_rejected_trials = unique([rejected_trials_RT rejected_trials]);
        cfg = []; 
        cfg.trials = 1:numel(data_to_clean.trialinfo);
        cfg.trials(tot_rejected_trials) = []; 
        cfg.channels = 1:64;
        cfg.channels(rejected_channels) = [];
        data = ft_selectdata(cfg, data_to_clean); 

        save('rejected_channels.mat', 'rejected_channels');
        save('rejected_trials_vEOG_hEOG_001.mat', 'rejected_trials');
%         save ('data_clean.mat', 'data_clean');
        clear data_clean filtered hand_cleaned raw_epoched raw 
    end

    %% PROPER INTERPOLATION
    %create yero channels where missing
    try
        load(neighbourslocation);
        [data badchan] = push_channels(data, labels);
        if length(badchan) > 0
            cfg               = [];
            cfg.badchannel    = data.label(badchan);
            cfg.neighbours = neighbours;
            data = ft_channelrepair(cfg, data);
        else
            disp("NO BAD CHANNEL SELECTED DURING BAD TRIAL SELECTION")
        end
    catch
    end

    %% REREF
    cfg = [];
    cfg.reref = 'yes';
    cfg.refchannel = 'all';
    cfg.refmethod = 'avg';
    data = ft_preprocessing(cfg, data);
    
    %% LPF
    cfg = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 45;
    cfg.lpfilttype = 'firws'; 
    data = ft_preprocessing(cfg, data);
    
    %% SAVE PREPROCESSED AND CLEANED TRIALS IN MAT FILE
    % save(strcat('preprocessed_final_', EOG_modality, 'march2024.mat'), 'data', '-v7.3');
    save(strcat('preprocessed_final_', EOG_modality, '.mat'), 'data', '-v7.3');
end

%% EXTRA FUNCTIONS - catch_RT_outliers
% This function will load the log file associated to the participant's ID
% and will determine which trials are outliers based on their response
% time.
function rejected_trials_RT = catch_RT_outliers(ID, events_list_EEG, data_to_clean)
    opts = delimitedTextImportOptions("NumVariables", 7);
    opts.DataLines = [4, Inf];
    opts.Delimiter = ["\t", " \t "];
    opts.VariableNames = ["Block", "Trial", "BlockType", "TrialType", "Stimulation", "Event_Name", "Time"];
    opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts = setvaropts(opts, ["Block", "Trial", "BlockType", "TrialType", "Stimulation", "Event_Name", "Time"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["Block", "Trial", "BlockType", "TrialType", "Stimulation", "Event_Name", "Time"], "EmptyFieldRule", "auto");
    Log = readmatrix(strcat("E:\Gian\GG_SensAtt_Prediction\02Data\ID", ID, "\00Behavioural\Experiment_Log_ID", ID, ".txt"), opts);
    for i = 1:length(Log)
        if contains(Log(i,1), "Start")
            Log(i,6) = Log(i,1);
            Log(i,7) = Log(i,2);
            Log(i,1) = Log(i+1,1);
            Log(i,2) = Log(i+1,2);
            Log(i,3) = Log(i+1,3);
            Log(i,4) = Log(i+1,4);
            Log(i,5) = Log(i+1,5);
        end
    end

    %first "correct" the event file with the Log file. This
    %makes sure that all events are in the correct order :)
    
    events_list_EEG_values = [events_list_EEG.value];

    events_list_EEG_values(events_list_EEG_values == 126) = [];
    events_list_EEG_values(events_list_EEG_values~=10 & events_list_EEG_values~=20 & events_list_EEG_values~=30 & events_list_EEG_values~=124 & events_list_EEG_values~=125) = []; 
    events_list_EEG_values = events_list_EEG_values'; 

    events_log_list = Log(:,6);
    index = startsWith(events_log_list,"Start_trial");
    events_log_list(index) = "10";
    events_log_list = strrep(events_log_list,'preGO','20');
    events_log_list = strrep(events_log_list,'GO','30');
    events_log_list = strrep(events_log_list,'Index_visual','125');
    events_log_list = str2double(events_log_list);
    events_log_list(Log(:,6) == "Index_visual" & Log(:,5) == "1",:) = 124; 

    %check for repeated ones in both files
    repeated_EEG = []; 
    for i = 2:length(events_list_EEG_values)
        if events_list_EEG_values(i) == events_list_EEG_values(i-1)
            repeated_EEG = [repeated_EEG i]; 
        end
    end
    events_list_EEG_values(repeated_EEG) = []; 

    repeated_LOG = []; 
    for i = 2:length(events_log_list)
        if events_log_list(i) == events_log_list(i-1)
            repeated_LOG = [repeated_LOG i]; 
        end
    end
    events_log_list(repeated_LOG) = []; 
    Log1 = Log; 
    Log1(repeated_LOG,:) = []; 

    i = 1;
    errors = [];
    events_list_new = events_list_EEG_values;
    while (i < length(events_log_list)) 
        if events_list_new(i) == events_log_list(i)
            i = i + 1;
        else
            idx = i;
            events_list_new = [events_list_new(1:length(events_list_new) < idx); events_log_list(i); events_list_new(1:length(events_list_new) >= idx)];
            errors = [errors i];
            i = i+2;
        end
    end
    if sum(events_list_new ~= events_log_list)>0 
        error('something went terribly wrong');
    end

    Log1(errors,:) = []; 
    list_conditions_nonCorr = Log1(Log1(:,6) == "Index_visual",[4 5 2 1 3]);

    list_conditions = list_conditions_nonCorr; 

    %check for correctness of data
    zeros_ones_EEG = abs(data_to_clean.trialinfo - 125);
    zeros_ones_Log = list_conditions_nonCorr(:,2); 
    if sum(zeros_ones_EEG ~= str2double(zeros_ones_Log)) > 0
        error('PORCAMADONNA IL LOG E` SBAGLIATO')
    else
        display('allesgut!')
    end

    for i = 1:size(events_list_EEG,2)
        if isempty(events_list_EEG(i).value)
            events_list_EEG(i).value = 0;
        end
    end
    movements = find(list_conditions(:,1) == "1" | list_conditions(:,1) == "3");
    touches = find([events_list_EEG.value] == 125 | [events_list_EEG.value] == 124 | [events_list_EEG.value] == 252);
    touches_move = touches(movements); 
    RTs_move = []; 
    for i = 1:length(touches_move)
        if events_list_EEG(touches_move(i)-1).value == 30
            RTs_move = [RTs_move events_list_EEG(touches_move(i)).sample - events_list_EEG(touches_move(i)-1).sample];
        end
    end

    RTs_all = []; 
    for i = 1:length(touches_move)
        if events_list_EEG(touches_move(i)-1).value == 30
            RTs_all = [RTs_all events_list_EEG(touches(i)).sample - events_list_EEG(touches(i)-1).sample];
        end
    end

    RTs_move = RTs_move/2048; %convert to seconds
    outliers_RT = find(isoutlier(RTs_move));
    rejected_trials_RT = []; 
    for i = 1:length(outliers_RT)
        rejected_trials_RT = [rejected_trials_RT find(touches == touches_move(outliers_RT(i)))];
    end

    false_starts = []; 
    false_starts = find(RTs_all < 0.050);
    rejected_trials_RT = unique([rejected_trials_RT false_starts]);
end


%% EXTRA FUNCTIONS - push channels
%This function will check missing channels, it will push new channels in
%the data filled with zeros (it will be needed for the interpolation
%procedure). 
function [restored badchanindx] = push_channels(data, label)
    [notmissing, dummy] = match_str(label, data.label);
    newtrial = cell(size(data.trial));
    for k = 1:numel(data.trial)
      newtrial{k} = zeros(numel(label), size(data.trial{k},2));
      newtrial{k}(notmissing,:) = data.trial{k};
    end
    goodchans   = false(numel(label),1);
    goodchans(notmissing) = true;
    badchanindx = find(goodchans==0);

    data.trial = newtrial; clear newtrial;
    data.label = label;
    restored = data;
end
