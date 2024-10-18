function [] = first_level_function(ID)
    % This function takes as input the ID (as character) of the participant
    % and fits a GLM through SPM to estimate B values as 3D images that are
    % then used at the second level. 
    % The preprocessed EEG file needs to be stored in the correct location
    % (see folder and subj_folder paths). 
    % The EEG epoched and preprocessed file is first divided into
    % conditions (12 conditions total) and saved into 4D images (space x
    % space x time x trials) that are then used to estimate the first level
    % GLM. 

    folder = 'E:/Gian/GG_SensAtt_Prediction/02Data';
    touch_events = [125 124 252];
    touch_labels = ["noTouch" "Touch"];
    
    move_events = [1, 3; 2, 4];
    move_labels = ["Move", "Stay"];
    
    prob_events = ["25/75", "50/50", "75/25"];
    prob_labels = ["LowProb", "EqualProb", "HighProb"];
    
    % CLEAR DATA? 
    CLEAR_FIRST = 1; 
    CLEAR_IMAGES = 0; 
    OVERWRITE = 0; 
    
    if CLEAR_IMAGES == 1
        subj_folder = strcat(folder, '/ID', ID, '/01EEG/');
        delete(strcat(subj_folder, 'spm/evEOG_thdMID', ID, '/*'))
    end
   
    ID = char(ID);
    
    %Stuff for SPM
    subj_folder = strcat(folder, '/ID', ID, '/01EEG/');

    %copy paste the loc file
    loc_dir = strcat(folder, '/ID', ID, '/00Behavioural/');
    bdf_file = strcat(folder, '/ID', ID, '/01EEG/*ID', ID, '*.bdf');
    bdf_file = dir(bdf_file); 
    bdf_file = strcat(bdf_file.folder, '\', bdf_file.name); 

    first_level_folder = strcat(folder, '/ID', ID, '/021stLevel/');
    if ~exist(first_level_folder, 'dir')
        mkdir(first_level_folder)
    end

    results_first_level = strcat(folder, '/ID', ID, '/021stLevel/00Results1stLevel');
    if CLEAR_FIRST == 1
        delete(strcat(results_first_level, '/*'))
    end
    if ~exist(results_first_level, 'dir')
        mkdir(results_first_level)
    end

    %Stuff for FT
    caplocation = 'E:/02Data/03Utils/biosemi64.lay';
    neighbourslocation = 'E:/02Data/03Utils/biosemi64_neighb.mat';

    cd(subj_folder) %jump in the right folder
    
    %% IMPORT RAW DATASET
    load('preprocessed_final_SPM_vEOG_hEOGApril2024_fully_manual.mat');
    load('rejected_trials_manual.mat'); 
    load('rejected_trials_RT.mat');

    %correct for missing samples (fieldtrip accepts non-homogenous trials,
    %SPM doesn't. Sometimes 1 datapoint is missed in fieldtrip for rounding
    %errors. Here we correct it so that SPM is happy. 
    for i = 1:length(data.sampleinfo)
        if data.sampleinfo(i,2) - data.sampleinfo(i,1) ~= 3072
            excess = 3072 - (data.sampleinfo(i,2) - data.sampleinfo(i,1));
            data.trial{i}(:,end+(1:excess)) = data.trial{i}(:,end-excess+1:end);
            data.sampleinfo(i,2) = data.sampleinfo(i,1)+3072; 
            data.time{i}(:,3073) = 3;
        end
    end

    cfg                         = [];
    cfg.dataset                 = bdf_file;
    cfg.trialfun                = 'ft_trialfun_general'; % this is the default
    cfg.trialdef.eventtype      = 'STATUS';
    cfg.trialdef.eventvalue     = touch_events; % the values of the stimulus trigger for the three conditions
    cfg.trialdef.prestim        = 3; % in seconds
    cfg.trialdef.poststim       = 3; % in seconds
    cfg = ft_definetrial(cfg);
    trl = cfg.trl;

    events_list_EEG = [cfg.event.value];
    samples_list_EEG = [cfg.event.sample]; 

    to_delete = []; 
    for i = 1:length(samples_list_EEG)
        if strcmp(cfg.event(i).type, 'STATUS')
            %all gut
        else
            to_delete = [to_delete i];
        end
    end
    samples_list_EEG(to_delete) = [];
    
    %% Import Log
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
    samples_list_EEG(events_list_EEG == 126) = [];
    events_list_EEG(events_list_EEG == 126) = [];
    samples_list_EEG(events_list_EEG~=10 & events_list_EEG~=20 & events_list_EEG~=30 & events_list_EEG~=124 & events_list_EEG~=125) = []; 
    events_list_EEG(events_list_EEG~=10 & events_list_EEG~=20 & events_list_EEG~=30 & events_list_EEG~=124 & events_list_EEG~=125) = []; 
    samples_list_EEG = samples_list_EEG'; 
    events_list_EEG = events_list_EEG'; 

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
    for i = 2:length(events_list_EEG)
        if events_list_EEG(i) == events_list_EEG(i-1)
            repeated_EEG = [repeated_EEG i]; 
        end
    end
    events_list_EEG(repeated_EEG) = []; 
    samples_list_EEG(repeated_EEG) = []; 

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
    events_list_new = events_list_EEG;
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
    list_conditions(unique([rejected_trials rejected_trials_RT]),:) = []; 
    
    %check for correctness of data
    zeros_ones_EEG = abs(data.trialinfo - 125); 
    zeros_ones_Log = list_conditions_nonCorr(:,2); 
    zeros_ones_Log(unique([rejected_trials rejected_trials_RT])) = []; 
    zerosss = [zeros_ones_EEG, zeros_ones_Log];
    if sum(zeros_ones_EEG ~= str2double(zeros_ones_Log)) > 0
        error('PORCAMADONNA IL LOG E` SBAGLIATO')
    else
        display('allesgut!')
    end

    %% Import participant's responses
    opts = delimitedTextImportOptions("NumVariables", 2);
    opts.DataLines = [3, Inf];
    opts.Delimiter = "\t";
    opts.VariableNames = ["Var1", "Time"];
    opts.SelectedVariableNames = "Time";
    opts.VariableTypes = ["string", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
    opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");
    Experimentresponse = readtable(strcat("E:\Gian\GG_SensAtt_Prediction\02Data\ID", ID, "\00Behavioural\Experiment_response_ID", ID, ".txt"), opts);
    Experimentresponse = table2array(Experimentresponse);
    
    %import TrialTable from experiment
    opts = delimitedTextImportOptions("NumVariables", 13);
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    opts.VariableNames = ["Block", "Block_num", "TrialInBlock", "Trial", "Block_type", "Trial_type", "Stimulation", "ITIs", "Speed", "Completed", "Attempts", "Skipped", "TrialTime"];
    opts.VariableTypes = ["double", "double", "double", "double", "categorical", "double", "double", "double", "double", "categorical", "double", "categorical", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts = setvaropts(opts, ["Block_type", "Completed", "Skipped"], "EmptyFieldRule", "auto");
    TrialTable = readtable(strcat("E:\Gian\GG_SensAtt_Prediction\02Data\ID", ID, "\00Behavioural\ID", ID, "TrialTable.csv"), opts);
    clear opts

    ExperimentCondition = string(TrialTable{TrialTable.TrialInBlock == 0, 5});
    ExperimentCondition(ExperimentCondition == "25/75") = "25"; 
    ExperimentCondition(ExperimentCondition == "50/50") = "50"; 
    ExperimentCondition(ExperimentCondition == "75/25") = "75"; 
    ExperimentCondition = str2double(ExperimentCondition); 

    Correctness = ExperimentCondition == Experimentresponse; 
    
    for i = 1:length(Correctness)
        list_conditions(list_conditions(:,4) == string(i-1),6) = string(Experimentresponse(i));
        list_conditions(list_conditions(:,4) == string(i-1),7) = string(double(Correctness(i)));
    end

    list_conditions(list_conditions(:,6) == "25",6) = "25/75"; 
    list_conditions(list_conditions(:,6) == "50",6) = "50/50"; 
    list_conditions(list_conditions(:,6) == "75",6) = "75/25"; 

    %% CALCULATE RT for each trial 
    touches = find(events_list_EEG == 124 | events_list_EEG == 125);
    samples_baseline(:,2) = samples_list_EEG(touches);
    samples_baseline(:,1) = samples_list_EEG(touches-1);
    samples_baseline(unique([rejected_trials rejected_trials_RT]),:) = []; 

    RTs = (samples_baseline(:,2) - samples_baseline(:,1))/2048;

    list_conditions(:,10) = string(RTs); 

    %% EPOCH and CONVERT TO IMAGES
    cd(first_level_folder)
    if isfile(strcat(subj_folder, 'spm/evEOG_thdMID', ID, '/condition_noTouch_Move_EqualProb.nii')) && OVERWRITE == 0
        disp('images already extracted for this condition')
    else
        trl = []; 
        trl = data.sampleinfo;
        trl(:,3) = repmat([-1536],length(trl),1);

        D_notpreprocessed = spm_eeg_load(strcat(subj_folder, 'spm/vEOG_thdMID', ID, '.mat'));

        S = []; 
        S.D = D_notpreprocessed; 
        S.bc = 0; 
        S.trl = trl; 
        S.conditionlabels = 'undefined';
        D_epoched = spm_eeg_epochs(S);
        %now prepare the full list of conditions
        conds = [];
        for i = 1:2 %touch and notouch
            for y = 1:2 %move and no move
                for n = 1:3 %prob
                    name = strcat(touch_labels(i), '_',  move_labels(y), '_', prob_labels(n)); 
                    %cfg.trials = find((str2double(list_conditions(:,3)) < 13) & (str2double(list_conditions(:,2)) == y-1));
                    %cfg.trials = find((list_conditions(:,5) == "25/75") & (str2double(list_conditions(:,1)) == move_events(i,1) | str2double(list_conditions(:,1)) == move_events(i,2)));
                    %cfg.trials = find((list_conditions(:,5) == "25/75") & (str2double(list_conditions(:,2)) == i-1));
                    trials = find((list_conditions(:,6) == prob_events(n)) & (str2double(list_conditions(:,2)) == i-1) & (str2double(list_conditions(:,1)) == move_events(y,1) | str2double(list_conditions(:,1)) == move_events(y,2)));
                    [conds{trials}] = deal(char(name));
                end
            end
        end

        if BASELINE == 1
            cfg = []; 
            cfg.demean = 'yes'; 
            cfg.baselinewindow = [-0.01 -0.005];
            data = ft_preprocessing(cfg, data);
        end

        D_epoched = conditions(D_epoched, ':', conds);
        D_epoched(1:64,:,:) = cat(3,data.trial{:});
        D_epoched.save(); 
        
        %and convert everything 
        S = [];
        S.D = D_epoched;
        S.timewin = [-50 500];
        S.mode = 'scalp x time';
        S.channels = {'EEG'}; 
        prova = spm_eeg_convert2images_jh(S);
    
       %  % SMOOTH!
       %  order_scans = ["noTouch_Stay_LowProb", "noTouch_Stay_EqualProb", "noTouch_Stay_HighProb",...
       % "noTouch_Move_LowProb", "noTouch_Move_EqualProb", "noTouch_Move_HighProb"...
       % "Touch_Stay_LowProb", "Touch_Stay_EqualProb", "Touch_Stay_HighProb",...
       % "Touch_Move_LowProb", "Touch_Move_EqualProb", "Touch_Move_HighProb"]; 
       %  scans_total = [];
       %  for sss = 1:length(order_scans)
       %      scans_total = [scans_total; cellstr(spm_select('expand', strcat(subj_folder, 'spm/evEOG_thdMID', ID, '/condition_', order_scans(sss),'_', folder_name_baseline, '.nii')))];
       %  end
       % 
       %  smoothing_kernel = [5 5 5];
       %  matlabbatch{1}.spm.spatial.smooth.data = scans_total;      
       %  matlabbatch{1}.spm.spatial.smooth.fwhm = smoothing_kernel;
       %  matlabbatch{1}.spm.spatial.smooth.dtype = 0;
       %  matlabbatch{1}.spm.spatial.smooth.im = 1;
       %  matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        % spm_jobman('run', matlabbatch);
        % clear matlabbatch
    end
    
    %% RUN COMPUTATIONAL MODEL
    % if COVARIATE == 1
    %     missing = []; 
    %     for i = 1:length(list_conditions_nonCorr)-1
    %         if str2double(list_conditions_nonCorr(i,3)) == 24
    %             %do nothing and skip
    %         elseif str2double(list_conditions_nonCorr(i,3)) == str2double(list_conditions_nonCorr(i+1,3))-1
    %             %all good
    %         elseif str2double(list_conditions_nonCorr(i,3)) ~= str2double(list_conditions_nonCorr(i+1,3))-1
    %             missing = [missing i];
    %             if str2double(list_conditions_nonCorr(i+1,3)) - str2double(list_conditions_nonCorr(i,3)) > 2 
    %                 missing = [missing i];
    %             end
    %         end
    %     end
    % 
    %     for i = 1:length(missing)
    %         miss_pos = missing(i) + (i-1); 
    %         string_to_insert = ["NaN", "NaN", string(str2double(list_conditions_nonCorr(miss_pos,3))+1), list_conditions_nonCorr(miss_pos,4), list_conditions_nonCorr(miss_pos,5)]; 
    %         list_conditions_nonCorr = [list_conditions_nonCorr(1:miss_pos,:); string_to_insert; list_conditions_nonCorr(miss_pos+1:end,:)]; 
    %     end
    % 
    %     %take out observations
    %     num_blocks = max(str2double(list_conditions_nonCorr(:,4))) + 1; 
    %     observations = []; 
    %     for i=1:num_blocks
    %         observations(i,:) = str2double(list_conditions_nonCorr(str2double(list_conditions_nonCorr(:,4)) == i-1, 2))';
    %     end
    % 
    %     Y_posterior = []; 
    %     Y_prior = []; 
    %     PE = []; 
    %     for i = 1:size(observations,1)
    %         a = 1; 
    %         b = 1; 
    %         theta = linspace(0,1,100); 
    %         binWidth = diff(theta);
    %         binWidth = [binWidth(end),binWidth];
    %         for y = 1:size(observations,2)
    %             if isnan(observations(i,y))
    %                 old_a = a;
    %                 old_b = b; 
    %                 Y_prior{i,y} = Y_prior{i,y-1};
    %                 Likelihood{i,y} = Likelihood{i,y-1};
    %                 Y_posterior{i,y} = Y_posterior{i,y-1};
    %                 predictive_surprise(i,y) = predictive_surprise(i,y-1);
    %                 bayesian_surprise(i,y) = bayesian_surprise(i,y-1);
    %                 confidence_surprise(i,y) = confidence_surprise(i,y-1);
    %                 PE(i,y) = PE(i,y-1);
    %                 wPE(i,y) = wPE(i,y-1);
    %             else
    %                 old_a = a;
    %                 old_b = b; 
    %                 Y_prior{i,y} = betapdf(theta, old_a, old_b);
    % 
    %                 Likelihood{i,y} = betapdf(theta, observations(i,y)+1, 1-observations(i,y)+1);
    % 
    %                 % Update posterior
    %                 a = old_a + observations(i,y); 
    %                 b = old_b + (1 - observations(i,y)); 
    %                 Y_posterior{i,y} = betapdf(theta, a, b); 
    % 
    %                 % Entropy
    %                 temp = Y_prior{i,y}.*log(Y_prior{i,y}./binWidth);
    %                 temp(isnan(temp)) = 0;
    %                 entropy(i,y) = -sum(temp);
    % 
    %                 % Update naive posterior
    %                 a_naive = 1 + observations(i,y); 
    %                 b_naive = 1 + (1 - observations(i,y)); 
    %                 Y_posterior_naive{i,y} = betapdf(theta, a_naive, b_naive);
    % 
    %                 % Predictive Surprise (Shannon surprise)
    %                 post_pred = [old_b/(old_a+old_b), old_a/(old_a+old_b)];
    %                 post_pred = post_pred(observations(i,y)+1);
    %                 predictive_surprise(i,y) = -log(post_pred);
    % 
    %                 % Bayesian surprise
    %                 %bayesian_surprise(i) = KLDiv(Y_prior{i}, Y_posterior{i}); 
    %                 bayesian_surprise(i,y) = KL_DIV([old_a old_b], [a b]); 
    % 
    %                 % Confidence-corrected
    %                 % KL(posterior, naive_posterior)
    %                 % naive posterior: flat prior (without any observation) updated with
    %                 % the new observation
    %                 %confidence_surprise(i) = KLDiv(Y_prior{i}, Y_posterior_naive{i}); 
    %                 confidence_surprise(i,y) = KL_DIV([old_a old_b], [a_naive b_naive]); 
    % 
    %                 %PE
    %                 mean_beta = old_a/(old_a + old_b); 
    %                 var_beta = (old_a*old_b)/((old_a+old_b)^2*(old_a+old_b+1));
    %                 SD_beta = sqrt((old_a*old_b)/((old_a+old_b)^2*(old_a+old_b+1)));
    %                 PE(i,y) = observations(i,y) - mean_beta; 
    %                 wPE(i,y) = (observations(i,y) - mean_beta)/SD_beta; 
    %             end
    %         end
    %     end
    % 
    %     % add parameters to our list_conditions_nonCorr matrix
    %     surprise = reshape(bayesian_surprise',1,[]);
    %     PredError = reshape(wPE', 1, []); 
    % 
    %     list_conditions_nonCorr(:,6) = string(normalize(surprise));
    %     list_conditions_nonCorr(:,7) = string(normalize(PredError));
    % 
    %     %prepare covariates for spm matrix
    %     surprise_ordered = [];
    %     PredError_ordered = [];
    %     for i = 1:2
    %         for y = 1:2
    %             for n = 1:3
    %                 indices = find((list_conditions_nonCorr(:,5) == prob_events(n)) & (str2double(list_conditions_nonCorr(:,2)) == i-1) & (str2double(list_conditions_nonCorr(:,1)) == move_events(y,1) | str2double(list_conditions_nonCorr(:,1)) == move_events(y,2)));
    % 
    %                 surprise_ordered = [surprise_ordered; double(list_conditions_nonCorr(indices,6))];
    %                 PredError_ordered = [PredError_ordered; double(list_conditions_nonCorr(indices,7))];
    % 
    %             end
    %         end
    %     end
    %     R = [];
    %     R(:,1) = surprise_ordered;
    %     %R(:,2) = PredError_ordered;
    %     R(rejected_trials,:) = [];
    %     names = {'Surprisal', 'PredError'};
    %     save('covs.mat', 'names', 'R')
    % end
        
    %% COVARIATE (for control analysis)
    move_events_new = [2,4;1,3]; 
    if COVARIATE == 1
        %prepare covariates for spm matrix
        RTs_ordered = [];
        for i = 1:2
            for y = 1:2
                for n = 1:3
                    trials = find((list_conditions(:,6) == prob_events(n)) & (str2double(list_conditions(:,2)) == i-1) & (str2double(list_conditions(:,1)) == move_events_new(y,1) | str2double(list_conditions(:,1)) == move_events_new(y,2)));
                    RTs_ordered = [RTs_ordered; double(list_conditions(trials,10))];
                end
            end
        end
        R = [];
        R(:,1) = RTs_ordered;
        %R(:,2) = PredError_ordered;
        names = {'RespTime'};
        save('covs.mat', 'names', 'R')
    end

    %% RUN 1ST LEVEL SPM
    order_scans = ["noTouch_Stay_LowProb", "noTouch_Stay_EqualProb", "noTouch_Stay_HighProb",...
                   "noTouch_Move_LowProb", "noTouch_Move_EqualProb", "noTouch_Move_HighProb"...
                   "Touch_Stay_LowProb", "Touch_Stay_EqualProb", "Touch_Stay_HighProb",...
                   "Touch_Move_LowProb", "Touch_Move_EqualProb", "Touch_Move_HighProb"]; 
    matlabbatch{1}.spm.stats.factorial_design.dir = {results_first_level};
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'Touch_noTouch';
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'Move_noMove';
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = 2;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).name = 'Prob';
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).levels = 3;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).dept = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).variance = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).ancova = 0;
    counter = 0; 
    for i = 1:2
        for y = 1:2
            for n = 1:3
                counter = counter+1; 
                matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(counter).levels = [i
                                                                                    y
                                                                                    n];
                dir_source = dir(strcat(subj_folder, 'spm/evEOG_thdMID', char(ID), '/', 'condition_', order_scans(counter) ,'.nii'));
                matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(counter).scans = cellstr(spm_select('expand', strcat(dir_source.folder, filesep, dir_source.name)));
            end
        end
    end
    matlabbatch{1}.spm.stats.factorial_design.des.fd.contrasts = 0;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    if COVARIATE == 0
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    elseif COVARIATE == 1
        matlabbatch{1}.spm.stats.factorial_design.multi_cov.files = {strcat('E:\Gian\GG_SensAtt_Prediction\02Data\ID', ID, '\021stLevel\covs.mat')};    
        matlabbatch{1}.spm.stats.factorial_design.multi_cov.iCFI = 1;
        matlabbatch{1}.spm.stats.factorial_design.multi_cov.iCC = 1;
    end
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;                            
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Touch_noTouch';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [-1 -1 -1 -1 -1 -1 1 1 1 1 1 1];
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Move_noMove';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 -1 -1 1 1 1 -1 -1 -1 1 1 1];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'TouchHigh_noTouchHigh';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 0 -1 0 0 -1 0 0 1 0 0 1];
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'TouchEqual_noTouchEqual';
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 -1 0 0 -1 0 0 1 0 0 1 0];
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'TouchLow_noTouchLow';
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [-1 0 0 -1 0 0 1 0 0 1 0 0];
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'Touch_noTouch_HighProb_Move';
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [0 0 0 ones(1,3)/(-3) 0 0 0 0 0 1];
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'Touch_noTouch_HighProb_noMove';
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [ones(1,3)/(-3) 0 0 0 0 0 1 0 0 0];
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'Touch_noTouch_EqualProb_Move';
    matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [0 0 0 ones(1,3)/(-3) 0 0 0 0 1 0];
    matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = 'Touch_noTouch_EqualProb_noMove';
    matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [ones(1,3)/(-3) 0 0 0 0 1 0 0 0 0];
    matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{10}.tcon.name = 'Touch_noTouch_LowProb_Move';
    matlabbatch{3}.spm.stats.con.consess{10}.tcon.weights = [0 0 0 ones(1,3)/(-3) 0 0 0 1 0 0];
    matlabbatch{3}.spm.stats.con.consess{10}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{11}.tcon.name = 'Touch_noTouch_LowProb_noMove';
    matlabbatch{3}.spm.stats.con.consess{11}.tcon.weights = [ones(1,3)/(-3) 0 0 0 1 0 0 0 0 0];
    matlabbatch{3}.spm.stats.con.consess{11}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 0;
    spm_jobman('run', matlabbatch);
    clear matlabbatch
end