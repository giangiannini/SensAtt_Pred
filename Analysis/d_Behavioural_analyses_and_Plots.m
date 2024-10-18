%% ANALYSE BEHAVOURAL DATA
% this script will go through log, EEG and response data to analyse the
% following behavioural data across participants:
% - Accuracy levels
% - Velocity across conditions

clear; clc; 
addpath('C:/Users/nnu02/Documents/MATLAB/fieldtrip-20220827');

subjects = ["04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "24" "25" "26" "27" "28" "29"];

%set up paths
folder = 'E:/Gian/GG_SensAtt_Prediction/02Data';
touch_events = [125 124 252];
touch_labels = ["noTouch" "Touch"];

move_events = [2, 4; 1, 3];
move_labels = ["Stay", "Move"];

prob_events = ["25/75", "50/50", "75/25"];
prob_labels = ["LowProb", "EqualProb", "HighProb"];

subj = [];
frequency_Unity = []; 

%% SET FOLDER TO EXPORT PLOTS 
img_folder = 'E:\Gian\GG_SensAtt_Prediction\02Data\behavioural_group_data';

%% IMPORT BEHAVIOURAL DATA
if ~exist(strcat(folder, '/behavioural_group_data/behavioural_data_all_removedtrials.mat'))
    for ID = subjects
    
        ID = char(ID);
    
        subj_folder = strcat(folder, '/ID', ID, '/01EEG/');
        bdf_file = strcat(folder, '/ID', ID, '/01EEG/*ID', ID, '*.bdf');
        bdf_file = dir(bdf_file); 
        bdf_file = strcat(bdf_file.folder, '\', bdf_file.name); 
    
        %% Load events file
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
    
        %% Load rejected trials
        %load(strcat(subj_folder, 'preprocessed_final_SPM_vEOG_hEOGlong_epoch_automated.mat'));
        load(strcat(subj_folder, 'rejected_trials_manual.mat'));
        load(strcat(subj_folder, 'rejected_trials_RT'));
        %rejected_trials = []; 
    
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
        
        events_list_EEG(events_list_EEG == 126) = [];
        events_list_EEG(events_list_EEG~=10 & events_list_EEG~=20 & events_list_EEG~=30 & events_list_EEG~=124 & events_list_EEG~=125) = []; 
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
        % zeros_ones_EEG = abs(data.trialinfo - 125);
        % zeros_ones_Log = list_conditions_nonCorr(:,2); 
        % zeros_ones_Log(unique([rejected_trials rejected_trials_RT])) = []; 
        % zerosss = [zeros_ones_EEG, zeros_ones_Log];
        % if sum(zeros_ones_EEG ~= str2double(zeros_ones_Log)) > 0
        %     error('PORCAMADONNA IL LOG E` SBAGLIATO')
        % else
        %     display('allesgut!')
        % end
    
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
    
        %% Save the list conditions matrix which contains accuracy levels
        accuracy = list_conditions(:,7);
        subj_response = list_conditions(:,6);
        touch_conditions = list_conditions(:, 2);
        touch_conditions(touch_conditions == "0") = "noTouch";
        touch_conditions(touch_conditions == "1") = "Touch"; 
        movement_conditions = list_conditions(:,1); 
        movement_conditions(movement_conditions == "1" | movement_conditions == "3") = "Move";
        movement_conditions(movement_conditions == "2" | movement_conditions == "4") = "Stay";
        prob_conditions = list_conditions(:,5);
    
        subj{find(subjects == string(ID))}.accuracy_table = table(repmat(ID, length(accuracy),1),touch_conditions,movement_conditions,prob_conditions,subj_response,accuracy);
    
        clear("touch_conditions", "movement_conditions", "prob_conditions")
    
        %% Also create simplified accuracy table
        subj{find(subjects == string(ID))}.accuracy_table_simplified = table(repmat(ID, length(Correctness),1),ExperimentCondition,Experimentresponse,Correctness);
    
        %% Now check RTs from Log file
        touches = find(Log1(:,6) == "Index_visual");
        RTs = []; 
        touch_conditions = Log1(touches, 5);
        touch_conditions(touch_conditions == "0") = "noTouch";
        touch_conditions(touch_conditions == "1") = "Touch"; 
        movement_conditions = Log1(touches,4); 
        movement_conditions(movement_conditions == "1" | movement_conditions == "3") = "Move";
        movement_conditions(movement_conditions == "2" | movement_conditions == "4") = "Stay";
        prob_conditions = Log1(touches,3);
        missing = [];
        for i = 1:length(touches)
            if Log1(touches(i)-1,6) == "GO"
                RTs = [RTs; double(Log1(touches(i),7)) - double(Log1(touches(i)-1,7))];
            else
                missing = [i];
            end
        end
        touch_conditions(missing) = []; 
        movement_conditions(missing) = []; 
        prob_conditions(missing) = [];  

        RTs_table = table(repmat(ID, length(touch_conditions),1),touch_conditions,movement_conditions,prob_conditions,RTs);
        if max(RTs_table.RTs) ~= max(RTs_table.RTs(rejected_trials_RT,:))
            error('Trial displacement!')
        else
            RTs_table(unique([rejected_trials rejected_trials_RT]),:) = []; 
        end    
        subj{find(subjects == string(ID))}.RTs_table = RTs_table;
    
        %% NOW CHECK VELOCITY FROM FILE
        %extract start and end of each condition
        touches = find(Log1(:,6) == "Index_visual");
        touches(unique([rejected_trials rejected_trials_RT])) = []; 
        start_end = [];
        for i = 1:length(touches)
            start_end = [start_end; double(Log1(touches(i)-1,7)), double(Log1(touches(i),7))];
        end
        velocity_table = [subj{find(subjects == string(ID))}.accuracy_table(:,1:4), table(start_end(:,1),'VariableNames',{'Start'}), table(start_end(:,2),'VariableNames',{'End'})];
        
        %FIRST FILL OUT THE MOVEMENT CONDITIONS
        %import hand movements
        opts = delimitedTextImportOptions("NumVariables", 7);
        opts.DataLines = [4, Inf];
        opts.Delimiter = "\t";
        opts.VariableNames = ["Timems", "Index_x", "Index_y", "Index_z", "Index_rotx", "Index_roty", "Index_rotz"];
        opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double"];
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        ExperimentHandPositions = readtable(strcat("E:\Gian\GG_SensAtt_Prediction\02Data\ID", ID, "\00Behavioural\Experiment_Hand_Positions_ID", ID, ".txt"), opts);
        clear opts

        frequency_Unity = [frequency_Unity mean(diff(ExperimentHandPositions.Timems))];
        
        movements = find(velocity_table.movement_conditions == "Move");
        for i = 1:length(movements)
            [value start_pos] = min(abs(ExperimentHandPositions.Timems - velocity_table{movements(i),5}));
            [value end_pos] = min(abs(ExperimentHandPositions.Timems - velocity_table{movements(i),6}));
            velocities = [];
            for j = 1:length(start_pos:end_pos)-1
                delta_space = sqrt((ExperimentHandPositions{start_pos+(j)-1,2}-ExperimentHandPositions{start_pos+(j),2})^2 + (ExperimentHandPositions{start_pos+(j)-1,4}-ExperimentHandPositions{start_pos+(j),4})^2);
                delta_time = (ExperimentHandPositions{start_pos+(j),1}-ExperimentHandPositions{start_pos+(j)-1,1})/1000; %in secs
                velocities = [velocities delta_space/delta_time];
            end
            velocities(velocities < 0.1) = []; 
            velocity_table(movements(i),7) = {mean(velocities)};
        end
    
        %THEN FILL OUT THE STAY CONDITIONS
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
    
        stays = find(velocity_table.movement_conditions == "Stay");
    
        counter = 0;
        for i = 1:length(stays)
            pos_events = find((TrialTable.Stimulation == (find(touch_labels == velocity_table{stays(i), 2})-1)) & ...
                (TrialTable.Trial_type == 2 | TrialTable.Trial_type == 4) & ...
                (TrialTable.Block_type == velocity_table{stays(i),4}));
            [value pos] = min(abs(pos_events - counter));
            TrialTable_corresp_trial = pos_events(pos); 
            counter = TrialTable_corresp_trial;
            velocity_table{stays(i),7} = TrialTable{TrialTable_corresp_trial,9};
        end
    
        subj{find(subjects == string(ID))}.velocity_table = velocity_table;         

        %% NOW DO MORE ACCURATE RT TABLE: include only RTs from proper movement onset
        %extract start and end of each condition
        RTs_table_accurate = subj{find(subjects == string(ID))}.RTs_table; 
        
        movements = find(RTs_table_accurate.movement_conditions == "Move");
        for i = 1:length(movements)
            [value start_pos] = min(abs(ExperimentHandPositions.Timems - velocity_table{movements(i),5}));
            [value end_pos] = min(abs(ExperimentHandPositions.Timems - velocity_table{movements(i),6}));
            velocities = [];
            for j = 1:length(start_pos:end_pos)-1
                delta_space = sqrt((ExperimentHandPositions{start_pos+(j)-1,2}-ExperimentHandPositions{start_pos+(j),2})^2 + (ExperimentHandPositions{start_pos+(j)-1,4}-ExperimentHandPositions{start_pos+(j),4})^2);
                delta_time = (ExperimentHandPositions{start_pos+(j),1}-ExperimentHandPositions{start_pos+(j)-1,1})/1000; %in secs
                velocities = [velocities delta_space/delta_time];
            end
            %plot also one trial
            % figure;
            % plot(ExperimentHandPositions{start_pos:end_pos-1,1},velocities)
            % xlabel('time during trial (ms)')
            % ylabel('velocity')
    
            to_delete = (velocities < 0.1); 
            RTs_real = ExperimentHandPositions{start_pos:end_pos,1};
            RTs_real(to_delete) = [];
            RTs_table_accurate(movements(i),5) = {max(RTs_real) - min(RTs_real)};
        end
    
        subj{find(subjects == string(ID))}.RTs_table_accurate = RTs_table_accurate; 
    
end
    save(strcat(folder, "/behavioural_group_data/behavioural_data_all_removedtrials.mat"), "subj")
else
    load(strcat(folder, '/behavioural_group_data/behavioural_data_all_removedtrials.mat'))
end

%% RUN STATS ON ACCURACY LEVELS
%put together the big table first
full_table_accuracy = []; 
for i = 1:length(subjects)
    full_table_accuracy = [full_table_accuracy; subj{i}.accuracy_table_simplified];
end
%plot data
full_table_accuracy.ExperimentCondition = categorical(full_table_accuracy.ExperimentCondition); 
full_table_accuracy.Experimentresponse = categorical(full_table_accuracy.Experimentresponse); 
full_table_accuracy.Var1 = cellstr(full_table_accuracy.Var1);
full_table_accuracy.Var1 = categorical(full_table_accuracy.Var1); 
full_table_accuracy.Correctness = double(full_table_accuracy.Correctness); 

full_table_accuracy = full_table_accuracy(full_table_accuracy.Var1 ~= '23',:);

avg_accuracies = groupsummary(full_table_accuracy, ["ExperimentCondition"], @(x)mean(x)*100,["Correctness"]);
error_accuracies = groupsummary(full_table_accuracy, ["ExperimentCondition"], @(x)std(x)/sqrt(8)*100,["Correctness"]);
avg_accuracies_subjects = groupsummary(full_table_accuracy, ["Var1","ExperimentCondition"], @(x)mean(x)*100,["Correctness"]);

figure;
bar(1, avg_accuracies.fun1_Correctness(1), 'FaceAlpha', 0.2, 'FaceColor', [1 0 1])
hold on
bar(2, avg_accuracies.fun1_Correctness(2), 'FaceAlpha', 0.5, 'FaceColor', [1 0 1])
bar(3, avg_accuracies.fun1_Correctness(3), 'FaceAlpha', 1, 'FaceColor', [1 0 1])
errorbar([1 2 3], avg_accuracies.fun1_Correctness, error_accuracies.fun1_Correctness,'k','LineWidth', 1.5,'linestyle','none','HandleVisibility','off'); 
for i = 1:25
    for j = 1:3
        r = 0.2 + (0.6-0.25) .* rand(1,1);
        r = j + r - 0.4;
        plot(r, avg_accuracies_subjects.fun1_Correctness(((i-1)*3)+j), 'ko');
    end
end
set(gca, 'xticklabel', cellstr(avg_accuracies.ExperimentCondition), 'ylim', [0 100], 'YTick', [0:10:100])
hold off
ylim([0 110])
set(gca,'XAxisLocation','top', 'box','off', 'XTick', [])
exportgraphics(gcf, strcat(img_folder, '/barPlot_accuracy.emf'))
exportgraphics(gcf, strcat(img_folder, '/barPlot_accuracy.png'))
close all


%% RUN STATS ON RTs
full_table_velocity = []; 
full_table_RTs_accurate = []; 
full_table_RTs = [];
for i = 1:length(subjects)
    full_table_velocity = [full_table_velocity; subj{i}.velocity_table];
    full_table_RTs_accurate = [full_table_RTs_accurate; subj{i}.RTs_table_accurate];
    full_table_RTs = [full_table_RTs; subj{i}.RTs_table];
end

full_table_velocity = full_table_velocity(categorical(cellstr(full_table_velocity.Var1)) ~= "23",:);
full_table_RTs_accurate = full_table_RTs_accurate(categorical(cellstr(full_table_RTs_accurate.Var1)) ~= '23',:);
full_table_RTs = full_table_RTs(categorical(cellstr(full_table_RTs.Var1)) ~= '23',:);

% full_table_velocity.touch_conditions = categorical(full_table_velocity.touch_conditions); 
% full_table_velocity.movement_conditions = categorical(full_table_velocity.movement_conditions); 
% full_table_velocity.movement_conditions = reordercats(full_table_velocity.movement_conditions, {'Stay', 'Move'}); 
% full_table_velocity.prob_conditions = categorical(full_table_velocity.prob_conditions); 
% full_table_velocity.Var1 = cellstr(full_table_velocity.Var1);
% full_table_velocity.Var1 = categorical(full_table_velocity.Var1);
% 

full_table_RTs_accurate.touch_conditions = categorical(full_table_RTs_accurate.touch_conditions); 
full_table_RTs_accurate.touch_conditions = reordercats(full_table_RTs_accurate.touch_conditions, {'noTouch', 'Touch'}); 
full_table_RTs_accurate.movement_conditions = categorical(full_table_RTs_accurate.movement_conditions); 
full_table_RTs_accurate.movement_conditions = reordercats(full_table_RTs_accurate.movement_conditions, {'Stay', 'Move'}); 
full_table_RTs_accurate.Var1 = cellstr(full_table_RTs_accurate.Var1);
full_table_RTs_accurate.Var1 = categorical(full_table_RTs_accurate.Var1); 

%
% do stats on RespT
%
model_RT = fitlme(full_table_RTs_accurate, 'RTs ~ movement_conditions*touch_conditions*prob_conditions + (1|Var1)');
anova(model_RT)

%% PLOT BIG GRAPH
avg_RespT = groupsummary(full_table_RTs_accurate, ["movement_conditions", "touch_conditions", "prob_conditions"], @(x)mean(x),["RTs"]);
avg_RespT_subj = groupsummary(full_table_RTs_accurate, ["Var1", "movement_conditions", "touch_conditions", "prob_conditions"], @(x)mean(x),["RTs"]);
error_RespT = groupsummary(avg_RespT_subj, ["movement_conditions", "touch_conditions", "prob_conditions"], @(x)std(x)/sqrt(length(x)),["fun1_RTs"]);

%barplot
colors_to_plot = [0.7, 0.7, 1; ... %very light blue
                0.4, 0.4, 1; ... %light blue
                0, 0, 1; ... %blue
                1, 0.7, 0.7;... %very light red
                1, 0.4, 0.4; ... %light red
                1, 0, 0]; %red
figure; 
hold on
    %means_barplot = mean(beta.(order_scans(i)).avg(index, times{cluster_num}));
aaa = bar(1, avg_RespT.fun1_RTs(1), 'FaceColor', colors_to_plot(1,:), 'FaceAlpha', 0.3);
aaa = bar(2, avg_RespT.fun1_RTs(2), 'FaceColor', colors_to_plot(2,:), 'FaceAlpha', 0.3);
aaa = bar(3, avg_RespT.fun1_RTs(3), 'FaceColor', colors_to_plot(3,:), 'FaceAlpha', 0.3);
aaa = bar(4, avg_RespT.fun1_RTs(4), 'FaceColor', colors_to_plot(4,:), 'FaceAlpha', 0.3);
aaa = bar(5, avg_RespT.fun1_RTs(5), 'FaceColor', colors_to_plot(5,:), 'FaceAlpha', 0.3);
aaa = bar(6, avg_RespT.fun1_RTs(6), 'FaceColor', colors_to_plot(6,:), 'FaceAlpha', 0.3);
aaa = bar(7, avg_RespT.fun1_RTs(7), 'FaceColor', colors_to_plot(1,:));
%hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(1,:))
aaa = bar(8, avg_RespT.fun1_RTs(8), 'FaceColor', colors_to_plot(2,:));
%hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(2,:))
aaa = bar(9, avg_RespT.fun1_RTs(9), 'FaceColor', colors_to_plot(3,:));
%hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(3,:))
aaa = bar(10, avg_RespT.fun1_RTs(10), 'FaceColor', colors_to_plot(4,:));
%hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(4,:))
aaa = bar(11, avg_RespT.fun1_RTs(11), 'FaceColor', colors_to_plot(5,:));
%hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(5,:))
aaa = bar(12, avg_RespT.fun1_RTs(12), 'FaceColor', colors_to_plot(6,:));
%hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(6,:))
title('RespTime across conditions')
errorbar([1], avg_RespT.fun1_RTs(1), error_RespT.fun1_fun1_RTs(1),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
errorbar([2], avg_RespT.fun1_RTs(2), error_RespT.fun1_fun1_RTs(2),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
errorbar([3], avg_RespT.fun1_RTs(3), error_RespT.fun1_fun1_RTs(3),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
errorbar([4], avg_RespT.fun1_RTs(4), error_RespT.fun1_fun1_RTs(4),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
errorbar([5], avg_RespT.fun1_RTs(5), error_RespT.fun1_fun1_RTs(5),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
errorbar([6], avg_RespT.fun1_RTs(6), error_RespT.fun1_fun1_RTs(6),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
errorbar([7], avg_RespT.fun1_RTs(7), error_RespT.fun1_fun1_RTs(7),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
errorbar([8], avg_RespT.fun1_RTs(8), error_RespT.fun1_fun1_RTs(8),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
errorbar([9], avg_RespT.fun1_RTs(9), error_RespT.fun1_fun1_RTs(9),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
errorbar([10], avg_RespT.fun1_RTs(10), error_RespT.fun1_fun1_RTs(10),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
errorbar([11], avg_RespT.fun1_RTs(11), error_RespT.fun1_fun1_RTs(11),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
errorbar([12], avg_RespT.fun1_RTs(12), error_RespT.fun1_fun1_RTs(12),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
set(gca,'xticklabel',{[]})
set(gca,'xtick',[])
ylabel('Response Times (ms)')
ylim([0 900])
xlim([0 18])
exportgraphics(gcf, strcat(img_folder, '/barPlot_RTs.emf'))
exportgraphics(gcf, strcat(img_folder, '/barPlot_RTs.png'))
close all