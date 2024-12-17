clc;clear;
addpath(genpath('C:/Users/nnu02/Documents/MATLAB/spm12'));
addpath('C:/Users/nnu02/Documents/MATLAB/fieldtrip-20220827');
addpath('E:/02Data/03Utils/Functions/');
ft_defaults

subjects = ["04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "24" "25" "26" "27" "28" "29"];

folder = 'E:/Gian/GG_SensAtt_Prediction/02Data';
touch_events = [125 124 252];
%touch_events = [20];
touch_labels = ["noTouch" "Touch"];

move_events = [1, 3; 2, 4];
move_labels = ["Move", "Stay"];

prob_events = ["25/75", "50/50", "75/25"];
prob_labels = ["LowProb", "EqualProb", "HighProb"];

diffs = []; 
diffs_total = []; 

%eog_modality = "SPM_vEOG_hEOGHPF05"; 
caplocation = 'E:/02Data/03Utils/biosemi64.lay';

colorscale = [45, 109, 224; 24, 115, 222; 0, 120, 221; 0, 124, 219; 0, 128, 218; 0, 132, 216; 0, 136, 214; 0, 140, 212; 0, 144, 209; 0, 148, 206; 0, 152, 203; 0, 156, 198; 0, 160, 194; 0, 164, 189; 0, 168, 184; 0, 171, 178; 0, 175, 173; 0, 180, 167; 0, 184, 160; 0, 188, 152; 0, 191, 142; 0, 193, 132; 44, 196, 121; 70, 197, 109; 92, 198, 96; 114, 198, 82; 134, 197, 67; 153, 196, 51; 172, 193, 32; 191, 191, 0]/255;


%% IMPORT DATA FROM SUBJECT
for sogg = 1:length(subjects)
    ID = subjects(sogg);
    
    verbose = 0; %specify whether you want each preproc step printed in a ps file (yes = 1; no = 0);
    notch = 0;
    
    %Stuff for SPM
    subj_folder = strcat(folder, '/ID', ID, '/01EEG/');

    %copy paste the loc file
    loc_dir = strcat(folder, '/ID', ID, '/00Behavioural/');
    bdf_file = strcat(folder, '/ID', ID, '/01EEG/*ID', ID, '*.bdf');
    bdf_file = dir(bdf_file); 
    bdf_file = strcat(bdf_file.folder, '\', bdf_file.name); 

    %Stuff for FT
    caplocation = 'E:/02Data/03Utils/biosemi64.lay';
    neighbourslocation = 'E:/02Data/03Utils/biosemi64_neighb.mat';
    
    cd(subj_folder) %jump in the right folder
    
    %% IMPORT RAW DATASET
    load('preprocessed_final_SPM_vEOG_hEOGApril2024_fully_manual.mat')
    load('rejected_trials_manual.mat')
    load('rejected_trials_RT.mat')
    load('list_conditions.mat')

    %% Filtering!
    % cfg = []; 
    % cfg.hpfilter = 'yes'; 
    % cfg.hpfreq = 0.5; 
    % data = ft_preprocessing(cfg, data); 

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

    %% Import Log
    
    %% EXTRACT EVENTS
    ERP = []; 
    for i = 1:2 %touch and notouch
        for y = 1:2 %move and no move
            for n = 1:3 %prob
                name = strcat(touch_labels(i), '_',  move_labels(y), '_', prob_labels(n)); 
                cfg = []; 
                cfg.trials = find((list_conditions(:,6) == prob_events(n)) & (str2double(list_conditions(:,2)) == i-1) & (str2double(list_conditions(:,1)) == move_events(y,1) | str2double(list_conditions(:,1)) == move_events(y,2)));
                trials_number = cfg.trials; 
                trials.(name){sogg} = ft_selectdata(cfg, data); 

                list_conditions(trials_number,9) = name; 

                % % % % % % % % % %HPF 
                % % % % % % % % % cfg = []; 
                % % % % % % % % % cfg.hpfilter = 'yes'; 
                % % % % % % % % % cfg.hpfreq = 1;
                % % % % % % % % % trials.(name){sogg} = ft_preprocessing(cfg, trials.(name){sogg});
                
                cfg = []; 
                cfg.demean = 'yes'; 
                cfg.baselinewindow = [-0.01 -0.005];
                trials.(name){sogg} = ft_preprocessing(cfg, trials.(name){sogg});
                
                %ERP
                cfg = []; 
                ERPs.(name){sogg} = ft_timelockanalysis(cfg, trials.(name){sogg}); 
                ERP.(name) = ERPs.(name){sogg};
            end
        end
    end
    save("ERP.mat", "ERP")
end

save("E:\Gian\GG_SensAtt_Prediction\02Data\EEG_group_data_old\ERPs.mat", "ERPs")

%% SET FOLDER TO EXPORT PLOTS 
img_folder = 'E:\Gian\GG_SensAtt_Prediction\02Data\EEG_group_data';

%% PLOT RTs HIST
tot_delays = []; 
for i = 1:length(subjects)
    ID = char(subjects(i)); 
    subj_folder = strcat(folder, '/ID', ID, '/01EEG/');
    cd(subj_folder); 
    load('list_conditions.mat')
    load('movement_delays.mat')
    tot_delays = [tot_delays; diff];
end

tot_delays = 0 - tot_delays./2048; 

figure('Renderer', 'painters', 'Position', [10 10 900 600])
histogram(tot_delays, 'NumBins', 100, 'FaceColor', 'k')
xlim([-2 1])
set(gca,'XTick',[])
set(gca,'YTick',[])
exportgraphics(gcf, strcat(img_folder, '/RTs_hist.emf'))
exportgraphics(gcf, strcat(img_folder, '/RTs_hist.png'))
close all


%% PLOT FULL IMAGE
cfg = []; 
noTouch_Stay_HighProb = ft_timelockgrandaverage(cfg, ERPs.noTouch_Stay_HighProb{:});
noTouch_Stay_LowProb = ft_timelockgrandaverage(cfg, ERPs.noTouch_Stay_LowProb{:});
noTouch_Stay_EqualProb = ft_timelockgrandaverage(cfg, ERPs.noTouch_Stay_EqualProb{:});
noTouch_Move_HighProb = ft_timelockgrandaverage(cfg, ERPs.noTouch_Move_HighProb{:});
noTouch_Move_LowProb = ft_timelockgrandaverage(cfg, ERPs.noTouch_Move_LowProb{:});
noTouch_Move_EqualProb = ft_timelockgrandaverage(cfg, ERPs.noTouch_Move_EqualProb{:});
Touch_Stay_HighProb = ft_timelockgrandaverage(cfg, ERPs.Touch_Stay_HighProb{:});
Touch_Stay_LowProb = ft_timelockgrandaverage(cfg, ERPs.Touch_Stay_LowProb{:});
Touch_Stay_EqualProb = ft_timelockgrandaverage(cfg, ERPs.Touch_Stay_EqualProb{:});
Touch_Move_HighProb = ft_timelockgrandaverage(cfg, ERPs.Touch_Move_HighProb{:});
Touch_Move_LowProb = ft_timelockgrandaverage(cfg, ERPs.Touch_Move_LowProb{:});
Touch_Move_EqualProb = ft_timelockgrandaverage(cfg, ERPs.Touch_Move_EqualProb{:});

%calculate error term
noTouch_Stay_HighProb.err = sqrt(noTouch_Stay_HighProb.var) ./ sqrt(noTouch_Stay_HighProb.dof); 
noTouch_Stay_LowProb.err = sqrt(noTouch_Stay_LowProb.var) ./ sqrt(noTouch_Stay_LowProb.dof); 
noTouch_Stay_EqualProb.err = sqrt(noTouch_Stay_EqualProb.var) ./ sqrt(noTouch_Stay_EqualProb.dof); 
noTouch_Move_HighProb.err = sqrt(noTouch_Move_HighProb.var) ./ sqrt(noTouch_Move_HighProb.dof); 
noTouch_Move_LowProb.err = sqrt(noTouch_Move_LowProb.var) ./ sqrt(noTouch_Move_LowProb.dof); 
noTouch_Move_EqualProb.err = sqrt(noTouch_Move_EqualProb.var) ./ sqrt(noTouch_Move_EqualProb.dof); 
Touch_Stay_HighProb.err = sqrt(Touch_Stay_HighProb.var) ./ sqrt(Touch_Stay_HighProb.dof); 
Touch_Stay_LowProb.err = sqrt(Touch_Stay_LowProb.var) ./ sqrt(Touch_Stay_LowProb.dof); 
Touch_Stay_EqualProb.err = sqrt(Touch_Stay_EqualProb.var) ./ sqrt(Touch_Stay_EqualProb.dof); 
Touch_Move_HighProb.err = sqrt(Touch_Move_HighProb.var) ./ sqrt(Touch_Move_HighProb.dof); 
Touch_Move_LowProb.err = sqrt(Touch_Move_LowProb.var) ./ sqrt(Touch_Move_LowProb.dof); 
Touch_Move_EqualProb.err = sqrt(Touch_Move_EqualProb.var) ./ sqrt(Touch_Move_EqualProb.dof); 

%% PLOT ALL DATA
cfg = [];
cfg.layout = caplocation;
cfg.colors = [245, 161, 95; ... %blue
                251, 193, 110; ... %light blue
                253, 224, 132; ... %very light blue
                186, 56, 75; ... %red
                226, 97, 82; ... %light red
                237, 130, 86;... %very light red
                138, 196, 131; ... %blue
                175, 216, 139; ... %light blue
                214, 235, 150; ... %very light blue
                18, 134, 108; ... %red
                65, 155, 116; ... %light red
                102, 176, 123]/255; %very light red
cfg.linestyle = {'-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-'};
%cfg.channel = [9 10 11 12 19 32 33 34 48 49 56];
cfg.xlim = [-2 1];
cfg.ylim = [-2 7];
%cfg.ylim = [-3.66 4.33];
%cfg.channel = [11 12 19 32 46 47 48 49 56];
cfg.channel = ["C1", "CP1", "CPz", "Cz"];
cfg.channel = ["FC3", "FC1", "C1", "C3", "CP3", "CP1", "P1", "P3", "Pz", ...
               "CPz", "FC2", "FCz", "Cz", "C2", "CP2", "P2"];
%cfg.channel = [14 15];
cfg.printlayers = 1;
cfg.output_printlayers = strcat(img_folder, '/big_epoch_groupERPs');
GIAN_plot_data(cfg, noTouch_Stay_HighProb, noTouch_Stay_EqualProb, noTouch_Stay_LowProb, Touch_Stay_HighProb, Touch_Stay_EqualProb, Touch_Stay_LowProb, ...
                noTouch_Move_HighProb, noTouch_Move_EqualProb, noTouch_Move_LowProb, Touch_Move_HighProb, Touch_Move_EqualProb, Touch_Move_LowProb);
exportgraphics(gcf, strcat(img_folder, '/big_epoch_groupERPs.emf'))
exportgraphics(gcf, strcat(img_folder, '/big_epoch_groupERPs.png'))
close all

cfg.output_printlayers = strcat(img_folder, '/small_epoch_groupERPs');
cfg.xlim = [-0.2 0.5];
cfg.ylim = [-1 7];
GIAN_plot_data(cfg, noTouch_Stay_HighProb, noTouch_Stay_EqualProb, noTouch_Stay_LowProb, Touch_Stay_HighProb, Touch_Stay_EqualProb, Touch_Stay_LowProb, ...
                noTouch_Move_HighProb, noTouch_Move_EqualProb, noTouch_Move_LowProb, Touch_Move_HighProb, Touch_Move_EqualProb, Touch_Move_LowProb);
exportgraphics(gcf, strcat(img_folder, '/small_epoch_groupERPs.emf'))
exportgraphics(gcf, strcat(img_folder, '/small_epoch_groupERPs.png'))
close all


%% plot layout
cfg = []; 
cfg.layout = caplocation; 
layout = ft_prepare_layout(cfg); 

cfg.channel = ["FC3", "FC1", "C1", "C3", "CP3", "CP1", "P1", "P3", "Pz", ...
               "CPz", "FC2", "FCz", "Cz", "C2", "CP2", "P2"];
index = arrayfun(@(k) find(strncmp(cfg.channel(k),data.label,3)), 1:length(cfg.channel));
figure('Renderer', 'painters', 'Position', [10 10 500 500])
plot(layout.pos(1:64,1), layout.pos(1:64,2), 'ko');
hold on
plot(layout.outline{1}(:,1), layout.outline{1}(:,2), 'k-');
plot(layout.pos(index,1), layout.pos(index,2), 'ko', 'MarkerFaceColor', [0 0 0], 'LineWidth', 5);
xlim([-0.6 0.6])
ylim([-0.6 0.6])
set(gcf, 'color', 'none')
exportgraphics(gcf, strcat(img_folder, '/layout.emf'))
exportgraphics(gcf, strcat(img_folder, '/layout.png'))
close all

%% PLOT EFF OF TOUCH
for i = 1:length(subjects)
    cfg = []; 
    cfg.operation = '(x1 + x2 + x3 + x4 + x5 + x6)/6';
    cfg.parameter = {'avg'}; 
    con_noTouch{i} = ft_math(cfg, ERPs.noTouch_Move_HighProb{i}, ERPs.noTouch_Move_EqualProb{i}, ERPs.noTouch_Move_LowProb{i}, ...
                             ERPs.noTouch_Stay_HighProb{i}, ERPs.noTouch_Stay_EqualProb{i}, ERPs.noTouch_Stay_LowProb{i});
    con_Touch{i} = ft_math(cfg, ERPs.Touch_Move_HighProb{i}, ERPs.Touch_Move_EqualProb{i}, ERPs.Touch_Move_LowProb{i}, ...
                           ERPs.Touch_Stay_HighProb{i}, ERPs.Touch_Stay_EqualProb{i}, ERPs.Touch_Stay_LowProb{i});

    cfg = []; 
    cfg.operation = 'x1 - x2'; 
    cfg.parameter = 'avg'; 
    con_Touch{i} = ft_math(cfg, con_Touch{i}, con_noTouch{i});
end

cfg = []; 
Touch = ft_timelockgrandaverage(cfg, con_Touch{:});
Touch.err = sqrt(Touch.var)./sqrt(Touch.dof); 

%plot butterfly plot manually
figure('Renderer', 'painters', 'Position', [10 10 900 600])
for i = 1:64
    plot(Touch.time, Touch.avg(i,:), "Color", [101 96 92]/255, 'LineStyle','-')
    % plot(Touch.time, Touch.avg(i,:), "Color", [1 0 0], 'LineStyle','-')
    % plot(noTouch.time, noTouch.avg(i,:), "Color", [0 0 1], 'LineStyle','-')
    hold on
end
xlim([-0.2 0.5])

xline(0.05)
xline(0.11)
%xline(0.20)
xline(0.25)
%xline(0.35)
exportgraphics(gcf, strcat(img_folder, '/butterfly.emf'))
exportgraphics(gcf, strcat(img_folder, '/butterfly.png'))
close all

%P50
cfg = []; 
cfg.layout = caplocation; 
cfg.xlim = [0.049 0.051];
cfg.comment = 'no';
cfg.zlim = [-2 2]; 
    cfg.marker = 'off';

cfg.colormap = colorscale; 
ft_topoplotER(cfg, Touch)
exportgraphics(gcf, strcat(img_folder, '/topo_50.emf'))
exportgraphics(gcf, strcat(img_folder, '/topo_50.png'))
close all


%P100
cfg = []; 
cfg.layout = caplocation; 
cfg.xlim = [0.110 0.112];
cfg.comment = 'no';
cfg.zlim = [-3 3];
cfg.colormap = colorscale; 
    cfg.marker = 'off';

ft_topoplotER(cfg, Touch)
exportgraphics(gcf, strcat(img_folder, '/topo_100.emf'))
exportgraphics(gcf, strcat(img_folder, '/topo_100.png'))
close all

% %P3b
% cfg = []; 
% cfg.layout = caplocation; 
% cfg.xlim = [0.249 0.251];
% cfg.zlim = [-6 6];
% cfg.comment = 'no';
% ft_topoplotER(cfg, Touch)
% exportgraphics(gcf, strcat(img_folder, '/topo_300.emf'))
% exportgraphics(gcf, strcat(img_folder, '/topo_300.png'))
% close all
