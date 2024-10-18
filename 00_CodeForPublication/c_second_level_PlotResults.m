%% This script will plot the EEG group-level results obtained from the "c_second_level_12regressors" script. 
clear;clc;
addpath(genpath('C:/Users/nnu02/Documents/MATLAB/spm12'));
addpath('C:/Users/nnu02/Documents/MATLAB/fieldtrip-20220827');
ft_defaults

subjects = ["04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "24" "25" "26" "27" "28" "29"];

folder = 'E:/Gian/GG_SensAtt_Prediction/02Data';
results_first_level = '00Results1stLevel';
results_second_level = strcat(folder, '\2nd_level');

%% PLOT IMGS!
beta = [];
for ID = subjects
    ID = char(ID);
    
    index_for_storing = find(strcmp(ID, subjects));
    
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

    first_level_folder = strcat(folder, '/ID', ID, '/021stLevel/');

    %% LOAD DATASET AND REJ TRIALS
    if strcmp(ID, '04')
        load(strcat(subj_folder, 'preprocessed_final_SPM_vEOG_hEOGApril2024_fully_manual.mat'))

        %create structure to store data from imgs
        cfg = []; 
        cfg.channel = 1:64;
        cfg.latency = [-0.05 0.5];
        data = ft_selectdata(cfg, data);
        cfg = []; 
        ERP = ft_timelockanalysis(cfg, data); 
    end

    %% IMPORT BETA IMAGES FOR PLOTTING
    for i = 1:12
        order_scans = ["noTouch_Stay_LowProb", "noTouch_Stay_EqualProb", "noTouch_Stay_HighProb",...
               "noTouch_Move_LowProb", "noTouch_Move_EqualProb", "noTouch_Move_HighProb"...
               "Touch_Stay_LowProb", "Touch_Stay_EqualProb", "Touch_Stay_HighProb",...
               "Touch_Move_LowProb", "Touch_Move_EqualProb", "Touch_Move_HighProb"];     
        name = order_scans(i); 
        beta.(name){index_for_storing} = ERP;
        if i < 10
            beta.(name){index_for_storing}.avg = gian_extract_SPM_images(strcat(first_level_folder, results_first_level, '/beta_000', num2str(i), '.nii'));
        elseif i >= 10
            beta.(name){index_for_storing}.avg = gian_extract_SPM_images(strcat(first_level_folder, results_first_level, '/beta_00', num2str(i), '.nii'));
        end
        beta.(name){index_for_storing};
    end

    %% IMPORT CON IMAGES FOR PLOTTING
    for i = 7:12
        order_scans = ["noTouch_Stay_LowProb", "noTouch_Stay_EqualProb", "noTouch_Stay_HighProb",...
               "noTouch_Move_LowProb", "noTouch_Move_EqualProb", "noTouch_Move_HighProb"...
               "Touch_Move_HighProb", "Touch_Stay_HighProb", "Touch_Move_EqualProb",...
               "Touch_Stay_EqualProb", "Touch_Move_LowProb", "Touch_Stay_LowProb"];     
        name = order_scans(i); 
        con.(name){index_for_storing} = ERP;
        if i < 11
            con.(name){index_for_storing}.avg = gian_extract_SPM_images(strcat(first_level_folder, results_first_level, '/con_000', num2str(i-1), '.nii'));
        elseif i >= 11
            con.(name){index_for_storing}.avg = gian_extract_SPM_images(strcat(first_level_folder, results_first_level, '/con_00', num2str(i-1), '.nii'));
        end
        con.(name){index_for_storing};
    end
end

%% beta imgs
cfg = []; 
noTouch_Stay_HighProb = ft_timelockgrandaverage(cfg, beta.noTouch_Stay_HighProb{:});
noTouch_Stay_LowProb = ft_timelockgrandaverage(cfg, beta.noTouch_Stay_LowProb{:});
noTouch_Stay_EqualProb = ft_timelockgrandaverage(cfg, beta.noTouch_Stay_EqualProb{:});
noTouch_Move_HighProb = ft_timelockgrandaverage(cfg, beta.noTouch_Move_HighProb{:});
noTouch_Move_LowProb = ft_timelockgrandaverage(cfg, beta.noTouch_Move_LowProb{:});
noTouch_Move_EqualProb = ft_timelockgrandaverage(cfg, beta.noTouch_Move_EqualProb{:});
Touch_Stay_HighProb = ft_timelockgrandaverage(cfg, beta.Touch_Stay_HighProb{:});
Touch_Stay_LowProb = ft_timelockgrandaverage(cfg, beta.Touch_Stay_LowProb{:});
Touch_Stay_EqualProb = ft_timelockgrandaverage(cfg, beta.Touch_Stay_EqualProb{:});
Touch_Move_HighProb = ft_timelockgrandaverage(cfg, beta.Touch_Move_HighProb{:});
Touch_Move_LowProb = ft_timelockgrandaverage(cfg, beta.Touch_Move_LowProb{:});
Touch_Move_EqualProb = ft_timelockgrandaverage(cfg, beta.Touch_Move_EqualProb{:});

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

%con imgs
cfg = []; 
con_Touch_Stay_HighProb = ft_timelockgrandaverage(cfg, con.Touch_Stay_HighProb{:});
con_Touch_Stay_EqualProb = ft_timelockgrandaverage(cfg, con.Touch_Stay_EqualProb{:});
con_Touch_Stay_LowProb = ft_timelockgrandaverage(cfg, con.Touch_Stay_LowProb{:});
con_Touch_Move_HighProb = ft_timelockgrandaverage(cfg, con.Touch_Move_HighProb{:});
con_Touch_Move_EqualProb = ft_timelockgrandaverage(cfg, con.Touch_Move_EqualProb{:});
con_Touch_Move_LowProb = ft_timelockgrandaverage(cfg, con.Touch_Move_LowProb{:});

%% Check that clusters imgs are extracted
load(strcat(results_second_level, '/SPM.mat'))
for i = 5:8
    matlabbatch{1}.spm.stats.results.spmmat = cellstr(strcat(results_second_level, '/SPM.mat'));
    matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
    matlabbatch{1}.spm.stats.results.conspec.contrasts = i;
    matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
    matlabbatch{1}.spm.stats.results.conspec.thresh = 0.001;
    matlabbatch{1}.spm.stats.results.conspec.extent = 0;
    matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
    matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
    matlabbatch{1}.spm.stats.results.units = 2;
    matlabbatch{1}.spm.stats.results.export{1}.nary.basename = SPM.xCon(i).name;
    spm_jobman('run', matlabbatch); 
    clear matlabbatch
end

%% PLOT INTERACTION TOUCH x PROB
%create folder for imgs
img_folder = strcat(results_second_level, '/00results_interaction_Touch_x_Prob');
if ~exist(img_folder)
    mkdir(img_folder)
end

%avg over conditions for plotting
cfg = []; 
TouchHighProb = ft_timelockgrandaverage(cfg, beta.Touch_Stay_HighProb{:}, beta.Touch_Move_HighProb{:});
noTouchHighProb = ft_timelockgrandaverage(cfg, beta.noTouch_Stay_HighProb{:}, beta.noTouch_Move_HighProb{:});
TouchLowProb = ft_timelockgrandaverage(cfg, beta.Touch_Stay_LowProb{:}, beta.Touch_Move_LowProb{:});
noTouchLowProb = ft_timelockgrandaverage(cfg, beta.noTouch_Stay_LowProb{:}, beta.noTouch_Move_LowProb{:});
TouchEqualProb = ft_timelockgrandaverage(cfg, beta.Touch_Stay_EqualProb{:}, beta.Touch_Move_EqualProb{:});
noTouchEqualProb = ft_timelockgrandaverage(cfg, beta.noTouch_Stay_EqualProb{:}, beta.noTouch_Move_EqualProb{:});
con_TouchHighProb = ft_timelockgrandaverage(cfg, con.Touch_Move_HighProb{:}, con.Touch_Stay_HighProb{:});
con_TouchEqualProb = ft_timelockgrandaverage(cfg, con.Touch_Move_EqualProb{:}, con.Touch_Stay_EqualProb{:});
con_TouchLowProb = ft_timelockgrandaverage(cfg, con.Touch_Move_LowProb{:}, con.Touch_Stay_LowProb{:});

TouchHighProb.err = sqrt(TouchHighProb.var)./sqrt(TouchHighProb.dof); 
noTouchHighProb.err = sqrt(noTouchHighProb.var)./sqrt(noTouchHighProb.dof); 
TouchLowProb.err = sqrt(TouchLowProb.var)./sqrt(TouchLowProb.dof); 
noTouchLowProb.err = sqrt(noTouchLowProb.var)./sqrt(noTouchLowProb.dof); 
TouchEqualProb.err = sqrt(TouchEqualProb.var)./sqrt(TouchEqualProb.dof); 
noTouchEqualProb.err = sqrt(noTouchEqualProb.var)./sqrt(noTouchEqualProb.dof); 
con_TouchHighProb.err = sqrt(con_TouchHighProb.var)./sqrt(con_TouchHighProb.dof); 
con_TouchEqualProb.err = sqrt(con_TouchEqualProb.var)./sqrt(con_TouchEqualProb.dof); 
con_TouchLowProb.err = sqrt(con_TouchLowProb.var)./sqrt(con_TouchLowProb.dof); 

[clust_chans times centroid time_centroid true_mask] = gian_extract_channels_from_cluster('spmF_0005_Interaction, factor Prob x Touch.nii', 'spmF_0005.nii');


cfg = []; 
cfg.layout = caplocation; 
layout = ft_prepare_layout(cfg);

for cluster_num = 1:length(clust_chans)
    % singleplot beta whole cluster
    cfg = []; 
    cfg.colors = [0, 0, 1; ... %blue
                    0.4, 0.4, 1; ... %light blue
                    0.7, 0.7, 1; ... %very light blue
                    1, 0, 0; ... %red
                    1, 0.4, 0.4 %light red
                    1, 0.7, 0.7]; %very light red
    cfg.linestyle = {'-', '-', '-', '-', '-', '-'};
    cfg.layout = caplocation;
    cfg.mask = true_mask{cluster_num};
    cfg.channel = clust_chans{cluster_num};
    %cfg.channel = "CP1"
    %cfg.baseline = [-0.05 -0.005];
    cfg.xlim = [-0.05 0.5];
    cfg.ylim = [-2 4.5];
    GIAN_plot_data(cfg, noTouchHighProb, noTouchEqualProb, noTouchLowProb, TouchHighProb, TouchEqualProb, TouchLowProb)
    exportgraphics(gcf, strcat(img_folder, '/singleplot_betas_cluster', int2str(cluster_num), '.emf'))
    close all

    % singleplot con whole cluster
    cfg = []; 
    cfg.colors = [1, 0, 1; ... %red
                    1, 0.25, 1 %light red
                    1, 0.5, 1]; %very light red
    cfg.linestyle = {'-', '-', '-'};
    cfg.layout = caplocation;
    cfg.channel = clust_chans{cluster_num};
    cfg.mask = true_mask{cluster_num};
    %cfg.baseline = [-0.05 -0.005];
    cfg.xlim = [-0.05 0.5];
    cfg.ylim = [-2 4.5];
    cfg.printlayers = 1;
    cfg.output_printlayers = strcat(img_folder, '/singleplot_cons_cluster', int2str(cluster_num));
    GIAN_plot_data(cfg, con_TouchHighProb, con_TouchEqualProb, con_TouchLowProb)
    exportgraphics(gcf, strcat(img_folder, '/singleplot_cons_cluster', int2str(cluster_num), '.emf'))
    close all

    % singleplot beta centroid
    cfg = []; 
    cfg.colors = [0, 0, 1; ... %blue
                    0.4, 0.4, 1; ... %light blue
                    0.7, 0.7, 1; ... %very light blue
                    1, 0, 0; ... %red
                    1, 0.4, 0.4 %light red
                    1, 0.7, 0.7]; %very light red
    cfg.linestyle = {'-', '-', '-', '-', '-', '-'};
    cfg.layout = caplocation;
    cfg.mask = true_mask{cluster_num};
    cfg.channel = centroid{cluster_num};
    %cfg.channel = "CP1"
    %cfg.baseline = [-0.05 -0.005];
    cfg.xlim = [-0.05 0.5];
    cfg.ylim = [-2 4.5];
    GIAN_plot_data(cfg, noTouchHighProb, noTouchEqualProb, noTouchLowProb, TouchHighProb, TouchEqualProb, TouchLowProb)
    exportgraphics(gcf, strcat(img_folder, '/singleplot_betas_centroid', int2str(cluster_num), '.emf'))
    close all

    % singleplot con centroid
    cfg = []; 
    cfg.colors = [1, 0, 1; ... %red
                    1, 0.25, 1 %light red
                    1, 0.5, 1]; %very light red
    cfg.linestyle = {'-', '-', '-'};
    cfg.layout = caplocation;
    cfg.channel = centroid{cluster_num};
    cfg.mask = true_mask{cluster_num};
    %cfg.baseline = [-0.05 -0.005];
    cfg.xlim = [-0.05 0.5];
    cfg.ylim = [-4 5.5];
    cfg.printlayers = 1;
    cfg.output_printlayers = strcat(img_folder, '/singleplot_cons_centroid', int2str(cluster_num));
    GIAN_plot_data(cfg, con_TouchHighProb, con_TouchEqualProb, con_TouchLowProb)
    exportgraphics(gcf, strcat(img_folder, '/singleplot_cons_centroid', int2str(cluster_num), '.emf'))
    close all

    %topoplot centroid
    cfg = [];
    cfg.operation = '(x1*(-1) + x2*(0) + x3*1) - (x4*(-1) + x5*(0) + x6*1)';
    cfg.parameter = 'avg';
    prova_prosciutto = ft_math(cfg, TouchHighProb, TouchEqualProb, TouchLowProb, noTouchHighProb, noTouchEqualProb, noTouchLowProb); 
    timelims = Touch_Move_HighProb.time(times{cluster_num}); 
    cfg = [];
    cfg.layout = caplocation;
    %cfg.baseline = [-0.05 -0.005];
    cfg.xlim = [Touch_Move_HighProb.time(time_centroid{cluster_num})-0.005 Touch_Move_HighProb.time(time_centroid{cluster_num})+0.005];
    cfg.zlim = [-1 1];
    cfg.colorbar = 'no'; 
    cfg.marker = 'off';
    ft_topoplotER(cfg, prova_prosciutto);
    hold on
    index = arrayfun(@(k) find(strncmp(centroid{cluster_num}{k},data.label,3)), 1:length(centroid{cluster_num}));
    for i = 1:length(index)
        plot(layout.pos(index,1),layout.pos(index,2),'o', 'Color', [0 0 0], 'MarkerSize', 30, 'MarkerFaceColor', [1 1 1], 'LineWidth', 3)
    end
    set(gcf,'units','normalized','outerpos',[0 0 1 1.2]);    
    exportgraphics(gcf, strcat(img_folder, '/topoplot_centroid', int2str(cluster_num), '.emf'))
    close all

    %topoplot cluster
    cfg = [];
    cfg.layout = caplocation;
    %cfg.baseline = [-0.05 -0.005];
    cfg.xlim = [Touch_Move_HighProb.time(times{cluster_num}(1)) Touch_Move_HighProb.time(times{cluster_num}(end))];
    cfg.zlim = [-1 1];
    cfg.colorbar = 'no';  
    cfg.marker = 'off';
    ft_topoplotER(cfg, prova_prosciutto)
    hold on
    index = arrayfun(@(k) find(strncmp(clust_chans{cluster_num}{k},data.label,3)), 1:length(clust_chans{cluster_num}));
    for i = 1:length(index)
        plot(layout.pos(index,1),layout.pos(index,2),'o', 'Color', [0 0 0], 'MarkerSize', 30, 'MarkerFaceColor', [1 1 1], 'LineWidth', 3)
    end
    set(gcf,'units','normalized','outerpos',[0 0 1 1.2]);    
    exportgraphics(gcf, strcat(img_folder, '/topoplot_cluster', int2str(cluster_num), '.emf'))
    close all

    %topoplot over time
    cfg = [];
    cfg.layout = caplocation;
    %cfg.baseline = [-0.05 -0.005];
    cfg.xlim = [0:0.025:0.5];
    cfg.zlim = [-1 1];
    cfg.marker = 'off';
    ft_topoplotER(cfg, prova_prosciutto)
    set(gcf,'units','normalized','outerpos',[0 0 1 1.2]);    
    exportgraphics(gcf, strcat(img_folder, '/topoplot_timedetail', int2str(cluster_num), '.emf'))
    close all

    %topoplot sequence, significant electrodes highligthed. 
    try
        prova_crudo = prova_prosciutto;
        prova_crudo.avg((true_mask{cluster_num}==0)) = nan;
        cfg = [];
        cfg.layout = caplocation;
        %cfg.baseline = [-0.05 -0.005];
        cfg.xlim = [0 min(timelims):((max(timelims)-min(timelims))/5):max(timelims) 0.5];
        cfg.zlim = [-1 1];
        cfg.colormap = '*RdBu';
        ft_topoplotER(cfg, prova_crudo)
        set(gcf,'units','normalized','outerpos',[0 0 1 1.2]);    
        exportgraphics(gcf, strcat(img_folder, '/topoplot_timedetail_cluster', int2str(cluster_num), '.emf'))
        close all
    end

    %multiplot betas
    figure;
    noTouchHighProb.mask = true_mask{cluster_num}; 
    cfg = []; 
    cfg.linecolor = [0, 0, 1; ... %blue
                    0.4, 0.4, 1; ... %light blue
                    0.7, 0.7, 1; ... %very light blue
                    1, 0, 0; ... %red
                    1, 0.4, 0.4 %light red
                    1, 0.7, 0.7]; %very light red
    cfg.layout = caplocation;
    %cfg.baseline = [-0.05 -0.005];
    cfg.maskparameter = 'mask';
    cfg.maskfacealpha = 0.2;
    set(gcf,'units','normalized','outerpos',[0 0 1 1.2]);    
    ft_multiplotER(cfg, noTouchHighProb, noTouchEqualProb, noTouchLowProb, TouchHighProb, TouchEqualProb, TouchLowProb)
    exportgraphics(gcf, strcat(img_folder, '/multiplot_beta_cluster', int2str(cluster_num), '.emf'))
    close all

    %multiplot cons
    figure;
    con_TouchHighProb.mask = true_mask{cluster_num}; 
    cfg = []; 
    cfg.linecolor = [1, 0, 0; ... %red
                    1, 0.4, 0.4 %light red
                    1, 0.7, 0.7]; %very light red
    cfg.layout = caplocation;
    %cfg.baseline = [-0.05 -0.005];
    cfg.maskparameter = 'mask';
    cfg.maskfacealpha = 0.2;
    set(gcf,'units','normalized','outerpos',[0 0 1 1.2]);    
    ft_multiplotER(cfg, con_TouchHighProb, con_TouchEqualProb, con_TouchLowProb)
    exportgraphics(gcf, strcat(img_folder, '/multiplot_cons_cluster', int2str(cluster_num), '.emf'))
    close all

    %barplot cluster
    index = arrayfun(@(k) find(strncmp(clust_chans{cluster_num}{k},data.label,3)), 1:length(clust_chans{cluster_num}));
    order_scans = ["noTouch_Stay_LowProb", "noTouch_Stay_EqualProb", "noTouch_Stay_HighProb",...
       "noTouch_Move_LowProb", "noTouch_Move_EqualProb", "noTouch_Move_HighProb"...
       "Touch_Stay_LowProb", "Touch_Stay_EqualProb", "Touch_Stay_HighProb",...
       "Touch_Move_LowProb", "Touch_Move_EqualProb", "Touch_Move_HighProb"]; 
    bar_values = []; 
    for i = 1:length(order_scans)
        values = []; 
        for y = 1:length(subjects)
            values = [values mean(beta.(order_scans(i)){y}.avg(index, times{cluster_num}), [1 2])];
        end
        bar_values.(order_scans(i)).mean = mean(values); 
        bar_values.(order_scans(i)).var = var(values); 
    end
 
    colors_to_plot = [0.7, 0.7, 1; ... %very light blue
                    0.4, 0.4, 1; ... %light blue
                    0, 0, 1; ... %blue
                    1, 0.7, 0.7;... %very light red
                    1, 0.4, 0.4; ... %light red
                    1, 0, 0]; %red
    figure; 
    hold on
        %means_barplot = mean(beta.(order_scans(i)).avg(index, times{cluster_num}));
    aaa = bar(1, bar_values.noTouch_Stay_LowProb.mean, 'FaceColor', colors_to_plot(1,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(2, bar_values.noTouch_Stay_EqualProb.mean, 'FaceColor', colors_to_plot(2,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(3, bar_values.noTouch_Stay_HighProb.mean, 'FaceColor', colors_to_plot(3,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(4, bar_values.Touch_Stay_LowProb.mean, 'FaceColor', colors_to_plot(4,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(5, bar_values.Touch_Stay_EqualProb.mean, 'FaceColor', colors_to_plot(5,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(6, bar_values.Touch_Stay_HighProb.mean, 'FaceColor', colors_to_plot(6,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(7, bar_values.noTouch_Move_LowProb.mean, 'FaceColor', colors_to_plot(1,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(1,:))
    aaa = bar(8, bar_values.noTouch_Move_EqualProb.mean, 'FaceColor', colors_to_plot(2,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(2,:))
    aaa = bar(9, bar_values.noTouch_Move_HighProb.mean, 'FaceColor', colors_to_plot(3,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(3,:))
    aaa = bar(10, bar_values.Touch_Move_LowProb.mean, 'FaceColor', colors_to_plot(4,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(4,:))
    aaa = bar(11, bar_values.Touch_Move_EqualProb.mean, 'FaceColor', colors_to_plot(5,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(5,:))
    aaa = bar(12, bar_values.Touch_Move_HighProb.mean, 'FaceColor', colors_to_plot(6,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(6,:))
    errorbar([1],  bar_values.noTouch_Stay_LowProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([2],  bar_values.noTouch_Stay_EqualProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([3],  bar_values.noTouch_Stay_HighProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([4],  bar_values.Touch_Stay_LowProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([5],  bar_values.Touch_Stay_EqualProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([6],  bar_values.Touch_Stay_HighProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([7],  bar_values.noTouch_Move_LowProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([8],  bar_values.noTouch_Move_EqualProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([9],  bar_values.noTouch_Move_HighProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([10],  bar_values.Touch_Move_LowProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([11],  bar_values.Touch_Move_EqualProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([12],  bar_values.Touch_Move_HighProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    set(gca,'xticklabel',{[]})
    set(gca,'xtick',[])
    ylim([-1 4.8])
    exportgraphics(gcf, strcat(img_folder, '/barplot_wholemask_cluster', int2str(cluster_num), '.emf'))
    close all


    %barplot centroid
    index = arrayfun(@(k) find(strncmp(centroid{cluster_num}{k},data.label,3)), 1:length(centroid{cluster_num}));
    order_scans = ["noTouch_Stay_LowProb", "noTouch_Stay_EqualProb", "noTouch_Stay_HighProb",...
       "noTouch_Move_LowProb", "noTouch_Move_EqualProb", "noTouch_Move_HighProb"...
       "Touch_Stay_LowProb", "Touch_Stay_EqualProb", "Touch_Stay_HighProb",...
       "Touch_Move_LowProb", "Touch_Move_EqualProb", "Touch_Move_HighProb"]; 
    bar_values = []; 
    for i = 1:length(order_scans)
        values = []; 
        for y = 1:length(subjects)
            values = [values mean(beta.(order_scans(i)){y}.avg(index, time_centroid{cluster_num}), [1 2])];
        end
        bar_values.(order_scans(i)).mean = mean(values); 
        bar_values.(order_scans(i)).var = var(values); 
    end
 
    colors_to_plot = [0.7, 0.7, 1; ... %very light blue
                    0.4, 0.4, 1; ... %light blue
                    0, 0, 1; ... %blue
                    1, 0.7, 0.7;... %very light red
                    1, 0.4, 0.4; ... %light red
                    1, 0, 0]; %red
    figure; 
    hold on
        %means_barplot = mean(beta.(order_scans(i)).avg(index, times{cluster_num}));
    aaa = bar(1, bar_values.noTouch_Stay_LowProb.mean, 'FaceColor', colors_to_plot(1,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(2, bar_values.noTouch_Stay_EqualProb.mean, 'FaceColor', colors_to_plot(2,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(3, bar_values.noTouch_Stay_HighProb.mean, 'FaceColor', colors_to_plot(3,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(4, bar_values.Touch_Stay_LowProb.mean, 'FaceColor', colors_to_plot(4,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(5, bar_values.Touch_Stay_EqualProb.mean, 'FaceColor', colors_to_plot(5,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(6, bar_values.Touch_Stay_HighProb.mean, 'FaceColor', colors_to_plot(6,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(7, bar_values.noTouch_Move_LowProb.mean, 'FaceColor', colors_to_plot(1,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(1,:))
    aaa = bar(8, bar_values.noTouch_Move_EqualProb.mean, 'FaceColor', colors_to_plot(2,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(2,:))
    aaa = bar(9, bar_values.noTouch_Move_HighProb.mean, 'FaceColor', colors_to_plot(3,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(3,:))
    aaa = bar(10, bar_values.Touch_Move_LowProb.mean, 'FaceColor', colors_to_plot(4,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(4,:))
    aaa = bar(11, bar_values.Touch_Move_EqualProb.mean, 'FaceColor', colors_to_plot(5,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(5,:))
    aaa = bar(12, bar_values.Touch_Move_HighProb.mean, 'FaceColor', colors_to_plot(6,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(6,:))
    errorbar([1],  bar_values.noTouch_Stay_LowProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([2],  bar_values.noTouch_Stay_EqualProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([3],  bar_values.noTouch_Stay_HighProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([4],  bar_values.Touch_Stay_LowProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([5],  bar_values.Touch_Stay_EqualProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([6],  bar_values.Touch_Stay_HighProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([7],  bar_values.noTouch_Move_LowProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([8],  bar_values.noTouch_Move_EqualProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([9],  bar_values.noTouch_Move_HighProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([10],  bar_values.Touch_Move_LowProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([11],  bar_values.Touch_Move_EqualProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([12],  bar_values.Touch_Move_HighProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    set(gca,'xticklabel',{[]})
    set(gca,'xtick',[])
    ylim([-1 6])
    exportgraphics(gcf, strcat(img_folder, '/barplot_wholemask_centroid', int2str(cluster_num), '.emf'))
    close all

end

%% PLOT INTERACTION TOUCH x MOVE
%create folder for imgs
img_folder = strcat(results_second_level, '/01results_interaction_Touch_x_Move');
if ~exist(img_folder)
    mkdir(img_folder)
end

%avg over conditions for plotting
cfg = []; 
cfg.operation = '(x1+x2+x3)/3';
cfg.parameter = {'avg', 'var', 'dof'};

cfg = []; 
noTouchMove = ft_timelockgrandaverage(cfg, beta.noTouch_Move_HighProb{:}, beta.noTouch_Move_EqualProb{:}, beta.noTouch_Move_LowProb{:});
noTouchStay = ft_timelockgrandaverage(cfg, beta.noTouch_Stay_HighProb{:}, beta.noTouch_Stay_EqualProb{:}, beta.noTouch_Stay_LowProb{:});
TouchMove = ft_timelockgrandaverage(cfg, beta.Touch_Move_HighProb{:}, beta.Touch_Move_EqualProb{:}, beta.Touch_Move_LowProb{:});
TouchStay = ft_timelockgrandaverage(cfg, beta.Touch_Stay_HighProb{:}, beta.Touch_Stay_EqualProb{:}, beta.Touch_Stay_LowProb{:});
con_TouchStay = ft_timelockgrandaverage(cfg, con.Touch_Stay_HighProb{:}, con.Touch_Stay_EqualProb{:}, con.Touch_Stay_LowProb{:});
con_TouchMove = ft_timelockgrandaverage(cfg, con.Touch_Move_HighProb{:}, con.Touch_Move_EqualProb{:}, con.Touch_Move_LowProb{:});

noTouchMove.err = sqrt(noTouchMove.var)./sqrt(noTouchMove.dof);
noTouchStay.err = sqrt(noTouchStay.var)./sqrt(noTouchStay.dof);
TouchMove.err = sqrt(TouchMove.var)./sqrt(TouchMove.dof);
TouchStay.err = sqrt(TouchStay.var)./sqrt(TouchStay.dof);
con_TouchStay.err = sqrt(con_TouchStay.var)./sqrt(con_TouchStay.dof);
con_TouchMove.err = sqrt(con_TouchMove.var)./sqrt(con_TouchMove.dof);

[clust_chans times centroid time_centroid true_mask] = gian_extract_channels_from_cluster('spmF_0006_Interaction, factor Movement x Touch.nii', 'spmF_0006.nii');

for cluster_num = 1:length(clust_chans)
    % singleplot beta whole cluster
    cfg = []; 
    cfg.colors = [0, 0, 1; ... %blue
                    1, 0, 0; ... %red
                    0, 0, 1; ... %blue
                    1, 0, 0]; %red
    cfg.linestyle = {'--', '--', '-', '-'};
    cfg.layout = caplocation;
    cfg.mask = true_mask{cluster_num};
    cfg.channel = clust_chans{cluster_num};
    %cfg.channel = "CP1"
    %cfg.baseline = [-0.05 -0.005];
    cfg.xlim = [-0.05 0.5];
    cfg.ylim = [-2 4.5];
    GIAN_plot_data(cfg, noTouchMove, TouchMove, noTouchStay, TouchStay)
    exportgraphics(gcf, strcat(img_folder, '/singleplot_betas_cluster', int2str(cluster_num), '.emf'))
    close all

    % singleplot con whole cluster
    cfg = []; 
    cfg.colors = [1, 0, 1; ... %red
                    1, 0, 1]; %red
    cfg.linestyle = {'-', '--'};
    cfg.layout = caplocation;
    cfg.channel = clust_chans{cluster_num};
    cfg.mask = true_mask{cluster_num};
    %cfg.baseline = [-0.05 -0.005];
    cfg.xlim = [-0.05 0.5];
    if cluster_num < 5 
        cfg.ylim = [-2 4.5];
    else
        cfg.ylim = [-4 2];
    end
    cfg.printlayers = 1;
    cfg.output_printlayers = strcat(img_folder, '/singleplot_cons_cluster', int2str(cluster_num));
    GIAN_plot_data(cfg, con_TouchStay, con_TouchMove)
    exportgraphics(gcf, strcat(img_folder, '/singleplot_cons_cluster', int2str(cluster_num), '.emf'))
    close all

    % singleplot beta centroid
    cfg = []; 
    cfg.colors = [0, 0, 1; ... %blue
                    1, 0, 0; ... %red
                    0, 0, 1; ... %blue
                    1, 0, 0]; %red
    cfg.linestyle = {'--', '--', '-', '-'};
    cfg.layout = caplocation;
    cfg.mask = true_mask{cluster_num};
    cfg.channel = centroid{cluster_num};
    %cfg.channel = "CP1"
    %cfg.baseline = [-0.05 -0.005];
    cfg.xlim = [-0.05 0.5];
    cfg.ylim = [-2 4.5];
    GIAN_plot_data(cfg, noTouchMove, TouchMove, noTouchStay, TouchStay)
    exportgraphics(gcf, strcat(img_folder, '/singleplot_betas_centroid', int2str(cluster_num), '.emf'))
    close all

    % singleplot con centroid
    cfg = []; 
    cfg.colors = [1, 0, 1; ... %red
                    1, 0, 1]; %red
    cfg.linestyle = {'-', '--'};
    cfg.layout = caplocation;
    cfg.channel = centroid{cluster_num};
    cfg.mask = true_mask{cluster_num};
    cfg.xlim = [-0.05 0.5];
    if cluster_num < 5 
        cfg.ylim = [-2 4.5];
    else
        cfg.ylim = [-4 2];
    end
    cfg.printlayers = 1;
    cfg.output_printlayers = strcat(img_folder, '/singleplot_cons_centroid', int2str(cluster_num));
    %cfg.baseline = [-0.05 -0.005];
    GIAN_plot_data(cfg, con_TouchStay, con_TouchMove)
    exportgraphics(gcf, strcat(img_folder, '/singleplot_cons_centroid', int2str(cluster_num), '.emf'))
    close all

    %topoplot centroid
    cfg = [];
    cfg.operation = '(x1*(-1) + x2*(1)) - (x3*(-1) + x4*(1))';
    cfg.parameter = 'avg';
    prova_prosciutto = ft_math(cfg, TouchStay, TouchMove, noTouchStay, noTouchMove); 
    %timelims = Touch_Move_HighProb.time(times{cluster_num}); 
    cfg = [];
    cfg.layout = caplocation;
    %cfg.baseline = [-0.05 -0.005];
    cfg.xlim = [Touch_Move_HighProb.time(time_centroid{cluster_num})-0.005 Touch_Move_HighProb.time(time_centroid{cluster_num})+0.005];
    cfg.zlim = [-1 1];
    cfg.marker = 'off';
    cfg.colorbar = 'no'; 
    ft_topoplotER(cfg, prova_prosciutto)
        hold on
    index = arrayfun(@(k) find(strncmp(centroid{cluster_num}{k},data.label,3)), 1:length(centroid{cluster_num}));
    for i = 1:length(index)
        plot(layout.pos(index,1),layout.pos(index,2),'o', 'Color', [0 0 0], 'MarkerSize', 30, 'MarkerFaceColor', [1 1 1], 'LineWidth', 3)
    end
    set(gcf,'units','normalized','outerpos',[0 0 1 1.2]);    
    exportgraphics(gcf, strcat(img_folder, '/topoplot_centroid', int2str(cluster_num), '.emf'))
    close all

    %topoplot cluster
    cfg = [];
    cfg.layout = caplocation;
    %cfg.baseline = [-0.05 -0.005];
    cfg.xlim = [Touch_Move_HighProb.time(times{cluster_num}(1)) Touch_Move_HighProb.time(times{cluster_num}(end))];
    cfg.zlim = [-1 1];
    cfg.colorbar = 'yes'; 
    cfg.marker = 'off';
    ft_topoplotER(cfg, prova_prosciutto)
    hold on
    index = arrayfun(@(k) find(strncmp(clust_chans{cluster_num}{k},data.label,3)), 1:length(clust_chans{cluster_num}));
    for i = 1:length(index)
        plot(layout.pos(index,1),layout.pos(index,2),'o', 'Color', [0 0 0], 'MarkerSize', 30, 'MarkerFaceColor', [1 1 1], 'LineWidth', 3)
    end
    set(gcf,'units','normalized','outerpos',[0 0 1 1.2]);    
    exportgraphics(gcf, strcat(img_folder, '/topoplot_cluster', int2str(cluster_num), '.emf'))
    close all

    %topoplot over time
    cfg = [];
    cfg.layout = caplocation;
    %cfg.baseline = [-0.05 -0.005];
    cfg.xlim = [0:0.025:0.5];
    cfg.zlim = [-1 1];
    cfg.marker = 'off';
    ft_topoplotER(cfg, prova_prosciutto)
    set(gcf,'units','normalized','outerpos',[0 0 1 1.2]);    
    exportgraphics(gcf, strcat(img_folder, '/topoplot_timedetail', int2str(cluster_num), '.emf'))
    close all

    %topoplot sequence, significant electrodes highligthed. 
    prova_crudo = prova_prosciutto;
    prova_crudo.avg((true_mask{cluster_num}==0)) = 0.0001;
    timelims = Touch_Move_HighProb.time(times{cluster_num}); 
    cfg = [];
    cfg.layout = caplocation;
    %cfg.baseline = [-0.05 -0.005];
    cfg.xlim = [0 min(timelims):((max(timelims)-min(timelims))/10):max(timelims) 0.5];
    cfg.zlim = [-1 1];
    cfg.colormap = '*RdBu';
    cfg.marker = 'off';
    ft_topoplotER(cfg, prova_crudo)
    set(gcf,'units','normalized','outerpos',[0 0 1 1.2]);    
    exportgraphics(gcf, strcat(img_folder, '/topoplot_timedetail_cluster', int2str(cluster_num), '.emf'))
    close all

    %multiplot betas
    figure;
    noTouchMove.mask = true_mask{cluster_num}; 
    cfg = []; 
    cfg.linecolor = [0, 0, 1; ... %blue
                    1, 0, 0; ... %red
                    0, 0, 1; ... %blue
                    1, 0, 0]; %red
    cfg.linestyle = {'--', '--', '-', '-'};
    cfg.layout = caplocation;
    cfg.maskparameter = 'mask';
    cfg.maskfacealpha = 0.2;
    %cfg.channel = cellstr(centroid{cluster_num});
    %cfg.baseline = [-0.05 -0.005];
    ft_multiplotER(cfg, noTouchMove, TouchMove, noTouchStay, TouchStay)
    set(gcf,'units','normalized','outerpos',[0 0 1 1.2]);    
    exportgraphics(gcf, strcat(img_folder, '/multiplot_beta_cluster', int2str(cluster_num), '.emf'))
    close all

    %multiplot cons
    figure;
    con_TouchStay.mask = true_mask{cluster_num};
    cfg = []; 
    cfg.linecolor = [1, 0, 0; ... %red
                    1, 0, 0]; %red
    cfg.linestyle = {'-', '--'};
    cfg.layout = caplocation;
    cfg.maskparameter = 'mask';
    cfg.maskfacealpha = 0.2;
    %cfg.channel = cellstr(centroid{cluster_num});
    % cfg.baseline = [-0.05 -0.005];
    ft_multiplotER(cfg, con_TouchStay, con_TouchMove)
    set(gcf,'units','normalized','outerpos',[0 0 1 1.2]);    
    exportgraphics(gcf, strcat(img_folder, '/multiplot_cons_cluster', int2str(cluster_num), '.emf'))
    close all

    %barplot cluster
    index = arrayfun(@(k) find(strncmp(clust_chans{cluster_num}{k},data.label,3)), 1:length(clust_chans{cluster_num}));
    order_scans = ["noTouch_Stay_LowProb", "noTouch_Stay_EqualProb", "noTouch_Stay_HighProb",...
       "noTouch_Move_LowProb", "noTouch_Move_EqualProb", "noTouch_Move_HighProb"...
       "Touch_Stay_LowProb", "Touch_Stay_EqualProb", "Touch_Stay_HighProb",...
       "Touch_Move_LowProb", "Touch_Move_EqualProb", "Touch_Move_HighProb"]; 
    bar_values = []; 
    for i = 1:length(order_scans)
        values = []; 
        for y = 1:length(subjects)
            values = [values mean(beta.(order_scans(i)){y}.avg(index, times{cluster_num}), [1 2])];
        end
        bar_values.(order_scans(i)).mean = mean(values); 
        bar_values.(order_scans(i)).var = var(values); 
    end
 
    colors_to_plot = [0.7, 0.7, 1; ... %very light blue
                    0.4, 0.4, 1; ... %light blue
                    0, 0, 1; ... %blue
                    1, 0.7, 0.7;... %very light red
                    1, 0.4, 0.4; ... %light red
                    1, 0, 0]; %red
    figure; 
    hold on
        %means_barplot = mean(beta.(order_scans(i)).avg(index, times{cluster_num}));
    aaa = bar(1, bar_values.noTouch_Stay_LowProb.mean, 'FaceColor', colors_to_plot(1,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(2, bar_values.noTouch_Stay_EqualProb.mean, 'FaceColor', colors_to_plot(2,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(3, bar_values.noTouch_Stay_HighProb.mean, 'FaceColor', colors_to_plot(3,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(4, bar_values.Touch_Stay_LowProb.mean, 'FaceColor', colors_to_plot(4,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(5, bar_values.Touch_Stay_EqualProb.mean, 'FaceColor', colors_to_plot(5,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(6, bar_values.Touch_Stay_HighProb.mean, 'FaceColor', colors_to_plot(6,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(7, bar_values.noTouch_Move_LowProb.mean, 'FaceColor', colors_to_plot(1,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(1,:))
    aaa = bar(8, bar_values.noTouch_Move_EqualProb.mean, 'FaceColor', colors_to_plot(2,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(2,:))
    aaa = bar(9, bar_values.noTouch_Move_HighProb.mean, 'FaceColor', colors_to_plot(3,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(3,:))
    aaa = bar(10, bar_values.Touch_Move_LowProb.mean, 'FaceColor', colors_to_plot(4,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(4,:))
    aaa = bar(11, bar_values.Touch_Move_EqualProb.mean, 'FaceColor', colors_to_plot(5,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(5,:))
    aaa = bar(12, bar_values.Touch_Move_HighProb.mean, 'FaceColor', colors_to_plot(6,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(6,:))
    errorbar([1],  bar_values.noTouch_Stay_LowProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([2],  bar_values.noTouch_Stay_EqualProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([3],  bar_values.noTouch_Stay_HighProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([4],  bar_values.Touch_Stay_LowProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([5],  bar_values.Touch_Stay_EqualProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([6],  bar_values.Touch_Stay_HighProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([7],  bar_values.noTouch_Move_LowProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([8],  bar_values.noTouch_Move_EqualProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([9],  bar_values.noTouch_Move_HighProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([10],  bar_values.Touch_Move_LowProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([11],  bar_values.Touch_Move_EqualProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([12],  bar_values.Touch_Move_HighProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    set(gca,'xticklabel',{[]})
    set(gca,'xtick',[])
    if cluster_num < 5 
        ylim([-1 4.8])
    else
        cfg.ylim = [-4.5 1];
    end
    exportgraphics(gcf, strcat(img_folder, '/barplot_wholemask_cluster', int2str(cluster_num), '.emf'))
    close all


    %barplot centroid
    index = arrayfun(@(k) find(strncmp(centroid{cluster_num}{k},data.label,3)), 1:length(centroid{cluster_num}));
    order_scans = ["noTouch_Stay_LowProb", "noTouch_Stay_EqualProb", "noTouch_Stay_HighProb",...
       "noTouch_Move_LowProb", "noTouch_Move_EqualProb", "noTouch_Move_HighProb"...
       "Touch_Stay_LowProb", "Touch_Stay_EqualProb", "Touch_Stay_HighProb",...
       "Touch_Move_LowProb", "Touch_Move_EqualProb", "Touch_Move_HighProb"]; 
    bar_values = []; 
    for i = 1:length(order_scans)
        values = []; 
        for y = 1:length(subjects)
            values = [values mean(beta.(order_scans(i)){y}.avg(index, time_centroid{cluster_num}), [1 2])];
        end
        bar_values.(order_scans(i)).mean = mean(values); 
        bar_values.(order_scans(i)).var = var(values); 
    end
 
    colors_to_plot = [0.7, 0.7, 1; ... %very light blue
                    0.4, 0.4, 1; ... %light blue
                    0, 0, 1; ... %blue
                    1, 0.7, 0.7;... %very light red
                    1, 0.4, 0.4; ... %light red
                    1, 0, 0]; %red
    figure; 
    hold on
        %means_barplot = mean(beta.(order_scans(i)).avg(index, times{cluster_num}));
    aaa = bar(1, bar_values.noTouch_Stay_LowProb.mean, 'FaceColor', colors_to_plot(1,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(2, bar_values.noTouch_Stay_EqualProb.mean, 'FaceColor', colors_to_plot(2,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(3, bar_values.noTouch_Stay_HighProb.mean, 'FaceColor', colors_to_plot(3,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(4, bar_values.Touch_Stay_LowProb.mean, 'FaceColor', colors_to_plot(4,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(5, bar_values.Touch_Stay_EqualProb.mean, 'FaceColor', colors_to_plot(5,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(6, bar_values.Touch_Stay_HighProb.mean, 'FaceColor', colors_to_plot(6,:), 'BarWidth',1, 'FaceAlpha', 0.3);
    aaa = bar(7, bar_values.noTouch_Move_LowProb.mean, 'FaceColor', colors_to_plot(1,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(1,:))
    aaa = bar(8, bar_values.noTouch_Move_EqualProb.mean, 'FaceColor', colors_to_plot(2,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(2,:))
    aaa = bar(9, bar_values.noTouch_Move_HighProb.mean, 'FaceColor', colors_to_plot(3,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(3,:))
    aaa = bar(10, bar_values.Touch_Move_LowProb.mean, 'FaceColor', colors_to_plot(4,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(4,:))
    aaa = bar(11, bar_values.Touch_Move_EqualProb.mean, 'FaceColor', colors_to_plot(5,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(5,:))
    aaa = bar(12, bar_values.Touch_Move_HighProb.mean, 'FaceColor', colors_to_plot(6,:), 'BarWidth',1);
    %hatchfill2(aaa(1), 'cross','HatchAngle',45,'hatchcolor',colors_to_plot(6,:))
    errorbar([1],  bar_values.noTouch_Stay_LowProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([2],  bar_values.noTouch_Stay_EqualProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([3],  bar_values.noTouch_Stay_HighProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([4],  bar_values.Touch_Stay_LowProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([5],  bar_values.Touch_Stay_EqualProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([6],  bar_values.Touch_Stay_HighProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([7],  bar_values.noTouch_Move_LowProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([8],  bar_values.noTouch_Move_EqualProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([9],  bar_values.noTouch_Move_HighProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([10],  bar_values.Touch_Move_LowProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([11],  bar_values.Touch_Move_EqualProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    errorbar([12],  bar_values.Touch_Move_HighProb.mean, sqrt(bar_values.noTouch_Stay_LowProb.var)/sqrt(length(subjects)),'k','LineWidth', 0.5,'linestyle','none','HandleVisibility','off'); 
    set(gca,'xticklabel',{[]})
    set(gca,'xtick',[])
    if cluster_num < 5 
        ylim([-1 4.8])
    else
        cfg.ylim = [-4.5 1];
    end
    exportgraphics(gcf, strcat(img_folder, '/barplot_wholemask_centroid', int2str(cluster_num), '.emf'))
    close all
end