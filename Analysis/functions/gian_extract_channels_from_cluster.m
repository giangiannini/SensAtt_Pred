function [chans times centroid times_centroid true_mask] = gian_extract_channels_from_cluster(path_to_cluster, path_to_map)
    addpath('E:/02Data/03Utils/Functions/');
    D = spm_eeg_load('E:\Gian\GG_SensAtt_Prediction\02Data\ID04\01EEG\spm\evEOG_thdMID04.mat');
    load('E:\Gian\GG_SensAtt_Prediction\02Data\ID04\021stLevel\00results_Baseline_April2024_fullyManual_participant_responses_noCovariate_avg_control_normal_Contrast\SPM.mat')
    load('E:\Gian\GG_SensAtt_Prediction\02Data\ID04\01EEG\preprocessed_final_SPM_vEOG_hEOGApril2024_fully_manual.mat');
    
    %read cluster positions
    clust = niftiread(path_to_cluster);

    %read cluster values
    map = niftiread(path_to_map); 
    
    for y = 1:5%max(max(max(clust)))
        mask = clust==y;
        [spacex,spacey,time] = ind2sub(size(clust),find(clust == y));
        coords = [spacex spacey];
        if size(coords,1) ~= 1
            chans{y} = []; 
            for i = 1:length(coords)
                chans{y} = [chans{y}, string(get_channels_at_spm_coords('VOX', D, coords(i,:), SPM))];
            end

            positions = []; 
            true_mask{y} = []; 
            %reconstruct a mask to apply to fieldtrip 
            true_mask{y} = zeros(size(data.label,1),size(map,3));
            for n = 1:length(chans{y})
                positions(n) = find(strcmp(chans{y}(n), data.label));
            end
            positions = positions';
            for i = 1:length(positions)
                true_mask{y}(positions(i), time(i)) = true; 
            end

            chans{y} = unique(chans{y}); 
            times{y} = unique(time); 
            

            %centroid
            [coords_centroid_x coords_centroid_y time_centroid] = ind2sub(size(map), find(map == max(map(mask))));
            centroid{y} = string(get_channels_at_spm_coords('VOX', D, [coords_centroid_x coords_centroid_y], SPM));
            times_centroid{y} = time_centroid;
        end
    end
end
