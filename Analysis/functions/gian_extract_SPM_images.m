function ft_data_matrix = gian_extract_SPM_images(image_path)
    D = spm_eeg_load('E:\Gian\GG_SensAtt_Prediction\02Data\ID04\01EEG\spm\hEOG_tvEOG_thdMID04.mat');
    chanind = D.selectchannels('all');
    n = 32;
    
    [Cel, x, y] = spm_eeg_locate_channels(D, n, chanind);

    image = niftiread(image_path);

    %check that the image is 32x32
    if size(image,1) == 32 && size(image,2) == 32
        disp('converting . . .')
    else
        error('image is not the right one!')
    end

    for i = 1:size(image,3)
        for y = 1:64
            ft_data_matrix(y,i) = image(Cel(y,1), Cel(y,2), i);
        end
    end
end