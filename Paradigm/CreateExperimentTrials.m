%% CREATE EXP RUNNER FILE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ID = "60";
trialperblock = 25;
numBlocks = 6; %must be divisible by 3!
numRuns = 4;
trial_pool = [3 4; 1 2];
ITIs_limits = [1750 2250];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% Setup general paths
%Path in which the table will be printed
% outputpath = strcat("C:/Users/CCNB-EEG1/Desktop/Gian/GG_Somato_Vision/02Data/ID", ID, "/00Behavioural/ID", ID, "TrialTable.csv");
folder = strcat("C:/Gian/GG_SensAtt_Prediction/02Data/ID", ID, "/00Behavioural");
if isfolder(folder) == 0
    mkdir(folder);
end
folder1 = strcat("C:/Gian/GG_SensAtt_Prediction/02Data/ID", ID, "/01EEG");
if isfolder(folder1) == 0
    mkdir(folder1);
end
outputpath = strcat(folder, "/ID", ID, "TrialTable.csv");

%% Check that training files are there
dir_struct = dir(fullfile(folder,'/Training_*'));
if numel(dir_struct) == 0
    disp("Training phase still to be run ... here's the file path for Unity")
    disp("C:/Gian/GG_SensAtt_Prediction/02Data/training_TrialTable.csv")
elseif numel(dir_struct) < 3
    error("something went wrong with the training phase");
elseif numel(dir_struct) == 3
    %Create some columns that will be concatenated in a table
    BlockOrder = repmat([0:1:2],1,numBlocks/3); %there are only 3 block types (25/75 /// 50/50 /// 75/25)
    BlockROrder = BlockOrder(randperm(length(BlockOrder)));
    
    RunOrder = 0:1:(numRuns-1);
    RunROrder = RunOrder(randperm(length(RunOrder)));
    
    defTable = table();
    Block_type = [];
    TrialInBlock = [];
    
    %This will loop for the number of trials
    Trial = 0:1:(trialperblock*numBlocks*numRuns - 1);
    Block = reshape(repmat([0:1:numBlocks*numRuns-1],trialperblock,1),[],1);
    Block_num = reshape(repmat([1:numBlocks*numRuns],trialperblock,1),[],1);
    Completed = repmat("false",1,trialperblock*numBlocks*numRuns);
    Attempts = zeros(1,trialperblock*numBlocks*numRuns);
    Skipped = repmat("false",1,trialperblock*numBlocks*numRuns);
    TrialTime = zeros(1,trialperblock*numBlocks*numRuns);
    TrialInBlock = 0:1:(trialperblock-1);
    TrialInBlock = repmat(TrialInBlock, 1, numBlocks*numRuns);
    ITIs = ITIs_limits(1) + (ITIs_limits(2) - ITIs_limits(1)).*rand(trialperblock*numBlocks*numRuns,1);
    
    for m = 1:double(ID)

        % Stuff for speed
        Speed = zeros(1,trialperblock*numBlocks*numRuns); 
        [avg_speed, std_speed] = TakeSpeedfromLog(strcat(folder, '/Training_Log_ID', ID, '.txt')); 
        Speed_vectors = [];
        for i = 1:numBlocks*numRuns
    %         r =  + (b-a).*rand(length(BlockROrder),1)
    
            r = (avg_speed-std_speed) + (avg_speed+std_speed-(avg_speed-std_speed)).*rand(25,1);
            Speed_vectors = [Speed_vectors; r];
        end
                
        % Stuff for random condition extraction
        def_trial_table = [];
        for i=1:length(RunROrder)
            for j = 1:length(BlockROrder)
                trial_table = [];
                trial_type = [];
                Block_type = [];
                trial_type(1) = randi(4,1);
                for k = 1:trialperblock-1
                    %extract pseudorandomised number
                    if trial_type(k) == 1 || trial_type(k) == 4
                        random_next_trial = trial_pool(1,round(rand(1)+1));
                    elseif trial_type(k) == 2 || trial_type(k) == 3
                        random_next_trial = trial_pool(2,round(rand(1)+1));
                    end
                    trial_type(k+1) = random_next_trial;
                end
        
                if BlockROrder(j) == 0 %%% Condition 25/75
                    Block_type = repmat("25/75", trialperblock,1);
                    random_trials = rand(trialperblock,1);
                    stimulation = zeros(trialperblock,1);
                    stimulation(random_trials < 0.25) = 1;
                elseif BlockROrder(j) == 1 %%% Condition 50/50
                    Block_type = repmat("50/50", trialperblock,1);
                    random_trials = rand(trialperblock,1);
                    stimulation = zeros(trialperblock,1);
                    stimulation(random_trials < 0.50) = 1;
                elseif BlockROrder(j) == 2 %%% Condition 75/25
                    Block_type = repmat("75/25", trialperblock,1);
                    random_trials = rand(trialperblock,1);
                    stimulation = zeros(trialperblock,1);
                    stimulation(random_trials < 0.75) = 1;
                end
        
                trial_table = table(Block_type, trial_type', stimulation); 
                def_trial_table = [def_trial_table; trial_table];
        
            end
        % %     for j=1:length(BlockROrder)
        % %         %This will loop for the number of blocks in every trial
        % %         if BlockROrder(j) == 0
        % %             Block_type = [Block_type repmat("TouchVision", 1, trialperblock)];
        % %             Block_num = [Block_num (repmat([j],1,trialperblock)+4*(i-1))];
        % %             Block = [Block (repmat((j-1+4*(i-1)),1,trialperblock))];
        % %         elseif BlockROrder(j) == 1
        % %             Block_type = [Block_type repmat("TouchnoVision", 1, trialperblock)];
        % %             Block_num = [Block_num (repmat([j],1,trialperblock)+4*(i-1))];
        % %             Block = [Block (repmat((j-1+4*(i-1)),1,trialperblock))];
        % %         elseif BlockROrder(j) == 2
        % %             Block_type = [Block_type repmat("noTouchVision", 1, trialperblock)];
        % %             Block_num = [Block_num (repmat([j],1,trialperblock)+4*(i-1))];
        % %             Block = [Block (repmat((j-1+4*(i-1)),1,trialperblock))];
        % %         elseif BlockROrder(j) == 3
        % %             Block_type = [Block_type repmat("noTouchnoVision", 1, trialperblock)];
        % %             Block_num = [Block_num (repmat([j],1,trialperblock)+4*(i-1))];
        % %             Block = [Block (repmat((j-1+4*(i-1)),1,trialperblock))];
        % %         end
        % %     end
        % %     expTable = table(Block', TrialInBlock', Trial', Block_type', Block_num', Vision', Completed', Attempts', Skipped', TrialTime');
        % %     defTable = [defTable; expTable];
        end
    end
    defTable = table(Block, Block_num, TrialInBlock', Trial', def_trial_table{:,1}, def_trial_table{:,2}, def_trial_table{:,3}, (ITIs/1000), Speed_vectors, Completed', Attempts', Skipped', TrialTime');
    
    defTable.Properties.VariableNames = {'Block' 'Block_num', 'TrialInBlock' 'Trial' 'Block_type' 'Trial_type' 'Stimulation' 'ITIs' 'Speed' 'Completed' 'Attempts' 'Skipped' 'TrialTime'};
    writetable(defTable, outputpath, 'Delimiter', ',');
    disp(outputpath)
end


function [avg_speed, std_speed] = TakeSpeedfromLog(filepath)
    opts = delimitedTextImportOptions("NumVariables", 7);

    % Specify range and delimiter
    opts.DataLines = [4, Inf];
    opts.Delimiter = "\t";
    % Specify column names and types
    opts.VariableNames = ["Block", "Trial", "BlockType", "TrialType", "Stimulation", "Event_Name", "Time"];
    opts.VariableTypes = ["double", "double", "categorical", "double", "double", "categorical", "double"];
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    % Specify variable properties
    opts = setvaropts(opts, ["BlockType", "Event_Name"], "EmptyFieldRule", "auto");
    % Import the data
    table_training = readtable(filepath, opts);
    
    touches = find(table_training.Event_Name == "Index_visual");
    touches_to_delete = [];
    for i = 1:length(touches)
        if table_training.Event_Name(touches(i)) == table_training.Event_Name(touches(i)-1)
            touches_to_delete = [touches_to_delete; touches(i)]; 
        end
    end

    touches(ismember(touches, touches_to_delete)) = []; 
    GOs = touches - 1; 
    types = table_training.TrialType;
    types = types(touches); 

    %calculate overall differences and select just the ones in which the
    %participant moved. 
    diff = table_training.Time(touches) - table_training.Time(GOs); 
    diff_MOVE = diff(find(types == 1 | types == 3));

    %the distance that the ball should move onto is 0.750 Unity Units (even
    %though it stops earlier!).
    %If set with a speed of 1, it takes 750 ms to go from one point to the
    %other. 

    %therefore, we need to calculate the different speeds based on this
    %value :)

    %take a look at the "NO MOVE" condition, to see where the subject
    %touches the ball on average (in terms of distance). 
    avg_distance_NOMOVE = mean(diff(find(types == 2 | types == 4)));

    

    avg_speed = mean(rmoutliers(avg_distance_NOMOVE ./ diff_MOVE));
    std_speed = std(rmoutliers(avg_distance_NOMOVE ./ diff_MOVE)); 

end


