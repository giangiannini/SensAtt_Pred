function D = spm_eeg_detect_eyeblinks_GIAN(S)
% Detects eyeblinks in spm continuous data file
% FORMAT  D = spm_eeg_detect_eyeblinks(S)
%
% S                    - input structure (optional)
% (optional) fields of S:
%          .stdthresh  - threshold to reject things that look like
%                         eye-blinks but probably aren't (default: 3)
%          .overwrite  - 1 - replace previous eybelink events (default)
%                        0 - append
% Output:
% D                 - MEEG object with added eyeblink events(also
%                     written on disk)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Laurence Hunt
% $Id: spm_eeg_detect_eyeblinks.m 4265 2011-03-28 13:18:31Z vladimir $

SVNrev = '$Rev: 4265 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','Eyeblink detect'); spm('Pointer','Watch');


%-Test for the presence of required Matlab toolbox - not needed anymore, 
% using ft_preproc_bandpassfilter instead
%--------------------------------------------------------------------------
if ~license('test','signal_toolbox')
%    error('Signal Processing Toolbox is required for eyeblink detection.');
end

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

D = spm_eeg_load(D);

if ~strcmp(D.type, 'continuous')
    error('Only continuous data is supported at the moment.');
end


% Get the indicex for EOG channel
%--------------------------------------------------------------------------
if ~(isfield(S, 'eogchan') && ~isempty(S.eogchan))
   eogchan = setdiff(D.eogchannels, D.badchannels);
   if length(eogchan)~=1
       [selection, ok]= listdlg('ListString', D.chanlabels(eogchan), 'SelectionMode', 'single' ,'Name', 'Select EOG channel' , 'ListSize', [400 300]);
       if ~ok
           return;
       end
       if ~isempty(eogchan)
           eogchan = eogchan(selection);
       else
           eogchan = selection;
       end
       S.eogchan = D.chanlabels(eogchan);
   end
elseif ~isnumeric(S.eogchan)
    eogchan = D.indchannel(S.eogchan);
else
    eogchan = S.eogchan;
end

try 
    stdthresh = S.stdthresh;
catch
    stdthresh = 4;
end

if ~isfield(S, 'overwrite')
    S.overwrite = spm_input('Overwrite previous?','+1','yes|no',[1 0], 1);
end

%% get EOG data
if length(eogchan)~=1
    error('More than one EOG channel - not currently supported')
end

eog_data = D(eogchan,:,:);

%% Convert back to fieldtrip and use fieldtrip function to determine hEOGs

prova = spm2fieldtrip(D);
cfg = []; 
cfg.channel = S.eogchan{1,1};
prova_EOG = ft_selectdata(cfg, prova); 

cfg = []; 
cfg.hpfreq = 2;

if strcmp(S.eogchan, 'hEOG')
    cfg.derivative = 'yes';
    cfg.hpfreq = 0.1;
end
cfg.hpfilter = 'yes';
cfg.hpfiltord = 3;
cfg.lpfilter = 'yes'; 
cfg.lpfreq = 10;
prova_EOG = ft_preprocessing(cfg, prova_EOG);

% cfg = []; 
% cfg.length = 10; 
% cfg.overlap = 0; 
% prova_EOG_epoched = ft_redefinetrial(cfg, prova_EOG); 

prova_EOG_epoched = prova_EOG;

%remove outliers


% cfg = []; 
% ft_databrowser(cfg, prova_EOG_epoched)

cfg            = [];
%cfg.trl        = trl;
cfg.continuous = 'yes';
cfg.artfctdef.zvalue.channel     = 1; %channel
cfg.artfctdef.zvalue.cutoff      = S.stdthresh; %threshold
cfg.artfctdef.zvalue.trlpadding  = 0;
cfg.artfctdef.zvalue.artpadding  = 2;
cfg.artfctdef.zvalue.fltpadding  = 0;
cfg.artfctdef.zvalue.bpfilter   = 'no';
cfg.artfctdef.zvalue.bpfilttype = 'but';
cfg.artfctdef.zvalue.bpfreq     = [2 15];
cfg.artfctdef.zvalue.bpfiltord  = 3;
cfg.artfctdef.zvalue.hilbert    = 'yes';
cfg.artfctdef.zvalue.interactive = 'no';
[cfg, artifact_eog] = ft_artifact_zvalue(cfg, prova_EOG_epoched);

cfg.artfctdef.reject          = 'nan';
prova_EOG_epoched = ft_rejectartifact(cfg, prova_EOG_epoched);

cfg = []; 
cfg.artfctdef.jump.channel = 1;
if strcmp(S.eogchan, 'vEOG')
    cfg.artfctdef.threshold.max = S.stdthresh; 
    cfg.artfctdef.threshold.min = -inf; 
    zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));
    prova_EOG_epoched.trial{1,1} = zscor_xnan(prova_EOG_epoched.trial{1,1}(:,:));
elseif strcmp(S.eogchan, 'hEOG')
    cfg.artfctdef.threshold.max = S.stdthresh; 
    cfg.artfctdef.threshold.min = -(S.stdthresh); 
end
cfg.artfctdef.threshold.bpfilter = 'no';
[cfg, artifact_eog] = ft_artifact_threshold(cfg, prova_EOG_epoched);

%% filter data at 1-15Hz (eyeblink duration typically 100-300ms) and demean
if strcmp(S.eogchan, 'vEOG')
    eog_filt = detrend(ft_preproc_bandpassfilter(eog_data, D.fsample, [2 10], 3, 'but'), 'constant');
    %eog_filt = zscore(eog_filt); 
elseif strcmp(S.eogchan, 'hEOG')
    eog_filt = detrend(ft_preproc_bandpassfilter(eog_data, D.fsample, [0.1 10], 3, 'but'), 'constant');
end

sd_eeg=(spm_percentile(eog_filt,85)-spm_percentile(eog_filt,15))/2; %robust estimate of standard deviation, suggested by Mark Woolrich

derivative_signal = prova_EOG_epoched.trial{1,1}(1,:);

if strcmp(S.eogchan, 'vEOG')
    eblength = round(D.fsample/6); %length of eyeblink(500 ms) in samples; %before it was
    %always positive, so no need to divide them
    spikes = []; 
    counter = 1;
    for i = 1:length(artifact_eog)-1
        [value pos] = max(eog_filt(artifact_eog(i,1):artifact_eog(i,2)));
        if artifact_eog(i,1)+pos > eblength && ... %not too close to beginning
                artifact_eog(i,1)+pos < length(derivative_signal)-eblength && ... %not too close to end
                    all(value >= eog_filt(artifact_eog(i,1)+pos-eblength:artifact_eog(i,1)+pos+eblength))
            spikes(counter) = artifact_eog(i,1)+pos;
            counter = counter+1;
        end
    end

    spikes = spikes';

elseif strcmp(S.eogchan, 'hEOG')
    eblength = round(D.fsample/2); %length of eyeblink(500 ms) in samples; %before it was
    %define the spikes as the midline between two consecutive drifts
    %first select only negative drifts
    for i = 1:length(artifact_eog)
        if mean(derivative_signal(artifact_eog(i,1):artifact_eog(i,2))) < 0
            artifact_eog(i,4) = 1; %negative peaks
        else
            artifact_eog(i,4) = 2; %positive peaks
        end
    end
    
    %now check whether in the surrounding 500 ms there is a positive drift.
    neg_peaks = artifact_eog(artifact_eog(:,4)==1,:);
    spikes = [0 0];
    counter = 1;
    for i = 1:length(neg_peaks)
        [value pos] = min(derivative_signal(neg_peaks(i,1):neg_peaks(i,2)));
        if neg_peaks(i,1)+pos > eblength && ...
                neg_peaks(i,1)+pos+eblength < length(derivative_signal) 
            [value_pos pos_pos] = max(derivative_signal(neg_peaks(i,1)+pos-eblength:neg_peaks(i,1)+pos+eblength));
            pos_pos_corr = pos_pos - eblength;
            if pos_pos_corr < 0
                start_drift = pos_pos_corr + neg_peaks(i,1);
                end_drift = pos+neg_peaks(i,1); 
            else
                start_drift = pos+neg_peaks(i,1); 
                end_drift = pos_pos_corr + neg_peaks(i,1);
            end

            if start_drift-eblength < 0 || end_drift + eblength > length(eog_filt)
                continue
            end

            % if eog_filt(round((start_drift+end_drift)/2)) > 0
            %     if ~all(eog_filt(start_drift:end_drift) > mean(eog_filt([start_drift-eblength:start_drift end_drift:end_drift+eblength])))
            %         continue
            %     end
            % elseif eog_filt(round((start_drift+end_drift)/2)) < 0
            %     if ~all(eog_filt(start_drift:end_drift) < mean(eog_filt([start_drift-eblength:start_drift end_drift:end_drift+eblength])))
            %         continue
            %     end
            % end


            %(all(eog_filt(start_drift:end_drift)>0) || all(eog_filt(start_drift:end_drift)<0))
            if value_pos > S.stdthresh && ...
                    all(abs(eog_filt(round((start_drift+end_drift)/2))) > abs(eog_filt(start_drift-eblength:start_drift))) && ...
                        all(abs(eog_filt(round((start_drift+end_drift)/2))) > abs(eog_filt(end_drift:end_drift+eblength))) && ...
                            end_drift - start_drift > 60
                pos_pos_corr = pos_pos - eblength;
                spikes(counter,1) = neg_peaks(i,1)+round(pos_pos_corr/2);

                %figure; plot(eog_filt(start_drift-200:end_drift+200)); hold on; xline(200); xline(end_drift-start_drift+200)

                %aaa = [];

                if pos_pos_corr < 0
                    spikes(counter,2) = 1; %left to right
                else
                    spikes(counter,2) = 2; %right to left
                end
                counter = counter+1; 
                % % % % figure; 
                % % % % plot(derivative_signal(neg_peaks(i,1)+pos-eblength:neg_peaks(i,1)+pos+eblength))
                % % % % hold on
                % % % % xline(eblength + pos_pos_corr/2)
            end
        end
    end
end

spikemat = zeros(eblength*2,length(spikes));
for i = 1:length(spikes)
    spikemat(:,i) = eog_filt(spikes(i)-eblength+1:spikes(i)+eblength);
end

% figure; 
% plot(spikemat)

% reject spikes whose peak is not within 1 s.d. of the mean (gets rid of most artefacts
%    etc. not removed by filtering):
mn_spike = mean(spikemat(eblength,:));
sd_spike = std(spikemat(eblength,:));
spikes((spikemat(eblength,:)>mn_spike+2*sd_spike | ...
       spikemat(eblength,:)<mn_spike-2*sd_spike),:) = [];
spikemat(:,find(spikemat(eblength,:)>mn_spike+2*sd_spike | ...
       spikemat(eblength,:)<mn_spike-2*sd_spike)) = [];

disp(['Number of putative eyeblinks detected: ' num2str(length(spikes))]);

% num_eb_per_min=60*length(spikes)/(D.time(end)-D.time(1));
% disp([num2str(num_eb_per_min) ' eye-blinks per minute '])
% if (num_eb_per_min<0.5)
%     error(['Only ' num2str(num_eb_per_min) ' eye-blinks per minute detected by algorithm. Try a lower threshold.'])
% end
% if (num_eb_per_min>60)
%     error(['As many as ' num2str(num_eb_per_min) ' eye-blinks per minute detected by algorithm. Try a higher threshold.'])
% end
% 
% plot
%----------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');
colormap(gray)
figure(Fgraph)
clf
subplot(2, 1 , 1)
plot(spikes,ones(length(spikes),1)*5*sd_eeg,'r.');
hold on;
plot(eog_filt);

subplot(2, 1 , 2)
hold on;
plot(spikemat);
%if (strcmp(S.eogchan, 'hEOG')) && S.direction_h_movement == 0
if (strcmp(S.eogchan, 'hEOG'))
    plot(mean(spikemat(:,spikes(:,2) == 1),2),'Color','k','LineWidth',4);
    plot(mean(spikemat(:,spikes(:,2) == 2),2),'Color','k','LineWidth',4);
else
    %plot(mean(spikemat,2),'Color','k','LineWidth',4);
end


%and save spikemat for future plotsss
save(strcat(string(S.eogchan{1}), "_spikemat.mat"), "spikemat")
%and save spikes for having artifacts positions
save(strcat(string(S.eogchan{1}), "_spikes.mat"), "spikes")

% Update the event structure
%----------------------------------------------------------------------

if S.direction_h_movement == 1 && strcmp(S.eogchan, 'hEOG') %left to right --> positive
    spikes = spikes(spikes(:,2) == 1);
elseif S.direction_h_movement == 2 && strcmp(S.eogchan, 'hEOG') %right to left --> negative
    spikes = spikes(spikes(:,2) == 2);
end
if ~isempty(spikes)  
    for n = 1:D.ntrials
        cspikes   = spikes(spikes(:,1)>(D.nsamples*(n-1)) & spikes(:,1)<(D.nsamples*n),1);
        ctime  = D.trialonset(n)+(cspikes - D.nsamples*(n-1)-1)/D.fsample;
        ctime  = num2cell(ctime);
        
        ev = events(D, n);
        
        if iscell(ev)
            ev = ev{1};
        end
        
        
        if ~isempty(ev) && S.overwrite
            ind1 = strmatch('artefact', {ev.type}, 'exact');
            if ~isempty(ind1)
                ind2 = strmatch('eyeblink', {ev(ind1).value}, 'exact');
                if ~isempty(ind2)
                    ev(ind1(ind2)) = [];
                end
            end
        end
        
        Nevents = numel(ev);
        for i=1:numel(ctime)
            ev(Nevents+i).type     = 'artefact';
            ev(Nevents+i).value    = 'eyeblink';
            ev(Nevents+i).duration = [];
            ev(Nevents+i).time     = ctime{i};
        end
        
        if ~isempty(ev)
            [tevent, I] = sort([ev.time]);
            ev = ev(I);
            D = events(D, n, ev);
        end
    end    
else
    warning(['No eye blinks events detected in the selected channel']);
end

%-Update history (not necessary; leave as call to spm_eeg_montage?) Save new evoked M/EEG dataset
%--------------------------------------------------------------------------
D = D.history(mfilename, S);

save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
spm('FigName','Eyeblink detect: done'); spm('Pointer','Arrow');




    
