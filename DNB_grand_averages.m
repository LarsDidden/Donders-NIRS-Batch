function DNB_grand_averages
% The grand average script averages the NIRS signal for every trial in every condition, similar to event-related potentials in EEG.
%
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO

isess=INFO.counter.iSess;

%%%%%%%% Skip this step if already done and INFO.overwrite == no.
if exist(INFO.file.grand_avg_name{isess});
    if strcmp(INFO.overwrite,'no')==1
        return
    end
end
%%%%%%%%

fprintf('## %s: running for subject %s ##\n',mfilename,INFO.dataselect.subjectnow);
load(INFO.file.pp_name);

INFO.gavg.window_size_samples = INFO.gavg.window_size_secs * INFO.conv.downfs; % window size in samples

% load onsets file
load(INFO.file.conditions_adj_ds_name{isess});

%% get grand average data
for iChannel = 1:INFO.conv.total_ch
    for iCond = 1:numel(onsets)
        condition_gavg(iCond,iChannel).gavg = NaN(numel(onsets{iCond}),INFO.gavg.window_size_samples);
        condition_gavg(iCond,iChannel).baseline = NaN(numel(onsets{iCond}),INFO.gavg.baseline_period_calc * INFO.conv.downfs);
        for iTrial = 1:numel(onsets{iCond})
            start_idx = onsets{iCond}(iTrial)-(INFO.gavg.baseline_period_plot*INFO.conv.downfs); % average window, defined relative to stimulus onset
            end_idx =   start_idx + INFO.gavg.window_size_samples-1; % find last sample for average
            if end_idx <= size(SPM_nirs.xX.KYsep{isess},1)
                condition_gavg(iCond,iChannel).gavg(iTrial,:) = SPM_nirs.xX.KYsep{isess}(start_idx:end_idx,iChannel)';
            else
                condition_gavg(iCond,iChannel).gavg(iTrial,1:size(SPM_nirs.xX.KYsep{isess},1)-start_idx+1) = SPM_nirs.xX.KYsep{isess}(start_idx:end,iChannel)';
            end
            bl_start_idx = onsets{iCond}(iTrial)-(INFO.gavg.baseline_period_calc*INFO.conv.downfs); % average window for baseline, defined relative to stimulus onset
            bl_end_idx = onsets{iCond}(iTrial)-1;
            if all(size(SPM_nirs.xX.KYsep{isess},1)>[bl_start_idx,bl_end_idx]);
                condition_gavg(iCond,iChannel).baseline(iTrial,:) = SPM_nirs.xX.KYsep{isess}(bl_start_idx:bl_end_idx,iChannel)';
            else
                fprintf('%s Warning: trial indices fall outside of sample range (subject %s)\n',mfilename,INFO.dataselect.subjectnow);
            end
            
            % correct for baseline
            condition_gavg(iCond,iChannel).gavg(iTrial,:) = condition_gavg(iCond,iChannel).gavg(iTrial,:) - mean(condition_gavg(iCond,iChannel).baseline(iTrial,:));
        end
        
    end
end
%% plot the grand averages

% plot per condition, over channels
for iCondition = 1:numel(onsets)
    cf=figure;
    x_axis_time = -INFO.gavg.baseline_period_plot:1/INFO.conv.downfs:INFO.gavg.window_size_secs-INFO.gavg.baseline_period_plot-1/INFO.conv.downfs; % time information on x-axis
    for iChannel = 1:INFO.conv.total_ch
        curr_data{iChannel} = nanmean(condition_gavg(iCondition,iChannel).gavg,1);
        plot(x_axis_time,curr_data{iChannel},'Color',[iChannel*(1/INFO.conv.total_ch) 0 0],'LineWidth',3);
        hold on
    end
    for iChan = 1:INFO.conv.total_ch
        if strcmp(INFO.SCI.check,'yes')==1 
            legendchan{iChan}= ['Channel ',num2str(INFO.SCI.sessremchannel{end}(iChan)),];
        else
            legendchan{iChan}= ['Channel ',num2str(iChan),];
        end
    end
    legend(legendchan,'Location','EastOutside');
    xlabel('Time relative to stimulus onset (s)','FontSize',20);
    ylabel(INFO.model.hb,'FontSize',20);
    
    % save the figures
    grand_avg_dir    = fullfile(INFO.file.figures{isess},INFO.file.name,'Grand Average plots'); %grand average plot map
    if ~exist(grand_avg_dir); mkdir(grand_avg_dir); end
    if INFO.sessions==1
        title(['Condition: ',names(iCondition),],'FontSize',20);
    else
        title(['Session: ',INFO.sess(isess).name,' Condition: ',names(iCondition),],'FontSize',20);
    end
    INFO.file.grand_avg_plot = fullfile(grand_avg_dir, ['grand_avg_',INFO.dataselect.subjectnow,'_',INFO.file.name,'_cond_',names{iCondition},'.',INFO.extension]);
    saveas(cf,INFO.file.grand_avg_plot);
    close(cf);
end

% plot per channel, over conditions
for iChannel = 1:INFO.conv.total_ch
    cf=figure;
    x_axis_time = -INFO.gavg.baseline_period_plot:1/INFO.conv.downfs:INFO.gavg.window_size_secs-INFO.gavg.baseline_period_plot-1/INFO.conv.downfs; % time information on x-axis
    for iCondition = 1:numel(onsets)
        curr_data{iCondition} = nanmean(condition_gavg(iCondition,iChannel).gavg,1);
        plot(x_axis_time,curr_data{iCondition},'Color',[iCondition*(1/numel(onsets)) 0 0],'LineWidth',3);
        hold on
    end
    legend(names,'Location','EastOutside');
    xlabel('Time relative to stimulus onset (s)','FontSize',20);
    ylabel(INFO.model.hb,'FontSize',20);
    
    % save the figures
    grand_avg_dir    = fullfile(INFO.file.figures{isess},INFO.file.name,'Grand Average plots'); %grand average plot map
    if ~exist(grand_avg_dir); mkdir(grand_avg_dir); end
    if INFO.sessions==1
        if strcmp(INFO.SCI.check,'yes')==1 
            title(['Channel: ',num2str(INFO.SCI.sessremchannel{end}(iChannel))],'FontSize',20);
            INFO.file.grand_avg_plot = fullfile(grand_avg_dir, ['grand_avg_',INFO.dataselect.subjectnow,'_',INFO.file.name,'_chan',num2str(INFO.SCI.sessremchannel{end}(iChannel)),'.',INFO.extension]);
        else
            title(['Channel: ',num2str(iChannel)],'FontSize',20);
            INFO.file.grand_avg_plot = fullfile(grand_avg_dir, ['grand_avg_',INFO.dataselect.subjectnow,'_',INFO.file.name,'_chan',num2str(iChannel),'.',INFO.extension]);
        end
    else
        if strcmp(INFO.SCI.check,'yes')==1
            title(['Session: ',INFO.sess(isess).name,' Channel: ',num2str(INFO.SCI.sessremchannel{end}(iChannel))],'FontSize',20);
            INFO.file.grand_avg_plot = fullfile(grand_avg_dir, ['grand_avg_',INFO.dataselect.subjectnow,'_',INFO.file.name,'_chan',num2str(INFO.SCI.sessremchannel{end}(iChannel)),'.',INFO.extension]);
        else
            title(['Session: ',INFO.sess(isess).name,' Channel: ',num2str(iChannel)],'FontSize',20);
            INFO.file.grand_avg_plot = fullfile(grand_avg_dir, ['grand_avg_',INFO.dataselect.subjectnow,'_',INFO.file.name,'_chan',num2str(iChannel),'.',INFO.extension]);
        end
    end
    
    saveas(cf,INFO.file.grand_avg_plot);
    close(cf);
end

% save grand averages to disk
save(INFO.file.grand_avg_name{isess},'condition_gavg');

