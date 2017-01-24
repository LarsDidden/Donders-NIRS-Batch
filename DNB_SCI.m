function DNB_SCI
% Run scalp coupling index as a quality measure over de dataset. Exclude
% channels/subjects that do not match the quality demands.
%
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO

fprintf('## %s: running for subject %s ##\n',mfilename,INFO.dataselect.subjectnow);

INFO.SCI.remchannel=[]; INFO.SCI.sessremchannel=[];
for isess=1:INFO.sessions
    
    bandpass_range = [0.5 2.5]; % typical frequency range for heart rate that excludes low-frequency NIRS activity.
    
    rawOD = oxysoft2matlab(INFO.file.Sess(isess).rawdata_dir{INFO.counter.iSubj},'rawOD');
    for iChan=1:numel(rawOD.label)
        corr_name{iChan} = ['Channel ',num2str(iChan),' HbO-HbR'];
    end
    ctr_channel = 1;
    
    for iChannel = 1:2:numel(rawOD.ODlabel)
        
        % filter the signal between 0.5-2.5 Hz
        dat = rawOD.OD(:,iChannel:iChannel+1)'; % fieldtrip wants channels x time
        filtered_dat = ft_preproc_bandpassfilter(dat,rawOD.Fs,bandpass_range)'; % store as timexchannels again
        % calculate SCI in sliding window average
        window_size = 5; % sliding window size, in seconds
        window_nsamples = window_size * rawOD.Fs;
        ctr = 1;
        for iLoc = 1:window_nsamples:size(dat,2)-window_nsamples+1
            curr_correlations = corr(filtered_dat(iLoc:iLoc+window_nsamples-1,:));
            SCI_data.channel_SCI(ctr_channel,ctr) = curr_correlations(2);
            ctr = ctr + 1;
        end
        % also get correlation over the whole timecourse
        curr_correlation = corr(filtered_dat);
        SCI_data.channel_SCI_total(ctr_channel) = curr_correlation(2);
        ctr_channel = ctr_channel+1;
    end
    % plot correlations
    if strcmp(INFO.plots,'yes')==1
        figure;plot(SCI_data.channel_SCI');
        legend(corr_name,'Location','SouthOutside');
        SCI_figdir=fullfile(INFO.file.figures{isess},INFO.file.name,'SCI figures');
        if ~exist(SCI_figdir); mkdir(SCI_figdir); end
        SCI_figname=fullfile(SCI_figdir,INFO.file.name);
        saveas(gcf,SCI_figname,INFO.extension);
        close(gcf);
    end
    
    % Select uncorrelated channel data
    for iChan=1:numel(rawOD.label)
        if SCI_data.channel_SCI_total(iChan)<INFO.SCI.corrrequired
            INFO.SCI.removech(iChan)=1;
        else
            INFO.SCI.removech(iChan)=0;
        end
    end
    % Make variable with approved channels
    delchan=find(INFO.SCI.removech==1);
    INFO.SCI.remchannel{isess}=[1:numel(rawOD.label)];
    INFO.SCI.remchannel{isess}(delchan)=[];
    remchancomb=1;
    if isess>2 %If using sessions, make sure each session has the same channels. So if a channel needs to be removed in one session, it is also taken out of other sessions.
        equalchanses=isequal(INFO.SCI.remchannel{isess},INFO.SCI.sessremchannel{isess-1});
    elseif isess>1
        equalchanses=isequal(INFO.SCI.remchannel{isess},INFO.SCI.remchannel{isess-1});
    else
        equalchanses=2;
    end
    if equalchanses==0
        for iChan=1:size(INFO.SCI.remchannel{isess},2)
            if isess>2
                NoRemidx=find(INFO.SCI.sessremchannel{isess-1}==INFO.SCI.remchannel{isess}(iChan));
            else
                NoRemidx=find(INFO.SCI.remchannel{isess-1}==INFO.SCI.remchannel{isess}(iChan));
            end
            if ~isempty(NoRemidx)
                if isess>2
                    INFO.SCI.sessremchannel{isess}(remchancomb)=INFO.SCI.sessremchannel{isess-1}(NoRemidx);
                else
                    INFO.SCI.sessremchannel{isess}(remchancomb)=INFO.SCI.remchannel{isess-1}(NoRemidx);
                end
                remchancomb=remchancomb+1;
            end
        end
    elseif equalchanses==1;
        INFO.SCI.sessremchannel{isess-1}=INFO.SCI.remchannel{isess};
    elseif equalchanses==2;
        INFO.SCI.sessremchannel{isess}=INFO.SCI.remchannel{isess};
    end
end

for isess=1:INFO.sessions
    %%%%%%%%%%%%%%%%%%%%%%
    % Remove unreliable datachannels based on SCI check.
    load(INFO.file.conv_name{isess});
    nirs_data.nchoriginal=nirs_data.nch;
    if isequal(size(INFO.SCI.sessremchannel{end},2),nirs_data.nch)==0
        nirs_data.oxyDataoriginal=nirs_data.oxyData; nirs_data.oxyData=[];
        nirs_data.dxyDataoriginal=nirs_data.dxyData; nirs_data.dxyData=[];
        nirs_data.oxyData=nirs_data.oxyDataoriginal(:,INFO.SCI.sessremchannel{end});
        nirs_data.dxyData=nirs_data.dxyDataoriginal(:,INFO.SCI.sessremchannel{end});
        nirs_data.nch=size(INFO.SCI.sessremchannel{end},2);
        save(INFO.file.conv_name{isess}, 'nirs_data');
    end
    INFO.conv.total_ch=nirs_data.nch;
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Report
    if isequal(nirs_data.nchoriginal,nirs_data.nch)
        if isempty(INFO.report.SCI)==1;
            INFO.report.SCI=[,INFO.dataselect.subjectnow,', Session: ',num2str(INFO.counter.iSess),', no channels deleted.'];
        else
            new_SCI=[,INFO.dataselect.subjectnow,', Session: ',num2str(INFO.counter.iSess),', no channels deleted.'];
            INFO.report.SCI=[INFO.report.SCI, {new_SCI}];
        end
    elseif nirs_data.nch<INFO.SCI.chanrequired==1
        if isempty(INFO.report.SCI)==1;
            INFO.report.SCI=[,INFO.dataselect.subjectnow,', Session: ',num2str(INFO.counter.iSess),', not analysed because less than ',num2str(INFO.SCI.chanrequired),' channels were reliable based on the SCI check.'];
        else
            new_SCI=[,INFO.dataselect.subjectnow,', Session: ',num2str(INFO.counter.iSess),', not analysed because less than ',num2str(INFO.SCI.chanrequired),' channels were reliable based on the SCI check.'];
            INFO.report.SCI=[INFO.report.SCI, {new_SCI}];
        end
    else
        if isempty(INFO.report.SCI)==1;
            INFO.report.SCI=[,INFO.dataselect.subjectnow,', Session: ',num2str(INFO.counter.iSess),', channels accepted after SCI check: ',num2str(INFO.SCI.sessremchannel{end}),];
        else
            new_SCI=[,INFO.dataselect.subjectnow,', Session: ',num2str(INFO.counter.iSess),', channels accepted after SCI check: ',num2str(INFO.SCI.sessremchannel{end}),];
            INFO.report.SCI=[INFO.report.SCI, {new_SCI}];
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%
% Stop analyses if there are more removable channels than INFO.SCI.chanrequired
if nirs_data.nch<INFO.SCI.chanrequired
    error(['More than halve of all channels (',num2str(size(INFO.SCI.sessremchannel{end},2)),') has unreliable data, based on a correlation check between HbO and HbR data. Therefore, analysis is aborted.']);
end

