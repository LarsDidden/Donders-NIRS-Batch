function DNB_runMARA
% Run MARA over the dataset
%
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO

fprintf('## %s: running for subject %s ##\n',mfilename,INFO.dataselect.subjectnow);

isess=INFO.counter.iSess;

%%%%%%%% Skip this step if already done and INFO.overwrite == no.
if exist(INFO.file.MARA_name{isess});
    if strcmp(INFO.overwrite,'no')==1
        return
    end
end
%%%%%%%%

try     % Load in NIRS-dataset
    fname_nirs=INFO.file.conv_ds_name{isess};
    load(fname_nirs);
catch
    fname_nirs=INFO.file.conv_name{isess};
    load(fname_nirs);
end
fs = INFO.conv.fs;
L  = round(INFO.MARA.L*fs);     % convert L from seconds to samples
T=INFO.MARA.T;
alpha = INFO.MARA.alpha;
if strcmp(INFO.filt.MARA,'interactive') == 1       % Give users the option to choose whether to use MARA or not using plots
    if strcmp(INFO.model.hb,'HbO') == 1
        x  = nirs_data.oxyData(:,1);
    elseif strcmp(INFO.model.hb,'HbR') == 1
        x = nirs_data.dxyData(:,1);
    end
    % Calculate moving standard deviation
    [y1] = spm_fnirs_MovStd(x,L);
    
    % Calculate threshold
    indx_n = find(isnan(y1) == 1); y2=y1; y2(indx_n) = [];
    mstd_y = mean(y2);
    Thres=mstd_y*T;
    
    if max(y1) < Thres || min(y1) > Thres
        figure; plot(x);
        disp(['No movement artifact was found in the first channel while using the threshold ',num2str(T),'']);
        inthresh1=input('Do you want to change the threshold? y for yes: ','s');
        if strcmp(inthresh1,'y')==1
            error('Change the threshold in the INFO file');
        end
        inthresh2=input('Do you want to continue running MARA? n for no: ','s');
        if strcmp(inthresh2,'n')==1
            return
        end
        close(gcf);
    end
    
    [y3,noMARA] = DNB_MARA(x,fs,T,L,alpha);            % present MARA of first channel and then give the choice of running MARA or not.
    figure;
    subplot(2,1,1); plot(x,'b'); title('NIRS oxydata channel 1 without MARA'); subplot(2,1,2); plot(x,'b'); hold on; plot(y3,'r'); title('NIRS oxydata channel 1 with MARA.');
    inMARA=input('Do you want to use the MARA script on this data? y for yes, n for no, s if you want to stop running to change the settings: ','s');
    close(gcf);
else
    inMARA='y';
end


if strcmp(inMARA,'y')==1
    % MARA script
    tic
    report_MARA=[];
    h_wait = waitbar(0, 'Running MARA over all channels, please wait...');
    if strcmp(INFO.model.hb,'HbO') == 1
        ox_MARA=nirs_data.oxyData; nirs_data.oxyData=[];
        for iChan=1:size(ox_MARA,2)
            waitbar(iChan/size(ox_MARA,2), h_wait);
            [nirs_data.oxyData(:,iChan),noMARA] = DNB_MARA(ox_MARA(:,iChan),fs,T,L,alpha);
            if noMARA==1
                disp(['NIRS oxychannel ',num2str(INFO.SCI.sessremchannel{end}(iChan)),' did not change after MARA.'])
                if strcmp(INFO.plots,'yes');
                    %Plotting the signal before and after filtering for every channel
                    DNB_plotchannel(ox_MARA(:,iChan),nirs_data.oxyData(:,iChan),'noMARA',iChan)
                end
            else
                report_MARA=[report_MARA, iChan];
                if strcmp(INFO.plots,'yes');
                    %Plotting the signal before and after filtering for every channel
                    DNB_plotchannel(ox_MARA(:,iChan),nirs_data.oxyData(:,iChan),'MARA',iChan)
                end
            end
        end
        
    elseif strcmp(INFO.model.hb,'HbR') == 1
        dx_MARA=nirs_data.dxyData; nirs_data.dxyData=[];
        for iChan=1:size(dx_MARA,2)
            waitbar(iChan/size(dx_MARA,2), h_wait);
            [nirs_data.dxyData(:,iChan),noMARA] = DNB_MARA(dx_MARA(:,iChan),fs,T,L,alpha);
            if noMARA==1
                disp(['NIRS oxychannel ',num2str(INFO.SCI.sessremchannel{end}(iChan)),' did not change after MARA.'])
                if strcmp(INFO.plots,'yes');
                    %Plotting the signal before and after filtering for every channel
                    DNB_plotchannel(dx_MARA(:,iChan),nirs_data.dxyData(:,iChan),'noMARA',iChan)
                end
            else
                report_MARA=[report_MARA, iChan];
                if strcmp(INFO.plots,'yes');
                    %Plotting the signal before and after filtering for every channel
                    DNB_plotchannel(dx_MARA(:,iChan),nirs_data.dxyData(:,iChan),'MARA',iChan)
                end
            end
        end
    end
    
    close(h_wait);
    
    %%%%%%%%%%%%%
    % MARA report
    if isempty(report_MARA)==1
        if isempty(INFO.report.MARA)==1;
            INFO.report.MARA=[,INFO.dataselect.subjectnow,', Session: ',num2str(INFO.counter.iSess),', no channels changed with MARA.'];
        else
            new_MARA=[,INFO.dataselect.subjectnow,', Session: ',num2str(INFO.counter.iSess),', no channels changed with MARA.'];
            INFO.report.MARA=[INFO.report.MARA, {new_MARA}];
        end
    else
        if isempty(INFO.report.MARA)==1;
            INFO.report.MARA=[,INFO.dataselect.subjectnow,', Session: ',num2str(INFO.counter.iSess),', channels changed with MARA: ',num2str(report_MARA),];
        else
            new_MARA=[,INFO.dataselect.subjectnow,', Session: ',num2str(INFO.counter.iSess),', channels changed with MARA: ',num2str(report_MARA),];
            INFO.report.MARA=[INFO.report.MARA, {new_MARA}];
        end
    end
    %%%%%%%%%%%%%
    
    % continue without running MARA
elseif strcmp(inMARA,'s')==1
    error('Change MARA settings in the INFO file.');
end


save(INFO.file.MARA_name{isess},'nirs_data');


