function DNB_extract_contrast_activation_group
% Calculates mean NIRS signal for the desired contrast vectors.
%
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO

fprintf('## %s: running for group %s ##\n',mfilename,INFO.dataselect.analysisname);

for iSubj=1:numel(INFO.sec.subjects)
    if INFO.sessions>1
        INFO.file.contrast_value_dir= fullfile(INFO.file.dir,INFO.sec.subjects{iSubj},INFO.dataselect.taskname,'Combined_Session','Contrast_value_files');
    else
        INFO.file.contrast_value_dir= fullfile(INFO.file.dir,INFO.sec.subjects{iSubj},INFO.dataselect.taskname,'Contrast_value_files');
    end
    sec_name=[INFO.sec.subjects{iSubj},'_',INFO.dataselect.taskname,'_',INFO.file.name,'_unit_',INFO.model.units,'.mat'];
    filename{iSubj}=fullfile(INFO.file.contrast_value_dir,sec_name);
end

iNoParMod=0; % Don't make barcharts of parametric modulations.
for iConvec=1:numel(INFO.stat.con)
    if strcmp(INFO.stat.con(iConvec).parmod,'yes')==0
        iNoParMod=iNoParMod+1;
        ConNoParMod(iNoParMod)=[iConvec];
    end
end

meanbetacon_group=zeros(iNoParMod,INFO.conv.total_ch); %pre-allocation
SCI_channel_count=zeros(iNoParMod,INFO.conv.total_ch);

for iSubj=1:numel(INFO.sec.subjects)    %Adding up meanbetacons of all subjects; taking SCI excluded channels into account
    load(filename{iSubj});
    for iConv=1:iNoParMod
        for iChan=1:size(meanbetacon,2)
            ChanLoc=numel(INFO.stat.con)+1; %Row where channelnrs are stored
            meanbetacon_group(iConv,meanbetacon(ChanLoc,iChan))=meanbetacon_group(iConv,meanbetacon(ChanLoc,iChan))+meanbetacon(iConv,iChan);
            SCI_channel_count(iConv,meanbetacon(ChanLoc,iChan))= SCI_channel_count(iConv,meanbetacon(ChanLoc,iChan))+1;
        end
    end
end

for iConv=1:iNoParMod
    for iChan=1:INFO.conv.total_ch      %Calculating mean of meanbetacons; taking SCI excluded channels into account
        meanbetacon_group(iConv,iChan)=meanbetacon_group(iConv,iChan)/SCI_channel_count(iConv,iChan);
    end
end

if strcmp(INFO.plots,'yes')==1
    figure;
    for iConvec=1:iNoParMod
        if INFO.groupplots==0
            subplot(1,iNoParMod,iConvec);
        end
        bar(meanbetacon_group(iConvec,:));
        title([,INFO.stat.con(ConNoParMod(iConvec)).name],'FontSize',10);
        xlabel('Channel','FontSize',12);  ylabel(['Contrast value (Experiment: ',INFO.dataselect.analysisname,', Task: ',INFO.dataselect.taskname,')'],'FontSize',12);
        if INFO.groupplots==1 % Save image separately
            imgname=fullfile(INFO.file.sec_config,[INFO.dataselect.taskname,'_',INFO.file.name,'_',INFO.stat.con(ConNoParMod(iConvec)).fname,'.',INFO.extension]);
            saveas(gcf,imgname);
            close(gcf);
        end
    end
    if INFO.groupplots==0 % Save image in one image
        imgname=fullfile(INFO.file.sec_config,[INFO.dataselect.taskname,'_',INFO.file.name,'.',INFO.extension]);
        saveas(gcf,imgname);
        close(gcf);
    end
end


filename=[INFO.dataselect.taskname,'_',INFO.file.name];;
INFO.file.sec_conname=fullfile(INFO.file.sec_confile,filename);
save(INFO.file.sec_conname,'meanbetacon_group');

