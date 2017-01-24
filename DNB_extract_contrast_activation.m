function DNB_extract_contrast_activation
% Calculates mean NIRS signal for the desired contrast vectors.
%
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO

%%%%%%%% Skip this step if already done and INFO.overwrite == no.
if exist(INFO.file.contrast_value_name);
    if strcmp(INFO.overwrite,'no')==1
        return
    end
end
%%%%%%%%

fprintf('## %s: running for subject %s ##\n',mfilename,INFO.dataselect.subjectnow);

load(INFO.file.pp_name);

iNoParMod=0; % Don't make barcharts of parametric modulations.
for iConvec=1:numel(INFO.stat.con)
    if strcmp(INFO.stat.con(iConvec).parmod,'yes')==0
        iNoParMod=iNoParMod+1;
        ConNoParMod(iNoParMod)=[iConvec];
    end
end

for iConvec=1:iNoParMod
    for iCond=1:size(SPM_nirs.nirs.beta,1)
        for iChan=1:size(SPM_nirs.nirs.beta,2)
            betacon{iConvec}(iCond,iChan)=INFO.stat.con(ConNoParMod(iConvec)).vec(iCond)*SPM_nirs.nirs.beta(iCond,iChan);
        end
    end
end

for iConvec=1:iNoParMod
    for iChan=1:size(SPM_nirs.nirs.beta,2)
        meanbetacon(iConvec,iChan)=sum(betacon{iConvec}(:,iChan));
    end
end

if strcmp(INFO.plots,'yes')==1
    figure;
    for iConvec=1:iNoParMod
        subplot(1,iNoParMod,iConvec); bar(meanbetacon(iConvec,:));
        title(['Contrast ',INFO.stat.con(ConNoParMod(iConvec)).name],'FontSize',10);
        xlabel('Channel','FontSize',12);
        if iConvec==1
            ylabel(['Contrast value (Subject: ',INFO.dataselect.subjectnow,', Task: ',INFO.dataselect.taskname,')'],'FontSize',10);
        end
        if strcmp(INFO.SCI.check,'yes')==1 
            set(gca,'XTickLabel',{INFO.SCI.sessremchannel{end}},'FontSize',1);
        else
            Chanidx=[1:INFO.conv.total_ch];
            set(gca,'XTickLabel',{Chanidx},'FontSize',1);
        end
    end
    
    % Save image
    plot_path = fullfile(INFO.file.dir,INFO.dataselect.subjectnow,INFO.dataselect.taskname,'Combined_Session','Figures',INFO.file.name,'Contrast value images');
    if ~exist(plot_path); mkdir(plot_path); end
    imgname=fullfile(plot_path,[INFO.file.name,'.',INFO.extension]);
    saveas(gcf,imgname);
    close(gcf);
end
% Save data for plotcontrastvalue_group
if strcmp(INFO.SCI.check,'yes')==1 
    meanbetacon((numel(INFO.stat.con)+1),:)=INFO.SCI.sessremchannel{end};
else
    meanbetacon((numel(INFO.stat.con)+1),:)=Chanidx;
end
save(INFO.file.contrast_value_name,'meanbetacon');


