function DNB_extract_activation
% Calculates mean (baseline corrected) NIRS signal for all conditions based on
% the model specification file (design matrix), plots the signal and saves the mean in a csv file.
%
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO

isess=INFO.counter.iSess;

%%%%%%%% Skip this step if already done and INFO.overwrite == no.
if exist(INFO.file.extract_name{isess});
    if strcmp(INFO.overwrite,'no')==1
        return
    end
end
%%%%%%%%

fprintf('## %s: running for subject %s ##\n',mfilename,INFO.dataselect.subjectnow);

baseline_cutoff = 0.6; % value under which total design matrix values should fall to be considered a baseline moment
fig_width=1000;
fig_height=600;

load(INFO.file.pp_name)
try     % Load in conditions-file
    load(INFO.file.conditions_adj_ds_name{isess});
catch
    try
        load(INFO.file.conditions_adj_name{isess});
    catch
        load(INFO.file.conditions_dir{isess});
    end
end

bf1find=strfind(SPM_nirs.xX.name,'bf(1)'); chanidx=cellfun(@isempty,bf1find); %exclude constants and hrf_type > 1
bf1findsess=strfind(SPM_nirs.xX.name, ['Sn(',num2str(isess),')']); chanidxsess=cellfun(@isempty,bf1findsess); %select session
chanidxfull=find(chanidx==0 & chanidxsess==0);
last_cond = size(chanidxfull,2)+1; % put implicit baseline at the end of the condition structure
condition(last_cond).name = 'implicit baseline';
sorted_baseline = sort(sum(SPM_nirs.xX.X(:,chanidxfull),2));
cutoffbase=sorted_baseline(round(numel(sorted_baseline)*0.1));
cutoffbase=min([cutoffbase 0.6]);
condition(last_cond).timepoints = find(sum(SPM_nirs.xX.X(:,chanidxfull),2)<baseline_cutoff);
condition(last_cond).avg_activation = mean(SPM_nirs.xX.KY(condition(last_cond).timepoints,:),1);
condition(last_cond).avg_corr_activation = zeros(size(condition(last_cond).avg_activation));
for iCond = 1:size(chanidxfull,2)
    condition(iCond).name = SPM_nirs.xX.name{chanidxfull(iCond)};
    sorted_regressor=sort(abs(SPM_nirs.xX.X(:,chanidxfull(iCond))));
    cutoff = sorted_regressor(round(numel(sorted_regressor)*0.9)); % take 90% value cutoff
    cutoff = max([cutoff 0.1]); % make sure the cutoff is at least 0.1
    condition(iCond).timepoints = find(abs(SPM_nirs.xX.X(:,chanidxfull(iCond))) > cutoff);
    condition(iCond).avg_activation = mean(SPM_nirs.xX.KY(condition(iCond).timepoints,:),1);
    condition(iCond).avg_corr_activation = condition(iCond).avg_activation-condition(last_cond).avg_activation;
end

if strcmp(INFO.plots,'yes')==1
    % plot activation for all conditions
    figure; set(gcf,'Position',[0 0 fig_width fig_height]);
    subplot(1,2,1);
    xdata = reshape([condition.avg_activation],numel(condition(1).avg_activation),numel(condition));
    bar(xdata');
    set(gca,'Xtick',1:size(xdata,2),'XTickLabel',{condition.name},'FontSize',1); xticklabel_rotate;
    if INFO.sessions==1
        ylabel([,INFO.dataselect.subjectnow,' ',INFO.dataselect.taskname,]);
    else
        ylabel([,INFO.dataselect.subjectnow,' ',INFO.dataselect.taskname,' Session: ',INFO.sess(isess).name,]);
    end
    title('Avg activation per condition');
    
    % plot baseline-corrected activation for all conditions
    subplot(1,2,2);
    xdata = reshape([condition(1:end-1).avg_corr_activation],numel(condition(1).avg_corr_activation),numel(condition)-1);
    bar(xdata');
    set(gca,'Xtick',1:size(xdata,2),'XTickLabel',{condition(1:end-1).name},'FontSize',1); xticklabel_rotate;
    title('Avg baseline-corrected activation');
        for iChan=1:size(INFO.SCI.sessremchannel{end},2);
            chanleg{iChan}=['Channel ',num2str(INFO.SCI.sessremchannel{end}(iChan)),];
        end
    legend(chanleg,'Location','NorthEastOutside')
    % save plot
    plot_path = fullfile(INFO.file.figures{isess},INFO.file.name,'Extract activation plots');
    if ~exist(plot_path); mkdir(plot_path); end
    imgname=fullfile(plot_path,[INFO.file.name,'.',INFO.extension]);
    saveas(gcf,imgname);
    close(gcf);
end

save(INFO.file.extract_name{isess},'condition');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write average activation to csv-file %
exactcsv_name=fullfile(INFO.file.extract_dir{isess},[,INFO.file.name,'.csv']);
ConChan{1,1}=['Subject']; ConChan{2,1}=[,INFO.dataselect.subjectnow,];

for iCon=1:size(condition,2)
        for iChan=1:size(INFO.SCI.sessremchannel{end},2)
            ConChan{1,(size(INFO.SCI.sessremchannel{end},2)*(iCon-1)+iChan)+1}=['Condition_',condition(iCon).name,'_Channel_',num2str(INFO.SCI.sessremchannel{end}(iChan)),];
            ConChan{2,(size(INFO.SCI.sessremchannel{end},2)*(iCon-1)+iChan)+1}=condition(iCon).avg_activation(iChan);
        end
end

ds=cell2dataset(ConChan);
export(ds,'file',exactcsv_name,'delimiter',',')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




