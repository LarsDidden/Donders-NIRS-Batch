function DNB_grand_averages_group
% The grand average group script averages the NIRS signal for every trial in every condition, for a group of subjects. Similar to event-related potentials in EEG.
%
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO

warning off MATLAB:legend:IgnoringExtraEntries

fprintf('## %s: running for group %s ##\n',mfilename,INFO.dataselect.analysisname);
% load onsets file
load(INFO.file.conditions_adj_ds_name);
load(INFO.file.stat_name);

iNoParMod=0; % Don't make grand averages plots of parametric modulations.
for iConvec=1:numel(INFO.stat.con)
    if strcmp(INFO.stat.con(iConvec).parmod,'yes')==0
        iNoParMod=iNoParMod+1;
        ConNoParMod(iNoParMod)=[iConvec];
    end
end

% Adjust names to exclude parmod names (parmod conditions are not used in these plots)
parvec=strfind(SPM_nirs.xX.name,'^');
cntr=1;
for iPVI=1:size(parvec,2)
    if isempty(parvec{iPVI})==1
        parvecidx(cntr)=iPVI;
        cntr=cntr+1;
    end
end
namesold=SPM_nirs.xX.name(parvecidx);

% Find session related vectors in contrast_vector
for isess=1:INFO.sessions
    sessvec{isess}=strfind(namesold,['Sn(',num2str(isess),')']);
    cntr=1;
    for iSVI=1:size(sessvec{isess},2)
        if isempty(sessvec{isess}{iSVI})==0
            sessvecidx{isess}(cntr)=iSVI;
            cntr=cntr+1;
        end
    end
    
    % Load grand avgs names
    for iSubj=1:numel(INFO.sec.subjects)
        if INFO.sessions>1
            INFO.file.grand_avg_dir = fullfile(INFO.file.dir,INFO.sec.subjects{iSubj},INFO.dataselect.taskname,INFO.sess(isess).name,'Grand_average_files');
        else
            INFO.file.grand_avg_dir = fullfile(INFO.file.dir,INFO.sec.subjects{iSubj},INFO.dataselect.taskname,'Grand_average_files');
        end
        sec_ga_name=[INFO.sec.subjects{iSubj},'_',INFO.dataselect.taskname,'_',INFO.file.name,'_unit_',INFO.model.units,'.mat'];
        ga_filename{iSubj}=fullfile(INFO.file.grand_avg_dir,sec_ga_name);
    end
    %% Take mean of group
    
    for iSubj = 1:numel(ga_filename);
        dat = load(ga_filename{iSubj});
        
        for iCond = 1:size(dat.condition_gavg,1);
            for iChan = 1:size(dat.condition_gavg,2);
                if size(dat.condition_gavg(iCond,iChan).gavg,1)<=1              % Calculate no mean or SEM if there is no or just one onset.
                    dat.condition_gavg(iCond,iChan).gavg(:)=[NaN];
                end
            end
        end
        
        nanidx{iSubj}=[];
        for iCond = 1:size(dat.condition_gavg,1);
            if isnan(dat.condition_gavg(iCond,1).gavg);
                nanidx{iSubj}=[nanidx{iSubj}, iCond];
            end
        end
        
        Condnr{iSubj}=[1:size(dat.condition_gavg,1)]; Condnr{iSubj}(nanidx{iSubj})=[];
        for iCond = 1:size(Condnr{iSubj},2);
            for iChan = 1:size(dat.condition_gavg,2);
                MEAN{(size(dat.condition_gavg,1)*(isess-1))+Condnr{iSubj}(iCond),iChan}(iSubj,:) = nanmean(dat.condition_gavg(Condnr{iSubj}(iCond),iChan).gavg,1);
            end
        end
    end
    
    for iSubj=1:numel(ga_filename);
        if ~isempty(nanidx{iSubj})
            for iChan=1:size(dat.condition_gavg,2)
                MEAN{(isess-1)*size(dat.condition_gavg,1)+nanidx{iSubj},iChan}(iSubj,:) = NaN;
            end
        end
    end
end

%Generate significant condition-names
for iConv=1:iNoParMod
    if numel(INFO.stat.con(ConNoParMod(iConv)).positive)>1
        signamehigh{iConv}=INFO.stat.con(ConNoParMod(iConv)).positive{1};
        for iPos=2:numel(INFO.stat.con(ConNoParMod(iConv)).positive)
            signamehigh{iConv}=[signamehigh{iConv},'+',INFO.stat.con(ConNoParMod(iConv)).positive{iPos}];
        end
    else
        signamehigh{iConv}=INFO.stat.con(ConNoParMod(iConv)).positive{1};
    end
    if numel(INFO.stat.con(ConNoParMod(iConv)).negative)>1
        signamelow{iConv}=INFO.stat.con(ConNoParMod(iConv)).negative{1};
        for iNeg=2:numel(INFO.stat.con(ConNoParMod(iConv)).negative)
            signamelow{iConv}=[signamelow{iConv},'+',INFO.stat.con(ConNoParMod(iConv)).negative{iNeg}];
        end
    else
        signamelow{iConv}=INFO.stat.con(ConNoParMod(iConv)).negative{1};
    end
end


if strcmp(INFO.plots,'yes')==1
    %Linecolors for boundedline-plot
    colorline{1}=['b'];colorline{2}=['r'];colorline{3}=['g'];colorline{4}=['m'];colorline{5}=['c'];colorline{6}=['k'];
    % plot per channel.
    for iChan=1:size(dat.condition_gavg,2);
        if INFO.groupplots==0
            cf=figure('Position',[100, 100, 1200, 1000]);
        end
        x_axis_time = -INFO.gavg.baseline_period_plot:1/INFO.conv.downfs:INFO.gavg.window_size_secs-INFO.gavg.baseline_period_plot-1/INFO.conv.downfs; % time information on x-axis
        for iCond = 1:iNoParMod;
            if INFO.groupplots==1
                cf=figure('Position',[100, 100, 1200, 1000]);
            end
            un_cv=unique(INFO.stat.con(ConNoParMod(iCond)).vec); x=find(un_cv==0); un_cv(x)=[];
            for iCV=1:numel(un_cv)
                for isess=1:INFO.sessions
                    idx{isess,iCV}=sessvecidx{isess}(find(INFO.stat.con(ConNoParMod(iCond)).vecold(sessvecidx{isess})==un_cv(iCV)));
                    nrsubj=numel(INFO.sec.subjects);
                    nrrows=nrsubj*size(idx,2);
                    for i=1:numel(idx{isess,iCV}) %Add all cond with the same contrast value in con_vector to one variable
                        MEANCONV{iCV,iCond}((((isess-1)*nrrows)+(i*nrsubj-(nrsubj-1)):(((isess-1)*nrrows)+i*nrsubj)),:)=MEAN{idx{isess,iCV}(i),iChan};
                    end
                end
                % Take mean and SEM of all subjecs, for all contrasts
                sig_mem{iCV,iCond} = nanmean(MEANCONV{iCV,iCond},1);
                dataline_sem{iCV,iCond} = nanstd(MEANCONV{iCV,iCond},1)/sqrt(size(MEANCONV{iCV,iCond},1));
                if INFO.groupplots==0 %Subplots for each contrast vector
                    subplot(2,iNoParMod,iCond);
                elseif INFO.groupplots==1 %Plots for each contrast vector
                    subplot(2,1,1);
                else
                    error('INFO.groupplots should be either 0 or 1.');
                end
                [l,p] = boundedline(x_axis_time,sig_mem{iCV,iCond},dataline_sem{iCV,iCond}, colorline{iCV},'alpha');
                if INFO.groupplots==0
                    if iCond==1 %Make sure the title is only on the most left subplot
                        title(['Experiment: ',INFO.dataselect.analysisname,' (task: ',INFO.dataselect.taskname,'; channel: ',num2str(iChan),')'],'FontSize',20);
                        ylabel(INFO.model.hb,'FontSize',10);
                    end
                else
                    title(['Experiment: ',INFO.dataselect.analysisname,' (task: ',INFO.dataselect.taskname,'; channel: ',num2str(iChan),')'],'FontSize',20);
                    ylabel(INFO.model.hb,'FontSize',10);
                end
                if INFO.groupplots==0 %Make sure the legend is only on the most left subplot
                    legendcond{1}=['SEM: ',signamelow{iCond}]; legendcond{2}=['Mean: ',signamelow{iCond}]; legendcond{3}=['SEM: ',signamehigh{iCond}]; legendcond{4}=['Mean: ',signamehigh{iCond}];
                    if iCond==1
                        legend(legendcond,'Location','NorthEast');
                    end
                else
                    legendcond{1}=['SEM: ',signamelow{iCond}]; legendcond{2}=['Mean: ',signamelow{iCond}]; legendcond{3}=['SEM: ',signamehigh{iCond}]; legendcond{4}=['Mean: ',signamehigh{iCond}];
                    legend(legendcond,'Location','NorthEast');
                end
                % Calculate and plot Condition A-Condition B
                MEANCONV2{iCV,iCond}=MEANCONV{iCV,iCond}.*un_cv(iCV);
                if iCV==1
                    ga_mean{iCond,iChan}= nanmean(MEANCONV2{iCV,iCond},1);
                else
                    ga_mean{iCond,iChan}=ga_mean{iCond,iChan}+nanmean(MEANCONV2{iCV,iCond},1);
                end
            end
            ga_SEM{iCond,iChan} = nanstd(ga_mean{iCond,iChan},1)/sqrt(size(ga_mean{iCond,iChan},1));
            if INFO.groupplots==0
                subplot(2,iNoParMod,(iCond+iNoParMod));
            else
                subplot(2,1,2);
            end
            [l,p] = boundedline(x_axis_time,ga_mean{iCond,iChan},ga_SEM{iCond,iChan}, colorline{1},'alpha');
            hold on; plot([-10 INFO.gavg.window_size_secs],[0 0],'-r');
            xlabel('Time relative to stimulus onset (s)','FontSize',10);
            yname=['Difference between ',signamehigh{iCond},' and ',signamelow{iCond}];
            ylabel(yname,'FontSize',10);
            if INFO.groupplots==1
                fig_path = fullfile(INFO.file.sec_gafig,[INFO.dataselect.taskname,'_',INFO.file.name,'_',INFO.stat.con(ConNoParMod(iCond)).fname,'_channel',num2str(iChan),'.jpg']);
                saveas(cf,fig_path);
                close(cf);
            end
        end
        
        if INFO.groupplots==0
            fig_path = fullfile(INFO.file.sec_gafig,[INFO.dataselect.taskname,'_',INFO.file.name,'_channel',num2str(iChan),'.jpg']);
            saveas(cf,fig_path);
            close(cf);
        end
        clear sig_mem x_axis_time xlabel ylabel MEANCONV MEANCONV2
        
    end
end

% Calculate AUC of GA
for iChan=1:size(dat.condition_gavg,2);
    for iCond=1:iNoParMod;
        SPM_nirs.nirs.AUC{iCond,iChan}=mean(ga_mean{iCond,iChan});
    end
end

%Save GA means and SEMs
filename=[INFO.dataselect.taskname,'_',INFO.file.name];;
INFO.file.sec_ganame=fullfile(INFO.file.sec_gafile,filename);
save(INFO.file.sec_ganame,'ga_mean','ga_SEM');

end
