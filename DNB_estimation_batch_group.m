function DNB_estimation_batch_group
% NIRS-SPM batch script for 'Estimation' routine, which estimates the GLM
% parameters and temporal correlation for groups of subjects (2nd level).
%
% Adjusted to the DNB by:
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO

fprintf('## %s: running for group %s ##\n',mfilename,INFO.dataselect.analysisname);



for iCon=1:numel(INFO.stat.con)
    for iView=1:numel(INFO.stat.spec_hemi)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Remove some columns fom con_vec because of INFO.model.hrf_type
        RemMat=[];  
        if INFO.model.hrf_type==1
            for i=1:((size(INFO.stat.con(iCon).vec,2)-1)/2)
                RemMat=[(RemMat),(i*2)];
            end
        elseif INFO.model.hrf_type==2
            for i=1:((size(INFO.stat.con(iCon).vec,2)-1)/3)
                RemMat=[(RemMat),(i*3-1),(i*3)];
            end
        end
        INFO.stat.con(iCon).vec(:,RemMat)=[];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Create map for group estimation
        
        % Assign pathways for all subject-files for group estimation into one variable
        for iSubj=1:numel(INFO.sec.subjects)
            beta_filename     = ['interp_beta_',INFO.stat.spec_hemi{iView},'_',INFO.model.hb,'_',INFO.stat.con(iCon).fname,'_',INFO.sec.subjects{iSubj},'.mat'];
            if INFO.sessions>1
                stat_dir = fullfile(INFO.file.dir,INFO.sec.subjects{iSubj},INFO.dataselect.taskname,'Combined_Session','Stat_files');                               % directory where datafile should be saved after Statistics/activation mapping
            else
                stat_dir = fullfile(INFO.file.dir,INFO.sec.subjects{iSubj},INFO.dataselect.taskname,'Stat_files');                               % directory where datafile should be saved after Statistics/activation mapping
            end
            fname_ginterp_betas{iSubj}=fullfile(stat_dir,beta_filename);
        end
        %Number of subjects and minimum number needed to overlap
        nsubj = numel(INFO.sec.subjects);
        min_subj = INFO.sec.overlapsubj;
        
        [gavg_beta, group_beta, xX, index_group, nsubj_mask, brain_view, fname_interp_cov] = nirs_spm_group(fname_ginterp_betas, min_subj);
        % gavg_beta: for group t-stat
        % group_beta: for group F-stat
        % common variable: index_group, nsubj_mask, brain_view, fname_interp_var
        dir_spm = INFO.file.sec_stat_dir{iCon}{iView};    %Dir where group statistics will be saved
        SPM_nirs.nirs.level = 'group';
        SPM_nirs.nirs.nsubj = nsubj;
        SPM_nirs.nirs.min_subj = min_subj;
        SPM_nirs.nirs.Hb=INFO.model.hb;
        SPM_nirs.xX = xX;
        SPM_nirs.nirs.brain_view = brain_view; %name,index,size
        SPM_nirs.nirs.index_group = index_group;
        SPM_nirs.nirs.nsubj_mask = nsubj_mask;
        SPM_nirs.nirs.fname_ginterp_cov = fname_interp_cov;
        SPM_nirs.nirs.fname_ginterp_beta = [dir_spm filesep 'ginterp_beta_' SPM_nirs.nirs.brain_view.name '_' SPM_nirs.nirs.Hb '.mat'];
        SPM_nirs.nirs.fname_ginterp_avg_beta = [dir_spm filesep 'ginterp_avgbeta_' SPM_nirs.nirs.brain_view.name '_' SPM_nirs.nirs.Hb '.mat'];
%         fname_others = cellstr(spm_select('FPList', dir_spm, ['^ginterp_T.*\' SPM_nirs.nirs.Hb '.mat']));
%         if strcmp(fname_others{1}, filesep) ~= 1
%             delete(fname_others{:});
%         end
%         fname_others = cellstr(spm_select('FPList', dir_spm, ['^ginterp_F.*\' SPM_nirs.nirs.Hb '.mat$']));
%         if strcmp(fname_others{1}, filesep) ~= 1
%             delete(fname_others{:});
%         end
     
        sec_SPM_name              = [INFO.dataselect.taskname,'_',INFO.file.name,'_',INFO.stat.con(iCon).fname,'_',INFO.stat.spec_hemi{iView},'.mat'];
        INFO.file.sec_estimationname               = fullfile(INFO.file.sec_estimationfile,sec_SPM_name);
        
        % save a SPM_nirs_group_brainview_HbX.mat in specified dir.
        save(INFO.file.sec_estimationname, 'SPM_nirs');
        % save a group-average beta
        save(SPM_nirs.nirs.fname_ginterp_avg_beta, 'gavg_beta');
        % save matrix containing individual interpolated betas
        save(SPM_nirs.nirs.fname_ginterp_beta, 'group_beta');
    end
end

