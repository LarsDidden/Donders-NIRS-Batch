function DNB_activation_map_group
% Calculation of the activation map over the threshold for groups of
% subjects (2nd level).
%
% Adjusted to the DNB by:
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO

fprintf('## %s: running for group %s ##\n',mfilename,INFO.dataselect.analysisname);
for iCon=1:numel(INFO.stat.con)
    for iView = 1:numel(INFO.stat.spec_hemi)
        
        %Load in contrast vector
        fnamecon = [INFO.file.name,'_',INFO.stat.con(iCon).fname];
        if INFO.sessions>1
            convecdir = fullfile(INFO.file.dir,INFO.sec.subjects{1},INFO.dataselect.taskname,'Combined_Session','Stat_files',fnamecon);                               % directory where datafile should be saved after Statistics/activation mapping
        else
            convecdir = fullfile(INFO.file.dir,INFO.sec.subjects{1},INFO.dataselect.taskname,'Stat_files',fnamecon);                               % directory where datafile should be saved after Statistics/activation mapping
        end
        load(convecdir); c=SPM_nirs.xCon.c;
        clear SPM_nirs
        
        %Load in estimation-group file
        sec_SPM_dir               = fullfile(INFO.file.sec_estimationfile,[INFO.dataselect.taskname,'_',INFO.file.name,'_',INFO.stat.con(iCon).fname,'_',INFO.stat.spec_hemi{iView},'.mat']);
        load(sec_SPM_dir);
        
        % specification of contrast vector
        %[Ic, xCon] = nirs_spm_conman(SPM_nirs, 'T&F', Inf, 'Select contrasts...', 'for conjunction', 1);
        if isfield(SPM_nirs, 'xCon') == 1
            xCon = SPM_nirs.xCon;
            Ic = size(SPM_nirs.xCon, 2) + 1;
        else
            Ic = 1;
        end
        
        xCon(1,Ic) = struct('name', INFO.stat.con(iCon).fname, 'STAT', INFO.stat.STAT, 'c', c, 'X0', [], 'iX0', [], 'X1o', [], 'eidf', [], 'Vcon', [], 'Vspm', []);
        SPM_nirs.xCon = xCon;
        
        %------------------------------------------------------------------
        % calculation of group-level T- or F-statistics
        %------------------------------------------------------------------
        stat = []; % T- or F-statistic matrix
        % loading a file (interpolated T- or F-statistics)
        filename = [INFO.file.sec_stat_dir{iCon}{iView} filesep 'ginterp_' xCon(Ic).STAT '_'  SPM_nirs.nirs.Hb '.mat'];
        %         fid = fopen(filename);
        %         if fid ~= -1 % if exists
        %             fclose(fid);
        %             try
        %                 load(filename);
        %             end
        %         end
        if isempty(stat) == 1 || isempty(index_mask) == 1
            index_group = SPM_nirs.nirs.index_group;
            nvox = length(index_group);
            nvox_brain = SPM_nirs.nirs.brain_view.size(1) * SPM_nirs.nirs.brain_view.size(2);
            nsubj = SPM_nirs.nirs.nsubj;
            L = SPM_nirs.nirs.nsubj_mask(index_group); % # of subject on each voxel
            
            % contrast * group-level interp beta
            load(SPM_nirs.nirs.fname_ginterp_avg_beta);
            cgavg_beta = xCon(Ic).c' * gavg_beta;
            load(SPM_nirs.nirs.fname_ginterp_beta);
            
            % group-level residuals
            Xg = ones(nsubj, 1);
            Yg = xCon(Ic).c' * [group_beta{:}];
            Yg = reshape(Yg', [nvox, nsubj])';
            res = Yg - Xg * cgavg_beta;
            % LKC calculation
            disp('Calculating Lipschitz-Killing curvatures ...');
            [L2] = calc_LKC(index_group, SPM_nirs.nirs.brain_view.size, res, SPM_nirs.nirs.level);
            disp('Completed.');
            r = sqrt(L2./pi);
            L1 = pi * r;
            L0 = 1;
            SPM_nirs.nirs.LKC{SPM_nirs.nirs.brain_view.index} = [L0 L1 L2];
            
            df = []; % degree of freedom
            switch xCon(Ic).STAT
                case 'T' % t-statistic calculation
                    var_s = zeros(1, nvox);
                    for kk = 1:nsubj
                        var_s = var_s + (xCon(Ic).c' * group_beta{kk} - cgavg_beta).^2;
                    end
                    var_s = var_s ./ (L-1);
                    gsum_var = zeros(1, nvox_brain);
                    term_df = zeros(1, nvox_brain);
                    
                    % calculation group-level variance
                    for kk = 1:nsubj
                        load(SPM_nirs.nirs.fname_ginterp_cov{kk});
                        indiv_cov = interp_var .* (xCon(Ic).c' * xCor * xCon(Ic).c);
                        gsum_var(index_mask) = gsum_var(index_mask) + indiv_cov;
                        term_df(index_mask) = term_df(index_mask) + (indiv_cov.^2)./df(2);
                    end
                    gsum_var = gsum_var(index_group);
                    term_df = term_df(index_group);
                    term = var_s .* L + gsum_var;
                    
                    stat = (xCon(Ic).c' * gavg_beta)./(sqrt(term)./L);
                    term2 = ((L.^2)./(L-1)).* (var_s.^2);
                    df(2) = sum((term.^2) ./ (term2 + term_df))./nvox; % degree of freedom
                case 'F' % F-statistic calculation
                    X1o = pinv(Xg');
                    [trMV, trMVMV] = spm_SpUtil('trMV', X1o, 1);
                    df(1) = trMV^2/trMVMV;
                    R = eye(nsubj) - Xg*pinv(Xg);
                    RVR = sum(res.^2)./trace(R);
                    MVM = (inv(sum(X1o.^2)).* (cgavg_beta.^2))./trMV;
                    stat = MVM./RVR;
                    df(2) = trace(R)^2./trace(R*R);
            end
            stat = stat(:);
            index_mask = index_group;
            save(filename, 'stat', 'df', 'index_mask');
        end
        side_hemi = SPM_nirs.nirs.brain_view.index;
        load([spm('dir') filesep 'rend' filesep 'render_single_subj.mat']);
        %%% modified (2008. 10. 13) %%%
        brain = rend{side_hemi}.ren;
        if issparse(brain),
            d = size(brain);
            B1 = spm_dctmtx(d(1),d(1));
            B2 = spm_dctmtx(d(2),d(2));
            brain = B1*brain*B2';
        end
        msk = find(brain>1);brain(msk)=1;
        msk = find(brain<0);brain(msk)=0;
        brain = brain(end:-1:1,:);
        brain = brain * 0.5;
        
        
        % min and max values of statistics for colormap control
        min_stat = min(stat);
        max_stat = max(stat);
        smin_stat = max_stat - ((max_stat - min_stat)./63) * 127;
        sbar = linspace(smin_stat, max_stat, 128);
        
        stat_brain = ((-sbar(1) + sbar(64))/(0.5)).*brain + sbar(1);
        stat_brain(index_mask) = stat;
        
        if strcmp(INFO.plots,'yes')==1
            figure;
            imagesc(stat_brain);
            load Split
            colormap(split)
            axis off
            axis image
            title([INFO.dataselect.taskname,' ',INFO.model.hb,' ' INFO.stat.STAT ' map: ',INFO.stat.con(iCon).name],'Fontsize',15);
            hc = colorbar;
            set(hc, 'YLim', [sbar(65) sbar(128)]);
            y_tick = linspace(sbar(65), sbar(128), 5)';
            set(hc, 'YTick', y_tick);
            set(hc, 'FontSize', 8);
        end
        
        fnamecon = [INFO.dataselect.taskname,'_',INFO.file.name,'_',INFO.stat.con(iCon).fname,'_',INFO.stat.spec_hemi{iView}];
        
        % Save T/P-map
        tp_path              = fullfile(INFO.file.sec_brainfig,[fnamecon,'_unthresholded']);
        saveas(gcf,tp_path,INFO.extension);
        close(gcf);
        % update & save the SPM_file
        save_dir = fullfile(INFO.file.sec_stat_dir{iCon}{iView},fnamecon);
        save(save_dir, 'SPM_nirs');
        % inference of brain activation
        
        % calculation of threshold value, depending on a correction method
        switch INFO.stat.correct_p
            case 'EC'
                threshold = DNB_calc_EC(SPM_nirs.nirs.LKC{side_hemi}, INFO.stat.p_value, SPM_nirs.xCon(Ic).STAT, df);
            case 'tube'
                threshold = calc_tube(SPM_nirs.nirs.kappa(Ic, side_hemi), INFO.stat.p_value);
            case 'none'
                threshold = spm_u(INFO.stat.p_value, df, SPM_nirs.xCon(Ic).STAT);
        end
        
        index_over = find(stat > threshold);
        load Split;
        if isempty(index_over) == 1
            disp('There is no significant voxel.');
            act_brain = brain;
            split = split(1:64,:);
            sbar = [];
        else
            index_mask = index_mask(index_over);
            stat = stat(index_over);
            % min and max values of statistics for colormap control
            min_stat = min(stat);
            max_stat = max(stat);
            smin_stat = max_stat - ((max_stat - min_stat)./63) * 127;
            sbar = linspace(smin_stat, max_stat, 128);
            act_brain = ((-sbar(1) + sbar(64))/(0.5)).*brain + sbar(1);
            act_brain(index_mask) = stat;
        end
        
        if strcmp(INFO.plots,'yes')==1
            figure;
            imagesc(act_brain);
            colormap(split)
            axis off
            axis image
            title([INFO.dataselect.taskname,' ',INFO.model.hb,' (p<' num2str(INFO.stat.p_value) ', ' INFO.stat.correct_p ' correction: ',num2str(numel(index_over)),' voxels): ',INFO.stat.con(iCon).name],'FontSize',15)
            try
                hc = colorbar;
                set(hc, 'YLim', [sbar(65) sbar(128)]);
                y_tick = linspace(sbar(65), sbar(128), 5)';
                set(hc, 'YTick', y_tick);
                set(hc, 'FontSize', 8);
            end
        end
        
        % Save activation-map
        activ_path              = fullfile(INFO.file.sec_brainfig,[fnamecon,'_thresholded']);
        saveas(gcf,activ_path,INFO.extension);
        close(gcf);
    end
end


