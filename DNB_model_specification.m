function DNB_model_specification
% NIRS-SPM batch script for 'specify 1st level' routine, which specifies
% the general linear model (GLM) such as the design matrix, temporal
% filtering, and temporal correlation estimation.
%
% Input variables
% 1. fname_nirs : the name of the nirs file.
% e.g., '...\NIRS_SPM_v3_2\Sample_data\converted_NIRS.mat';
% 2. hb : specific hemoglobin, e.g., 'hbo' or 'hbr', or 'hbt'
% 3. HPF : detrending method (wavelet MDL or DCT based method), e.g.,
% 1) 'wavelet' : Wavelet-MDL based
% 2) 'DCT, 128' : DCT-based detrending with cut-off 128 (sec).
% 4. LPF : low pass filtering (hrf or Gaussian smoothing), e.g.,
% 1) 'hrf' : hemodynamic response function smoothing or
% 2) 'gaussian, 4' : Gaussian smoothing with FWHM 4 (sec)
% 5. method_cor : specific method for estimating the temporal correlation
% (precoloring or prewhitening)
% e.g., 0 : precoloring or 1 : prewhitening
% 6. dir_save : the directory to save the variable 'SPM_nirs'
% e.g., '...\NIRS_SPM_v3_2\Sample_data\categorical_indiv_HbO\';
% 7. flag_window : 0 : do not show the design matrix, 1 : show the design
% matrix
% 8. hrf_type : the basis function to model the hemodynamic response
% e.g., 0 - 'hrf', 1 - hrf (time der.), 2 - hrf( time & dispersion
% der.)
% 9. units : the unit for design, e.g., 0 - scan, 1 - secs
% 10. names : name for condition, e.g. 'right finger tpping'
% Note that the input variable for multiple conditions should be cell
% array,
% e.g., name{1} = 'finger tapping', name{2} = 'n-back test';
% 11. onsets : a vector of onset times, e.g. [42 93 144 195]
% 12. durations : the task durations, e.g. [21 21 21 21]
%
% example usage
% 1.
% >> fname_nirs = 'C:\NIRS_SPM_v3_2\Sample_data\converted_NIRS.mat';
% >> hb = 'hbo';
% >> HPF = 'wavelet';
% >> LPF = 'hrf';
% >> method_cor = 0;
% >> dir_save = 'C:\NIRS_SPM_v3_2\Sample_data\categorical_indiv_HbO\';
% >> flag_window = 1;
% >> hrf_type = 2;
% >> units = 1;
% >> names{1} = 'right finger tapping';
% >> onsets{1} = 42:51:501;
% >> durations{1} = 21 * ones(10,1);
% >> [SPM_nirs] = specification_batch(fname_nirs, hb, HPF, LPF, method_cor, dir_save, flag_window, hrf_type, units,  names, onsets, durations);
% 2. In the NIRS data from Hitachi ETG-4000 and ISS Imagent system, there
% is marker column that contains the vector of onsets and durations. In
% that case, this function automatically read the vector of onsets and
% durations from the data.
% >> fname_nirs = 'C:\NIRS_SPM_v3_2\Sample_data\converted_NIRS.mat';
% >> hb = 'hbo';
% >> HPF = 'wavelet';
% >> LPF = 'hrf';
% >> method_cor = 0;
% >> dir_save = 'C:\NIRS_SPM_v3_2\Sample_data\categorical_indiv_HbO\';
% >> flag_window = 1;
% >> hrf_type = 2;
% if you want to specify the name of condition,
% >> names{1} = 'right finger tapping';
% >> [SPM_nirs] = specification_batch(fname_nirs, hb, HPF, LPF, method_cor, dir_save, flag_window, hrf_type, names);
% if not,
% >> [SPM_nirs] = specification_batch(fname_nirs, hb, HPF, LPF, method_cor, dir_save, flag_window, hrf_type);
% 3. Simultaneous entry of multiple condition names, onsets, and durations
% using *.mat file is allowed. This option can be used to load all the
% required information (e.g. condition names, onset, and durations) in
% one-go. You will first need to create *.mat file containing the relevant
% information. This *.mat file must include the following cell arrays (each
% 1 x n) : name, onsets, and durations. Please refer to the sample file;
% e.g.,¡¦\Sample_data\NIRS_data_file\sample_multiple_condition.mat.
% >> fname_nirs = 'C:\NIRS_SPM_v3_2\Sample_data\converted_NIRS.mat';
% >> hb = 'hbo';
% >> HPF = 'wavelet';
% >> LPF = 'hrf';
% >> method_cor = 0;
% >> dir_save = 'C:\NIRS_SPM_v3_2\Sample_data\categorical_indiv_HbO\';
% >> flag_window = 1;
% >> hrf_type = 2;
% >> units = 1;
% >> names =
% 'C:\NIRS_SPM_v3_2\\Sample_data\NIRS_data_file\sample_multiple_condition.m
% at';
% >> [SPM_nirs] = specification_batch(fname_nirs, hb, HPF, LPF, method_cor,
% dir_save, flag_window, hrf_type, units,  names);
%
% Adjusted to the DNB by:
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO
%%%%%%%% Skip this step if already done and INFO.overwrite == no.
if exist(INFO.file.spec_name);
    if strcmp(INFO.overwrite,'no')==1
        return
    end
end
%%%%%%%%

fprintf('## %s: running for subject %s ##\n',mfilename,INFO.dataselect.subjectnow);

if isfield(INFO.model,'hb')==0 || isfield(INFO.model,'hrf_type')==0 || isfield(INFO.model,'units')==0 || isfield(INFO.filt,'HPF')==0 || isfield(INFO.filt,'LPF')==0
    error('Not all INFO components necessary to run model specification are available. Add the missing components to the INFO-file to continue.');
else
    hb                      =INFO.model.hb;
    flag_window             =INFO.plots;
    hrf_type                =INFO.model.hrf_type;
    units                   =INFO.model.units;
    
    HPF                     =INFO.filt.HPF;
    LPF                     =INFO.filt.LPF;
    if isfield(INFO.model,'method_cor')==1
        method_cor              =INFO.filt.method_cor;
    else
        method_cor=1;
    end
end

for s=1:INFO.sessions
    file_name                        = ['',INFO.dataselect.subjectnow,'_',INFO.dataselect.taskname,'_',INFO.file.name,'_unit_',INFO.model.units,'.mat'];
    conv_name                        = ['',INFO.dataselect.subjectnow,'_',INFO.dataselect.taskname,'_unit_',INFO.model.units,'.mat'];
    
    MARA_file       = fullfile(INFO.file.general_dir{s},'MARA_files',file_name);
    Conv_ds_file= fullfile(INFO.file.general_dir{s},'Converted_downsampled_files',file_name);
    Conv_file= fullfile(INFO.file.general_dir{s},'Converted_files',file_name);
    Conditions_adj_ds_file= fullfile(INFO.file.general_dir{s},'Adjusted_Downsampled_Conditions_files',sprintf('conditions_%s',file_name));
    Conditions_adj_file= fullfile(INFO.file.general_dir{s},'Adjusted_Conditions_files',sprintf('condition_%s',conv_name));
    Conditions_file= INFO.file.Sess(s).conditions_dir{INFO.counter.iSubj};
    
    try     % Load in NIRS-dataset
        fname_nirs{s}=MARA_file;
        load(fname_nirs{s});
    catch
        try
            fname_nirs{s}=Conv_ds_file;
            load(fname_nirs{s});
        catch
            fname_nirs{s}=Conv_file;
            load(fname_nirs{s});
        end
    end
    
    try     % Load in conditions-file
        load(Conditions_adj_ds_file);
    catch
        try
            load(Conditions_adj_file);
        catch
            load(Conditions_file);
        end
    end
    
    SPM.nscan(s) = size(nirs_data.oxyData,1);
    
    if isempty(units)==1 & isempty(names)==1 & isempty(onsets)==1 & isempty(durations)==1   % vector_onset
        vector_onset = nirs_data.vector_onset;
        SPM.xBF.UNITS = 'scans';
    elseif isempty(units)==1 & isempty(onsets)==1 & isempty(durations)==1                   % vector onset and name
        vector_onset = nirs_data.vector_onset;
        SPM.xBF.UNITS = 'scans';
        for kk = 1:size(names,2)
            SPM.Sess(s).U(kk).name = names(kk);
        end
    elseif isempty(onsets)==1 & isempty(durations)==1                                       % multiple condition
        if strcmp(units,'scans') == 1
            SPM.xBF.UNITS = 'scans';
        elseif strcmp(units,'secs') == 1
            SPM.xBF.UNITS = 'secs';
        end
        for kk = 1:size(names, 2)
            SPM.Sess(s).U(kk).name = names(kk);
            SPM.Sess(s).U(kk).ons = onsets{kk};
            SPM.Sess(s).U(kk).dur = durations{kk};
            if exist('P')
                SPM.Sess(s).U(kk).P = P;
            end
        end
    else
        if strcmp(units,'scans') == 1
            SPM.xBF.UNITS = 'scans';
        elseif strcmp(units,'secs') == 1
            SPM.xBF.UNITS = 'secs';
        end
        for kk = 1:size(names, 2)
            SPM.Sess(s).U(kk).name = names(kk);
            SPM.Sess(s).U(kk).ons = onsets{kk};
            SPM.Sess(s).U(kk).dur = durations{kk};
            if exist('P')
                SPM.Sess(s).U(kk).P = P(kk);
            end
            
        end
    end
end

SPM.xY.RT = 1/nirs_data.fs;
SPM.xBF.T = 10;
SPM.xBF.T0 = 1;
SPM.xBF.dt = SPM.xY.RT/SPM.xBF.T;


rep     = 0;
switch hrf_type
    case 0
        SPM.xBF.name = 'hrf';
    case 1
        SPM.xBF.name = 'hrf (with time derivative)';
    case 2
        SPM.xBF.name = 'hrf (with time and dispersion derivatives)';
end

SPM.xBF = nirs_spm_get_bf(SPM.xBF);

switch hb
    case 'HbO'
        bf = SPM.xBF.bf;
    case 'HbR'
        bf = SPM.xBF.bf * (-1);
end

V = 1;
SPM.xBF.Volterra = V; % model interactions (Volterra) : no

Xx    = [];
Xb    = [];
Xname = {};
Bname = {};
for s = 1:length(SPM.nscan)
    if (s == 1) | ~rep
        k   = SPM.nscan(s);
        if isempty(units)==1 & isempty(names)==1 & isempty(onsets)==1 & isempty(durations)==1 | isempty(names)==1 & isempty(onsets)==1 & isempty(durations)==1
            tSPM = SPM;
            tSPM.vector_onset = vector_onset;
            U = DNB_nirs_spm_get_ons_batch(tSPM, s, 2);
        else
            U = DNB_nirs_spm_get_ons_batch(SPM, s, 1);
        end
        [X,Xn,Fc] = spm_Volterra(U,bf,V);
        try
            X = X([0:(k - 1)]*SPM.xBF.T + SPM.xBF.T0 + 32,:);
        end
        
        for i = 1:length(Fc)
            X(:,Fc(i).i) = spm_orth(X(:,Fc(i).i));
        end
        
        try
            C     = SPM.Sess(s).C.C;
            Cname = SPM.Sess(s).C.name;
        catch
            %             str   = sprintf('Session %d',s);
            %             spm_input('Other regressors',1,'d',str)
            C     = [];
            %             c     = spm_input('user specified','+1','w1',0);
            c = 0;
            while size(C,2) < c
                str = sprintf('regressor %i',size(C,2) + 1);
                C  = [C spm_input(str,2,'e',[],[k Inf])];
            end
            
            Cname = {};
            for i = 1:size(C,2)
                str      = sprintf('regressor %i',i);
                Cname{i} = spm_input('name of','+0','s',str);
            end
        end
        
        X      = [X spm_detrend(C)];
        Xn     = {Xn{:}   Cname{:}};
        
        B      = ones(k,1);
        Bn{1}  = sprintf('constant');
        
    end
    
    SPM.Sess(s).U      = U;
    SPM.Sess(s).C.C    = C;
    SPM.Sess(s).C.name = Cname;
    SPM.Sess(s).row    = size(Xx,1) + [1:k];
    SPM.Sess(s).col    = size(Xx,2) + [1:size(X,2)];
    SPM.Sess(s).Fc     = Fc;
    
    % Append names
    %---------------------------------------------------------------
    for i = 1:length(Xn)
        Xname{end + 1} = [sprintf('Sn(%i) ',s) Xn{i}];
    end
    for i = 1:length(Bn)
        Bname{end + 1} = [sprintf('Sn(%i) ',s) Bn{i}];
    end
    
    % append into Xx and Xb
    %===============================================================converted_data
    Xx    = blkdiag(Xx,X);
    Xb    = blkdiag(Xb,B);
    
end %- for s

% finished
%-----------------------------------------------------------------------
SPM.xX.X      = [Xx Xb];
SPM.xX.iH     = [];
SPM.xX.iC     = [1:size(Xx,2)];
SPM.xX.iB     = [1:size(Xb,2)] + size(Xx,2);
SPM.xX.iG     = [];
SPM.xX.name   = {Xname{:} Bname{:}};
% end

nscan = SPM.nscan;
nsess = length(nscan);

%%% updated for wavelet-MDL detrending 2009-03-19
str = 'Detrending?';

for s = 1:length(SPM.nscan)
    if isempty(strfind(HPF, 'wavelet')) == 0 % wavelet-MDL
        SPM.xX.K(s).HParam.type = 'Wavelet-MDL';
    elseif isempty(strfind(HPF, 'DCT')) == 0 || isempty(strfind(HPF, 'detrend')) == 0 % DCT
        index_cutoff = find(HPF == ',');
        if isempty(index_cutoff) == 1
            cutoff = 128;
        else
            cutoff = str2num(HPF(index_cutoff+1:end));
        end
        SPM.xX.K(s).HParam.type = 'DCT';
        SPM.xX.K(s).HParam.M = cutoff;
    end
    
    if isempty(strfind(LPF, 'hrf')) == 0 % hrf smoothing
        SPM.xX.K(s).LParam.type = 'hrf';
    elseif isempty(strfind(LPF, 'gaussian')) == 0 % Gaussian smoothing
        index_FWHM = find(LPF == ',');
        if isempty(index_FWHM) == 1
            FWHM = 4;
        else
            FWHM = str2num(LPF(index_FWHM+1:end));
        end
        SPM.xX.K(s).LParam.FWHM = FWHM;
        SPM.xX.K(s).LParam.type = 'Gaussian';
    else
        SPM.xX.K(s).LParam.type = 'none';
    end
end

for i=1:nsess
    K(i) = struct( 'HParam', SPM.xX.K(i).HParam,...
        'row', SPM.Sess(i).row,...
        'RT', SPM.xY.RT,...
        'LParam', SPM.xX.K(i).LParam);
end
    SPM.xX.K = DNB_spm_filter_HPF_LPF_WMDL(K);

% related spm m-file : spm_fmri_spm_ui.m
if method_cor == 0
    cVi = 'none';
elseif method_cor == 1
    cVi = 'AR(1)';
end

if ~ischar(cVi)	% AR coeficient[s] specified
    SPM.xVi.Vi = spm_Ce(nscan,cVi(1:3));
    cVi        = ['AR( ' sprintf('%0.1f ',cVi) ')'];
    
else
    switch lower(cVi)
        case 'none'		%  xVi.V is i.i.d
            %---------------------------------------------------------------
            SPM.xVi.V  = speye(sum(nscan));
            cVi        = 'i.i.d';
        otherwise		% otherwise assume AR(0.2) in xVi.Vi
            %---------------------------------------------------------------
            SPM.xVi.Vi = spm_Ce(nscan,0.2);
            cVi        = 'AR(0.2)';
    end
end
SPM.xVi.form = cVi;
SPM.xsDes = struct('Basis_functions', SPM.xBF.name, 'Sampling_period_sec', num2str(SPM.xY.RT), 'Total_number_of_samples', num2str(SPM.nscan));
if strcmp(flag_window,'yes') == 1
    spm_DesRep('DesMtx',SPM.xX,[],SPM.xsDes)
end

SPM.nirs.step = 'specification';
SPM.nirs.fname = fname_nirs;

SPM_nirs = SPM;
switch hb
    case 'HbO'
        SPM_nirs.nirs.Hb = 'HbO';
        SPM_nirs.nirs.level = 'individual';
    case 'HbR'
        SPM_nirs.nirs.Hb = 'HbR';
        SPM_nirs.nirs.level = 'individual';
end

if strcmp(flag_window,'yes')==1
    % save design figure
    fg = spm_figure('FindWin','Graphics');
    spec_design_dir        = fullfile(INFO.file.spec_dir,'Model_estimation');
    if ~exist(spec_design_dir); mkdir(spec_design_dir); end
    spec_design_name       = fullfile(spec_design_dir,['',INFO.dataselect.subjectnow,'_',INFO.dataselect.taskname,'_',INFO.file.name,'_design.pdf']);
    saveas(fg,spec_design_name);
    close(fg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if contrasts given in INFO.stat.con exist in SPM structure
if INFO.counter.iSubj==1 %Only check/change contrasts once
for iCon=1:numel(INFO.stat.con)
    for iConparts=1:numel(INFO.stat.con(iCon).positive)
        conname=INFO.stat.con(iCon).positive{iConparts};
        forbiddenchar1=regexpi(conname,'\('); %Search if '(' is in contrast-name. If so, place a backward slash before it, so regexpi still works.
        if isempty(forbiddenchar1)==0; conname2=[conname(1:(forbiddenchar1-1)),'\',conname(forbiddenchar1:end)]; else; conname2=conname; end
        forbiddenchar2=regexpi(conname2,'\)'); %Search if ')' is in contrast-name. If so, place a backward slash before it, so regexpi still works.
        if isempty(forbiddenchar2)==0; conname3=[conname2(1:(forbiddenchar2-1)),'\',conname2(forbiddenchar2:end)]; else; conname3=conname2; end
        IndexP = regexpi(SPM_nirs.xX.name, conname3);
        IndexPos{iConparts} = find(not(cellfun('isempty', IndexP)));
        IsNone = strcmp('none', conname3);
        if isempty(IndexPos{iConparts})==1 & IsNone==0 % Error when contrast not found.
            error(['Contrast ',INFO.stat.con(iCon).positive{iConparts},' was not found in the SPM_nirs structure. Please use contrasts that agree with the names given in the conditionsfile.']);
        end
        clear conname conname2 conname3 forbiddenchar1 forbiddenchar2    
    end
    for iConparts=1:numel(INFO.stat.con(iCon).negative)
        conname=INFO.stat.con(iCon).negative{iConparts};
        forbiddenchar1=regexpi(conname,'\('); %Search if '(' is in contrast-name. If so, place a backward slash before it, so regexpi still works.
        if isempty(forbiddenchar1)==0; conname2=[conname(1:(forbiddenchar1-1)),'\',conname(forbiddenchar1:end)]; else; conname2=conname; end
        forbiddenchar2=regexpi(conname2,'\)'); %Search if ')' is in contrast-name. If so, place a backward slash before it, so regexpi still works.
        if isempty(forbiddenchar2)==0; conname3=[conname2(1:(forbiddenchar2-1)),'\',conname2(forbiddenchar2:end)]; else; conname3=conname2; end
        IndexN = regexpi(SPM_nirs.xX.name, conname3); %Search for names given in INFO
        IndexNeg{iConparts} = find(not(cellfun('isempty', IndexN)));
        IsNone = strcmp('none', conname3);
        if isempty(IndexNeg{iConparts})==1 & IsNone==0 % Error when contrast not found.
            error(['Contrast ',INFO.stat.con(iCon).negative{iConparts},' was not found in the SPM_nirs structure. Please use contrasts that agree with the names given in the conditionsfile.']);
        end
        clear conname conname2 conname3 forbiddenchar1 forbiddenchar2
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save file
save(INFO.file.spec_name, 'SPM_nirs');



