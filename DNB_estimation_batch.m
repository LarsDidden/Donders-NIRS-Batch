function DNB_estimation_batch
% NIRS-SPM batch script for 'Estimation' routine, which estimates the GLM
% parameters and temporal correlation.
% fname_SPM : the name of file which results from model specifiction
% e.g.,'...\NIRS_SPM_v3_3\Sample_data\categorical_indiv_HbO\SPM_indiv_HbO.mat';
% fname_nirs : the name of the nirs file
% e.g.,'...\NIRS_SPM_v3_3\Sample_data\converted_NIRS.mat';
% example usage
% >> fname_SPM =
% 'C:\NIRS_SPM_v3_3\Sample_data\categorical_indiv_HbO\SPM_indiv_HbO.mat';
% >> fname_nirs = 'C:\NIRS_SPM_v3_3\Sample_data\converted_NIRS.mat';
% >> [SPM_nirs] = estimation_batch(fname_SPM, fname_nirs);
%
% Adjusted to the DNB by:
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO;

warning off images:initSize:adjustingMag

fprintf('## %s: running for subject %s ##\n',mfilename,INFO.dataselect.subjectnow);

for isess=1:INFO.sessions
    INFO.counter.iSess=isess;
    
    %%%%%%%% Skip this step if already done and INFO.overwrite == no.
    if exist(INFO.file.pp_name);
        if strcmp(INFO.overwrite,'no')==1
            return
        end
    end
    %%%%%%%%
    
    if strcmp(INFO.filt.LPF(1:3),'lpf')==1; %If lpf is selected, set LPF_SMP to 'none' to prevent other lpf options to run.                                     % When simple LPF is selected, the SPM_LPF settings should be set to 'none'.
        INFO.filt.LPF_SPM='none';
    else
        INFO.filt.LPF_SPM=INFO.filt.LPF;
    end
    
    
    if isfield(INFO.file,'spec_name')==0
        error('Not all INFO components necessary to run model estimation are available. Add the missing components to the INFO-file to continue.');
    else
        fname_SPM       = INFO.file.spec_name;
        
        pp_name         = INFO.file.pp_name;
    end
    
    try
        [pathn, name, ext] = fileparts(fname_SPM);
        pathn = [pathn filesep];
    catch
        index = find(fname_SPM == filesep);
        pathn = fname_SPM(1:pathn(end));
    end
    load(fname_SPM);
    
    try     % Load in NIRS-dataset
        fname_nirs=INFO.file.MARA_name{INFO.counter.iSess};
        load(fname_nirs);
    catch
        try
            fname_nirs=INFO.file.conv_ds_name{INFO.counter.iSess};
            load(fname_nirs);
        catch
            fname_nirs=INFO.file.conv_name{INFO.counter.iSess};
            load(fname_nirs);
        end
    end
    switch SPM_nirs.nirs.step
        case 'estimation'
            disp('The process of estimation has been already done. Please do the step of model specification, again.');
            SPM_nirs = [];
            return;
    end
    
    disp('Model parameter estimation starts...');
    switch SPM_nirs.nirs.Hb
        case 'HbO'
            Y{INFO.counter.iSess} = nirs_data.oxyData;
        case 'HbR'
            Y{INFO.counter.iSess} = nirs_data.dxyData;
        case 'HbT'
            tf = isfield(nirs_data, 'tHbData');
            if tf == 0
                Y{INFO.counter.iSess} = nirs_data.oxyData + nirs_data.dxyData;
            elseif tf == 1
                Y{INFO.counter.iSess} = nirs_data.tHbData;
            end
    end
end

% estimation of GLM parameters using either percoloring or prewhitening
if isfield(SPM_nirs.xVi, 'V') == 1 % precoloring method
    SPM_nirs = rmfield(SPM_nirs, 'xVi');
    [SPM_nirs] = DNB_precoloring(SPM_nirs, Y);
elseif isfield(SPM_nirs.xVi, 'V') == 0
    [SPM_nirs] = DNB_prewhitening(SPM_nirs, Y, pathn);
end
disp('Completed.');

% delete precalculated files (e.g. interpolated beta, its
% covariance and t- or F-statistics)
% fname_others = cellstr(spm_select('FPList', pathn, ['^interp.*\' SPM_nirs.nirs.Hb '.mat$']));
% if strcmp(fname_others{1}, filesep) ~= 1
%     delete(fname_others{:});
% end
% fname_others = cellstr(spm_select('FPList', pathn, '^interp_matrix.*\.mat$'));
% if strcmp(fname_others{1}, filesep) ~= 1
%     delete(fname_others{:});
% end
disp('Estimation of model parameters has been completed.');

% % % % % % % Check correlations of HbO vs HbR before/after filtering
data_bf=load(INFO.file.conv_ds_name{isess});
for iChan=1:size(SPM_nirs.xX.KYsep{1,isess},2) 
    INFO.corrfilt(iChan,1)=corr(data_bf.nirs_data.oxyData(:,iChan),SPM_nirs.xX.KYsep{1,isess}(:,iChan));
    INFO.corrfilt(iChan,2)=corr(data_bf.nirs_data.dxyData(:,iChan),SPM_nirs.xX.KYsep{1,isess}(:,iChan));
    switch SPM_nirs.nirs.Hb
        case 'HbO'
            if INFO.corrfilt(iChan,1)<=INFO.corrfilt(iChan,2)
                warning(['Correlation between HbO data before filtering is higher with HbR after filtering than with HbO after filtering in channel ',num2str(iChan),'.']);
            end
        case 'HbR'
            if INFO.corrfilt(iChan,2)<=INFO.corrfilt(iChan,1)
                warning(['Correlation between HbR data before filtering is higher with HbO after filtering than with HbR after filtering in channel ',num2str(iChan),'.']);
            end
    end
end

% Save data
save(pp_name,'SPM_nirs');

% Plot/save before/after preprocessing figure
for isess=1:INFO.sessions
    INFO.counter.iSess=isess;
    
    if strcmp(INFO.plots,'yes');
        try
            load(INFO.file.MARA_name{isess});
        catch
            try
                load(INFO.file.conv_ds_name{isess});
            catch
                load(INFO.file.conv_name{isess});
            end
        end
        if strcmp(INFO.model.hb,'HbO');
            DNB_plotchannel(nirs_data.oxyData,SPM_nirs.xX.KYsep{isess},'all');
        else
            DNB_plotchannel(nirs_data.dxyData,SPM_nirs.xX.KYsep{isess},'all');
        end
    end
    
    % Make PDF-figures of all preprocessing steps
    % Place INFO.filt.pdfplot number of plots on one PDF. Afterwards, place
    % the remnants (remPDF) in another PDF.
    if strcmp(INFO.plots,'yes')
        remPDF=rem(nirs_data.nch,3);
        nrPDF=(nirs_data.nch-remPDF)/3;
        
        % Last PDF has less than 3 img
        if remPDF ~=0
            for iChan2=1:nrPDF+1
                if iChan2 == nrPDF+1; Chan1page = remPDF; else; Chan1page = 3; end
                for iChan=1:Chan1page
                    Chan=((iChan2-1)*3)+iChan;
                    if strcmp(INFO.SCI.check,'yes')==1
                        Chan=num2str(INFO.SCI.sessremchannel{end}(Chan));
                    else
                        Chan=num2str(Chan);
                    end
                    imgchanname=['channel ',Chan,'.jpg'];
                    % Read in images with 3 preprocessing steps and before/after image.
                    if strcmp(INFO.filt.MARA,'yes')==1 || strcmp(INFO.filt.MARA,'interactive')==1
                        MARAimg{iChan}=imread(fullfile(INFO.file.figures{isess},INFO.file.name,'MARA_figures',imgchanname));
                    end
                    LPFimg{iChan}=imread(fullfile(INFO.file.figures{isess},INFO.file.name,'PP_afterLPF',imgchanname));
                    HPFimg{iChan}=imread(fullfile(INFO.file.figures{isess},INFO.file.name,'PP_afterHPF',imgchanname));
                    filtimg{iChan}=imread(fullfile(INFO.file.figures{isess},INFO.file.name,'PP_before_after_filtering',imgchanname));
                end
                % Title master-pdf files
                if strcmp(INFO.SCI.check,'yes')==1
                    Chan=num2str(INFO.SCI.sessremchannel{end}(((iChan2-1)*3)+1):(((iChan2-1)*3)+Chan1page)); Channame=['Channels ',Chan];
                else
                    Chan=num2str((((iChan2-1)*3)+1):(((iChan2-1)*3)+Chan1page)); Channame=['Channels ',Chan];
                end
                % Saving directory
                PPfullimg_dir=fullfile(INFO.file.figures{isess},INFO.file.name,'PP_PDF_allfiltering');
                if ~exist(PPfullimg_dir); mkdir(PPfullimg_dir); end
                PPfullimg_name=fullfile(PPfullimg_dir,Channame);
                % Plot master-pdf and save
                if strcmp(INFO.filt.MARA,'yes')==1 || strcmp(INFO.filt.MARA,'interactive')==1
                    if Chan1page == 3
                        masterimg=[MARAimg{1}, LPFimg{1}, HPFimg{1}, filtimg{1}; MARAimg{2}, LPFimg{2}, HPFimg{2}, filtimg{2}; MARAimg{3}, LPFimg{3}, HPFimg{3}, filtimg{3}];
                    elseif Chan1page == 2
                        masterimg=[MARAimg{1}, LPFimg{1}, HPFimg{1}, filtimg{1}; MARAimg{2}, LPFimg{2}, HPFimg{2}, filtimg{2}];
                        
                    else
                        masterimg=[MARAimg{1}, LPFimg{1}, HPFimg{1}, filtimg{1}];
                    end
                else
                    if Chan1page == 3
                        masterimg=[LPFimg{1}, HPFimg{1}, filtimg{1}; LPFimg{2}, HPFimg{2}, filtimg{2}; LPFimg{3}, HPFimg{3}, filtimg{3}];
                    elseif Chan1page == 2
                        masterimg=[LPFimg{1}, HPFimg{1}, filtimg{1}; LPFimg{2}, HPFimg{2}, filtimg{2}];
                    else
                        masterimg=[LPFimg{1}, HPFimg{1}, filtimg{1}];
                    end
                    
                end
                imshow(masterimg);
                saveas(gcf,PPfullimg_name,'pdf');
                close(gcf);
            end
            % Last PDF has 3 img
        else
            for iChan2=1:nrPDF; Chan1page = 3;
                for iChan=1:Chan1page
                    Chan=((iChan2-1)*3)+iChan;
                    if strcmp(INFO.SCI.check,'yes')==1
                        Chan=num2str(INFO.SCI.sessremchannel{end}(Chan));
                    else
                        Chan=num2str(Chan);
                    end
                    imgchanname=['channel ',Chan,'.jpg'];
                    % Read in images with 3 preprocessing steps and before/after image.
                    if strcmp(INFO.filt.MARA,'yes')==1 || strcmp(INFO.filt.MARA,'interactive')==1
                        MARAimg{iChan}=imread(fullfile(INFO.file.figures{isess},INFO.file.name,'MARA_figures',imgchanname));
                    end
                    LPFimg{iChan}=imread(fullfile(INFO.file.figures{isess},INFO.file.name,'PP_afterLPF',imgchanname));
                    HPFimg{iChan}=imread(fullfile(INFO.file.figures{isess},INFO.file.name,'PP_afterHPF',imgchanname));
                    filtimg{iChan}=imread(fullfile(INFO.file.figures{isess},INFO.file.name,'PP_before_after_filtering',imgchanname));
                end
                % Title master-pdf files
                if strcmp(INFO.SCI.check,'yes')==1
                    Chan=num2str(INFO.SCI.sessremchannel{end}(((iChan2-1)*3)+1):(((iChan2-1)*3)+Chan1page)); Channame=['Channels ',Chan];
                else
                    Chan=num2str((((iChan2-1)*3)+1):(((iChan2-1)*3)+Chan1page)); Channame=['Channels ',Chan];
                end
                PPfullimg_dir=fullfile(INFO.file.figures{isess},INFO.file.name,'PP_PDF_allfiltering');
                if ~exist(PPfullimg_dir); mkdir(PPfullimg_dir); end
                PPfullimg_name=fullfile(PPfullimg_dir,Channame);
                % Plot master-pdf and save
                if strcmp(INFO.filt.MARA,'yes')==1 || strcmp(INFO.filt.MARA,'interactive')==1
                    if Chan1page == 3
                        masterimg=[MARAimg{1}, LPFimg{1}, HPFimg{1}, filtimg{1}; MARAimg{2}, LPFimg{2}, HPFimg{2}, filtimg{2}; MARAimg{3}, LPFimg{3}, HPFimg{3}, filtimg{3}];
                    elseif Chan1page == 2
                        masterimg=[MARAimg{1}, LPFimg{1}, HPFimg{1}, filtimg{1}; MARAimg{2}, LPFimg{2}, HPFimg{2}, filtimg{2}];
                    else
                        masterimg=[MARAimg{1}, LPFimg{1}, HPFimg{1}, filtimg{1}];
                    end
                else
                    if Chan1page == 3
                        masterimg=[LPFimg{1}, HPFimg{1}, filtimg{1}; LPFimg{2}, HPFimg{2}, filtimg{2}; LPFimg{3}, HPFimg{3}, filtimg{3}];
                    elseif Chan1page == 2
                        masterimg=[LPFimg{1}, HPFimg{1}, filtimg{1}; LPFimg{2}, HPFimg{2}, filtimg{2}];
                    else
                        masterimg=[LPFimg{1}, HPFimg{1}, filtimg{1}];
                    end
                end
                imshow(masterimg);
                saveas(gcf,PPfullimg_name,'pdf');
                close(gcf);
            end
        end
    end
end

