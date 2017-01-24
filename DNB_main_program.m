function DNB_main_program
% Main program of the Donders NIRS batch. This batch preprocesses NIRS data using various filters (see manual or INFO-file for options)
% to clean NIRS data. Subsequently the DNB processes the data using the General Linear Model and plots various figures showing results.

% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

clear all;
close all;
warning on verbose
warning off MATLAB:dispatcher:nameConflict

global INFO;
dbstop if error;

INFO=[];
INFO.infofile='';                                                          % Enter the INFO-file used for this analysis.
run(INFO.infofile);                                                        % Get INFO structure

if isempty(whos('INFO')); error('Please provide INFO structure'); end      % Check whether we have a valid INFO structure

DNB_pre_loopscript                                                         % 1.1./2.1. Constructs INFO.file_name; changes contrast vector and contrast names to the version required; add defaults.
if strcmp(INFO.firstsecond,'second')==0                                    % The following steps are not used in group analysis (2nd level).
    INFO.report.failed_subjects = []; INFO.report.failed_error=[]; INFO.report.SCI=[]; INFO.report.MARA=[]; % Pre-allocate filenames for report
    for iSubj=1:numel(INFO.dataselect.subjects)                            % Run several subjects at the same time
        INFO.dataselect.subjectnow=INFO.dataselect.subjects{iSubj}; INFO.counter.iSubj=iSubj; % Add current subject to the INFO structure.
        for iSess=1:INFO.sessions
            INFO.counter.iSess=iSess;
            try
                DNB_name;                                                  % 1.2. Construct names for maps/files/figures produced in the batch.
                DNB_data_conversion;                                       % 1.3. Read in raw data file and convert to mat file (NIRS_SPM format).
            catch err                                                      % When an error occurs: either skip the subject or show the error (base on INFO.stopwhenerror).
                if strcmp(INFO.stopwhenerror,'yes')==1
                    rethrow(err);
                else
                    fprintf('## Subject %s failed to run properly. Continuing with the next subject.. ##\n',INFO.dataselect.subjectnow);
                    INFO.report.failed_subjects = [INFO.report.failed_subjects, {INFO.dataselect.subjects{iSubj}}];
                    INFO.report.failed_error = [INFO.report.failed_error, {err.message}];
                end
            end
        end
        try
            if strcmp(INFO.SCI.check,'yes') ==1
                DNB_SCI;                                                   % 1.4. Quality check data using the scalp coupling index.
            end
        catch err                                                          % When an error occurs: either skip the subject or show the error (base on INFO.stopwhenerror).
            if strcmp(INFO.stopwhenerror,'yes')==1
                rethrow(err);
            else
                fprintf('## Subject %s failed to run properly. Continuing with the next subject.. ##\n',INFO.dataselect.subjectnow);
                INFO.report.failed_subjects = [INFO.report.failed_subjects, {INFO.dataselect.subjects{iSubj}}];
                INFO.report.failed_error = [INFO.report.failed_error, {err.message}];
            end
        end
        for iSess=1:INFO.sessions
            INFO.counter.iSess=iSess;
            try
                eval(['DNB_adjust_onsets_',INFO.dataselect.taskname]);     % 1.5. Adjusting onsets from onset file to match presentation onsets.
                if strcmp(INFO.conv.downfs,'none') == 0
                    DNB_downsampling;                                      % 1.6. Downsampling the dataset to avoid high memory loads.
                end
                if strcmp(INFO.filt.MARA,'yes')==1 || strcmp(INFO.filt.MARA,'interactive')==1
                    DNB_runMARA                                            % 1.7. Delete movement artifacts from the data.
                end
            catch err                                                      % When an error occurs: either skip the subject or show the error (base on INFO.stopwhenerror).
                if strcmp(INFO.stopwhenerror,'yes')==1
                    rethrow(err);
                else
                    fprintf('## Subject %s failed to run properly. Continuing with the next subject.. ##\n',INFO.dataselect.subjectnow);
                    INFO.report.failed_subjects = [INFO.report.failed_subjects, {INFO.dataselect.subjects{iSubj}}];
                    INFO.report.failed_error = [INFO.report.failed_error, {err.message}];
                end
            end
        end
        try
            DNB_model_specification;                                       % 1.8. Specifies the GLM including design matrix, temporal filtering and temporal correlation estimation.
            DNB_estimation_batch;                                          % 1.9. Estimates the GLM parameters and temporal correlation.
            if INFO.counter.iSubj==1
                DNB_construct_contrast_vector;
            end
            for iSess=1:INFO.sessions
                INFO.counter.iSess=iSess;
                DNB_extract_activation;                                    % 1.10. Calculates mean (baseline corrected) NIRS signal for all conditions.
                DNB_grand_averages;                                        % 1.11. Calculates average activation per 'block'.
            end
            DNB_MNI_coord2render;                                          % 1.12. Calculates coordinates of channels on the head from MNI-locations.
            DNB_extract_contrast_activation;                               % 1.13. Calculates mean NIRS signal for the desired contrast vectors.
            DNB_activation_map;                                            % 1.14. Calculates statistics (activation map).
            if strcmp(INFO.delsteps,'yes')==1;
                DNB_delete_steps;                                          % 1.15. Delete mat files made before the final file.
            end
        catch err                                                          % When an error occurs: either skip the subject or show the error (base on INFO.stopwhenerror).
            if strcmp(INFO.stopwhenerror,'yes')==1
                rethrow(err);
            else
                fprintf('## Subject %s failed to run properly. Continuing with the next subject.. ##\n',INFO.dataselect.subjectnow);
                INFO.report.failed_subjects = [INFO.report.failed_subjects, {INFO.dataselect.subjects{iSubj}}];
                INFO.report.failed_error = [INFO.report.failed_error, {err.message}];
            end
        end
    end
end


if strcmp(INFO.firstsecond,'first')==0                                     % Following steps are not used in single subject analysis (1st level).
    DNB_name_group;                                                        % 2.2. Construct names for maps/files/figures produced in the group analysis.
    DNB_estimation_batch_group;                                            % 2.3. Estimates the GLM parameters and temporal correlation for group analysis.
    DNB_activation_map_group;                                              % 2.4. Calculates statistics (activation map) for group statistics for group analysis.
    DNB_grand_averages_group;                                              % 2.5. Calculates average activation per 'block' for group analysis.
    DNB_extract_contrast_activation_group;                                 % 2.6. Calculates mean NIRS signal for the desired contrast vectors for group analysis.
end

DNB_report;                                                                % 1.16./2.8. Constructing a report file giving specific information about the analysis done (e.g. SCI/MARA results).




