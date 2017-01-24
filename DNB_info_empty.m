function DNB_info
% INFO-file for the DNB. All analysis/filtering options can be changed here.
%
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO

%% 1. DNB: Data selection
INFO.dataselect.analysisname    = '';                                      % Name given to analysis. Important for report and second level analysis.
INFO.dataselect.subjects        = {};                                      % Name(s) of the subjects analyzed.
INFO.dataselect.taskname        = '';                                      % Name of the task performed

INFO.sessions                   = 1;                                       % Number of sessions measured for one subject. Enter 1 if there were no sessions used (INFO.sess.name is not neccessary to enter if no different sessions were used).
INFO.sess(1).name               = '';
INFO.sess(2).name               = '';

%% 2. DNB: General options
% DNB skips parts of the DNB_script if the results of these steps are
% already available. When INFO.overwrite = yes, DNB runs the scripts again
% if the results are available.
INFO.overwrite                  = 'yes';                                   % If 'no', before every step, the batch searches if this step has been done. If the step has already been taken, the batch will use the resulting database instead of running the step again. Options: 'yes' or 'no'.
INFO.plots                      = 'yes';                                   % Make plots of every step. Options: 'yes' or 'no'.
INFO.extension                  = 'jpg';                                   % Extension of all plots. Options: 'jpg' or 'fig'.
INFO.stopwhenerror              = 'no';                                    % Stops the script if an error occurs. If 'no' is selected, the batch will skip the subject with the error and move on to the next one.
INFO.delsteps                   = 'no';                                    % Delete all between mat files after generating the final mat file. Options: 'yes' or 'no'.
INFO.firstsecond                = 'first';                                 % First/second level analysis or both. 'first' : 1st level (individual) analysis; 'second' : 2nd level (group) analysis. 2nd level analysis can only be done when there are some 1st level analysis files available. 'both': First do 1st level analysis (over INFO.dataselect.subjects), followed by 2nd level over selected PP (in INFO.sec.subjects)).
INFO.SCI.check                  = 'yes';                                   % Make use of SCI script if 'yes'. Options: 'yes' or 'no'.
INFO.SCI.corrrequired           = 0.9;                                     % Correlation per channel between HbO and HbR is calculated (in heart-beat frequency). If the correlation is lower than INFO.SCI.correquired, and INFO.SCI.deletechan is 'yes', the channel-data will not be analysed.
INFO.SCI.chanrequired           = 5;                                       % Number of channels required after deleting channels with SCI to continue analyzing. Only used if INFO.SCI.deletechan is 'yes'.
INFO.report                     = 'yes';                                   % Construct a report after every analysis containing details about the analysis done. Options: 'yes' or 'no'.

%% 3. DNB: Add data directories
INFO.file.dir                   = '';                                      % Directory where raw datafile and processed files are stored.
INFO.file.extension_type        = 'oxy3';                                  % Extension type of raw_datafile. This should be synced with INFO.conv.system, because the NIRS-system used determines the extension (eg. Oxymon uses oxy3 extension).

% The following loop automatically searches for the right files because the file-names contain subjectname and taskname. Adapt all filenames to this formula if you want
% to use this algorithm. In case you don't want to use the algorithm: delete (or comment) the loop part and put all desired filenames manually
% in the INFO.file.rawdata_name and INFO.file.conditions_name matrices.
for iSubj = 1:numel(INFO.dataselect.subjects)
    % Names of files to load.
    INFO.file.rawdata_name{iSubj} = ['',INFO.dataselect.taskname,'_',INFO.dataselect.subjects{iSubj},'.',INFO.file.extension_type,'']; % Name of raw datafile including extension
    INFO.file.conditions_name{iSubj} = ['conditions_',INFO.dataselect.taskname,'_',INFO.dataselect.subjects{iSubj},'']; % Name of conditionsfile. The conditions_file is often given by the NIRS-system and should consist of a 'names', 'onsets' and 'durations' variable.
    %Directories of files to load.
    if INFO.sessions == 1
        INFO.file.Sess(1).rawdata_dir{iSubj} = fullfile(INFO.file.dir,INFO.dataselect.subjects{iSubj},INFO.dataselect.taskname,'Raw_datafiles',INFO.file.rawdata_name{iSubj});
        INFO.file.Sess(1).conditions_dir{iSubj} = fullfile(INFO.file.dir,INFO.dataselect.subjects{iSubj},INFO.dataselect.taskname,'Raw_datafiles',INFO.file.conditions_name{iSubj});
    else
        for iSess=1:INFO.sessions
            INFO.file.Sess(iSess).rawdata_dir{iSubj} = fullfile(INFO.file.dir,INFO.dataselect.subjects{iSubj},INFO.dataselect.taskname,INFO.sess(iSess).name,'Raw_datafiles',INFO.file.rawdata_name{iSubj});
            INFO.file.Sess(iSess).conditions_dir{iSubj} = fullfile(INFO.file.dir,INFO.dataselect.subjects{iSubj},INFO.dataselect.taskname,INFO.sess(iSess).name,'Raw_datafiles',INFO.file.conditions_name{iSubj});
        end
    end
end

INFO.file.MNI_file_name         = fullfile(INFO.file.dir,'Channel_information','MNI','channel_NIRS_WP4.mat'); % The MNI file gives information on the location of the NIRS cap and is required to produce brain images. 

%% 4. DNB: Add script directories
INFO.file.scriptdir              = '';                                     % Directory where scripts/toolboxes are stored
% Fieldtrip scripts
INFO.file.ft_dir='';                                                       % Some fieldtrip scripts are used in this batch. Fieldtrip can be downloaded for free at http://www.fieldtriptoolbox.org/
addpath(genpath(INFO.file.ft_dir));
% Batch scripts
INFO.file.DNB_dir=fullfile(INFO.file.scriptdir,'NIRS_analysis_scripts','DNB scripts'); 
addpath(genpath(INFO.file.DNB_dir));
% Other scripts
INFO.file.tb_dir=fullfile(INFO.file.scriptdir,'NIRS_analysis_scripts_linux','NIRS_analysis_scripts','Toolboxes');
addpath(genpath(INFO.file.tb_dir));

%% 5. DNB: Data conversion
INFO.conv.system                = 'oxymon';                                % Specific NIR system, e.g., 'oxymon'.
INFO.conv.fs                    = 250;                                     % Sampling frequency(Hz). Will be filled in automatically for the Oxymon.
INFO.conv.total_ch              = 8;                                       % Total number of channels, default in OXYMON: 24
INFO.conv.dist                  = [3.5,3.5,3.5,3.5,3.5,3.5];               % Distance between source and detector, enter in (cm). Will be filled in automatically for the Oxymon.
INFO.conv.wavelength            = [765 858];                               % Wavelength of the light source (nm). Will be filled in automatically for the Oxymon.
INFO.conv.DPF                   = 4;                                       % Differential pathlength factor.  Will be filled in automatically for the Oxymon.
INFO.conv.flag_DPF_corr         = 1;                                       % 0: without DPF correction, 1: with DPF correction. Will be filled in automatically for the Oxymon.

% *Only NIRX Dynot/Manual OD system need manual input of ext_coef/ch_config, for other systems [] should be entered,because the values are already given by the system.
INFO.conv.ext_coef              = [];                                      % Extinction coefficient (e.g. [1.4866 3.8437 2.2314 1.7917])*
INFO.conv.ch_config             = [];                                      % Name of channel configuration file*

% Memory issues solving
INFO.conv.downfs_method         = 'ft_resampling'                          % Downsampling the dataset using this method. Method choices: fieldtrip downsampling (ft_downsampling), fieltrip resampling (ft_resampling). If empty, default ft_resampling is used.
INFO.conv.downfs                = 5;                                       % Dataset is downsampled for downfs samples per second to avoid memory issues in large datasets.

%% 6. DNB: Model specification
INFO.model.hb                   = 'HbO';                                   % Specific hemoglobin, e.g., 'HbO' or 'HbR'.
INFO.model.hrf_type             = 0;                                       % Basis function to model the hemodynamic response, 0: 'hrf', 1: hrf (time der.), 2: hrf(time & dispersion der.).
INFO.model.units                = 'scans';                                  % Unit desired from model design. Options: 'scans' or 'secs'.

%% 7. DNB: Filtering
INFO.filt.HPF                   = 'DCT';                                   % High pass method (wavelet, DCT, detrend), enter cut-off value after DCT, default: 'DCT, 128'.
INFO.filt.LPF                   = 'hrf';                                   % Low pass method (hrf, gaussian, lpf), enter FWHM value after gaussian, default: 'gaussian, 4' and cut-off frequency after lpf, default: 'lpf, 0.5').
INFO.filt.MARA                  = 'yes';                                   % Use MARA on the dataset or not. Options are: 'yes', 'no', 'interactive'. When the option 'interactive' is used, the user can choose to use MARA or not while running the script (do not use when running overnight).
INFO.MARA.T                     = 3;                                       % Threshold for artifact detection in MARA, used in the moving standard deviation signal, and expressed in T * mean of the moving standard deviation signal.
INFO.MARA.L                     = 1;                                       % Length of the moving-window to calculated the moving standard deviation (in sec.).
INFO.MARA.alpha                 = 5;                                       % Parameter that defines how much high-frequency information should be preserved by the removal of the artifact (i.e., it corresponds to the length of the LOESS smoothing window)
% *MARA: Movement Artifact Removal Algorithm. Should only be used when
% large spike-like artifacts are present in the data.

%% 8. DNB: Grand average calculation
INFO.gavg.baseline_period_calc  = 1;                                       % Time (s) before onsets to calculate baseline on
INFO.gavg.baseline_period_plot  = 1;                                       % Time (s) before onsets to plot before onsets
INFO.gavg.window_size_secs      = 15;                                      % Time (s) window to calculate grand average over after stimulus onset

%% 9. DNB: Single subject statistics

INFO.stat.con(1).name        = '2>O';                                      % Names of the contrast vectors. Multiple contrast can be analyzed at the same time.
INFO.stat.con(1).positive    = {'2B'};                                     % Names of the positive contrasts in one contrast vector. Make sure the names agree with names given in the conditionsfile (not case-sensitive).
INFO.stat.con(1).negative    = {'OB'};                                     % Names of the negative contrasts in one contrast vector. Make sure the names agree with names given in the conditionsfile (not case-sensitive).
INFO.stat.con(2).name        = 'O>2';
INFO.stat.con(2).positive    = {'OB'};
INFO.stat.con(2).negative    = {'2B'};

INFO.stat.STAT                  = 'T';                                     % Either T- or F- statistics, e.g., 'T'
INFO.stat.spec_hemi             = {'frontal'};                             % Specific view of the rendered brain, including: ventral, dorsal, right, left, frontal, occipital
INFO.stat.p_value               = 0.05;                                    % P-value, e.g., 0.05
INFO.stat.correct_p             = 'EC';                                    % P-value correction method: 'EC': Lipschitz-Killing curvature based expected Euler characteristics, 'tube': Sun's tube formula, 'none': uncorrected
INFO.stat.disp_fig              = 1;                                       % 1 : show interpolated t- or F-map and activation map over the specific threshold. 0 : do not show the result.
% * If multiple contrast vectors should be investigated, add multiple names
% to con_name and add vectors in a new row of con_vec.

%% 10. DNB: Group statistics
INFO.sec.subjects               = {};                                      % Subjects used for the group analysis. These options are only applicable if INFO.firstsecond = second or both.                  
INFO.sec.overlapsubj            = 2;                                       % Required overlapping subjects. should be over 1; not more than number of subjects indicated in INFO.sec.subjects.
INFO.groupplots                 = 1;                                       % All grand average lines and extract activations will be plotted in 1 figure or each separately for the 2nd level plots. Options: 0: all in 1 figure, 1: all figures separately
