function DNB_name
%Constructing file names and directories to use when writing away files.
%
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO

fprintf('## %s: running for subject %s ##\n',mfilename,INFO.dataselect.subjectnow);

isess=INFO.counter.iSess;

%% Construct name files/figures and pathways to save
if INFO.sessions==1
    INFO.file.general_dir{isess}=fullfile(INFO.file.dir,INFO.dataselect.subjectnow,INFO.dataselect.taskname);
else
    INFO.file.general_dir{isess}=fullfile(INFO.file.dir,INFO.dataselect.subjectnow,INFO.dataselect.taskname,INFO.sess(isess).name);
end

%Construct save-directories per session
conv_dir      = fullfile(INFO.file.general_dir{isess},'Converted_files');                        % directory where datafile should be saved after conversion
if ~exist(conv_dir); mkdir(conv_dir); end
conv_ds_dir          = fullfile(INFO.file.general_dir{isess},'Converted_downsampled_files');                        % directory where datafile should be saved after conversion
if ~exist(conv_ds_dir); mkdir(conv_ds_dir); end
INFO.file.MARA_dir{isess}  = fullfile(INFO.file.general_dir{isess},'MARA_files');                        % directory where datafile should be saved after MARA
if ~exist(INFO.file.MARA_dir{isess}); mkdir(INFO.file.MARA_dir{isess}); end
INFO.file.figures{isess}     = fullfile(INFO.file.general_dir{isess},'Figures');       % directory where datafile should be saved after SPM modeling
if ~exist(INFO.file.figures{isess}); mkdir(INFO.file.figures{isess}); end
INFO.file.channel{isess} = fullfile(INFO.file.dir,'Channel_information');
if ~exist(INFO.file.channel{isess}); mkdir(INFO.file.channel{isess}); end
INFO.file.extract_dir{isess} = fullfile(INFO.file.general_dir{isess},'Extract_activation_files');
if ~exist(INFO.file.extract_dir{isess}); mkdir(INFO.file.extract_dir{isess}); end
INFO.file.grand_avg_dir{isess} = fullfile(INFO.file.general_dir{isess},'Grand_average_files');
if ~exist(INFO.file.grand_avg_dir{isess}); mkdir(INFO.file.grand_avg_dir{isess}); end

% Construct save-directories (only one necessary over all sessions)
if INFO.sessions==1
    INFO.file.spec_dir      = fullfile(INFO.file.general_dir{isess},'Model_specification_files');                       % directory where datafile should be saved after model specification
    pp_dir                  = fullfile(INFO.file.general_dir{isess},'Estimation_batch_files');                               % directory where datafile should be saved after SPM modeling
    INFO.file.coor_dir      = fullfile(INFO.file.general_dir{isess},'Channel_Coordinates_files');
    INFO.file.stat_dir      = fullfile(INFO.file.general_dir{isess},'Stat_files');                               % directory where datafile should be saved after Statistics/activation mapping
    INFO.file.info_dir      = fullfile(INFO.file.general_dir{isess},'INFO_files');                               % directory where datafile should be saved after Statistics/activation mapping
    INFO.file.contrast_value_dir = fullfile(INFO.file.general_dir{isess},'Contrast_value_files');
else
    INFO.file.spec_dir      = fullfile(INFO.file.dir,INFO.dataselect.subjectnow,INFO.dataselect.taskname,'Combined_Session','Model_specification_files');                       % directory where datafile should be saved after model specification
    pp_dir                  = fullfile(INFO.file.dir,INFO.dataselect.subjectnow,INFO.dataselect.taskname,'Combined_Session','Estimation_batch_files');                               % directory where datafile should be saved after SPM modeling
    INFO.file.coor_dir      = fullfile(INFO.file.dir,INFO.dataselect.subjectnow,INFO.dataselect.taskname,'Combined_Session','Channel_Coordinates_files');
    INFO.file.stat_dir      = fullfile(INFO.file.dir,INFO.dataselect.subjectnow,INFO.dataselect.taskname,'Combined_Session','Stat_files');                               % directory where datafile should be saved after Statistics/activation mapping
    INFO.file.info_dir      = fullfile(INFO.file.dir,INFO.dataselect.subjectnow,INFO.dataselect.taskname,'Combined_Session','INFO_files');                               % directory where datafile should be saved after Statistics/activation mapping
    INFO.file.contrast_value_dir = fullfile(INFO.file.dir,INFO.dataselect.subjectnow,INFO.dataselect.taskname,'Combined_Session','Contrast_value_files');
end

if ~exist(INFO.file.spec_dir); mkdir(INFO.file.spec_dir); end
if ~exist(pp_dir); mkdir(pp_dir); end
if ~exist(INFO.file.coor_dir); mkdir(INFO.file.coor_dir); end
if ~exist(INFO.file.stat_dir); mkdir(INFO.file.stat_dir); end
if ~exist(INFO.file.info_dir); mkdir(INFO.file.info_dir); end
if ~exist(INFO.file.contrast_value_dir); mkdir(INFO.file.contrast_value_dir); end

INFO.file.report_dir    = fullfile(INFO.file.dir,'Reports');
if ~exist(INFO.file.report_dir); mkdir(INFO.file.report_dir); end

conv_name                                   = ['',INFO.dataselect.subjectnow,'_',INFO.dataselect.taskname,'_unit_',INFO.model.units,'.mat'];
file_name                                   = ['',INFO.dataselect.subjectnow,'_',INFO.dataselect.taskname,'_',INFO.file.name,'_unit_',INFO.model.units,'.mat'];
fig_name                                    = ['',INFO.dataselect.subjectnow,'_',INFO.dataselect.taskname,'_',INFO.file.name,'.',INFO.extension,];
INFO.file.conv_name{isess}                  = fullfile(conv_dir,conv_name);
INFO.file.conv_ds_name{isess}               = fullfile(conv_ds_dir,file_name);
INFO.file.conditions_logfile{isess}         = fullfile(INFO.file.general_dir{isess},'Raw_datafiles',['conditions_',INFO.dataselect.taskname,'_',INFO.dataselect.subjectnow,'_log.txt']);
INFO.file.MARA_name{isess}                  = fullfile(INFO.file.MARA_dir{isess},file_name);
INFO.file.extract_name{isess}               = fullfile(INFO.file.extract_dir{isess},file_name);
INFO.file.grand_avg_name{isess}             = fullfile(INFO.file.grand_avg_dir{isess},file_name);

INFO.file.contrast_value_name               = fullfile(INFO.file.contrast_value_dir,file_name);
INFO.file.spec_name                         = fullfile(INFO.file.spec_dir,file_name);
INFO.file.pp_name                           = fullfile(pp_dir,file_name);
INFO.file.coor_name                         = fullfile(INFO.file.coor_dir,file_name);

%Construct condition save-files
INFO.file.conditions_adj_dir{isess}         = fullfile(INFO.file.general_dir{isess},'Adjusted_Conditions_files');
if ~exist(INFO.file.conditions_adj_dir{isess}); mkdir(INFO.file.conditions_adj_dir{isess}); end
INFO.file.conditions_adj_name{isess}        = fullfile(INFO.file.conditions_adj_dir{isess},sprintf('conditions_%s',conv_name));
INFO.file.conditions_adj_ds_dir{isess}      = fullfile(INFO.file.general_dir{isess},'Adjusted_Downsampled_Conditions_files');
if ~exist(INFO.file.conditions_adj_ds_dir{isess}); mkdir(INFO.file.conditions_adj_ds_dir{isess}); end
INFO.file.conditions_adj_ds_name{isess}     = fullfile(INFO.file.conditions_adj_ds_dir{isess},sprintf('conditions_%s',file_name));
end


