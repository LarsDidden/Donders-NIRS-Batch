function DNB_name_group
%Constructing file names and directories to use when writing away files during the group statistics.
%
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO

fprintf('## %s: running for group %s ##\n',mfilename,INFO.dataselect.analysisname);

INFO.file.sec_dir=fullfile(INFO.file.dir,'Group Analysis');
if ~exist(INFO.file.sec_dir); mkdir(INFO.file.sec_dir); end

%Figure maps
INFO.file.sec_brainfig=fullfile(INFO.file.sec_dir,'Figures',INFO.dataselect.analysisname,'Brain maps');
if ~exist(INFO.file.sec_brainfig); mkdir(INFO.file.sec_brainfig); end
INFO.file.sec_config=fullfile(INFO.file.sec_dir,'Figures',INFO.dataselect.analysisname,'Contrast value images');
if ~exist(INFO.file.sec_config); mkdir(INFO.file.sec_config); end
INFO.file.sec_gafig=fullfile(INFO.file.sec_dir,'Figures',INFO.dataselect.analysisname,'Grand Average plots');
if ~exist(INFO.file.sec_gafig); mkdir(INFO.file.sec_gafig); end


%File maps
INFO.file.sec_confile=fullfile(INFO.file.sec_dir,'Contrast_value_files',INFO.dataselect.analysisname);
if ~exist(INFO.file.sec_confile); mkdir(INFO.file.sec_confile); end
INFO.file.sec_gafile=fullfile(INFO.file.sec_dir,'Grand_average_files',INFO.dataselect.analysisname);
if ~exist(INFO.file.sec_gafile); mkdir(INFO.file.sec_gafile); end
INFO.file.sec_estimationfile = fullfile(INFO.file.sec_dir,'Estimation_batch_files',INFO.dataselect.analysisname);
if ~exist(INFO.file.sec_estimationfile); mkdir(INFO.file.sec_estimationfile); end

% Existing files used to gain information on the analysis
file_name=['',INFO.sec.subjects{1},'_',INFO.dataselect.taskname,'_',INFO.file.name,'_unit_',INFO.model.units,'.mat'];
if INFO.sessions>1
    INFO.file.conditions_adj_ds_name = fullfile(INFO.file.dir,INFO.sec.subjects{1},INFO.dataselect.taskname,INFO.sess(1).name,'Adjusted_Downsampled_Conditions_files',sprintf('conditions_%s',file_name));
else
    INFO.file.conditions_adj_ds_name = fullfile(INFO.file.dir,INFO.sec.subjects{1},INFO.dataselect.taskname,'Adjusted_Downsampled_Conditions_files',sprintf('conditions_%s',file_name));
end
% report file location
INFO.file.report_dir=fullfile(INFO.file.dir,'Reports');
if ~exist(INFO.file.report_dir); mkdir(INFO.file.report_dir); end

DNB_construct_contrast_vector

file_name_stat=[INFO.file.name,'_',INFO.stat.con(1).fname];
if INFO.sessions>1
    INFO.file.stat_name                  = fullfile(INFO.file.dir,INFO.sec.subjects{1},INFO.dataselect.taskname,'Combined_Session','Stat_files',file_name_stat);                               % directory where datafile should be saved after SPM modeling
else
    INFO.file.stat_name                  = fullfile(INFO.file.dir,INFO.sec.subjects{1},INFO.dataselect.taskname,'Stat_files',file_name_stat);                               % directory where datafile should be saved after SPM modeling
end
for iCon=1:numel(INFO.stat.con)
    for iView=1:numel(INFO.stat.spec_hemi)
        % Create map for group estimation
        INFO.file.sec_stat_name{iCon}{iView}     = [INFO.dataselect.taskname,'_',INFO.file.name,'_',INFO.stat.con(iCon).fname,'_',INFO.stat.spec_hemi{iView}];
        INFO.file.sec_stat_dir{iCon}{iView}      = fullfile(INFO.file.sec_dir,'Stat_files',INFO.dataselect.analysisname,INFO.file.sec_stat_name{iCon}{iView});
        if ~exist(INFO.file.sec_stat_dir{iCon}{iView}); mkdir(INFO.file.sec_stat_dir{iCon}{iView}); end
    end
end

