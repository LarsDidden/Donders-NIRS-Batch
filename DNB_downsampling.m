function DNB_downsampling
% Downsampling the dataset. Three methods can be used avg_downsampling, ft_downsampling
% and ft_resampling) selected in the INFO-file.
%
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO;

isess=INFO.counter.iSess;

% set a default choice
if ~isfield(INFO.conv,'method')
    INFO.conv.method = 'ft_resample';
end

%Downsampling dataset/onsets/samplefrequency to avoid memory issues.

%%%%%%%% Skip this step if already done and INFO.overwrite == no.
if exist(INFO.file.conditions_adj_ds_name{isess}) & exist(INFO.file.conv_ds_name{isess});
    if strcmp(INFO.overwrite,'no')==1
        return
    end
end
%%%%%%%%

fprintf('## %s: running for subject %s ##\n',mfilename,INFO.dataselect.subjectnow);

load(INFO.file.conv_name{isess});
try
    load(INFO.file.conditions_adj_name{isess});
catch
    load(INFO.file.Sess(isess).conditions_dir{INFO.counter.iSubj});
end

switch INFO.conv.downfs_method
    case 'ft_downsampling'
        ds_data=ft_downsampling(nirs_data);
    case 'ft_resampling'
        ds_data=ft_resampling(nirs_data);
    otherwise
        fprintf('No existing downsampling method was selected; the default option (ft_resampling) was used');
        ds_data=ft_resampling(nirs_data);
end
nirs_data=[]; nirs_data=ds_data;

%% Downsample onsets/durations
Dsfactor=INFO.conv.fs/INFO.conv.downfs;

if strcmp(INFO.dataselect.taskname,'manikin')==1
    onsets=onsets'
end

for iOns=1:size(onsets,2)
    if strcmp(INFO.model.units,'scans')==1 %Change onsets/durations from scans to ds scans
        onsets{iOns}=onsets{iOns}/Dsfactor;
        if isempty(durations{iOns})==0 & durations{iOns}>0 %only change dur if the matrix is not empy and >0
            durations{iOns}=durations{iOns}/Dsfactor';
        end
    end
end

if strcmp(INFO.dataselect.taskname,'manikin')==1
    durations=[];
    for iOns=1:size(onsets,2)
        durations{iOns}=0;
    end
end

% update sampling frequency
if isess==INFO.sessions
    nirs_data.fs=INFO.conv.downfs; INFO.conv.fs=nirs_data.fs;
end

%Save dataset and onsets/durations
save(INFO.file.conv_ds_name{isess}, 'nirs_data');
if exist('P')
    save(INFO.file.conditions_adj_ds_name{isess}, 'onsets', 'durations', 'names','P');
else
    save(INFO.file.conditions_adj_ds_name{isess}, 'onsets', 'durations', 'names');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% downsampling functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nirs_data = ft_downsampling(nirs_data)
%
global INFO

nirs_data.oxyData = ft_preproc_resample(nirs_data.oxyData', INFO.conv.fs, INFO.conv.downfs, 'downsample')';
nirs_data.dxyData = ft_preproc_resample(nirs_data.dxyData', INFO.conv.fs, INFO.conv.downfs, 'downsample')';


function nirs_data = ft_resampling(nirs_data)

global INFO

nirs_data.oxyData = ft_preproc_resample(nirs_data.oxyData', INFO.conv.fs, INFO.conv.downfs, 'resample')';
nirs_data.dxyData = ft_preproc_resample(nirs_data.dxyData', INFO.conv.fs, INFO.conv.downfs, 'resample')';

