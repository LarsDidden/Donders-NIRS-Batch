function DNB_pre_loopscripts
% Changes some variables supplied in the INFO-file to match with the
% required format. Also constructing some variables required for the batch.
%
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Construct name for saving files %
hrftype =       num2str(INFO.model.hrf_type);
index_HPF = find(INFO.filt.HPF == ','); if isempty(index_HPF) ==1; HPFname=INFO.filt.HPF; else; HPFname= INFO.filt.HPF(1:(index_HPF-1)); end
index_LPF = find(INFO.filt.LPF == ','); if isempty(index_LPF) ==1; LPFname=INFO.filt.LPF; else; LPFname= INFO.filt.LPF(1:(index_LPF-1)); end
if strcmp(INFO.filt.MARA,'yes')==1; MARA_name='MARA'; else; MARA_name='noMARA'; end
if strcmp(INFO.conv.downfs_method,'ft_downsampling')==1; dsname='Ft_ds';
elseif strcmp(INFO.conv.downfs_method,'ft_resampling')==1; dsname='Ft_rs';
end
dsname2=[,dsname,'(',num2str(INFO.conv.downfs),')'];%Downsampling method + rate name
INFO.file.name                          = ['',INFO.model.hb,'_',hrftype,'der_',dsname2,'_',HPFname,'_',LPFname,'_',MARA_name,'',];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add default value to downfs if necessary
if isempty(INFO.conv.downfs) == 1
    INFO.conv.downfs='ft_resampling';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
