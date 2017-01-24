function DNB_MNI_coord2render
%Converts the MNI coordinates of the optodes given by the user into real
%coordinates, which can be used to construct brain maps.
%
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO

%%%%%%%% Skip this step if already done and INFO.overwrite == no.
if exist(INFO.file.coor_name);
    if strcmp(INFO.overwrite,'no')==1
        return
    end
end
%%%%%%%%

fprintf('## %s: running for subject %s ##\n',mfilename,INFO.dataselect.subjectnow);
template_info = spm_vol(fullfile(spm('dir'),'templates','T1.nii'));
% Fill in location where channel information was stored.
load(INFO.file.MNI_file_name);
% Remove unreliable datachannels based on SCI check.

if strcmp(INFO.SCI.check,'yes')==1 
    mni_locations=mni_locations(:,INFO.SCI.sessremchannel{end});
end

ch_MNI_mm = [mni_locations]; % fill in coordinates
ch_MNI_mm = [ch_MNI_mm;ones(1,size(ch_MNI_mm,2))];
ch_MNI_vx = inv(template_info.mat) * ch_MNI_mm;
[rend,rendered_MNI] = render_MNI_coordinates(ch_MNI_vx,template_info);

% Store information in a file shaped like the original NIRS_SPM channel_file
for iView=1:6
    preproc_info.rend_ch_pos{iView}.rchn=rendered_MNI{iView}.rchn;
    preproc_info.rend_ch_pos{iView}.cchn=rendered_MNI{iView}.cchn;
    preproc_info.rend_ch_pos{iView}.ren=rend{iView}.ren;
end

preproc_info.ch_set.name=['Set #1'];
if strcmp(INFO.SCI.check,'yes')==1
    preproc_info.ch_set.ch=[(INFO.SCI.sessremchannel{end})];
else
    preproc_info.ch_set.ch=[1:size(mni_locations,2)];
end

save(INFO.file.coor_name,'preproc_info');