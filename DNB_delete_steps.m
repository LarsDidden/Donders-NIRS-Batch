function DNB_delete_steps
% Deleting all mat-files made before the last step of the DNB.
%
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO

if strcmp(INFO.firstsecond,'second')==0
    
    for isess=1:INFO.sessions
        delete(INFO.file.conv_name{isess});            %data conversion files
        delete(INFO.file.conditions_adj_name{isess});  %data conversion conditions files
        delete(INFO.file.conv_ds_name{isess});         %downsampling files
        delete(INFO.file.MARA_name{isess});            %MARA files
    end
    
    delete(INFO.file.spec_name);            %model specification files
    delete(INFO.file.pp_name);              %Model estimation files
end
