function DNB_report

% Writing some information to a report-file after the analysis so the user
% may check details of the analysis.
%
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO

fprintf('## %s: running. ##\n',mfilename);

% Call in reportfiles
fnINFO=fullfile(INFO.file.DNB_dir,INFO.infofile);
report_name=[,INFO.dataselect.analysisname,'_report.txt'];
INFOreport_name=[,INFO.dataselect.analysisname,'_report_INFO.txt'];
fn_report=fullfile(INFO.file.report_dir,report_name);
fn_report_INFO=fullfile(INFO.file.report_dir,INFOreport_name);

if strcmp(INFO.firstsecond,'second')==0
    
    fid = fopen(fn_report,'wt');
    
    %% Subjects run
    title='Subjects run:';
    fprintf(fid,'%s',title);  %string
    fprintf(fid,'\n');  %new line
    for iSubj=1:size(INFO.dataselect.subjects,2);
        fprintf(fid,'%s',INFO.dataselect.subjects{iSubj});  %string
        fprintf(fid,'\n');  %new line
    end
    fprintf(fid,'\v');   %new alinea
    
    %% Subjects failed + error
    title='Subjects failed:';
    fprintf(fid,'%s',title);  %string
    fprintf(fid,'\n');  %new line
    if isempty(INFO.report.failed_subjects)==1;
        strfailed='The analysis was succesful for all subjects.';
        fprintf(fid,'%s',strfailed); %string
        fprintf(fid,'\n');  %new line
    elseif iscell(INFO.report.failed_subjects)==0; %Only 1 subject
        fprintf(fid,'%s',INFO.report.failed_subjects)
        fprintf(fid,'\n');  %new line
        fprintf(fid,'%s',INFO.report.failed_error);  %string
        fprintf(fid,'\n');  %new line
    else
        for iSubj=1:size(INFO.report.failed_subjects,2);
            fprintf(fid,'%s',INFO.report.failed_subjects{iSubj});  %string
            fprintf(fid,'\n');  %new line
            fprintf(fid,'%s',INFO.report.failed_error{iSubj});  %string
            fprintf(fid,'\n');  %new line
        end
    end
    fprintf(fid,'\v');   %new alinea
    
    %% Channels deleted beause of SCI
    title='Channels deleted after SCI check:';
    fprintf(fid,'%s',title);  %string
    fprintf(fid,'\n');  %new line
    if iscell(INFO.report.SCI)==0 %Only 1 subject
        fprintf(fid,'%s',INFO.report.SCI);
        fprintf(fid,'\n');  %new line
    else
        for iSCI=1:size(INFO.report.SCI,2)
            fprintf(fid,'%s',INFO.report.SCI{iSCI});  %string
            fprintf(fid,'\n');  %new line
        end
    end
    fprintf(fid,'\v');   %new alinea
    
    %% Channels changed after MARA
    title='Channels changed after MARA:';
    fprintf(fid,'%s',title);  %string
    fprintf(fid,'\n');  %new line
    if iscell(INFO.report.MARA)==0 %Just 1 subject is analysed
        fprintf(fid,'%s',INFO.report.MARA);
        fprintf(fid,'\n');  %new line
    else
        for iMARA=1:size(INFO.report.MARA,2)
            fprintf(fid,'%s',INFO.report.MARA{iMARA});  %string
            fprintf(fid,'\n');  %new line
        end
    end
    fprintf(fid,'\v');   %new alinea
    
end

%% INFO_file
copyfile(fnINFO,fn_report_INFO);


