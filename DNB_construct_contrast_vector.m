function DNB_construct_contrast_vector

global INFO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change illegal characters in contrastnames        %
for iCon=1:numel(INFO.stat.con)
    if isempty(strfind(INFO.stat.con(iCon).name,'>'))==0;
        xname=strfind(INFO.stat.con(iCon).name,'>');
    elseif isempty(strfind(INFO.stat.con(iCon).name,'<'))==0;
        xname=strfind(INFO.stat.con(iCon).name,'<');
    else
        xname=[];
    end
    if isempty(strfind(INFO.stat.con(iCon).name,'_'))==0;
        yname=strfind(INFO.stat.con(iCon).name,'_');
    else
        yname=[];
    end
    INFO.stat.con(iCon).fname=INFO.stat.con(iCon).name;
    if isempty(xname)==0; INFO.stat.con(iCon).fname(xname)='_'; end
    if isempty(yname)==0; INFO.stat.con(iCon).name(yname)=' '; end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(INFO.firstsecond,'second')==0
    load(INFO.file.pp_name);
else
    file_name                        = [INFO.file.name,'_',INFO.stat.con(iCon).fname,'.mat'];
    if INFO.sessions>1
        pp_dir=fullfile(INFO.file.dir,INFO.sec.subjects{1},INFO.dataselect.taskname,'Combined_Session','Stat_files',file_name);
    else
        pp_dir=fullfile(INFO.file.dir,INFO.sec.subjects{1},INFO.dataselect.taskname,'Stat_files',file_name);
    end
    load(pp_dir);
end

% Load in conditions-file
if strcmp(INFO.firstsecond,'second')==0
    load(INFO.file.conditions_adj_ds_name{1});
else
    load(INFO.file.conditions_adj_ds_name);
end
% Check if a parametric modulation is analyzed.
for iCon=1:size(INFO.stat.con,2)
    if exist('P','var')==1
        for iP=1:numel(P)
            for iConparts=1:size(INFO.stat.con(iCon).positive,2)
                if isempty(P(iP).name)==0
                    Parmodp(iCon).parmodname{iP,iConparts} = strfind(INFO.stat.con(iCon).positive{iConparts},P(iP).name);
                else
                    Parmodp(iCon).parmodname{iP,iConparts} = [];
                end
            end
            for iConparts=1:size(INFO.stat.con(iCon).negative,2)
                if isempty(P(iP).name)==0
                    Parmodn(iCon).parmodname{iP,iConparts} = strfind(INFO.stat.con(iCon).negative{iConparts},P(iP).name);
                else
                    Parmodn(iCon).parmodname{iP,iConparts} = [];
                end
            end
        end
        
        IndexEmptyp{iCon}= find(not(cellfun('isempty', Parmodp(iCon).parmodname)));
        IndexEmptyn{iCon}= find(not(cellfun('isempty', Parmodn(iCon).parmodname)));
        
        if ~isempty(IndexEmptyp{iCon})==1 || ~isempty(IndexEmptyn{iCon})==1
            INFO.stat.con(iCon).parmod='yes';
        else
            INFO.stat.con(iCon).parmod='no';
        end
    else
        INFO.stat.con(iCon).parmod='no';
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructing a contrast vector based on the contrast names given in the INFO structure

for iCon=1:size(INFO.stat.con,2)
    Conpartspos=size(INFO.stat.con(iCon).positive,2);
    for iConparts=1:Conpartspos
        conname=INFO.stat.con(iCon).positive{iConparts};
        forbiddenchar1=regexpi(conname,'\('); %Search if '(' is in contrast-name. If so, place a backward slash before it, so regexpi still works.
        if isempty(forbiddenchar1)==0; conname2=[conname(1:(forbiddenchar1-1)),'\',conname(forbiddenchar1:end)]; else; conname2=conname; end
        forbiddenchar2=regexpi(conname2,'\)'); %Search if ')' is in contrast-name. If so, place a backward slash before it, so regexpi still works.
        if isempty(forbiddenchar2)==0; conname3=[conname2(1:(forbiddenchar2-1)),'\',conname2(forbiddenchar2:end)]; else; conname3=conname2; end
        INFO.stat.con(iCon).positive{iConparts}=conname3; clear conname conname2 conname3 forbiddenchar1 forbiddenchar2
        IndexP = regexpi(SPM_nirs.xX.name, INFO.stat.con(iCon).positive{iConparts}); %Search for names given in INFO
        IndexPos{iConparts} = find(not(cellfun('isempty', IndexP)));
        IndexPbas = regexpi(SPM_nirs.xX.name(IndexPos{1}),'\*bf\(1\)');              %Dont use second/third basis function.
        IndexPosbas{iConparts} = find(cellfun('isempty', IndexPbas));
        IndexPos{iConparts}(IndexPosbas{1})=[];
        IndexPosPar=[]; j=1;
        if strcmp(INFO.stat.con(iCon).parmod,'no')==1
            for i=1:size(IndexPos{iConparts},2)
                if strfind(SPM_nirs.xX.name{IndexPos{iConparts}(i)},'^')>0;              % Dont use parametric modulation contrast.
                    IndexPosPar(j)=i;
                    j=1+1;
                end
            end
            
            IndexPos{iConparts}(IndexPosPar)=[];
        end
    end
    
    Conpartsneg=size(INFO.stat.con(iCon).negative,2);
    for iConparts=1:Conpartsneg
        conname=INFO.stat.con(iCon).negative{iConparts};
        forbiddenchar1=regexpi(conname,'\('); %Search if '(' is in contrast-name. If so, place a backward slash before it, so regexpi still works.
        if isempty(forbiddenchar1)==0; conname2=[conname(1:(forbiddenchar1-1)),'\',conname(forbiddenchar1:end)]; else; conname2=conname; end
        forbiddenchar2=regexpi(conname2,'\)'); %Search if ')' is in contrast-name. If so, place a backward slash before it, so regexpi still works.
        if isempty(forbiddenchar2)==0; conname3=[conname2(1:(forbiddenchar2-1)),'\',conname2(forbiddenchar2:end)]; else; conname3=conname2; end
        INFO.stat.con(iCon).negative{iConparts}=conname3; clear conname conname2 conname3 forbiddenchar1 forbiddenchar2
        IndexN = regexpi(SPM_nirs.xX.name, INFO.stat.con(iCon).negative{iConparts}); %Search for names given in INFO (first basis function only).
        IndexNeg{iConparts} = find(not(cellfun('isempty', IndexN)));
        IndexNbas = regexpi(SPM_nirs.xX.name(IndexNeg{1}),'\*bf\(1\)');              %Dont use second/third basis function.
        IndexNegbas{iConparts} = find(cellfun('isempty', IndexNbas));
        IndexNeg{iConparts}(IndexNegbas{1})=[];
        IndexNegPar=[]; j=1;
        if strcmp(INFO.stat.con(iCon).parmod,'no')==1
            for i=1:size(IndexNeg{iConparts},2)
                if strfind(SPM_nirs.xX.name{IndexNeg{iConparts}(i)},'^')>0; % Dont use parametric modulation contrast.
                    IndexNegPar(j)=i;
                    j=1+1;
                end
            end
            IndexNeg{iConparts}(IndexNegPar)=[];
        end
    end
    
    INFO.stat.con(iCon).vec=zeros(size(SPM_nirs.xX.name,2),1);
    for iConparts=1:size(INFO.stat.con(iCon).positive,2)
        INFO.stat.con(iCon).vec(IndexPos{iConparts})=INFO.stat.con(iCon).vec(IndexPos{iConparts})+Conpartsneg;
    end
    for iConparts=1:size(INFO.stat.con(iCon).negative,2)
        INFO.stat.con(iCon).vec(IndexNeg{iConparts})=INFO.stat.con(iCon).vec(IndexNeg{iConparts})-Conpartspos;
    end
    
    if sum(INFO.stat.con(iCon).vec)~=0 & (strcmp(INFO.stat.con(iCon).positive,'none')==0 || strcmp(INFO.stat.con(iCon).negative,'none')==0)
        warning('The sum of the contrast vector is not equal to zero. Make sure all contrast vectors are entered correctly.')
    end
    
    % Make contrast vector without zeroes of extra basis functions and parametric modulations.
    INFO.stat.con(iCon).vecold=INFO.stat.con(iCon).vec;
    if INFO.model.hrf_type==1
        INFO.stat.con(iCon).vecold(2:2:size(INFO.stat.con(iCon).vec))=[];
    elseif INFO.model.hrf_type==2
        INFO.stat.con(iCon).vecold(2:3:size(INFO.stat.con(iCon).vecold))=[];
        INFO.stat.con(iCon).vecold(2:2:size(INFO.stat.con(iCon).vecold))=[];
    end
    if strcmp(INFO.stat.con(iCon).parmod,'no')==1
        cntr=1; cntr2=0;
        try
            for iP=1:numel(P)
                if isempty(P(iP).P)==0
                    if INFO.sessions==1
                        INFO.stat.con(iCon).vecold(iP+1)=[];
                    else
                        for iSess=1:INFO.sessions
                            remcon(cntr)=(iP+iSess)+((iSess-1)*numel(P))+cntr2;
                            cntr=cntr+1;
                        end
                    end
                    cntr2=cntr2+1;
                end
            end
            if INFO.sessions~=1
                INFO.stat.con(iCon).vecold(remcon)=[];
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


