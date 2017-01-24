function [SPM] = DNB_prewhitening(SPM, Y, save_directory)
% this function is used for estimation of GLM parameters.
% especially, 'prewhitening' method is used for estimating temporal
% correlation.
% last update : March 8th, 2010.

global INFO;

xX            = SPM.xX;
[nScan nBeta] = size(xX.X);

% update (April 28th, 2009)
if ~isfield(xX,'K')
    xX.K  = 1;
else
    if strcmp(xX.K(1).HParam.type, 'Wavelet-MDL')
        xX.K.X = xX.X;
    end
end

xVi = SPM.xVi;

try
    V = xVi.V;
catch
    V = speye(nScan, nScan);
end

try
    %-If W is specified, use it
    %----------------------------------------------------------------------
    W     = xX.W;
catch
    if isfield(xVi,'V')
        % otherwise make W a whitening filter W*W' = inv(V)
        %------------------------------------------------------------------
        [u s] = spm_svd(xVi.V);
        s     = spdiags(1./sqrt(diag(s)),0,length(s),length(s));
        W     = u*s*u';
        W     = W.*(abs(W) > 1e-6);
        xX.W  = sparse(W);
    else
        % unless xVi.V has not been estimated - requiring 2 passes
        %------------------------------------------------------------------
        W     = speye(nScan,nScan);
    end
end

INFO.counter.iFilt=2; %Counter.iFilt decides in spm_filter script if Y is returned per session or as a whole.
% Design space and projector matrix [pseudoinverse] for WLS
%==================================================================

switch xX.K(1).HParam.type
    case 'Wavelet-MDL'
        INFO.plots2='no' %Set plots to no because plots would be incomplete at this point; they should only be plotted later on in the script.
        tmp_K = xX.K;
        tmp_K.HParam.type = '';
        xX.xKXs = spm_sp('Set', DNB_spm_filter_HPF_LPF_WMDL(tmp_K, W*xX.X)); % KWX
        if strcmp(INFO.plots,'yes')==1; INFO.plots2='yes'; end
    case 'DCT'
        INFO.plots2='no' %Set plots to no because plots would be incomplete at this point; they should only be plotted later on in the script.
        xX.xKXs = spm_sp('Set', DNB_spm_filter_HPF_LPF_WMDL(xX.K, W*xX.X)); % KWX
        if strcmp(INFO.plots,'yes')==1; INFO.plots2='yes'; end
end

xX.xKXs.X = full(xX.xKXs.X);
xX.pKX = spm_sp('x-', xX.xKXs); % projector

global defaults

%-If xVi.V is not defined compute Hsqr and F-threshold under i.i.d.
%--------------------------------------------------------------------------
if ~isfield(xVi,'V')
    Fcname = 'effects of interest';
    iX0    = [SPM.xX.iB SPM.xX.iG];
    xCon   = spm_FcUtil('Set',Fcname,'F','iX0',iX0,xX.xKXs);
    X1o    = spm_FcUtil('X1o', xCon(1),xX.xKXs);
    Hsqr   = spm_FcUtil('Hsqr',xCon(1),xX.xKXs);
    trRV   = spm_SpUtil('trRV',xX.xKXs);
    trMV   = spm_SpUtil('trMV',X1o);
    clear X1o
    
    % Threshold for voxels entering non-sphericity esimtates
    %----------------------------------------------------------------------
    try
        UFp = eval(['defaults.stats.' lower(defaults.modality) '.ufp']);
    catch
        UFp = 0.001;
    end
    UF     = spm_invFcdf(1 - UFp,[trMV,trRV]);
end

% fit model & write parameter images
% initialize variables used in the loop

S = 0; % Volumes (voxels)
s = 0; % Volumes (voxels > UF)
Cy = 0; % <Y*Y'> spatially whitened
CY = 0; % <Y*Y'> for ReML
EY = 0; % <Y> for ReML
nCh = size(Y,2);

Yn=[]; KWYn=[];
INFO.counter.iFilt=1;
for isess=1:size(xX.K,2)
    INFO.counter.iSess=isess;
    load(SPM.nirs.fname{isess});

    Yn=[Yn; Y{isess}];
    %% Actual filtering
    if strcmp(INFO.filt.LPF(1:3),'lpf')==1; %If lpf is selected, set LPF_SMP to 'none' to prevent other lpf options to run.                                     % When simple LPF is selected, the SPM_LPF settings should be set to 'none'.
        INFO.filt.LPF_SPM='none';
    else
        INFO.filt.LPF_SPM=INFO.filt.LPF;
    end
    
    % Simple low pass filtering
    if strcmp(INFO.filt.LPF(1:3),'lpf')==1
        y=DNB_simplelpf(Y{isess});
        % Plot and save before/after plots of every channel
        if strcmp(INFO.plots,'yes')==1
            DNB_plotchannel(Y{isess},y,'LPF')
        end
        Y{isess}=y; clear y;
    end
    W=speye(size(Y{isess},1));
    % SPM_filtering: hrf or gaussian, wavelet or DCT
    try
        if strcmp(nirs_data.cL.type, 'none') == 0 & strcmp(nirs_data.cH.type, 'none') == 0
            KWY{isess} = W*Y{isess};
        elseif strcmp(nirs_data.cL.type, 'none') == 0 & strcmp(nirs_data.cH.type, 'none') == 1
            K_tmp = xX.K(isess);
            K_tmp.LParam.type = 'none';
            KWY{isess} = DNB_spm_filter_HPF_LPF_WMDL(K_tmp, W*Y{isess});
        elseif strcmp(nirs_data.cL.type, 'none') == 1 & strcmp(nirs_data.cH.type, 'none') == 0
            K_tmp = xX.K(isess);
            K_tmp.HParam.type = 'none';
            KWY{isess} = DNB_spm_filter_HPF_LPF_WMDL(K_tmp, W*Y{isess});
        end
        clear nirs_data;
    catch
        KWY{isess} = DNB_spm_filter_HPF_LPF_WMDL(xX.K(isess), W*Y{isess});
    end
    clear nirs_data;
    
    % Simple detrending
    if strcmp(INFO.filt.HPF,'detrend')==1
        y=DNB_simpledetrend(KWY{isess});
        % Plot and save before/after plots of every channel
        if strcmp(INFO.plots,'yes')==1
            DNB_plotchannel(Y{isess},y,'HPF')
        end
        KWY{isess}=y; clear y
    end
    xX.KYsep{isess}=KWY{isess};
    KWYn=[KWYn; KWY{isess}];
end
xX.KY=KWYn; %store filtered data in vertical allignment
beta = xX.pKX * KWYn; % Parameter estimates
res = spm_sp('r', xX.xKXs, KWYn); % Residuals
ResSS = sum(res.^2); % Residual SSQ

clear KWY KWYn; % clear to save memory

% if ~isfield(xVi, 'V')
%     save('tmp_prewhiten.mat', 'Y','save_directory');
% end


%- if ReML hyperparameters are needed for xVi.V
if ~isfield(xVi, 'V')
    
    %-F-threshold & accumulate spatially whitened Y*Y'
    %--------------------------------------------------
    j = sum((Hsqr*beta).^2,1)/trMV > UF*ResSS/trRV;
    j = find(j);
    if length(j)
        q = size(j,2);
        s = s+q;
        q = spdiags(sqrt(trRV./ResSS(j)'),0,q,q);
        Yn  = Yn(:,j)*q;
        %~original calc Cy = Cy + Yn*Yn';
        temp=Yn*Yn';
        Cy=bsxfun(@plus,Cy,temp);
        clear temp
    end
end

% if we are saving the WLS parameters
%---------------------------------------------
if isfield(xX,'W')
    % sample covariance and mean of Y (all voxels)
    CY = CY + Yn*Yn';
    EY = EY + sum(Yn,2);
end

clear Y Yn; % clear to save memory

%-------------------------
% Post estimation cleanup
%-------------------------

% average sample covariance and mean of Y
CY = CY/nCh;
EY = EY/nCh;
CY = CY - EY*EY';

% if not defined, compute non-sphericity V using ReML Hyperparemters
if ~ isfield(xVi, 'V')
    
    %-REML estimate of residual correlations through hyperparameters (h)
    
    Cy = Cy/s;
    
    % ReML for separable designs and covariance components
    if isstruct(xX.K)
        m = length(xVi.Vi);
        h = zeros(m,1);
        V = sparse(nScan, nScan);
        for i = 1:length(xX.K)
            % extract blocks from bases
            q = xX.K(i).row;
            p = [];
            Qp= {};
            for j = 1:m
                if nnz(xVi.Vi{j}(q,q)) % nnz : number of nonzero matrix element
                    Qp{end + 1} = xVi.Vi{j}(q,q);
                    p = [p j];
                end
            end
            
            % design space for ReML (with confounds in filter)
            Xp = xX.X(q,:);
            try
                Xp = [Xp xX.K(i).X0];
            catch
            end
            
            % ReML
            [Vp, hp] = DNB_spm_reml(Cy(q,q), Xp, Qp);
            %save([save_directory '/temp2.mat'],'Vp','hp');
            %~original calc~ V(q,q) = V(q,q) + Vp;
            V(q,q)=bsxfun(@plus,V(q,q),Vp);
            h(p) = hp;
        end
    end
    clear Qp Vp Xp
    
    % normalize non-sphericity and save hyperparemters
    V = V*nScan/trace(V);
    xVi.h = h;
    xVi.V = V;
    xVi.Cy = Cy; % spatially whitened <Y*V'> (used by ReML to estimate h)
    clear Cy
    SPM.xVi = xVi; % non-sphericity structure
    
    %         if ~isfield(xX,'W')
    %             if spm_matlab_version_chk('7') >= 0
    %                 save([save_directory 'SPM.mat'],'SPM','-V6');
    %             else
    %                 save([save_directory 'SPM.mat'],'SPM');
    %             end
    %             clear
    %             load('tmp_prewhiten.mat');
    %             load([save_directory 'SPM.mat']);
    %             delete('tmp_prewhiten.mat');
    %             delete([save_directory 'SPM.mat']);
    %             [SPM]= nirs_spm_indiv(SPM, Y, save_directory);
    %             return
    %         end
end

SPM.xVi = [];
xVi.CY = CY;   %-<(Y - <Y>)*(Y - <Y>)'> (used by spm_spm_Bayes)
clear xVi;

INFO.counter.iFilt=2;
W=speye(size(xX.X,1));
% ~original calc~ W*V*W'
WVW=sparse(size(W,1),size(W,1));
temp=bsxfun(@times,W,V);
WVW=bsxfun(@times,temp,W');
clear temp W V

%V : V matrix (K*W*Vi*W'*K') = correlations after K*W is applied
switch xX.K(1).HParam.type
    case 'DCT'
        INFO.plots2='no' %Set plots to no because no filtering is done, V is made.
        % ~original calc~ V = spm_filter_HPF_LPF_WMDL(xX.K, spm_filter_HPF_LPF_WMDL(xX.K, W*V*W')'); %KWVW'K'
        Vtemp=DNB_spm_filter_HPF_LPF_WMDL(xX.K, WVW);
        clear WVW
        V = DNB_spm_filter_HPF_LPF_WMDL(xX.K,Vtemp');  %KWVW'K'
        clear Vtemp
        if strcmp(INFO.plots,'yes')==1; INFO.plots2='yes'; end
    case 'Wavelet-MDL'
        INFO.plots2='no' %Set plots to no because plots would be incomplete at this point; they should only be plotted later on in the script.
        K_tmp = xX.K;
        K_tmp.HParam.type = 'none';
        V = DNB_spm_filter_HPF_LPF_WMDL(K_tmp, DNB_spm_filter_HPF_LPF_WMDL(K_tmp, WVW)'); %KWVW'K'
        clear K_tmp;
        if strcmp(INFO.plots,'yes')==1; INFO.plots2='yes'; end
end


[trRV trRVRV] = spm_SpUtil('trRV', xX.xKXs, V); % trRV for X
xX.trRV = trRV; % <R'*y'*y*R>
xX.trRVRV = trRVRV; %- Satterthwaite
xX.erdf = trRV^2/trRVRV;
xX.Bcov = xX.pKX*V*xX.pKX';
clear V;

% compute scaled design matrix for display purpose
xX.nKX = spm_DesMtx('sca', xX.xKXs.X, xX.name);
SPM.nirs.beta  = beta;
SPM.nirs.res   = res;
SPM.nirs.ResSS = res' * res; % Residual SSQ; %%~original calc~ SPM.nirs.ResSS=ResSS. Changed because activation map needs a 'nrchannel X nrchannel' ResSS;

% update (April 28th, 2009)
try
    K = xX.K;
    K = rmfield(K, 'X');
    xX.K = K;
    clear K;
end

SPM.xX         = xX;                %-design structure
%SPM.xCon       = struct([]);        %-contrast structure
SPM.nirs.step = 'estimation';
disp('Model parameter estimation has been completed');

