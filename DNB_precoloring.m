function [SPM] = DNB_precoloring(SPM, Y)
% this function is used for estimation of GLM parameters.
% especially, 'precoloring' method is used for estimating temporal
% correlation.
% 'out of memory' problem is sufficiently solved.
% last update : March 8th, 2010.
%
% Adjusted to the DNB by:
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO;


[nScan nBeta] = size(SPM.xX.X);

% Design space and projector matrix [pseudoinverse] for WLS
%==================================================================

switch SPM.xX.K(1).HParam.type
    case 'Wavelet-MDL'
        INFO.plots2='no' %Set plots to no because plots would be incomplete at this point; they should only be plotted later on in the script.
        tmp_K = SPM.xX.K;
        tmp_K.HParam.type = '';
        SPM.xX.xKXs = spm_sp('Set', DNB_spm_filter_HPF_LPF_WMDL(tmp_K, SPM.xX.X)); % KX
        SPM.xX.K.X = SPM.xX.X;
        clear tmp_K;
        if strcmp(INFO.plots,'yes')==1; INFO.plots2='yes'; end
    case 'DCT'
        INFO.plots2='no';
        SPM.xX.xKXs = spm_sp('Set', DNB_spm_filter_HPF_LPF_WMDL(SPM.xX.K, SPM.xX.X)); % KX
        if strcmp(INFO.plots,'yes')==1; INFO.plots2='yes'; end
end

SPM.xX.xKXs.X = full(SPM.xX.xKXs.X);
SPM.xX.pKX = spm_sp('x-', SPM.xX.xKXs); % projector

KYn=[];
for isess=1:size(SPM.xX.K,2)
    
    load(SPM.nirs.fname{isess});
    
    
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
    
    % SPM_filtering: hrf or gaussian, wavelet or DCT
    try
        if strcmp(nirs_data.cL.type, 'none') == 0 & strcmp(nirs_data.cH.type, 'none') == 0
            KY{isess} = Y{isess};
        elseif strcmp(nirs_data.cL.type, 'none') == 0 & strcmp(nirs_data.cH.type, 'none') == 1
            K_tmp = SPM.xX.K(isess);
            K_tmp.LParam.type = 'none';
            KY{isess} = DNB_spm_filter_HPF_LPF_WMDL(K_tmp, Y{isess});
        elseif strcmp(nirs_data.cL.type, 'none') == 1 & strcmp(nirs_data.cH.type, 'none') == 0
            K_tmp = SPM.xX.K(isess);
            K_tmp.HParam.type = 'none';
            KY{isess} = DNB_spm_filter_HPF_LPF_WMDL(K_tmp, Y{isess});
        end
        clear K_tmp;
    catch
        KY{isess} = DNB_spm_filter_HPF_LPF_WMDL(SPM.xX.K(isess), Y{isess});
    end
    clear nirs_data;
    
    % Simple detrending
    if strcmp(INFO.filt.HPF,'detrend')==1
        y=DNB_simpledetrend(KY{isess});
        % Plot and save before/after plots of every channel
        if strcmp(INFO.plots,'yes')==1
            DNB_plotchannel(Y{isess},y,'HPF')
        end
        KY{isess}=y; clear y
    end
    % end
    
    KYn=[KYn; KY{isess}];
end

SPM.nirs.beta = SPM.xX.pKX * KYn; % beta : least square estimate
res = spm_sp('r', SPM.xX.xKXs, KYn); % Residuals

% update for calculating the channel-wise least-square residual correlation
% date: Aug 10, 2011

ResSS = (KYn' * KYn) - (KYn' * SPM.xX.xKXs.X) * SPM.xX.pKX * KYn;
SPM.nirs.ResSS = ResSS;
SPM.nirs.res = res;
% end of update

% for ichan=1:size(KY,2)      %op 0 leggen signaal
%     verschil=KY(1,ichan);
%     for idata=1:size(KY,1)
%         KY(idata,ichan)=KY(idata,ichan)-verschil;
%     end
% end

SPM.xX.KY = KYn; % store filtered data

clear KY KYn;
clear Y;
for isess=1:size(SPM.xX.K,2)
    switch SPM.xX.K(isess).LParam.type
        case {'hrf', 'Gaussian'}
            S = SPM.xX.K.KL;
        case 'none'
            S = speye(nScan);
    end
    
    switch SPM.xX.K(isess).HParam.type
        case 'DCT'
            try
                S = S - SPM.xX.K(isess).X0 * (SPM.xX.K(isess).X0' * S);
            catch % in case that 'out of memory' appears
                % writing S matrix as text file
                str = [];
                for kk = 1:nScan
                    str = [str '%g '];
                end
                str = [str '\n'];
                fid = fopen('S.txt','wt');
                fprintf(fid, str, full(S));
                fclose(fid);
                clear fid;
                
                var1 = SPM.xX.K(isess).X0 * (SPM.xX.K(isess).X0' * S);
                fid = fopen('S2.txt', 'wt');
                fprintf(fid, str, full(var1));
                fclose(fid);
                clear var1;
                clear fid;
                
                S = [];
                fid1 = fopen('S.txt', 'rt');
                fid2 = fopen('S2.txt', 'rt');
                
                %h_wait = waitbar(0, 'Generating the filtering matrix... Please wait.');
                disp('Generating a filtering matrix... Please wait.');
                for kk = 1:nScan
                    %   waitbar(kk/nScan, h_wait);
                    S = [S; str2num(fgetl(fid1)) - str2num(fgetl(fid2))];
                end
                disp('Completed.');
                %close(h_wait);
                fclose(fid1);
                fclose(fid2);
                S = S';
                delete('S.txt');
                delete('S2.txt');
            end
    end
    
    % if the 'out of memory' happens,
    % new code for calculating the trace(RV) & trace(RVRV)
    trRV = 0;
    trRVRV = 0;
    disp('Calculating the trace(RV) & trace(RVRV)... Please wait.');
    h_wait = waitbar(0, 'Calculating the trace(RV) & trace(RVRV)... Please wait.');
    
    % Calculate traces
    tic
    if isempty(INFO.conv.stepsizeTr)==1;
        
        % Determine stepsize
        stepflag = true;
        stepA=[1,2,4,8,16,32];
        iStep=0;
        
        while stepflag                              %DNB tries to perform this calculation with the biggest stepsize possible. If not possible, it automatically takes a smaller stepsize.
            try
                iStep=iStep+1;
                stepsize=round(nScan/stepA(iStep));
                for kk = 1:stepsize:nScan
                    waitbar(kk/nScan, h_wait);
                    if kk+stepsize-1 <= size(S,2)
                        var1 = S * S(kk:kk+stepsize-1,:)';
                    else
                        var1 = S * S(kk:end,:)';
                    end
                    var2 = var1 - SPM.xX.xKXs.X * (SPM.xX.pKX * var1); % RV * e(kk)
                    var3 = S * (S' * var2);
                    var4 = var3 - SPM.xX.xKXs.X * (SPM.xX.pKX * var3);
                    if kk+stepsize-1 <= size(S,2)
                        idx = sub2ind(size(var2),kk:kk+stepsize-1,1:stepsize);
                    else
                        idx = sub2ind(size(var2),kk:size(var2,1),1:(size(var2,1)-kk+1));
                    end
                    trRV = trRV + sum(var2(idx));
                    trRVRV = trRVRV + sum(var4(idx));
                end
                stepflag=false;
            end
        end
        close(h_wait);
        disp('Completed.');
    else
        stepsize = round(INFO.conv.stepsizeTr);       %DNB performs the calculation with a stepsize determined in the INFO-file.
        for kk = 1:stepsize:nScan
            waitbar(kk/nScan, h_wait);
            if kk+stepsize-1 <= size(S,2)
                var1 = S * S(kk:kk+stepsize-1,:)';
            else
                var1 = S * S(kk:end,:)';
            end
            keyboard
            var2 = var1 - SPM.xX.xKXs.X * (SPM.xX.pKX * var1); % RV * e(kk)
            var3 = S * (S' * var2);
            var4 = var3 - SPM.xX.xKXs.X * (SPM.xX.pKX * var3);
            if kk+stepsize-1 <= size(S,2)
                idx = sub2ind(size(var2),kk:kk+stepsize-1,1:stepsize);
            else
                idx = sub2ind(size(var2),kk:size(var2,1),1:(size(var2,1)-kk+1));
            end
            trRV = trRV + sum(var2(idx));
            trRVRV = trRVRV + sum(var4(idx));
        end
        close(h_wait);
        disp('Completed.');
    end
    toc
end

SPM.xX.trRV = trRV; % <R'*y'*y*R>
SPM.xX.trRVRV = trRVRV; %- Satterthwaite
SPM.xX.erdf = trRV^2/trRVRV; % effective degrees of freedom
% SPM.xX.Bcov = SPM.xX.pKX*V*SPM.xX.pKX';
SPM.xX.Bcov = (SPM.xX.pKX * S);
SPM.xX.Bcov = SPM.xX.Bcov * SPM.xX.Bcov';
SPM.nirs.step = 'estimation';

try
    K = SPM.xX.K;
    K = rmfield(K, 'X');
    SPM.xX.K = K;
    clear K;
end

disp('Model parameter estimation has been completed');

