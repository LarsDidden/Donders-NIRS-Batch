function [argout] = DNB_spm_filter_HPF_LPF_WMDL(K,Y)
% Removes low frequency confounds X0
% FORMAT [Y] = spm_filter(K,Y)
% FORMAT [K] = spm_filter(K)
%
% K           - filter matrix or:
% K(s)        - struct array containing partition-specific specifications
%
% K(s).RT     - observation interval in seconds
% K(s).row    - row of Y constituting block/partition s
% K(s).HParam - cut-off period in seconds
% K(s).LPF_type - The shape of LPF filter which can be either Gaussian or
%                 hrf. Difference between these two are slight but hrf may
%                 provide a better sensitivity for event-related data
%                 modelled using a hrf-basis function.
% K(s).LParam - FWHM of Gaussian
% K(s).X0     - low frequencies to be removed (DCT)
%
% Y           - data matrix
%
% K           - filter structure
% Y           - filtered data
%___________________________________________________________________________
%
% spm_filter implements high-pass filtering in an efficient way by
% using the residual forming matrix of X0 - low frequency confounds
%.spm_filter also configures the filter structure in accord with the
% specification fields if called with one argument
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_filter.m 184 2005-05-31 13:23:32Z john $
% updated 09-04-28
% this code has been updated for NIRS-SPM and wavelet-MDL detrending
%
% Adjusted to the DNB by:
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO;

% set or apply
%---------------------------------------------------------------------------
if nargin == 1 && isstruct(K)
    
    % set K.X0
    %-------------------------------------------------------------------
    for s = 1:length(K)
        k = length(K(s).row);
        % low pass filter
        switch K(s).LParam.type
            case 'hrf'
                h = spm_hrf(K(s).RT);
                h = [h; zeros(size(h))];
                g = abs(fft(h));
                h = real(ifft(g));
                h = fftshift(h)';
                n = length(h);
                d = [1:n] - n/2 -1;
                K(s).KL = spdiags(ones(k,1)*h, d, k,k);
                K(s).KL = spdiags(1./sum(K(s).KL')',0,k,k)*K(s).KL;
            case 'Gaussian'
                sigma   = K(s).LParam.FWHM/K(s).RT;
                h       = round(4*sigma);
                h       = exp(-[-h:h].^2/(2*sigma^2));
                n       = length(h);
                d       = [1:n] - (n + 1)/2;
                if      n == 1, h = 1; end
                K(s).KL = spdiags(ones(k,1)*h, d, k,k);
                K(s).KL = spdiags(1./sum(K(s).KL')',0,k,k)*K(s).KL;
        end
        % high pass filter
        switch K(s).HParam.type
            case 'DCT'
                k       = length(K(s).row);
                n       = fix(2*(k*K(s).RT)/K(s).HParam.M + 1);
                if strcmp(INFO.filt.HPF,'detrend')==1               % If HPF_type=detrend, run 'DCT', but nullify it by changing X0 to zeros matrix.
                    X0  = zeros(k,n);
                else
                    X0  = spm_dctmtx(k,n);
                end
                K(s).X0 = X0(:,2:end);
        end
    end
    
    % return structure
    %-------------------------------------------------------------------
    argout = K;
    
else
    % apply filters (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
    %-------------------------------------------------------------------
    if isstruct(K)
        
        for s = 1:length(K)
            
            % select data
            %---------------------------------------------------
            if INFO.counter.iFilt==2
                y = Y(K(s).row,:);
            elseif INFO.counter.iFilt==1;
                y=Y;
            end
            
            switch K(s).LParam.type
                case {'hrf', 'Gaussian'}
                    y = K(s).KL*y;
                    % This function filters the design matrix and the data itself, but only the data should be plotted.
                    % Plot and save before/after plots of every channel
                    if strcmp(INFO.plots2,'yes')==1
                        DNB_plotchannel(Y,y,'LPF')
                    end
            end
            
            y_afterLParam=y;
            
            % apply high pass filter
            %---------------------------------------------------
            switch K(s).HParam.type
                case 'Wavelet-MDL'
                    [biasM] = detrend_wavelet_MDL(full(y), K(s).X(:,1:end-1));
                    % ~original calc~ y = y - biasM;
                    y = bsxfun(@minus,y,biasM);
                    % This function filters the design matrix and the data itself, but only the data should be plotted.
                    % Plot and save before/after plots of every channel
                    if strcmp(INFO.plots2,'yes')==1
                        DNB_plotchannel(y_afterLParam,y,'HPF')
                    end
                case 'DCT'
                    % ~original calc~ y = y - K(s).X0*(K(s).X0'*y);
                    temp=K(s).X0*(K(s).X0'*y);
                    y=bsxfun(@minus,y,temp);
                    clear temp
                    
                    % This function filters the design matrix and the data itself, but only the data should be plotted.
                    % Plot and save before/after plots of every channel
                    if strcmp(INFO.plots2,'yes')==1 & strcmp(INFO.filt.HPF,'DCT')==1
                        DNB_plotchannel(y_afterLParam,y,'HPF')
                    end
            end
            
            clear y_afterLParam;
            
            % reset filtered data in Y
            %---------------------------------------------------
            tic
            tmpsave=tempname; save(tmpsave,'-v7.3');
            rowidx=K(s).row;
            var=whos;
            varnames={var.name};
            varnamesnew=[];
            keepvar={'Y','y','tmpsave','rowidx','varnamesnew'};
            cntr=1;
            for i=1:size(keepvar,2)
                idxVar{i}=strfind(varnames,keepvar{i});
                for j=1:size(idxVar{i},2)
                    if ~isempty(idxVar{i}{j})==1
                        idxVar2(cntr)=j;
                        cntr=cntr+1;
                    end
                end
            end

            for i=1:size(idxVar2,2)
                     varnames{idxVar2(i)}=[];
            end
            varnamesnew = varnames(~cellfun(@isempty, varnames));
            Y2(rowidx,size(y,2))=zeros; % [JW] als je 'Y2=sparse(rowidx,size(y,2));' gebruikt heb je meteen een sparse matrix met nullen (scheelt een stap, die even heel veel geheugen kost)
            Y2=sparse(Y2);
            
            if INFO.counter.iFilt==2
                keep Y y tmpsave rowidx varnamesnew; % [JW] Waarom Y bewaren? Die gebruik je hieronder niet
               %% Methode 1: stapjesmethode: iets minder memory nodig, veel meer tijd
                stepsize=5;
                r4=round(size(rowidx,2)/stepsize)+1;
                for ir=1:(stepsize-1)
                    Y2(rowidx(1,r4*(ir-1)+1:(r4*ir)),:)=y(r4*(ir-1)+1:(r4*ir),:);
                end
                Y2(rowidx(1,r4*(stepsize-1)+1:end),:)=y(r4*(stepsize-1)+1:end,:);
                clear y
                
                %% Methode 2: alles in 1 keer overzetten: heel veel memory nodig, beste tijd
                 Y2(rowidx,:)=y;
                 clear y
                 
                load(tmpsave,varnamesnew{:}); % [JW] Ik zou de temp file na gebruik direct deleten, anders kan dit een systeem flink doen vollopen
                toc
            end
        end
        
        if INFO.counter.iFilt==1
            clear Y;
            Y=y;
            clear y;
        else
            clear Y;
            Y=Y2;
            clear Y2;
        end
        % K is simply a filter matrix
        %-------------------------------------------------------------------
    else
        Y = K*Y;
    end
    
    % return filtered data
    %-------------------------------------------------------------------
    %if any(~finite(Y)), warning('Found non-finite values in Y (could be the data).'); end;
    argout = Y;
end

