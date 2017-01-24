function [V,h,Ph,F,Fa,Fc] = spm_reml(YY,X,Q,N);
% ReML estimation of [improper] covariance components from y*y'
% FORMAT [C,h,Ph,F,Fa,Fc] = spm_reml(YY,X,Q,N);
%
% YY  - (m x m) sample covariance matrix Y*Y'  {Y = (m x N) data matrix}
% X   - (m x p) design matrix
% Q   - {1 x q} covariance components
% N   - number of samples
%
% C   - (m x m) estimated errors = h(1)*Q{1} + h(2)*Q{2} + ...
% h   - (q x 1) ReML hyperparameters h
% Ph  - (q x q) conditional precision of h
%
% F   - [-ve] free energy F = log evidence = p(Y|X,Q) = ReML objective
%
% Fa  - accuracy
% Fc  - complexity (F = Fa - Fc)
%
% Performs a Fisher-Scoring ascent on F to find ReML variance parameter
% estimates.
%
% see also: spm_reml_sc for the equivalent scheme using log-normal
% hyperpriors
%__________________________________________________________________________
%
% SPM ReML routines:
%
%      spm_reml:    no positivity constraints on covariance parameters
%      spm_reml_sc: positivity constraints on covariance parameters
%      spm_sp_reml: for sparse patterns (c.f., ARD)
%
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner & Karl Friston
% $Id: spm_reml.m 1165 2008-02-22 12:45:16Z guillaume $

% assume a single sample if not specified
%--------------------------------------------------------------------------
try, N; catch, N  = 1;  end

% default number of iterations
%--------------------------------------------------------------------------
try, K; catch, K  = 128; end

% catch NaNs
%--------------------------------------------------------------------------
W     = Q;
q     = find(all(isfinite(YY)));
YY    = YY(q,q);
for i = 1:length(Q)
    Q{i} = Q{i}(q,q);
end

% initialise h
%--------------------------------------------------------------------------
n     = length(Q{1});
m     = length(Q);
h     = zeros(m,1);
dh    = zeros(m,1);
dFdh  = zeros(m,1);
dFdhh = zeros(m,m);

% ortho-normalise X
%--------------------------------------------------------------------------
if isempty(X)
    X = sparse(n,0);
else
    X = orth(full(X(q,:)));
end

% initialise and specify hyperpriors
%==========================================================================
for i = 1:m
    h(i) = any(diag(Q{i}));
end
hE  = sparse(m,1);
hP  = speye(m,m)/exp(32);


% ReML (EM/VB)
%--------------------------------------------------------------------------
dF    = Inf;
t     = 2;

for k = 1:K
    
    % compute current estimate of covariance
    %----------------------------------------------------------------------
    C     = sparse(n,n);
    for i = 1:m
        %~original calc~ C = C + Q{i}*h(i);
        temp    = bsxfun(@times,Q{i},h(i));
        C       = bsxfun(@plus,C,temp);
        clear temp
    end
    %~original calc~ iC    = inv(C + speye(n,n)/exp(32));
    temp    = speye(n,n)/exp(32);
    temp    = bsxfun(@plus,C,temp);
    iC      = inv(temp);
    clear temp
    
    % E-step: conditional covariance cov(B|y) {Cq}
    %======================================================================
    iCX    = iC*X;
    
    if length(X)
        Cq = inv(X'*iCX);
    else
        Cq = sparse(0);
    end
    
    % M-step: ReML estimate of hyperparameters
    %======================================================================
    
    % Gradient dF/dh (first derivatives)
    %----------------------------------------------------------------------
    
    P     = iC - iCX*Cq*iCX';
    
    %~original calc~ U     = speye(n) - P*YY/N;
    temp   = bsxfun(@times,P,YY);
    temp   = temp/N;
    U      = bsxfun(@minus,speye(n),temp);
    clear temp
    
    for i = 1:m
        
        % dF/dh = -trace(dF/diC*iC*Q{i}*iC)
        %------------------------------------------------------------------
        %~original calc~ PQ{i}   = P*Q{i};
        PQ{i}   = bsxfun(@times,P,Q{i});
        %~original calc~ dFdh(i) = -trace(PQ{i}*U)*N/2;
        temp    = bsxfun(@times,PQ{i},U);
        dFdh(i) = -trace(temp)*N/2;
        clear temp
        
    end
    
    % Expected curvature E{dF/dhh} (second derivatives)
    %----------------------------------------------------------------------
    for i = 1:m
        for j = i:m
            
            % dF/dhh = -trace{P*Q{i}*P*Q{j}}
            %--------------------------------------------------------------
            %~original calc~ dFdhh(i,j) = -trace(PQ{i}*PQ{j})*N/2;
            temp        = bsxfun(@times,PQ{i},PQ{j});
            dFdhh(i,j)  = -trace(temp)*N/2;
            clear temp
            
            dFdhh(j,i)  =  dFdhh(i,j);
        end
    end
    
    % add hyperpriors
    %----------------------------------------------------------------------
    e     = h     - hE;
    dFdh  = dFdh  - hP*e;
    dFdhh = dFdhh - hP;
    
    % Fisher scoring: update dh = -inv(ddF/dhh)*dF/dh
    %----------------------------------------------------------------------
    dh    = spm_dx(dFdhh,dFdh)/log(k + 2);
    h     = h + dh;
    
    % Convergence (1% change in log-evidence)
    %======================================================================
    
    % update regulariser
    %----------------------------------------------------------------------
    dF    = dFdh'*dh;
    fprintf('%-30s: %i %30s%e\n','  ReML Iteration',k,'...',full(dF));
    
    % final estimate of covariance (with missing data points)
    %----------------------------------------------------------------------
    if dF < 1e-1
        V     = 0;
        for i = 1:m
            if i==1
                V = V + W{i}*h(i);
            else
                %~original calc~ V = V + W{i}*h(i);
                temp=W{i}*h(i);
                V=bsxfun(@plus,V,temp);
                clear temp
            end
        end
        break
    end
    
end

% log evidence = ln p(y|X,Q) = ReML objective = F = trace(R'*iC*R*YY)/2 ...
%--------------------------------------------------------------------------
Ph    = -dFdhh;
if nargout > 3
    
    % tr(hP*inv(Ph)) - nh + tr(pP*inv(Pp)) - np (pP = 0)
    %----------------------------------------------------------------------
    Ft = trace(hP*inv(Ph)) - length(Ph) - length(Cq);
    
    % complexity - KL(Ph,hP)
    %----------------------------------------------------------------------
    Fc = Ft/2 + e'*hP*e/2 + spm_logdet(Ph*inv(hP))/2 - N*spm_logdet(Cq)/2;
    
    % Accuracy - ln p(Y|h)
    %----------------------------------------------------------------------
    Fa = Ft/2 - trace(C*P*YY*P)/2 - N*n*log(2*pi)/2 - N*spm_logdet(C)/2;
    
    % Free-energy
    %----------------------------------------------------------------------
    F  = Fa - Fc;
    
end
