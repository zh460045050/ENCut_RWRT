function [YNew,gammaVals,xiVals,costJ,costLtau,timing] = IRRQiteration(W,k,Y,tau,Pi,e,options,defaultopts)
% IRRQiteration: compute minimum of majorizing function

% author: Nathan D. Cahill
% email: nathan.cahill@rit.edu
% date: 21 February 2018

% start timer
timerHandle = tic;

% find indices of nonzero weights
[i,j,wVals] = find(W);

% compute weights for numerator
gammaVals = zeros(numel(i),1);
if tau<=2
    for m = 1:numel(i)
        Yi = Y(i(m),:);
        Yj = Y(j(m),:);
        %gammaVals(m) = 1/sqrt(wVals(m).*sum((Yi-Yj).^2) + e.^2);
        gammaVals(m) = (wVals(m).*sum((Yi-Yj).^2) + (e.n).^2)^((tau-2)/2);
    end
else
    for m = 1:numel(i)
        Yi = Y(i(m),:);
        Yj = Y(j(m),:);
        gammaVals(m) = (wVals(m).*sum((Yi-Yj).^2))^((tau-2)/2);
    end
end
gammaVals = gammaVals./mean(gammaVals);

% multiply by graph weights
WCap = sparse(i,j,wVals.*gammaVals,size(W,1),size(W,2));
WCap = (WCap + WCap')/2;

% degree and laplacian matrices
DCap = diag(sum(WCap,2));
LCap = DCap - WCap;
LCap = (LCap + LCap')/2;

% compute eigenvectors
opts.p = min(20*(k+1),size(W,1)-1);
opts.issym = 1;
constr = ccbcutGet(options,'IRRQConstraint',defaultopts,'fast');
xiVals = [];
switch constr
    case 'soft' % Voronin
        lambda = ccbcutGet(options,'IRRQLambda',defaultopts,'fast');
        d = sum(W,2);
        [YNew,~] = eigs(4*lambda*LCap + d*d',DCap,k+1,'SA',opts);
        % discard eigenvector corresponding to zero eigenvalue
        YNew = YNew(:,2:k+1);
        
    case 'none' % no orthogonality constraint
        deltaTinv = max(diag(LCap))/500;
        YNew = (deltaTinv*speye(size(LCap,1)) + LCap)\Y;
        YNew = YNew*(norm(Y)/norm(YNew));
        
    case 'hard' % our proposed technique
        
        % degree/size vector
        d = diag(Pi);

        % D^(-1/2)
        DNegHalf = spdiags(1./sqrt(d),0,size(W,1),size(W,2));
        
        % q vectors
        dhalf = sqrt(d);
        q = dhalf/norm(dhalf);
        q1 = q(1);
        qcap = q(2:end);

        % projection into subspace orthogonal to q
        M = [qcap';-q1*speye(size(W,1)-1)];        
        Mp = M';
        
        % solution to generalized eigenvector problem
        eigfun = @(x) geneigfun(x,LCap,DNegHalf,q1,qcap,M,Mp);
        %[H,~] = eigs(eigfun,size(W,1)-1,k,'sa',opts);
        try
            [H,~,~]=lobpcg(randn(size(W,1)-1,k),@(x) eigfun(x),[],[],1e-5,100*size(W,1));
        catch
            [H,~] = eigs(eigfun,size(W,1)-1,k,'sm');
        end
        
        % transform
        YNew = DNegHalf*M*((1/q1)*H - qcap*((qcap'*H)/(q1*(1+q1))));
        
end

% stop timer
timing = toc(timerHandle);

% compute cost function value
switch constr
    case 'soft' % Voronin
        costJ = trace(YNew'*(4*lambda*LCap + d*d')*YNew);
        costLtau = trace((YNew'*d)*(d'*YNew)) + 4*lambda*computeLtauCost(W,YNew,tau);
    case {'none','hard'} % Daubechies + ours
        costJ = trace(YNew'*LCap*YNew)/2;
        costLtau = computeLtauCost(W,YNew,tau);
end

end

function y = geneigfun(x,LCap,DNegHalf,q1,qcap,M,Mp)

% first stage
z1 = x/q1 - qcap*(qcap'*x)/(q1*(1+q1));

% second stage
z2 = Mp*(DNegHalf*(LCap*(DNegHalf*(M*z1))));

% third stage
y = z2/q1 - qcap*(qcap'*z2)/(q1*(1+q1));

end
