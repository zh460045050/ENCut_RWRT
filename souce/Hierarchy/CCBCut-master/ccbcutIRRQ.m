function [Y,cost,t,gamma,e] = ccbcutIRRQ(W,k,Y,options,defaultopts)
% ccbcutIRRQ: solve for CCBCut embeddings using IRRQ algorithm

% author: Nathan D. Cahill
% email: nathan.cahill@rit.edu
% date: 21 February 2018

% construct degree matrix
d = sum(W,2);
D = diag(d);

% get option settings
displayFlag = ccbcutGet(options,'Display',defaultopts,'fast');
balanceType = ccbcutGet(options,'BalanceType',defaultopts,'fast');
tau = ccbcutGet(options,'Tau',defaultopts,'fast');
minIters = ccbcutGet(options,'MinIRRQIters',defaultopts,'fast');
maxIters = ccbcutGet(options,'MaxIRRQIters',defaultopts,'fast');
thresh.X = ccbcutGet(options,'TolIRRQX',defaultopts,'fast');
thresh.F = ccbcutGet(options,'TolIRRQF',defaultopts,'fast');
thresh.L1 = ccbcutGet(options,'TolIRRQL1',defaultopts,'fast');
thresh.L1Curv = ccbcutGet(options,'TolIRRQL1Curv',defaultopts,'fast');
e.n = ccbcutGet(options,'Epsilon',defaultopts,'fast');
e.d = ccbcutGet(options,'EpsilonDenom',defaultopts,'fast');
alpha = ccbcutGet(options,'IRRQAlpha',defaultopts,'fast');
constr = ccbcutGet(options,'IRRQConstraint',defaultopts,'fast');

% construct normalization matrix
switch balanceType
    case 'normalized'
        Pi = D;
    case 'ratio'
        Pi = speye(size(W));
end

% find indices of nonzero weights in original adjacency matrix
[i,j,wVals] = find(W);
nzW = numel(i);
ut = i<j;
i = i(ut);
j = j(ut);
wVals = wVals(ut);

% initialize cost vector and timing array
cost = nan(maxIters+1,2);
t = nan(maxIters,1);

% initialize gamma
gamma = ones(nzW,1);

% compute initial costs
costJ = trace(Y'*(D-W)*Y)/2;
costLtau = computeLtauCost(W,Y,tau);
cost(1,:) = [costJ costLtau];

% IRRQ Algorithm
currentIter = 1;

% if maximum iterations is zero, stop here
if ccbcutGet(options,'MaxIRRQIters',defaultopts,'fast')==0
    return;
end

doItAgain = true;
firstIter = true;
while doItAgain
    
    doItAgain = false;
    
    % do one iteration of IRRQ function
    [YNew,gammaNew,~,costJNew,costLtauNew,t(currentIter)] = IRRQiteration(W,k,Y,tau,Pi,e,options,defaultopts);
    if length(YNew) == 2
        continue;
    end
    cost(currentIter+1,:) = [costJNew costLtauNew];
    
    % compute some tolerances
    tol.X = norm(YNew-Y);
%     tol.F = abs(costJNew - costJ)/(1 + abs(costJ));
%     tol.L1 = abs(costLtauNew - costLtau)/(1 + abs(costLtau));
    tol.F = abs(costJNew - costJ)/abs(costJ);
    tol.L1 = abs(costLtauNew - costLtau)/abs(costLtau);
    switch displayFlag
        case {'iter','notify'}
            fprintf('Iteration: %g\n\tEmbedding Difference: %g\n\tCost (J): %g\tRel. Error (J): %g\n\tCost (Ltau): %g\tRel. Error (Ltau): %g\n',currentIter,tol.X,costJNew,tol.F,costLtauNew,tol.L1);
    end
    
    if currentIter>minIters 
        if costLtauNew>costLtau
            switch displayFlag
                case {'iter','notify','final'}
                    fprintf('\nLtau cost increase, so stopping at previous iteration.');
                    fprintf('\nTotal Iterations: %g\n\tCost (J): %g\n\tCost (Ltau): %g\n',currentIter-1,costJ,costLtau);
            end
            cost(currentIter+1,:) = [NaN NaN];
            break;
        elseif (tol.F < thresh.F) || (tol.L1 < thresh.L1)
            switch displayFlag
                case {'iter','notify','final'}
                    fprintf('\nRelative difference too small, so stopping at previous iteration.');
                    fprintf('\nTotal Iterations: %g\n\tCost (J): %g\n\tCost (Ltau): %g\n',currentIter-1,costJ,costLtau);
            end
            cost(currentIter+1,:) = [NaN NaN];
            break;
        elseif (cost(currentIter+1,2) - 2*cost(currentIter,2) + cost(currentIter-1,2)) > thresh.L1Curv*cost(1,2)
            switch displayFlag
                case {'iter','notify','final'}
                    fprintf('\nCurvature too great, so stopping at previous iteration.');
                    fprintf('\nTotal Iterations: %g\n\tCost (J): %g\n\tCost (Ltau): %g\n',currentIter-1,costJ,costLtau);
            end
            cost(currentIter+1,:) = [NaN NaN];
            break;
        end
    end
    if (tol.X < thresh.X) || (currentIter == maxIters)
        switch displayFlag
            case {'iter','notify','final'}
                fprintf('\nTotal Iterations: %g\n\tCost (J): %g\tRel. Error (J): %g\n\tCost (Ltau): %g\tRel. Error (Ltau): %g\n',currentIter,costJNew,tol.F,costLtauNew,tol.L1);
        end
    else
        doItAgain = true;
        currentIter = currentIter + 1;
    end
    
    % update value of epsilon
    if tau<=2
        switch constr
            case 'soft' % follow scheme of newer Voronin paper
                
                e.n = min(e.n,sqrt(norm(YNew-Y) + alpha^(currentIter-1)));
                if e.n <= eps
                    break
                end
                %             sprintf('min val %.2f',sqrt(norm(YNew-Y) + alpha^(currentIter-1)))
                
            case {'none','hard'} % follow scheme of original Daubechies paper
                
                if firstIter % estimate value of K to get good convergence rate
                    
                    % cluster points using current embedding
                    idx = kmeans(YNew,k+1);
                    
                    % compute number of edges in original adjacency matrix that connect
                    % different clusters
                    numEdges = 0;
                    for ii = 1:numel(i)
                        numEdges = numEdges + (idx(i(ii))~=idx(j(ii)));
                    end
                    
                    % estimate K to be halfway between this number of edges and total
                    % edges in complete graph over all vertices
                    KRatio.n = ccbcutGet(options,'KRatio',defaultopts,'fast');
                    KRatio.d = ccbcutGet(options,'KRatioDenom',defaultopts,'fast');
                    K.n = ceil(KRatio.n.*(numEdges + numel(i)));
                    %K.d = ceil(KRatio.d.*(numEdges + numel(i)));
                    
                end
                
                Z = zeros(numel(i),k);
                for m = 1:numel(i)
                    Yi = YNew(i(m),:);
                    Yj = YNew(j(m),:);
                    % Are we doing the right thing here by multiplying by
                    % gamma? Probably not!
                    % Z(m,:) = wVals(m).*gamma(m).*(Yi-Yj);
                    Z(m,:) = wVals(m).*(Yi-Yj);
                end
                %rVec = sort(sum(abs(Z),2),1,'descend');
                rVec = sort(sqrt(sum(Z.^2,2)),1,'descend');
                r = rVec(K.n+1);
                %r = rVec(1)
                %e = min(e,r/numel(i))
                e.n = min(e.n,r);
                if e.n <= eps
                    fprintf('Epsilon parameter has reached zero. Terminating.\n');
                    break
                end
                
        end
        
    end
    Y = YNew;
    gamma = gammaNew;
    
    costJ = costJNew;
    costLtau = costLtauNew;

    firstIter = false;
    
end

% change this later
e = e.n;

end

