function [clusters,cuts] = computeMultiPartitioning_ccbcut(W,optsCCBCut,k,verbosity)
% Computes a multipartitioning of the data given by a similarity matrix W 
% by recursively computing bipartitions using solutions to the relaxed
% CCBCut problem.
%
% Usage: [clusters,cuts] 
%           = computeMultiPartitioning(W,optsCCBCut,k,verbosity);
%
% Input:
%   W: Sparse weight matrix. Has to be symmetric.
%   optsCCBCut: ccbcut options structure. If empty, defaults will be used.
%   k: number of clusters
%   verbosity: Controls how much information is displayed. Levels 0-3,
%   default is 2.
%
% Output:
%   clusters: mx(k-1) matrix containing in each column the computed 
%   clustering for each partitioning step.
%   cuts: (k-1)x1 vector containing the CCBCut values after each
%   partitioning step.
%
% author: Nathan D. Cahill, based on modifications of code from:
%
% (C)2010-11 Thomas Buehler and Matthias Hein
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%

% set default verbosity
if nargin<4
    verbosity = 2;
end

% set default options for CCBCut
defaultopts = struct('Display','iter','Algorithm','irrq','Tau',1,...
    'BalanceType','normalized','RoundingType','kmeans',...
    'KmeansNumStartPoints',10,'KmeansDistance','sqeuclidean',...
    'MaxSOCIters',100,'MaxL1ProbIters',100,'TolBregSOCX',1e-3,...
    'TolBregSOCF',1e-3,'TolBregL1ProbX',1e-3,'TolBregL1ProbF',1e-3,...
    'BregLambda',1000,'R',100,'MinIRRQIters',5,'MaxIRRQIters',100,'TolIRRQX',1e-3,'TolIRRQF',1e-3,...
    'TolIRRQL1',1e-3,'TolIRRQL1Curv',1e-1,'Epsilon',1,'KRatio',0.5,'IRRQConstraint','none',...
    'IRRQLambda',1,'IRRQq',1,'IRRQAlpha',0.5,'EpsilonDenom',1,...
    'KRatioDenom',0.5); 

% set up empty structure if nothing is passed in
if isempty(optsCCBCut)
    optsCCBCut = ccbcutSet('Display','iter');
end

% get tau
tau = ccbcutGet(optsCCBCut,'Tau',defaultopts,'fast');

num=size(W,1);

assert(k>=2,'Wrong usage. Number of clusters has to be at least 2.');
assert(k<=num, 'Wrong usage. Number of clusters is larger than size of the graph.');
assert(isnumeric(W) && issparse(W),'Wrong usage. W should be sparse and numeric.');
assert(sum(sum(W~=W'))==0,'Wrong usage. W should be symmetric.');
assert(isempty(find(diag(W)~=0,1)),'Wrong usage. Graph contains self loops. W has to have zero diagonal.');

threshold_type = -1;

clusters=zeros(num,k-1);
cuts=zeros(1,k-1);
cutParts=zeros(1,k);

deg=full(sum(W,2)); %% will be needed in createSubClusters2 (also in unnormalized case)
normalized = isequal(ccbcutGet(optsCCBCut,'BalanceType',defaultopts,'fast'),'normalized');
if normalized
    deg2=deg;%full(sum(W));
else
    deg2=ones(num,1);
end

cut=inf;

% Check if graph is connected
[comp,connected,sizes]=connectedComponents(W);
if(~connected)
    if(verbosity>=1), disp('WARNING! GRAPH IS NOT CONNECTED!');end
    if(verbosity>=2), disp('Optimal Cut achieved by separating connected components.');end
    allClusters = balanceConnectedComponents(comp,sizes,W,normalized);
    [cut,cutPart1,cutPart2]=deal(0);
else
    if(verbosity>=2), disp('Computing partitioning.'); end
    vmin=computeEigenvectorGeneral_ccbcut(W,optsCCBCut);
    [allClusters_temp, cut_temp, cutPart1_temp, cutPart2_temp] =  createClustersGeneral_ccbcut(vmin,W,optsCCBCut,threshold_type,deg2);
    
    % Display current objective
    if(verbosity>=2), displayCurrentObjective(cut_temp,normalized); end
    
    % Check if we're better
    if (cut_temp<cut)
        [allClusters, cut, cutPart1, cutPart2] = deal(allClusters_temp, cut_temp, cutPart1_temp, cutPart2_temp);
    end
end

allClusters=allClusters+1;
clusters(:,1)=allClusters;

cuts(:,1)=cut;

cutParts(1)=cutPart1;
cutParts(2)=cutPart2;

subCutParts=zeros(k,2);
subClusters=cell(1,k);

if(verbosity>=1)
    fprintf('Finished Clustering into 2 parts.\n');
    displayCurrentObjective(cut,normalized);
    fprintf('\n');
end

%Perform the 2nd to (k-1)th partitioning step
for l=3:k
    bestCut=inf;
    % in each step consider each of the current l-1 clusters
    for m=1:l-1
        
        index_m=find(allClusters==m);
        
        % if we have already solved this subproblem
        if (~isempty(subClusters{m}))
            allClustersInCluster_m = subClusters{m};
            cutPart1 = subCutParts(m,1);
            cutPart2 = subCutParts(m,2);
            % if the current cluster has size 1 it cannot be further divided
        elseif(length(index_m)==1)
            allClustersInCluster_m=[];
            cutPart1=inf;
            cutPart2=inf;
            subClusters{m}=allClustersInCluster_m;
            subCutParts(m,1:2)=[cutPart1 cutPart2];
        elseif(length(index_m)==2)
            allClustersInCluster_m=[0;1];
            cutPart1=deg(index_m(1))-W(index_m(1),index_m(1));
            cutPart2=deg(index_m(2))-W(index_m(2),index_m(2));
            if normalized
                if cutPart1>0
                    cutPart1=cutPart1/((2*deg(index_m(1)))^(tau/2));
                end
                if cutPart2>0
                    cutPart2=cutPart2/((2*deg(index_m(2)))^(tau/2));
                end
            end
            subClusters{m}=allClustersInCluster_m;
            subCutParts(m,1:2)=[cutPart1 cutPart2];
            % otherwise we have to compute the partition
        else
            if(verbosity>=2) fprintf('Computing partitioning of subgraph %d.\n',m);end
            
            % extract subgraph and its connected components
            Wm=W(index_m,index_m);
            size_m=size(Wm,1);
            [comp,connected,sizes]=connectedComponents(Wm);
            
            if(verbosity>=2 && ~connected) disp('...Subgraph is not connected.'); end
            cutPart1=inf;
            cutPart2=inf;
            
            % if subgraph is not connected
            if (~connected)
                
                %check partition which has connected component as one
                %cluster and rest as second cluster
                for m1=1:length(sizes)
                    %if(true)
                    if(cutPart1+cutPart2>0)
                        if(verbosity>=2) fprintf('...Checking partition found by isolating connected component %d of %d.\n',m1,length(sizes)); end
                        
                        allClustersInCluster_m_temp = double(comp==m1);
                        
                        cluster_m2=zeros(size(allClusters,1),1);
                        cluster_m2(index_m)=allClustersInCluster_m_temp;
                        cutPart2_temp = computeCutValue_ccbcut(cluster_m2,W,optsCCBCut,deg2);
                        
                        cluster_m1=zeros(size(allClusters,1),1);
                        cluster_m1(index_m)=(allClustersInCluster_m_temp==0);
                        cutPart1_temp = computeCutValue_ccbcut(cluster_m1,W,optsCCBCut,deg2);
                        
                        % Display current objective
                        if(verbosity>=2)
                            displayCurrentObjective(cutPart1_temp+cutPart2_temp,normalized);
                        end
                        
                        %Check if we're better
                        if (cutPart1_temp+cutPart2_temp<cutPart1+cutPart2)
                            [cutPart1,cutPart2,allClustersInCluster_m]=deal(cutPart1_temp,cutPart2_temp,allClustersInCluster_m_temp);
                            assert(length(allClustersInCluster_m)==length(index_m));
                        end
                        assert(logical(exist('allClustersInCluster_m','var')));
                    end
                end
            end
            %if(true)
            if(cutPart1+cutPart2>0)
                for m1=1:length(sizes)
                    index_comp=find(comp==m1);
                    % if the size of the current connected component is larger than 1, try to partition it
                    if (length(index_comp)>1)
                        Wm_comp=sparse(Wm(index_comp,index_comp));
                        Wm_comp2=Wm_comp;
                        if(2*max(sum(Wm_comp.^2))<eps)
                            Wm_comp2=Wm_comp/max(max(Wm_comp));
                        end
                        if(~connected && verbosity>=2) fprintf('...Computing partitioning of connected component %d of %d of subgraph %d.\n',m1,length(sizes), m); end
                        if (~connected)
                            index_rest=find(comp~=m1); % all other components in the current cluster
                            cut_rest=sum(sum(W(index_m(index_rest),setdiff(1:num,index_m(index_rest)))));
                            size_rest=sum(deg2(index_m(index_rest)));
                        else
                            cut_rest=0;
                            size_rest=0;
                        end
                        
                        if(verbosity>=2) fprintf('...Computing solution of relaxed CCBCut problem.\n'); end
                        vmin_comp=computeEigenvectorGeneral_ccbcut(Wm_comp2,optsCCBCut);
                        [allClustersInCluster_m_temp, cutPart1_temp, cutPart2_temp] =  createSubClusters2(vmin_comp,Wm_comp,optsCCBCut,normalized,deg,index_comp,index_m,cut_rest,size_rest,size_m);
                        
                        % Display current objective
                        if(verbosity>=2)
                            displayCurrentObjective(cutPart1_temp+cutPart2_temp,normalized);
                        end
                        
                        %Check if we're better
                        if (cutPart1_temp+cutPart2_temp<cutPart1+cutPart2)
                            [cutPart1,cutPart2,allClustersInCluster_m]=deal(cutPart1_temp,cutPart2_temp,allClustersInCluster_m_temp);
                            assert(length(allClustersInCluster_m)==length(index_m));
                        end
                        
                    end
                    
                end
            end
            % store current best partition
            subClusters{m}=allClustersInCluster_m;
            subCutParts(m,1:2)=[cutPart1 cutPart2];
            
        end
        
        % print out best cut possible by partitioning of current subgraph
        cut=computeCutCheeger(cutParts,cutPart1,cutPart2,m,l);
        if(verbosity>=2)
            fprintf('Best result achievable by partitioning of subgraph %d:\n',m);
            displayCurrentObjective(cut,normalized);
            fprintf('\n');
        end
        
        % check if partitoning of the current subgraph gives better cut
        if (cut<bestCut)
            [bestCut,bestCutPart1,bestCutPart2,best_m]= deal(cut,cutPart1,cutPart2,m);
            clusters_new=allClusters;
            clusters_new(index_m)=(l-m)*allClustersInCluster_m+clusters_new(index_m);
            
            assert(bestCut>=0);
        end
        
        % if we have already found a partition with cut 0, we don't
        % need to consider the other subgraphs
        if bestCut==0
            break;
        end
    end
    
    % Update
    allClusters=clusters_new;
    cuts(1,l-1)=bestCut;
    clusters(:,l-1)=allClusters;
    
    cutParts(best_m)=bestCutPart1;
    cutParts(l)=bestCutPart2;
    
    % Check that we have the right number of clusters
    assert(length(unique(allClusters))==l);
    
    % Reset subcutparts and subclusters;
    subCutParts(best_m,:)=0;
    subClusters{best_m}= [];
    subCutParts(l,:)=0;
    subClusters{l}= [];
    
    % Print out current objective
    if(verbosity>=1)
        fprintf('Decided to partition subgraph %d. Finished Clustering into %d parts.\n',best_m,l);
        displayCurrentObjective(bestCut,normalized);
        fprintf('\n');
    end
    
end


end

% Computes Rcut/Ncut and Cheeger Cut values
function [cut,cheeger]=computeCutCheeger(cutParts,cutPart1,cutPart2,m,l)

cut= sum(cutParts)-cutParts(m)+cutPart1+cutPart2;
cheeger=max([cutParts((1:l-1)~=m) cutPart1 cutPart2]);

end

% Displays the current objective value
function displayCurrentObjective(cut_temp,normalized)

if (normalized)
    fprintf('...CCNCut: %g\n',cut_temp);
else
    fprintf('...CCRCut: %g\n',cut_temp);
end

end


% Creates two clusters by thresholding the vector vmin_comp obtained on a
% connected component of a subgraph. Given the two clusters on the
% connected component, there are two ways of constructing the final clusters
% on the subgraph, as we can keep each of the clusters on the connected
% component as cluster and merge the other one with the remaining connected
% components. The method takes the one yielding the lower CCBCut.
function [allClustersInClusterM, cutPart1,cutPart2] =  createSubClusters2(vmin_comp,W_comp,opts,normalized,deg,index_comp,index_m,cut_rest,size_rest,size_m)

% input parameter deg has to be the degree vector (also in unnormalised case)
%deg=full(sum(W));
%Make deg a row vector;
if (size(deg,1)>1)
    deg=deg';
end

% get tau value
tau = ccbcutGet(opts,'Tau');

[vminM_sorted, index]=sort(vmin_comp);
[vminU,indexU]=unique(vminM_sorted);

W_sorted=W_comp(index,index);

% calculate cuts
deg_comp=deg(index_m(index_comp));
volumes_threshold=cumsum(deg_comp(index));
triup=triu(W_sorted);
tempcuts_threshold=volumes_threshold - 2*cumsum(full(sum(triup)));
tempcuts_threshold2=(volumes_threshold(end)-volumes_threshold) - (sum(sum(W_sorted))-2*cumsum(full(sum(triup,2)))');

% it may happen that (due to numerical imprecision) the tempcuts
% are a small factor of epsilon below zero.
tempcuts_threshold(tempcuts_threshold<0)=0;
tempcuts_threshold2(tempcuts_threshold2<0)=0;

tempcuts_threshold=tempcuts_threshold(indexU);
tempcuts_threshold2=tempcuts_threshold2(indexU);
volumes_threshold=volumes_threshold(indexU);

% divide by size/volume
if(normalized)
    cutparts1_threshold=(tempcuts_threshold(1:end-1)+cut_rest)./((2*(volumes_threshold(1:end-1)+size_rest)).^(tau/2));
    cutparts1_threshold(isnan(cutparts1_threshold))=0;
    cutparts2_threshold=tempcuts_threshold2(1:end-1)./((2*(volumes_threshold(end)-volumes_threshold(1:end-1))).^(tau/2));
    cutparts2_threshold(isnan(cutparts2_threshold))=0;
    
    cutparts1b_threshold=tempcuts_threshold(1:end-1)./((2*volumes_threshold(1:end-1)).^(tau/2));
    cutparts1b_threshold(isnan(cutparts1b_threshold))=0;
    cutparts2b_threshold=(tempcuts_threshold2(1:end-1)+cut_rest)./((2*((volumes_threshold(end)-volumes_threshold(1:end-1))+size_rest)).^(tau/2));
    cutparts2b_threshold(isnan(cutparts2b_threshold))=0;
else
    sizes_threshold=cumsum(ones(1,size(vmin_comp,1)-1));
    sizes_threshold=sizes_threshold(indexU(1:end-1));
    cutparts1_threshold=(tempcuts_threshold(1:end-1)+cut_rest)./((2*(sizes_threshold+size_rest)).^(tau/2));
    cutparts2_threshold=tempcuts_threshold2(1:end-1)./((2*(size(vmin_comp,1)-sizes_threshold)).^(tau/2));
    
    cutparts1b_threshold=tempcuts_threshold(1:end-1)./((2*sizes_threshold).^(tau/2));
    cutparts2b_threshold=(tempcuts_threshold2(1:end-1)+cut_rest)./((2*((size(vmin_comp,1)-sizes_threshold)+size_rest)).^(tau/2));
end

%calculate total cuts
cuts_threshold=cutparts1_threshold+cutparts2_threshold;
[cut1,threshold_index]=min(cuts_threshold);

cuts_threshold_b=cutparts1b_threshold+cutparts2b_threshold;
[cut1b,threshold_index_b]=min(cuts_threshold_b);

comp_case=1;
if (cut1b<cut1)
    comp_case=2;
end

if(comp_case==1)
    cutPart1=cutparts1_threshold(threshold_index);
    cutPart2=cutparts2_threshold(threshold_index);
    
    allClustersInClusterM_comp= (vmin_comp>vminU(threshold_index));
    
    allClustersInClusterM= zeros(size_m,1);
    allClustersInClusterM(index_comp)=allClustersInClusterM_comp;
else
    cutPart1=cutparts1b_threshold(threshold_index_b);
    cutPart2=cutparts2b_threshold(threshold_index_b);
    
    allClustersInClusterM_comp= (vmin_comp>vminU(threshold_index_b));
    
    allClustersInClusterM= ones(size_m,1);
    allClustersInClusterM(index_comp)=allClustersInClusterM_comp;
end

end


% Tries to separate the connected components into two clusters which have
% roughly the same cardinality/volume
function comp2 = balanceConnectedComponents(comp,sizes,W,normalized)

% for normalized variant, compute the volume for every connected
% component
if(normalized)
    deg=sum(W);
    volumes=zeros(length(sizes),1);
    for l=1:length(sizes)
        volumes(l)=sum(deg(comp==l));
    end
    sizes=volumes;
end

% fill up clusters, trying to balance the size
[sizes_sort,ind]=sort(sizes,'descend');
size_a=0;
size_b=0;
ind_a=[];
for l=1:length(sizes_sort)
    if(size_a<=size_b)
        size_a=size_a+sizes_sort(l);
        ind_a=[ind_a ind(l)];
    else size_b=size_b+sizes_sort(l);
    end
end
comp2=double(ismember(comp,ind_a));
end
