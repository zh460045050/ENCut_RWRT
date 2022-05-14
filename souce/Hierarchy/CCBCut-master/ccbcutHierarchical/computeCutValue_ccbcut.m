function [cutpart1,cutpart2] = computeCutValue_ccbcut(clusters,W,opts,deg)
% Computes the components in the CCBCut expression. 
%
% Usage: [cutpart1,cutpart2] = computeCutValue_ccbcut(clusters,W,opts)
%
% One then has CCBCut = cutpart1 + cutpart2
%
% author: Nathan D. Cahill, based on modifications of code from:
%
% (C)2010-11 Thomas Buehler and Matthias Hein
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de

tau = ccbcutGet(opts,'Tau');

W3= sum(W(clusters==1,:),1);
cut=full(sum(W3(clusters~=1),2));

if(cut==0)
    cutpart1=0;
    cutpart2=0;
else
    if isequal(ccbcutGet(opts,'BalanceType'),'ratio')
        sizeA = sum(clusters==1);
        sizeB = size(clusters,1)-sizeA;
        
        cutpart1=cut/((2*sizeA)^(tau/2));
        cutpart2=cut/((2*sizeB)^(tau/2));
    else
        degA=deg(clusters==1);
        volA=sum(degA);
        
        degB=deg(clusters~=1);
        volB=sum(degB);
        
        cutpart1=cut/((2*volA)^(tau/2));
        cutpart2=cut/((2*volB)^(tau/2));
    end
end

end