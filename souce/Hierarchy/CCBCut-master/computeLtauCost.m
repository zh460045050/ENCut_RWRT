function costLtau = computeLtauCost(W,Y,tau)
% computeLtauCost: compute the Ltau cost function 
%   sum_{i,j} W(i,j)*Ltaunorm(Y_{i}-y_{j})^tau
%

% author: Nathan D. Cahill
% email: nathan.cahill@rit.edu
% date: 21 February 2018

% find indices of nonzero weights
[i,j,wVals] = find(W);
ut = i<j;
i = i(ut);
j = j(ut);
wVals = wVals(ut);

% compute errors
costLtau = 0;
for m = 1:numel(i)
    costLtau = costLtau + wVals(m).*sum(abs(Y(i(m),:)-Y(j(m),:)).^tau);
end

end