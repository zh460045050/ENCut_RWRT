function [fold,FctValOuter]=computeEigenvectorGeneral_ccbcut(W,opts)
% Computes the solution to relaxed CCBCuts for k=2
%
% Usage:
%   [fold,FctValOuter] = computeEigenvectorGeneral_ccbcut(W,opts)
%
% Input:
%   W: Sparse symmetric weight matrix.
%   opts: CCBCut options structure.
%
% Output:
%	fold: Final eigenvector.
%   FctValOuter: Values of the functional in each iteration.
%
% author: Nathan D. Cahill, based on modifications of code from:
%
% (C)2010-11 Matthias Hein and Thomas Buehler
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de

assert(isnumeric(W) && issparse(W),'Wrong usage. W should be sparse and numeric.');

% compute solution to relaxed ccbcut problem with k=2
[~,fold] = ccbcut(W,2,[],opts);

FctValOuter = [];

end