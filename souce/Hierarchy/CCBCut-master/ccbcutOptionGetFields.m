function ccbcutFields = ccbcutOptionGetFields
% ccbcutOptionGetFields: fieldnames of options in ccbcut Toolbox
% usage: ccbcutFields = ccbcutOptionGetFields;
%
% Note: This is a helper function for ccbcutSet and ccbcutGet
%

% author: Nathan D. Cahill
% email: nathan.cahill@rit.edu
% date: 21 February 2018

ccbcutFields =   {...
    'Display'; ...
    'Algorithm'; ...
    'Tau'; ...
    'BalanceType'; ...
    'RoundingType'; ...
    'KmeansNumStartPoints'; ...
    'KmeansDistance'; ...
    'GraclusPath'; ...
    'GraclusTempDir'; ...
    'MaxSOCIters';  ...
    'MaxL1ProbIters';  ...
    'TolBregSOCX';...
    'TolBregSOCF';...
    'TolBregL1ProbX';...
    'TolBregL1ProbF'; ...
    'BregLambda'; ...
    'R';...
    'MinIRRQIters'; ...
    'MaxIRRQIters'; ...
    'TolIRRQX';...
    'TolIRRQF';...
    'TolIRRQL1'; ...
    'TolIRRQL1Curv'; ...
    'Epsilon';...
    'KRatio';...
    'IRRQConstraint';...
    'IRRQLambda';...
    'IRRQq';...
    'IRRQAlpha';...
    'EpsilonDenom';...
    'KRatioDenom'};
