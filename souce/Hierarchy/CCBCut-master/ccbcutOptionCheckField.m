function [validvalue, errmsg, errid, validfield] = ccbcutOptionCheckField(field,value)
% ccbcutOptionCheckField: check validity of structure field contents
% usage: [VALIDVALUE, ERRMSG, ERRID, VALIDFIELD] = ccbcutOptionCheckField('field',V)
%
% Note: This is a helper function for ccbcutSet and ccbcutGet

% author: Nathan D. Cahill
% email: nathan.cahill@rit.edu
% date: 21 February 2018

% empty matrix is always valid
if isempty(value)
    validvalue = true;
    errmsg = '';
    errid = '';
    validfield = true;
    return
end

% Some fields are checked in ccbcutSet/checkfield: Display.
validfield = true;
switch field
    case {'Tau','TolBregSOCX','TolBregSOCF','TolBregL1ProbX','TolBregL1ProbF','BregLambda','R','TolIRRQX','TolIRRQF','TolIRRQL1','TolIRRQL1Curv','Epsilon','IRRQLambda','EpsilonDenom'}
        % real positive scalar
        [validvalue, errmsg, errid] = posReal(field,value);
    case {'KRatio','KRatioDenom'}
        [validvalue, errmsg, errid] = nonNegReal(field,value);
    case {'IRRQConstraint'}
        [validvalue, errmsg, errid] = stringsType(field,value,{'none','soft','hard'});
    case {'IRRQq'}
        [validvalue, errmsg, errid] = boundedRealClosed(field,value,[1,2]);
    case {'IRRQAlpha'}
        [validvalue, errmsg, errid] = boundedRealOpen(field,value,[0,1]);        
    case {'MaxSOCIters','MaxL1ProbIters'}
        [validvalue, errmsg, errid] = posInteger(field,value,{'100000*numberofvariables'});
    case {'MinIRRQIters','MaxIRRQIters','KmeansNumStartPoints'}
        [validvalue, errmsg, errid] = nonNegInteger(field,value,{'100000*numberofvariables'});
    case {'Algorithm'}
        [validvalue, errmsg, errid] = stringsType(field,value,{'bregman','irrq'});
    case {'BalanceType'}
        [validvalue, errmsg, errid] = stringsType(field,value,{'normalized','ratio'});
    case {'RoundingType'}
        [validvalue, errmsg, errid] = stringsType(field,value,{'kmeans','graclus','none'});
    case {'KmeansDistance'}
        [validvalue, errmsg, errid] = stringsType(field,value,{'sqeuclidean','cityblock','cosine','correlation','hamming'});
    otherwise
        validfield = false;
        validvalue = false;
        % No need to set an error. If the field isn't valid for MATLAB or ccbcut,
        % will have already errored in ccbcutSet. If field is valid for MATLAB,
        % then the error will be an invalid value for MATLAB.
        errid = '';
        errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = posReal(field,value,string)
% Any nonnegative real scalar or sometimes a special string
valid =  isreal(value) && isscalar(value) && (value > 0) ;
if nargin > 2
    valid = valid || isequal(value,string);
end
if ~valid
    if ischar(value)
        errid = 'ccbcutLib:ccbcutOptionCheckField:PosReal:nonpositiveNum';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real positive scalar (not a string).',field);
    else
        errid = 'ccbcutLib:ccbcutOptionCheckField:PosReal:nonpositiveNum';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real positive scalar.',field);
    end
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = nonNegReal(field,value,string)
% Any nonnegative real scalar or sometimes a special string
valid =  isreal(value) && isscalar(value) && (value >= 0) ;
if nargin > 2
    valid = valid || isequal(value,string);
end
if ~valid
    if ischar(value)
        errid = 'ccbcutLib:ccbcutOptionCheckField:NonNegReal:negativeNum';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative scalar (not a string).',field);
    else
        errid = 'ccbcutLib:ccbcutOptionCheckField:NonNegReal:negativeNum';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative scalar.',field);
    end
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = nonPosReal(field,value,string)
% Any nonpositive real scalar or sometimes a special string
valid =  isreal(value) && isscalar(value) && (value <= 0) ;
if nargin > 2
    valid = valid || isequal(value,string);
end
if ~valid
    if ischar(value)
        errid = 'ccbcutLib:ccbcutOptionCheckField:NonPosReal:positiveNum';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-positive scalar (not a string).',field);
    else
        errid = 'ccbcutLib:ccbcutOptionCheckField:NonPosReal:positiveNum';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-positive scalar.',field);
    end
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = posInteger(field,value,strings)
% Any positive real integer scalar or sometimes a special string
valid =  isreal(value) && isscalar(value) && (value > 0) && value == floor(value) ;
if nargin > 2
    valid = valid || any(strcmp(value,strings));
end
if ~valid
    if ischar(value)
        errid = 'ccbcutLib:ccbcutOptionCheckField:notAPosInteger';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real positive scalar (not a string).',field);
    else
        errid = 'ccbcutLib:ccbcutOptionCheckField:notAPosInteger';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real positive scalar.',field);
    end
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = nonNegInteger(field,value,strings)
% Any nonnegative real integer scalar or sometimes a special string
valid =  isreal(value) && isscalar(value) && (value >= 0) && value == floor(value) ;
if nargin > 2
    valid = valid || any(strcmp(value,strings));
end
if ~valid
    if ischar(value)
        errid = 'ccbcutLib:ccbcutOptionCheckField:notANonNegInteger';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative scalar (not a string).',field);
    else
        errid = 'ccbcutLib:ccbcutOptionCheckField:notANonNegInteger';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative scalar.',field);
    end
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = matrixType(field,value,strings)
% Any matrix
valid =  isa(value,'double');
if nargin > 2
    valid = valid || any(strcmp(value,strings));
end
if ~valid
    if ischar(value)
        errid = 'ccbcutLib:ccbcutOptionCheckField:notANonNegInteger';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a matrix (not a string).',field);
    else
        errid = 'ccbcutLib:ccbcutOptionCheckField:posMatrixType:notAPosMatrix';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a matrix.',field);
    end
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = posMatrixType(field,value)
% Any positive scalar or all positive vector
valid =  isa(value,'double') && all(value > 0) && isvector(value);
if ~valid
    errid = 'ccbcutLib:ccbcutOptionCheckField:posMatrixType:notAPosMatrix';
    errmsg = sprintf('Invalid value for OPTIONS parameter %s: \n must be a positive scalar or a vector with positive entries.',field);
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = functionType(field,value)
% Any function handle or string (we do not test if the string is a function name)
valid =  ischar(value) || isa(value, 'function_handle');
if ~valid
    errid = 'ccbcutLib:ccbcutOptionCheckField:functionType:notAFunction';
    errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a function handle.',field);
else
    errid = '';
    errmsg = '';
end
%-----------------------------------------------------------------------------------------
function [valid, errmsg, errid] = stringsType(field,value,strings)
% One of the strings in cell array strings
valid =  ischar(value) && any(strcmp(value,strings));

% To print out the error message beautifully, need to get the commas and "or"s
% in all the correct places while building up the string of possible string values.
if ~valid
    allstrings = ['''',strings{1},''''];
    for index = 2:(length(strings)-1)
        % add comma and a space after all but the last string
        allstrings = [allstrings, ', ''', strings{index},''''];
    end
    if length(strings) > 2
        allstrings = [allstrings,', or ''',strings{end},''''];
    elseif length(strings) == 2
        allstrings = [allstrings,' or ''',strings{end},''''];
    end
    errid = 'ccbcutLib:ccbcutOptionCheckField:stringsType:notAStringsType';
    errmsg = sprintf('Invalid value for OPTIONS parameter %s:\n must be %s.',field, allstrings);
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------
function [valid, errmsg, errid] = boundedRealClosed(field,value,bounds)
% Scalar in the bounds
valid =  isa(value,'double') && isscalar(value) && ...
    (value >= bounds(1)) && (value <= bounds(2));
if ~valid
    errid = 'ccbcutLib:ccbcutOptionCheckField:boundedRealClosed:notAboundedReal';
    errmsg = sprintf('Invalid value for OPTIONS parameter %s: \n must be a scalar in the range [%6.3g, %6.3g].', ...
        field, bounds(1), bounds(2));
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------
function [valid, errmsg, errid] = boundedRealOpen(field,value,bounds)
% Scalar in the bounds
valid =  isa(value,'double') && isscalar(value) && ...
    (value > bounds(1)) && (value < bounds(2));
if ~valid
    errid = 'ccbcutLib:ccbcutOptionCheckField:boundedRealOpen:notAboundedReal';
    errmsg = sprintf('Invalid value for OPTIONS parameter %s: \n must be a scalar in the range (%6.3g, %6.3g).', ...
        field, bounds(1), bounds(2));
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------
function [valid, errmsg, errid] = logicalScalar(field,value)
% Scalar in the bounds
valid =  isa(value,'logical') && isscalar(value);
if ~valid
    errid = 'ccbcutLib:ccbcutOptionCheckField:logicalScalar:notALogicalScalar';
    errmsg = sprintf('Invalid value for OPTIONS parameter %s: \n must be a logical scalar.', ...
        field);
else
    errid = '';
    errmsg = '';
end

