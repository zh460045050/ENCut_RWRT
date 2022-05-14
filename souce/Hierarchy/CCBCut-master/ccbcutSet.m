function options = ccbcutSet(varargin)
% ccbcutSet: Create/alter ccbcut OPTIONS structure.
%   OPTIONS = ccbcutSet('PARAM1',VALUE1,'PARAM2',VALUE2,...) creates a
%   ccbcut options structure OPTIONS in which the named parameters have the 
%   specified values.  Any unspecified parameters are set to [] 
%   (parameters with value [] indicate to use the default value for that 
%   parameter when OPTIONS is passed to the ccbcut function). It is sufficient
%   to type only the leading characters that uniquely identify the 
%   parameter.  Case is ignored for parameter names. NOTE: For values that
%   are strings, the complete string is required.
%
%   OPTIONS = ccbcutSet(OLDOPTS,'PARAM1',VALUE1,...) creates a copy of 
%   OLDOPTS with the named parameters altered with the specified values.
%
%   OPTIONS = ccbcutSet(OLDOPTS,NEWOPTS) combines an existing options 
%   structure OLDOPTS with a new options structure NEWOPTS.  Any parameters
%   in NEWOPTS with non-empty values overwrite the corresponding old 
%   parameters in OLDOPTS.
%
%   ccbcutSet with no input arguments and no output arguments displays 
%   all parameter names and their possible values, with defaults shown in {}
%   when the default is the same for all functions that use that option -- use
%   ccbcutSet(ccbcutFunction) to see options for a specific function.).
%
%   OPTIONS = ccbcutSet (with no input arguments) creates an options 
%   structure OPTIONS where all the fields are set to [].
%
%   OPTIONS = ccbcutSet(ccbcutFUNCTION) creates an options structure
%   with all the parameter names and default values relevant to the npable
%   registration function named in ccbcutFunction. For example,
%           ccbcutSet('ccbcut2D')
%   or
%           ccbcutSet(@ccbcut2D)
%   returns an options structure containing all the parameter names and
%   default values relevant to the function 'ccbcut'.
%
%ccbcutSet PARAMETERS for MATLAB
%Display - Level of display [ off | iter | notify | {final} ]
%Algorithm - Algorithm for compute CCBCut embedding
%           [ bregman | {irrq} ]
%Tau - Order of CCBCut objective function
%           [ positive scalar {1} ]
%BalanceType - Choice of balance type in CCBCut cost
%           [ {normalized} | ratio ]
%RoundingType - Choice of method for rounding embedding to approximate
%   solution to discrete problem
%           [ {kmeans} | graclus | none]
%KmeansNumStartPoints - number of times to run kmeans algorithm on 
%   computed embedding from random starting points; only used if
%   RoundingType == 'kmeans'
%           [ positive scalar {10} ]
%KmeansDistance - distance measure to use in kmeans algorithm; only used if
%   RoundingType == 'kmeans'
%           [ {sqeuclidean} | cityblock | cosine | correlation | hamming ]
%
% Options specific to Bregman algorithm:
%MaxSOCIters - Maximum number of iterations of the Splitting Orthogonality 
%       Constraint algorithm (outer loop of Bregman iteration)  
%           [ positive scalar {100} ]
%MaxL1ProbIters - Maximum number of iterations of the L1 Problem algorithm 
%       (inner loop of Bregman iteration)  
%           [ positive scalar {100} ]
%TolBregSOCX - SOC iteration termination tolerance on the L2 norm of the
%       difference between successive estimates of embedding coordinates
%           [ positive scalar {1e-3} ]
%TolBregSOCF - SOC iteration termination tolerance on the relative
%       difference between successive estimates of the objective function
%           [ positive scalar {1e-3} ]
%TolBregL1ProbX - L1 subproblem iteration termination tolerance on the L2 
%       norm of the difference between successive estimates of embedding
%       coordinates
%           [ positive scalar {1e-3} ]
%TolBregL1ProbF - L1 subproblem iteration termination tolerance on the
%       relative difference between successive estimates of the objective
%       function
%           [ positive scalar {1e-3} ]
%BregLambda - Nonnegative parameter for use in first two steps of L1 Problem
%       iteration
%           [ nonnegative scalar {1000} ]
%R - Nonnegative parameter for use in first step of L1 Problem iteration
%           [ nonnegative scalar {100} ]
%
% Options specific to IRRQ algorithm:
%MinIRRQIters - Minimum number of iterations of the Majorize/minimize algorithm  
%           [ nonnegative scalar {5} ]
%MaxIRRQIters - Maximum number of iterations of the Majorize/minimize algorithm  
%           [ nonnegative scalar {100} ]
%TolIRRQX - Majorize/minimize iteration termination tolerance on the L2 norm 
%       of the difference between successive estimates of embedding
%       coordinates
%           [ positive scalar {1e-3} ]
%TolIRRQF - Majorize/minimize iteration termination tolerance on the relative
%       difference between successive estimates of the (majorizing) 
%       objective function
%           [ positive scalar {1e-3} ]
%TolIRRQL1 - Majorize/minimize iteration termination tolerance on the
%       relative difference between successive estimates of the (original
%       weighted L1) objective function
%           [ positive scalar {1e-3} ]
%TolIRRQL1Curve - Termination tolerance on the curvature of the Ltau 
%       objective function relative to the Ltau cost at the first 
%       iteration. CCBCut will be rolled back one iteration if the
%       curvature becomes greater than TolIRRQL1Curve. 
%           [ positive scalar {1e-1} ]
%Epsilon - Starting value of weight regularization parameter computed in IRRQ
%       iteration.
%           [ positive scalar {1} ]
%KRatio - Fraction of the distance from k (assuming the solution will be
%       k-sparse) to the total number of nonzero weights. k will be
%       automatically estimated in the first iteration. KRatio will be used
%       to compute K (see Daubechies et al.)
%           [ nonnegative scalar {0.5} ]
%IRRQConstraint - choice of how to incorporate the balance constraint
%           [ {none} | soft | hard ]
%IRRQq - power used to define l^q norm in numerator/denominator of Rayleigh
%       quotient. Must satisfy 1<= IRRQq <= 2.
%           [ positive scalar in [1,2] {1} ]
%
% Options specific to IRRQ algorithm with IRRQConstraint == 'soft':
%IRRQLambda - Nonnegative parameter trading off L2 penalty on balance
%       constraint with Rayleigh quotient
%           [ positive scalar {1} ]
%IRRQAlpha - parameter in (0,1) used in reweighting calculation
%           [ positive scalar in (0,1) {0.5} ]
%
%   Note: To see ccbcutSet parameters for the ccbcut Toolbox
%         (if you have the ccbcut Toolbox installed), type
%             help ccbcutOptions
%
% Options specific to IRRQ algorithm with IRRQConstraint == 'hard':
%EpsilonDenom - Starting value of weight regularization parameter computed 
%       in denominator of IRRQ iteration.
%           [ positive scalar {1} ]
%KRatioDenom - KRatio for denominator term in IRRQ iteration
%           [ nonnegative scalar {0.5} ]
%
%   Note: To see ccbcutSet parameters for the ccbcut Toolbox
%         (if you have the ccbcut Toolbox installed), type
%             help ccbcutOptions
%   Examples
%     To create options with the default options for ccbcut
%       options = ccbcutSet('ccbcut');
%     To create an options structure with TolSOC equal 
%     to 1e-2
%       options = ccbcutSet('TolSOC',1e-2);
%     To change the Display value of options to 'iter'
%       options = ccbcutSet(options,'Display','iter');
%
%   See also ccbcutGet, ccbcut.

% author: Nathan D. Cahill
% email: nathan.cahill@rit.edu
% date: 21 February 2018

% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    fprintf('                  Display: [ off | iter | notify | {final} ]\n');
    fprintf('                Algorithm: [ bregman | {irrq} ]\n');
    fprintf('                      Tau: [ positive scalar ]\n');
    fprintf('              BalanceType: [ {normalized} | ratio ]\n');
    fprintf('             RoundingType: [ {kmeans} | graclus | none ]\n');
    fprintf('              MaxSOCIters: [ positive integer ]\n');
    fprintf('           MaxL1ProbIters: [ positive integer ]\n');
    fprintf('              TolBregSOCX: [ positive scalar ]\n');
    fprintf('              TolBregSOCF: [ positive scalar ]\n');
    fprintf('           TolBregL1ProbX: [ positive scalar ]\n');
    fprintf('           TolBregL1ProbF: [ positive scalar ]\n');
    fprintf('               BregLambda: [ positive scalar ]\n');
    fprintf('                        R: [ positive scalar ]\n');
    fprintf('             MinIRRQIters: [ nonnegative integer ]\n');
    fprintf('             MaxIRRQIters: [ nonnegative integer ]\n');
    fprintf('                 TolIRRQX: [ positive scalar ]\n');
    fprintf('                 TolIRRQF: [ positive scalar ]\n');
    fprintf('                TolIRRQL1: [ positive scalar ]\n');
    fprintf('            TolIRRQL1Curv: [ positive scalar ]\n');
    fprintf('                  Epsilon: [ positive scalar ]\n');
    fprintf('                   KRatio: [ positive scalar ]\n');
    fprintf('           IRRQConstraint: [ {none} | soft | hard ]\n');
    fprintf('               IRRQLambda: [ positive scalar ]\n');
    fprintf('                    IRRQq: [ positive scalar in [1,2] ]\n');
    fprintf('                IRRQAlpha: [ positive scalar in (0,1) ]\n');
    fprintf('             EpsilonDenom: [ positive scalar ]\n');
    fprintf('              KRatioDenom: [ positive scalar ]\n');
        
    try
        ccbcutoptions;
    catch
        lasterrstruct = lasterror;
        if strcmp(lasterrstruct.identifier, 'MATLAB:UndefinedFunction')
            % Function ccbcutOPTIONS not found, so we assume ccbcut Toolbox not on path
            %   and print the "MATLAB only" fields.
            % clean up last error
            lasterr('');
        else
            rethrow(lasterror);
        end
    end

    fprintf('\n');
    return;
end

% Create a struct of all the fields with all values set to 
allfields = {'Display';'Algorithm';'Tau';'BalanceType';'RoundingType';...
    'KmeansNumStartPoints';'KmeansDistance';'MaxSOCIters';'MaxL1ProbIters';...
    'TolBregSOCX';'TolBregSOCF';'TolBregL1ProbX';'TolBregL1ProbF';...
    'BregLambda';'R';'MinIRRQIters';'MaxIRRQIters';'TolIRRQX';'TolIRRQF';'TolIRRQL1';...
    'TolIRRQL1Curv';'Epsilon';...
    'KRatio';'IRRQConstraint';'IRRQLambda';'IRRQq';'IRRQAlpha';'EpsilonDenom';...
    'KRatioDenom'};
try 
    % assume we have the ccbcut Toolbox
    ccbcuttbx = true;
    ccbcutfields = ccbcutOptionGetFields;  
    allfields = [allfields; optimfields];
catch
    lasterrstruct = lasterror;
    if strcmp(lasterrstruct.identifier, 'MATLAB:UndefinedFunction')
        % Function ccbcutOPTIONGETFIELDS not found, so we assume ccbcut Toolbox not on path
        %   and use the "MATLAB only" struct.
        ccbcuttbx = false;
        % clean up last error
        lasterr('');
    else
        rethrow(lasterror);
    end
end
% create cell array
structinput = cell(2,length(allfields));
% fields go in first row
structinput(1,:) = allfields';
% []'s go in second row
structinput(2,:) = {[]};
% turn it into correctly ordered comma separated list and call struct
options = struct(structinput{:});

numberargs = nargin; % we might change this value, so assign it
% If we pass in a function name then return the defaults.
if (numberargs==1) && (ischar(varargin{1}) || isa(varargin{1},'function_handle') )
    if ischar(varargin{1})
        funcname = lower(varargin{1});
        if ~exist(funcname)
            msg = sprintf(...
                'No default options available: the function ''%s'' does not exist on the path.',funcname);
            error('MATLAB:ccbcutSet:FcnNotFoundOnPath', msg)
        end
    elseif isa(varargin{1},'function_handle')
        funcname = func2str(varargin{1});
    end
    try 
        optionsfcn = feval(varargin{1},'defaults');
    catch
        msg = sprintf(...
            'No default options available for the function ''%s''.',funcname);
        error('MATLAB:ccbcutSet:NoDefaultsForFcn', msg)
    end
    % The defaults from the optim functions don't include all the fields
    % typically, so run the rest of ccbcutset as if called with
    % ccbcutSet(options,optionsfcn)
    % to get all the fields.
    varargin{1} = options;
    varargin{2} = optionsfcn;
    numberargs = 2;
end

Names = allfields;
m = size(Names,1);
names = lower(Names);

i = 1;
while i <= numberargs
    arg = varargin{i};
    if ischar(arg)                         % arg is an option name
        break;
    end
    if ~isempty(arg)                      % [] is a valid options argument
        if ~isa(arg,'struct')
            error('MATLAB:ccbcutSet:NoParamNameOrStruct',...
                ['Expected argument %d to be a string parameter name ' ...
                'or an options structure\ncreated with ccbcutSet.'], i);
        end
        for j = 1:m
            if any(strcmp(fieldnames(arg),[Names{j,:}]))
                val = arg.(Names{j,:});
            else
                val = [];
            end
            if ~isempty(val)
                if ischar(val)
                    val = lower(deblank(val));
                end
                checkfield(Names{j,:},val,ccbcuttbx);
                options.(Names{j,:}) = val;
            end
        end
    end
    i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(numberargs-i+1,2) ~= 0
    error('MATLAB:ccbcutSet:ArgNameValueMismatch',...
        'Arguments must occur in name-value pairs.');
end
expectval = 0;                          % start expecting a name, not a value
while i <= numberargs
    arg = varargin{i};

    if ~expectval
        if ~ischar(arg)
            error('MATLAB:ccbcutSet:InvalidParamName',...
                'Expected argument %d to be a string parameter name.', i);
        end

        lowArg = lower(arg);
        j = strmatch(lowArg,names);
        if isempty(j)                       % if no matches
            [wasinmatlab, optionname] = checkccbcutonlylist(lowArg);
            if ~wasinmatlab
                error('MATLAB:ccbcutSet:InvalidParamName',...
                    'Unrecognized parameter name ''%s''.', arg);
            else
                warning('MATLAB:ccbcutSet:InvalidParamName',...
                    ['The option ''%s'' is a ccbcut Toolbox option and is not\n', ...
                     'used by any MATLAB functions. This option will be ignored and not included\n', ...
                     'in the options returned by ccbcutSet. Please change your code to not use \n', ...
                     'this option as it will error in a future release.'], ...
                     optionname);
                i = i + 2; % skip this parameter and its value; go to next parameter
                continue; % skip the rest of this loop
            end
        elseif length(j) > 1                % if more than one match
            % Check for any exact matches (in case any names are subsets of others)
            k = strmatch(lowArg,names,'exact');
            if length(k) == 1
                j = k;
            else
                msg = sprintf('Ambiguous parameter name ''%s'' ', arg);
                msg = [msg '(' Names{j(1),:}];
                for k = j(2:length(j))'
                    msg = [msg ', ' Names{k,:}];
                end
                msg = sprintf('%s).', msg);
                error('MATLAB:ccbcutSet:AmbiguousParamName', msg);
            end
        end
        expectval = 1;                      % we expect a value next

    else
        if ischar(arg)
            arg = lower(deblank(arg));
        end
        checkfield(Names{j},arg,ccbcuttbx);
        options.(Names{j}) = arg;
        expectval = 0;
    end
    i = i + 1;
end

if expectval
    error('MATLAB:ccbcutSet:NoValueForParam',...
        'Expected value for parameter ''%s''.', arg);
end

%-------------------------------------------------
function checkfield(field,value,ccbcuttbx)
%CHECKFIELD Check validity of structure field contents.
%   CHECKFIELD('field',V,ccbcutTBX) checks the contents of the specified
%   value V to be valid for the field 'field'. ccbcutTBX indicates if 
%   the ccbcut Toolbox is on the path.
%

% empty matrix is always valid
if isempty(value)
    return
end

% See if it is one of the valid MATLAB fields.  It may be both a ccbcut
% and MATLAB field, e.g. MaxFunEvals, in which case the MATLAB valid
% test may fail and the ccbcut one may pass.
validfield = true;
switch field
    case {'Display'} % off,none,iter,final,notify,testing
        [validvalue, errmsg, errid] = displayType(field,value);
    otherwise
        validfield = false;  
        validvalue = false;
        errmsg = sprintf('Unrecognized parameter name ''%s''.', field);
        errid = 'MATLAB:ccbcutSet:checkfield:InvalidParamName';
end

if validvalue 
    return;
elseif ~ccbcuttbx && validfield  
    % Throw the MATLAB invalid value error
    error(errid, errmsg);
else % Check if valid for ccbcut Tbx
    [optvalidvalue, opterrmsg, opterrid, optvalidfield] = ccbcutOptionCheckField(field,value);
    if optvalidvalue
        return;
    elseif optvalidfield
        % Throw the ccbcut invalid value error
        error(opterrid, opterrmsg)
    else % Neither field nor value is valid for ccbcut
        % Throw the MATLAB invalid value error (can't be invalid field for
        % MATLAB & ccbcut or would have errored already in ccbcutSet).
        error(errid, errmsg)
    end
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
        errid = 'MATLAB:funfun:ccbcutset:NonNegReal:negativeNum';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative scalar (not a string).',field);
    else
        errid = 'MATLAB:funfun:ccbcutset:NonNegReal:negativeNum';
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
        errid = 'MATLAB:funfun:ccbcutset:NonPosReal:positiveNum';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-positive scalar (not a string).',field);
    else
        errid = 'MATLAB:funfun:ccbcutset:NonPosReal:positiveNum';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-positive scalar.',field);
    end
else
    errid = '';
    errmsg = '';
end
%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = nonNegInteger(field,value,string)
% Any nonnegative real integer scalar or sometimes a special string
valid =  isreal(value) && isscalar(value) && (value >= 0) && value == floor(value) ;
if nargin > 2
    valid = valid || isequal(value,string);
end
if ~valid
    if ischar(value)
        errid = 'MATLAB:funfun:ccbcutSet:nonNegInteger:notANonNegInteger';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative integer (not a string).',field);
    else
        errid = 'MATLAB:funfun:ccbcutSet:nonNegInteger:notANonNegInteger';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative integer.',field);
    end
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = displayType(field,value)
% One of these strings: on, off, none, iter, final, notify
valid =  ischar(value) && any(strcmp(value,{'on';'off';'none';'iter';'final';'notify';'testing'}));
if ~valid
    errid = 'MATLAB:funfun:ccbcutSet:displayType:notADisplayType';
    errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''off'',''on'',''iter'',''notify'', or ''final''.',field);
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = onOffType(field,value)
% One of these strings: on, off
valid =  ischar(value) && any(strcmp(value,{'on';'off'}));
if ~valid
    errid = 'MATLAB:funfun:ccbcutSet:onOffType:notOnOffType';
    errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''off'' or ''on''.',field);
else
    errid = '';
    errmsg = '';
end

%-----------------------------------------------------------------------------------------

function [valid, errmsg, errid] = functionType(field,value)
% Any function handle or string (we do not test if the string is a function name)
valid =  ischar(value) || isa(value, 'function_handle');
if ~valid
    errid = 'MATLAB:funfun:ccbcutSet:functionType:notAFunction';
    errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a function.',field);
else
    errid = '';
    errmsg = '';
end

%--------------------------------------------------------------------------------

function [wasinmatlab, optionname] = checkccbcutonlylist(lowArg)
% Check if the user is trying to set an option that is only used by
% ccbcut Toolbox functions -- this used to have no effect.
% Now it will warn. In a future release, it will error.  
names =  {};
lowernames = lower(names);
k = strmatch(lowArg,lowernames);
wasinmatlab = ~isempty(k);
if wasinmatlab
    optionname = names{k};
else
    optionname = '';
end
