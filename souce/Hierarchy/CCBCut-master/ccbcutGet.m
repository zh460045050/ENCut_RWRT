function o = ccbcutGet(options,name,default,flag)
% ccbcutGet: get ccbcut OPTIONS parameters.
%   val = ccbcutGet(OPTIONS,'NAME') extracts the value of the named 
%   parameter from ccbcut options structure OPTIONS, returning an empty
%   matrix if the parameter value is not specified in OPTIONS.  It is 
%   sufficient to type only the leading characters that uniquely identify 
%   the parameter.  Case is ignored for parameter names. [] is a valid 
%   OPTIONS argument.
%
%   VAL = ccbcutGet(OPTIONS,'NAME',DEFAULT) extracts the named parameter
%   as above, but returns DEFAULT if the named parameter is not specified 
%   (is []) in OPTIONS.  For example
%
%     val = ccbcutGet(opts,'tolSOC',1e-4);
%
%   returns val = 1e-4 if the tolSOC property is not specified in opts.
%
%   See also ccbcutSet.

% author: Nathan D. Cahill
% email: nathan.cahill@rit.edu
% date: 21 February 2018

% undocumented usage for fast access with no error checking
if (nargin == 4) && isequal(flag,'fast')
    o = ccbcutGetFast(options,name,default);
    return
end

if nargin < 2
    error('MATLAB:ccbcutGet:NotEnoughInputs', 'Not enough input arguments.');
end
if nargin < 3
    default = [];
end

if ~isempty(options) && ~isa(options,'struct')
    error('MATLAB:ccbcutGet:Arg1NotStruct',...
        'First argument must be an options structure created with ccbcutSET.');
end

if isempty(options)
    o = default;
    return;
end

%allfields = {'Display'};
try
    allfields = ccbcutOptionGetFields;  
    %allfields = [allfields; ccbcutFields];
catch
    lasterrstruct = lasterror;
    if strcmp(lasterrstruct.identifier, 'MATLAB:UndefinedFunction')
        % Function ccbcutOptionGetFields not found, so we assume ccbcut Toolbox not on path
        %   and use the "MATLAB only" struct.
        % clean up last error
        lasterr('');
    else
        rethrow(lasterror);
    end
end

Names = allfields;

name = deblank(name(:)'); % force this to be a row vector
j = find(strncmpi(name,Names,length(name)));
if isempty(j)               % if no matches
    error('MATLAB:ccbcutGet:InvalidPropName',...
        ['Unrecognized property name ''%s''.  ' ...
        'See ccbcutSet for possibilities.'], name);
elseif length(j) > 1            % if more than one match
    % Check for any exact matches (in case any names are subsets of others)
    k = find(strcmpi(name,Names));
    if length(k) == 1
        j = k;
    else
        msg = sprintf('Ambiguous property name ''%s'' ', name);
        msg = [msg '(' Names{j(1),:}];
        for k = j(2:length(j))'
            msg = [msg ', ' Names{k,:}];
        end
        msg = sprintf('%s).', msg);
        error('MATLAB:ccbcutGet:AmbiguousPropName', msg);
    end
end

if any(strcmp(Names,Names{j,:}))
    o = options.(Names{j,:});
    if isempty(o)
        o = default;
    end
else
    o = default;
end

%------------------------------------------------------------------
function value = ccbcutGetFast(options,name,defaultopt)
%ccbcutGetFast Get ccbcut OPTIONS parameter with no error checking.
%   VAL = ccbcutGetFast(OPTIONS,FIELDNAME,DEFAULTOPTIONS) will get the
%   value of the FIELDNAME from OPTIONS with no error checking or
%   fieldname completion. If the value is [], it gets the value of the
%   FIELDNAME from DEFAULTOPTIONS, another OPTIONS structure which is
%   probably a subset of the options in OPTIONS.
%

if isempty(options)
     value = defaultopt.(name);
     return;
end
% We need to know if name is a valid field of options, but it is faster to use 
% a try-catch than to test if the field exists and if the field name is
% correct. If the options structure is from an older version of the
% toolbox, it could be missing a newer field.
try
    value = options.(name);
catch
    value = [];
    lasterr('');  % clean up last error
end

if isempty(value)
    value = defaultopt.(name);
end


