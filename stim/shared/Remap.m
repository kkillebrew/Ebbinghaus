function newMatrix = Remap(x, oldVals, newVals)

% Remap.m
%    [newMatrix] = Remap(x, oldVals, newVals)
%
% DESCRIPTION
%    Given matrix x, create a new matrix by remapping any occurances of the
%    values defined in oldVals with the index-corresponding value supplied
%    by newVals.
%    
%    x = [1 2 1 2;
%         3 4 4 3];
%
%    Remap(x,[1 4],[10 14]) = [10 2 10 2;
%                              3 14 14 3;]
%
%
% ARGUMENTS
%    'x' is m x n input matrix.
%    'oldVals' is a vector of values that (might) occur in x to be remapped
%    'newVals' is a vector of the same length as oldVals that defines the
%              re-mapped value of items in x.
%
%    Currently, all data is expected to be numerical (i.e., no strings) and
%    in standard matrix format (no cell arrays)
%
% RETURN
%    newMatrix, same size as x
%
% SEE ALSO
%    REMAPCELLSTR

% 08.07.08  rehbm  Wrote it.
% 09.29.08  rehbm  Algorithm failed on logical x. Fix by changing datatype.
% 07.29.09  rehbm  Added check that input is all numeric


%% validate arguments
if nargin ~= 3
    error('Usage Remap(x,oldVals,newVals).')
end

% oldVals and newVals must be the same length
if (length(oldVals) ~= length(newVals))
    error('oldVals and newVals must be the same length.')
end

if (~isnumeric(x) & ~islogical(x)) || ~isnumeric(oldVals) || ~isnumeric(newVals)
    error('all inputs must be numeric')
end


%% datatype fix
if islogical(x)
    x = 1*x; % so it isn't logical
end


%% do the work
newMatrix = x;

% special cases, nothing to do
if isempty(x) || isempty(oldVals)
    return
end

for i = 1:length(oldVals)
    newMatrix(x==oldVals(i)) = newVals(i);
end
