function ll = listreplace(ll, repind, replist)
% Convenience function that replaces a specified element of a list with a given
% replacement.
%
% Syntax:
%   ll = listreplace(ll0, repind, replist) returns a version of the input list
%       ll0, but with the repind-th element replaced by replist
%
% Inputs:
%   ll - {n} size, [list] type, input list to have one of its elements replaces
%   repind - [1] size, [int] type, index of the element in ll to be replaced
%   replist - arbitrary size/type, replacement for ll{repind}
%
% Outputs:
%   ll - {n} size, [list] type, list identical to the input ll, but with
%       repind{repind}=replist

ll{repind} = replist;
