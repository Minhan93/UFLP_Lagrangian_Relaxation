function [X, TC] = UFLP_Feas(k, C, y)
% The Feasible Solution Generated From UFLP-LR Subproblem Solution.
% 
% [X, TC] = UFLP_Feas(k, C, y)
% 
% INPUT ARGUMENTS:
%     k = n-element fixed cost vector, where k(i) is cost of NF at Site i
%         (if scalar, then same fixed cost)
%     C = n x m variable cost matrix,
%         where C(i,j) is the cost of serving EF j from NF i
%     y = NF site index vector from UFLP-LR subproblem sulotion
% 
% OUTPUT ARGUMENTS:
%    TC = total cost
%       = sum(k(y)) + sum(sum(C(X)))
%     X = n x m logical matrix, where X(i,j) = 1 if EF j allocated to NF i

[n,m] = size(C);
if isscalar(k), k = repmat(k,n,1); else k = k(:); end

% Calculate X:
% initialize X
X = zeros(size(C)); 

% find the facility that is nearest to each costumers in the set of facilities are open
C_1 = repmat(y, 1, m) .* C;
C_1(find(C_1 == 0)) = Inf;
X(sub2ind(size(X), argmin(C_1, 1), [1:m])) = 1;

% Calculate TC:
TC = sum(k .* y) + sum(sum(C .* X));

end