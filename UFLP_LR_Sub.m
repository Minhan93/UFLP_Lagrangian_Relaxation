function [X, y, TC, beta] = UFLP_LR_Sub(k, C, lambda)
% The Subproblem of UFLP Lagrangian Relaxation Algorithm.
% 
% [X, y, TC, beta] = UFLP_LR_Sub(k, C, lambda)
% 
% INPUT ARGUMENTS:
%     k = n-element fixed cost vector, where k(i) is cost of NF at Site i
%         (if scalar, then same fixed cost)
%     C = n x m variable cost matrix,
%         where C(i,j) is the cost of serving EF j from NF i
%lambda = m-element Lagrange mutiplier vector, 
%         where lambda(j) is the penalty if EF j cannot be allocated to 
%         any NF
% 
% OUTPUT ARGUMENTS:
%     y = NF site index vector
%    TC = total cost
%       = sum(k(y)) + sum(sum(C(X))) + penalty
%     X = n x m logical matrix, where X(i,j) = 1 if EF j allocated to NF i
%  beta = m-element benefit vector, 
%         where beta(i) is benefit of building NF i

[n,m] = size(C);
if isscalar(k)
    k = repmat(k,n,1); 
else
    k = k(:);
end

% Calculate beta:
beta = sum(min(0, C - repmat(lambda, n, 1)), 2);

% Calculate y:
y = beta + k < 0;

% Calculate X:
X = C - repmat(lambda, n, 1) < 0 & repmat(y, 1, m) == 1;

% Calculate TC:
TC = sum(k .* y) + sum(sum((C - repmat(lambda, n, 1)) .* X)) + sum(lambda);

% [n,m] = size(C);
% if isscalar(k), k = repmat(k,n,1); else k = k(:); end
% beta = sum(min(0, C - repmat(lambda, size(C, 1), 1)), 2);
% y = beta + k < 0;
% x = C - repmat(lambda, size(C, 1), 1) < 0 & repmat(y, 1, size(C, 2)) == 1;
% z = sum(k .* y) + sum(sum((C - repmat(lambda, size(C, 1), 1)) .* x)) + sum(lambda);
end