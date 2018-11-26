function [y,TC,X] = dropcheck(k,C,y,p)
%UFLP adding node heuristic

C = C';  % Make column-based to speed up minimization
[m,n] = size(C);
if isscalar(k)
    k = repmat(k,1,n);
else
    k = k(:)';
end

if nargin < 3
    y = [];
else
    y = y(:)';
end
if nargin < 4
    p = [];
end


TC = sum(k(y)) + sum(min(C(:,y),[],2));

TC1 = Inf;
if length(y) > 1, done = false; else done = true; end
while ~done
   [c1,idx] = min(C(:,y),[],2);
   k1 = sum(k(y));
   for i = 1:length(y)
      is = i == idx;
      TCi = k1 - k(y(i)) + sum(c1(~is)) + ...
         sum(min(C(is,y([1:i-1 i+1:end])),[],2));
      if TCi < TC1
         [i1,TC1] = deal(i,TCi);
      end
   end
   if (isempty(p) && TC1 < TC) || (~isempty(p) && length(y) > p)
      y(i1) = [];
      TC = TC1;
      if ~isempty(p), TC1 = Inf; end
   else
      done = true;
   end
end

y = sort(y(:)');
X = logical(sparse(y(argmin(C(:,y),2)),1:m,1,n,m));
