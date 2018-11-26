function [y,TC,X] = addcheck(k,C,y,p)
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

if isempty(y)
   [TC,y] = min(sum(C,1) + k);
else
   TC = sum(k(y)) + sum(min(C(:,y),[],2));
end

ny = 1:n; ny(y) = [];
TC1 = Inf;
done = false;
while ~done
   c1 = min(C(:,y),[],2);
   k1 = sum(k(y));
   for i = 1:length(ny)
      TCi = k1 + k(ny(i)) + sum(min(c1,C(:,ny(i))));
      if TCi < TC1
         [i1,TC1] = deal(i,TCi);
      end
   end
   if (isempty(p) && TC1 < TC) || (~isempty(p) && length(y) < p)
      y = [y ny(i1)];
      ny(i1) = [];
      TC = TC1;
      if ~isempty(p), TC1 = Inf; end
   else
      done = true;
   end
end

y = sort(y);
X = logical(sparse(y(argmin(C(:,y),2)),1:m,1,n,m));

if isinf(TC)
    y = [];
    X = [];
end
