function [y,TC,X] = exchangecheck(k,C,y)
%UFLXCHG Exchange best improvement procedure for uncap. facility location.
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
   ny = 1:n; ny(y) = [];
   for i = 1:length(y)
      is = i == idx;
      c1i = c1; c1i(is) = min(C(is,y([1:i-1 i+1:end])),[],2);
      for j = 1:length(ny)
         TCij = k1 - k(y(i)) + k(ny(j)) + sum(min(c1i,C(:,ny(j))));
         if TCij < TC1
            [i1,j1,TC1] = deal(i,j,TCij);
         end
      end
   end
   if TC1 < TC
      [ny(j1),y(i1)] = deal(y(i1),ny(j1));
      TC = TC1;
   else
      done = true;
   end
end

y = sort(y);
X = logical(sparse(y(argmin(C(:,y),2)),1:m,1,n,m));
