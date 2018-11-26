function [y,TC,X] = getupper(k,C)
%UFLP Heuristically get a upperbound.
[n,m] = size(C);
[y1,TC1] = addcheck(k,C); 
fprintf('  Add: %f\n',TC1)
done = false;
while ~done
   [y,TC] = exchangecheck(k,C,y1);
   fprintf(' Xchg: %f\n',TC)
   if ~isequal(y,y1)
      [y1,TC1] = addcheck(k,C,y);
      fprintf('  Add: %f\n',TC1);
      [y2,TC2] = dropcheck(k,C,y);
      fprintf(' Drop: %f\n',TC2)
      if TC2 < TC1
          TC1 = TC2;
          y1 = y2;
      end
      if TC1 >= TC
          done = true;
      end
   else
      done = true;
   end
end
fprintf('Final: %f\n',TC)

y = sort(y);
X = logical(sparse(y(argmin(C(y,:),1)),1:m,1,n,m));
