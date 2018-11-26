function mp = MILP_model_of_UFLP(C,k)
% Create MILP model of UFLP
clear mp
mp = Milp('UFLP');

[n,m] = size(C);
kn = iff(isscalar(k),repmat(k,1,n),k(:)');  % expand if k is constant value
mp.addobj('min',kn,C)  % min sum_i(ki*yi) + sum_i(sum_j(cij*xij))

for j = 1:m
   mp.addcstr(0,{':',j},'=',1)   % sum_i(xij) = 1
end
for i = 1:n
   mp.addcstr({m,{i}},'>=',{i,':'})  % m*yi >= sum_j(xij)  (weak form.)
end
mp.addub(1,1)
mp.addctype('B','C')         % only k are integer (binary)
end

