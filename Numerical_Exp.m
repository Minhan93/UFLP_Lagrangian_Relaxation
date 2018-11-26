% 1st Numerical Example (Small Scale)
clear

k = [100 100 100 200 200]';
D = gallery('integerdata', 1000, [5, 8], 1);
h = [5 10 5 10 5 10 5 10];
C = repmat(h, size(D, 1), 1).*D;
lambda = gallery('normaldata', [1, 8], 10) * 1000 + 2000;
alpha = 2;

[X, y, zLR, beta] = UFLP_LR_Sub(k, C, lambda);
[X_feas, z] = UFLP_Feas(k, C, y);
[lambda_new, delta] = UFLP_Upda_Mult(lambda, X, zLR, z, alpha);

IGDM = 'SubGradient'; % Improved Gradient Descent method: SubGradient
st = 120; % determin maximum running time length (unit: sec)

% run Lagrangian Relaxation Algorithm of UFLP
run_UFLP_LR

UB = z_opt; % current minimum upper bound (i.e. current minmum feasible solution)
LB = zLR; % current maximum lower bound (i.e. solution of UFLP-LR subproblem)
gap = Gap(end); % gap bt UB and LB


% 2nd Numerial Example (Large Scale)
clear
% UFLPP-LR
load data

lambda = gallery('uniformdata', [1, size(D_us, 2)], 10) * 1 * k + 3 * k;
k = repmat(k, size(D_us, 1), 1);
C = repmat(f_us, size(D_us, 1), 1).*D_us;
alpha = 2;

% Improved Gradient Descent methods
IGDMS = {'SubGradient', 'Momentum', 'Nesterov Momentum', 'AdaGrad', 'RMSProp', 'Adam'};  

GapC = {}; % cell of Gap
Z_optC = {}; % cell of Z_opt
ZLRC = {}; % cell of ZLR

for i = 1:length(IGDMS)
    IGDM = IGDMS{i}; % Assign Improved Gradient Descent method
    st = 10; % determin maximum running time length (unit: sec)
    
    % run Lagrangian Relaxation Algorithm of UFLP
    run_UFLP_LR    
    
    GapC{end + 1} = Gap;
    Z_optC{end + 1} = Z_opt;
    ZLRC{end + 1} = ZLR;
end
% save result_120s GapC Z_optC ZLRC

% Heuristic result
load data

r = 1;
C_us = r*(f_us(:)'.*D_us);

tic, [yufl, TCufl, Xufl] = getupper(k,C_us);toc

% Gurobi
load data

r = 1;
C_us = r*(f_us(:)'.*D_us);

% Create MILP model of UFLP
mp = MILP_model_of_UFLP(C_us,k);

% Solve using Gurobi
clear model params
model = mp.milp2gb;
params.outputflag = 1;
result = gurobi(model,params);
x = result.x;
TC = result.objval;