function [lambda_new, delta] = UFLP_Upda_Mult(lambda, X, LB, UB, alpha, varargin)
% Update Lagrange multipliers by applying different Improved Gradient Descent method.
% 
% [lambda_new, delta] = UFLP_Upda_Mult(lambda, X, LB, UB, alpha) % using SubGradient by default
%                     = UFLP_Upda_Mult(lambda, X, LB, UB, alpha, 'SubGradient', delta0) % using SubGradient
%                     = UFLP_Upda_Mult(lambda, X, LB, UB, alpha, 'Momentum', delta0) % using Momentum
%                     = UFLP_Upda_Mult(lambda, X, LB, UB, alpha, 'Nesterov Momentum', delta0) % using Nesterov Momentum
%                     = UFLP_Upda_Mult(lambda, X, LB, UB, alpha, 'AdaGrad', delta0) % using AdaGrad
%                     = UFLP_Upda_Mult(lambda, X, LB, UB, alpha, 'RMSProp', delta0) % using RMSProp
%                     = UFLP_Upda_Mult(lambda, X, LB, UB, alpha, 'Adam', delta0) % using Adam
% 
% INPUT ARGUMENTS:
%     X = n x m logical matrix, where X(i,j) = 1 if EF j allocated to NF i
%lambda = m-element Lagrange mutiplier vector, 
%         where lambda(j) is the penalty if EF j cannot be allocated to 
%         any NF
%    LB = lower bound (optimal objective value of UFLP-LR) at current 
%         iteration
%    UB = best-known upper bound
% alpha = constant used in step-size calculation
%delta0 = initial step size
% 
% OUTPUT ARGUMENTS:
%lambda_ new = new Lagrange multipliers
%      delta = the step size at thie iteration
persistent  grad_squared;
if nargin==5
    delta = alpha * (UB - LB)/sum((1 - sum(X, 1)).^2);
    lambda_new = lambda + delta * (1 - sum(X, 1));
else
    switch varargin{1}
        case 'SubGradient'
            delta = alpha * (UB - LB)/sum((1 - sum(X, 1)).^2);
            lambda_new = lambda + delta * (1 - sum(X, 1));
        case 'Momentum'
            rho = 0.9; % friction rate
            delta = (1 - sum(X, 1))+rho*varargin{2};
            lambda_new = lambda + delta * alpha;
        case 'Nesterov Momentum'
            rho = 0.9;
            old_delta=varargin{2};
            delta=rho* old_delta- alpha*(1 - sum(X, 1));
            lambda_new=lambda+rho*old_delta-(1+rho)*delta;
        case 'AdaGrad'
            if isempty(grad_squared)
                grad_squared=0;
            end
            grad_squared=grad_squared+(1 - sum(X, 1)).^2;
            lambda_new=lambda+alpha*(1 - sum(X, 1))./((grad_squared.^0.5)+(1e-7));
            delta=[];
        case 'RMSProp'
            if isempty(grad_squared)
                grad_squared=0;
            end
            decay_rate=0.5;
            grad_squared=decay_rate*grad_squared+(1-decay_rate)*(1 - sum(X, 1)).^2;
            lambda_new=lambda+alpha*(1 - sum(X, 1))./((grad_squared.^0.5)+(1e-7));
            delta=[];
        case 'Adam'
            beta1=0.9;beta2=0.999;
            first_moment=0; second_moment=0;
            for t=1:3
                first_moment=beta1* first_moment+(1-beta1)*(1 - sum(X, 1));
                second_moment=beta2*second_moment+(1-beta2)*(1 - sum(X, 1)).^2;
                first_unbias=first_moment/(1-beta1.^t);
                second_unbias=second_moment/(1-beta2.^t);
            end
            lambda_new=lambda+alpha*first_unbias./(second_unbias.^0.5+1e-7);
            delta=[];
        otherwise
            error('wrong input')
    end
end
end

