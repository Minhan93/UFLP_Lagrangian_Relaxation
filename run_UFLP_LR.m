tic
% initialize
[X, y, zLR, beta] = UFLP_LR_Sub(k, C, lambda);

[X_feas, z] = UFLP_Feas(k, C, y);
z_opt = z;
X_opt = X_feas;
y_opt = y;

Gap = (z_opt - zLR) / z_opt;
Z_opt = z_opt;
ZLR = zLR;

while Gap(end) > 0.005 && toc < st
    % update lambda
    [lambda1, delta] = UFLP_Upda_Mult(lambda, X, zLR, z_opt, alpha, IGDM, 0);
    
    % update X, y and zLR
    [X1, y1, zLR1, beta1] = UFLP_LR_Sub(k, C, lambda1);
    
    % update z
    [X_feas1, z1] = UFLP_Feas(k, C, y1);
    if z1 < z_opt
        z_opt = z1;
        X_opt = X_feas1;
        y_opt = y1;
    end
    
    count = 0;
    while zLR1 <= zLR
        if count == 20
            count = 0;
            alpha = alpha / 2; % update alpha
        else
            count = count + 1;
        end
        [lambda1, delta] = UFLP_Upda_Mult(lambda1, X1, zLR1, z_opt, alpha, IGDM, delta); % update lambda1
        [X1, y1, zLR1, beta1] = UFLP_LR_Sub(k, C, lambda1); % update X1, y1, and zLR1
        [X_feas1, z1] = UFLP_Feas(k, C, y1);% update z_opt
        if z1 < z_opt
            z_opt = z1;
            X_opt = X_feas1;
            y_opt = y1;
        end
    end
    
    zLR = zLR1; % if zLR1 > zLR, then update zLR
    lambda = lambda1; % if zLR1 > zLR, then update lambda
    Gap(end + 1) = (z_opt - zLR) / z_opt; % update gap
    Z_opt(end + 1) = z_opt;
    ZLR(end + 1) = zLR;
end