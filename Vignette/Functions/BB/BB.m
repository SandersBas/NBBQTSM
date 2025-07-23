function [theta_baseline, theta_mat_bb] = BB(df, psi, I, theta_init, B)
% Inputs are: 
% 1. df: dataframe in which each row consists of an observation X and
%        indicators for the units involved;
% 2. psi: moment function, takes in observation X and parameter theta and
%         outputs realized moments;
% 3. I: polyadic order (dyadic has I=2, triadic has I=3, etc);
% 4. theta_init: initial values for theta, dimension K x 1;
% 5. B: number of bootstrap draws.
%
% Outputs are:
% 1. theta_baseline: GMM point estimate for theta;
% 2. theta_mat_bb: B x K matrix of Bayesian bootstrap estimates;
 
% Obtain K
K = length(theta_init);
% Obtain data matrix
X = df(:,1:size(df,2)-I);
% Back out number of moment conditions
L = length(psi(X(1,:)',theta_init));
% Gather all involved units
units = sort(unique(df(:,size(X,2)+1:end)));
% Find number of exchangeable units
N = length(units);

%Initialize output
theta_mat_bb = zeros(B,K);

if K==L % just-identified case
    %% Just-identified case, baseline
    options = optimset('Display', 'off');
    justIdObj_func = @(theta) justIdObj(1/size(X,1).*ones(size(X,1),1),psi,X,theta);
    theta_baseline = fsolve(justIdObj_func, theta_init, options);

    %% Just-identified case, Bayesian bootstrap
    parfor b=1:B
        V = exprnd(1,N,1);
        W = V./sum(V);       
        
        weight = zeros(size(X,1),1);
        for o = 1:size(X,1)
            weight_o = 1;
            for i = 1:I
                unit = df(o,size(X,2)+i);
                weight_o = weight_o * W(units==unit);
            end      
             weight(o) = weight_o;
        end   
        normalization  = sum(weight);     
        weight = weight./normalization;

        justIdObj_func_b = @(theta) justIdObj(weight,psi,X,theta);
        theta_mat_bb(b,:) = fsolve(justIdObj_func_b, theta_baseline, options);
        
        if mod(b,100)==0
            b/B
        end
    end  
elseif L>K % over-identified case
    %% Over-identified case, baseline
    options = optimset('Display', 'off');
    
    % First step
    oneStepObj_func = @(theta) gmmObj(1/size(X,1).*ones(size(X,1),1),psi,X,theta,eye(L));
    theta_baseline_os = fminsearch(oneStepObj_func, theta_init, options);
    
    % Weight matrix
    Weight_mat = inv(varMat(1/size(X,1).*ones(size(X,1),1),psi,X,theta_baseline_os));
    
    % Second step
    twoStepObj_func = @(theta) gmmObj(1/size(X,1).*ones(size(X,1),1),psi,X,theta,Weight_mat);
    theta_baseline = fminsearch(twoStepObj_func, theta_init, options);

    %% Over-identified case, Bayesian bootstrap
    parfor b=1:B
        V = exprnd(1,N,1);
        W = V./sum(V);       
        
        weight = zeros(size(X,1),1);
        for o = 1:size(X,1)
            weight_o = 1;
            for i = 1:I
                unit = df(o,size(X,2)+i);
                weight_o = weight_o * W(units==unit);
            end      
             weight(o) = weight_o;
        end   
        normalization  = sum(weight);
        weight = weight./normalization;
        
        % First step
        oneStepObj_func_b = @(theta) gmmObj(weight,psi,X,theta,eye(L));
        theta_os_b = fminsearch(oneStepObj_func_b, theta_baseline_os, options);
        
        % Weight matrix
        Weight_mat_b = inv(varMat(weight,psi,X,theta_os_b));
        
        % Second step
        twoStepObj_func_b = @(theta) gmmObj(weight,psi,X,theta,Weight_mat_b);
        theta_mat_bb(b,:) = fminsearch(twoStepObj_func_b, theta_baseline, options);
    
        if mod(b,100)==0
            b/B
        end
    end
else
    error('L should be greater than or equal to K')
end
end