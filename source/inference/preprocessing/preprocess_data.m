function [s,Cov,it] = preprocess_data(Y,t,lambda,params)
% Preprocesses the antibunching data by computing pixel-wise the corresponding
% PSF order convolution images
% Y --> cell with raw data including 0 photon counts
% no_pulses --> number of illumination pulses
% lambda --> estimated background counts
% params - Parameters:
%       md --> number of detectors used for measurements
%       mD --> number of channels used in reconstruction (<= md)
%       mQ --> number of Q's used in approximation
%       mS --> number of S's used in approximation formula for Q
%       ms --> number of s's to be determined
%       mp --> number of p's or lambda's used in formula for Q
%       NOTE: meaninful choices should satisfy mQ,mp <= mS <= ms
%       maxit_preprocessing --> maximal Number of Newton iterations used for preprocessing
%       tol_preprocessing --> Relative tolerance for Newton's method for preprocessing
% s - convolution images
% Cov - estimated covariance matrix of the output
%
% (c) Frank Werner, IMS Uni Göttingen, 25.02.2019 or later

% Check consistency
if params.mD>params.md
    error('One cannot use more channels for reconstruction than actually measured!');
end

% keep Newton as default
if ~isfield(params, 'preprocessing_inversion_method')
    params.preprocessing_inversion_method = 'newton';
end

% Determine highest order needed
params.max = max([params.mQ params.mS params.ms params.mp]);
% Precompute factorials (to save time)
params.f = factorial(0:params.max);  % f(i) = (i-1)!

A = setup_transition_matrix(params.mQ, params.md);
A = A(1:params.mD+1,:);

% Now proceed pixelwise: Normalize data, estimate covariance matrix,
% preprocess, transform covariance matrix

s = cell(params.ms,1);
it = zeros(size(Y{1}));
Cov = cell(size(Y{1})); % Each pixel has its own estimated covariance matrix
for i=1:size(Y{1},1)
    for j=1:size(Y{1},2)
        % Estimate probabilities on this pixel
        P = zeros(params.mD+1,1);
        for k=1:params.mD+1
            P(k) = Y{k}(i,j)/t;
        end
        % Estimate also covariance matrix on this pixel
        Sigma = diag(P.*(1-P))/t;
        for k=1:params.mD+1
            for l=1:params.mD+1
                if l~=k
                    Sigma(k,l) = -P(k) * P(l) / t;
                end
            end
        end
        % Preprocess pixel
        switch params.preprocessing_inversion_method
            case 'newton'
                [aux,Cov{i,j}] = preprocess_pixel_newton(P,Sigma,lambda(i,j)/t,A,params);
            case 'lsqnonlin'
                [aux,Cov{i,j}] = preprocess_pixel_lsqnonlin(P, Sigma, lambda(i,j)/t, A, params);
            case 'fmincon'
                [aux,Cov{i,j}] = preprocess_pixel_fmincon(P, Sigma, lambda(i,j)/t, A, params);
            case 'gauss_newton'
                [aux,Cov{i,j},it(i,j)] = preprocess_pixel_gauss_newton(P,Sigma,lambda(i,j)/t,A,params);
            otherwise
                error('Unknown pixel preprocessing method!');
        end
        
        for k=1:params.ms
            s{k}(i,j) = aux(k);
        end
    end
end

end

function [s,Cov] = preprocess_pixel_fmincon(P,Sigma,lambda,A,params)
% stripped down variant of preprocess_pixel_lsqnonlin

% internal order
o = params.mQ;

% scale s to bring to same value
% k = 0.05.^(1:o);
k = 0.1.^(1:o);

% initial conditions
s0 = zeros(1, o);

% solver options (global variable to save a bit of time)
persistent opt;
if isempty(opt)
    opt = optimoptions('fmincon', 'display', 'off', 'algorithm', 'interior-point');
end
% fmincon optimizer (LSE)
minimizer = @(s) sum((A*g(s.*k,lambda,params)-P).^2);
s = fmincon(minimizer, s0, [], [], [], [], [], [], [], opt);
s = s .* k;

% Now we have to update Cov. We basically solve
% dgds^T(s)*A^T (A*g(s) - P) = 0, hence by the Delta method (similar to MLE
% routine):
tmp = A*dgds(s,lambda, params);
tmp = tmp(:, 1:params.ms);
Cov = (tmp'*tmp)\(tmp'*Sigma*tmp)/((tmp'*tmp)');

% truncate as to desired length
s = s(1:params.ms);
end

function [s,Cov] = preprocess_pixel_lsqnonlin(P,Sigma,lambda,A,params)
% internal order
o = params.mQ;

% scale s to bring to same value
if params.preprocessing_inversion.scale_s
    k = 0.05.^(1:o);
else
    k = ones(1, o);
end

% initial conditions
switch params.preprocessing_inversion.initial_s
    case 'zeros'
        s0 = zeros(1, o);
    case 'one'
        s0 = 17*0.013.^(1:o) ./ k; % random start point n=17,p=0.013
    case 'many'
        s0 = [17*0.013.^(1:o) ./ k; 3*0.048.^(1:o) ./ k; 87*0.032.^(1:o) ./ k; 31*0.009.^(1:o) ./ k]; % four random start points
    otherwise
        error('unknown initial s');
end

% bounds
if params.preprocessing_inversion.bounded
    lb = zeros(1, o); % s are positive
    ub = 200 * 0.05.^(1:o) ./ k; % not more than n=200,p=0.05
else
    lb = [];
    ub = [];
end

% solver
switch params.preprocessing_inversion.solver
    case 'fmincon'
        opt = optimoptions('fmincon', 'display', 'off', 'algorithm', 'interior-point');
        % fmincon optimizer (LSE)
        minimizer = @(s) sum((A*g(s.*k,lambda,params)-P).^2);
        fs = Inf;
        % loop over all start values
        for i = 1 : size(s0, 1)
            [x, fval] = fmincon(minimizer, s0(i, :), [], [], [], [], lb, ub, [], opt);
            % keep the best one
            if fval < fs
                fs = fval;
                s = x;
            end
        end
        s = s .* k;
    case 'lsqnonlin'
        assert(params.preprocessing_inversion.bounded == false); % otherwise you get: The Levenberg-Marquardt algorithm does not handle bound constraints and the trust-region-reflective algorithm requires at least as many equations as variables; aborting.
        opt = optimoptions('lsqnonlin', 'display', 'off', 'algorithm', 'levenberg-marquardt');
        fs = Inf;
        % loop over all start values
        for i = 1 : size(s0, 1)
            [x, resnorm] = lsqnonlin(@(s) A*g(s.*k,lambda,params)-P, s0(i, :), lb, ub, opt);
            % keep the best one
            if resnorm < fs
                fs = resnorm;
                s = x;
            end
        end
        s = s .* k;
    otherwise
        error('unknown solver');
end

% Now we have to update Cov. We basically solve
% dgds^T(s)*A^T (A*g(s) - P) = 0, hence by the Delta method (similar to MLE
% routine):
tmp = A*dgds(s,lambda, params);
Cov = (tmp'*tmp)\(tmp'*Sigma*tmp)/((tmp'*tmp)');

% truncate as to desired length
s = s(1:params.ms);
end

function [s,Cov] = preprocess_pixel_newton(P,Sigma,lambda,A,params)

% First apply A^{-1} to the data
Q = zeros(params.max+1,1);
Q(1:params.mQ+1) = A\P;

Cov = zeros(params.max+1,params.max+1);
Cov(1:params.mQ+1,1:params.mQ+1)=A\Sigma/A';

% Now we have to invert the function (s_1,...,s_{params.ms}) -> g(s_1,...,s_{params.ms},lambda)
% on the Q's.
% Here we use Newton's method
s0 = zeros(params.max,1);

% optimization loop
s = s0;
res = g(s, lambda,params) - Q;
it = 0;
res0 = norm(res,2);
while (norm(res,2)>params.tol_preprocessing*res0 && it<params.maxit_preprocessing)
    s = s - dgds(s,lambda, params)\res;
    res = g(s,lambda, params) - Q;
    it = it+1;
end

if it == params.maxit_preprocessing
    % Newton's method did not converge. Use initial guess.
    s = s0;
end

% Now we have to update Cov. Delta method:
% Cov <- (D(g^(-1))) Cov (D(g^(-1)))^T
% Inverse function theorem:
% D(g^(-1)) = (Dg)^(-1)
% Hence
tmp = dgds(s,lambda, params);
Cov = tmp\Cov/tmp';

%Truncate as to desired length
s = s(1:params.ms);
end

function [s,Cov,it] = preprocess_pixel_gauss_newton(P,Sigma,lambda,A,params)
% We have to invert the function (s_1,...,s_{params.ms}) -> A*g(s_1,...,s_{params.ms},lambda)
% on the P's.
% We relax this problem by minimizing ||A*g(s) - P||^2 via the Gauss-Newton method
% Further we scale s as s*k with the vector k = 0.05.^(1:params.max);
k = ones(params.max,1); %0.05.^(1:params.max);

% Starting value
s0 = zeros(params.max,1);

% optimization loop
s = s0;
res = A*g(s.*k, lambda,params) - P;
it = 0;
update = 1;
while (norm(update,2)>params.tol_preprocessing && it<params.maxit_preprocessing)
    dg = A*dgds(s.*k,lambda,params)*diag(k);
    update = (dg'*dg)\(dg'*res);
    s = s - update;
    res = A*g(s.*k,lambda, params) - P;
    it = it+1;
end
s = s.*k;

% Now we have to update Cov. We basically solve
% dgds^T(s)*A^T (A*g(s) - P) = 0, hence by the Delta method (similar to MLE
% routine):
tmp = A*dgds(s,lambda, params);
Cov = (tmp'*tmp)\(tmp'*Sigma*tmp)/((tmp'*tmp)');

%Truncate as to desired length
s = s(1:params.ms);
end