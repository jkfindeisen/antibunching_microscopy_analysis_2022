function Bnp = multiscale_test(Y, Psi, dist_params, alpha)
% Performs the multiscale test to infer on the support of the underlying mean
%
% Y     data, generated from the sample
% Psi   the test itself
% dist_params   struct containing the data distribution parameters
% alpha         the level, i.e. we ensure (asymptotically) that in
%               100*(1-alpha)% there are no false positives

assert(nargin == 4, 'Not enough arguments');

%% Preparations

if ~isfield(Psi,'realizations')
    Psi = simulate_gaussian_approximation(Psi.filename, 1e3);
end

if ~isfield(Psi,'Phi_FFTd')
    Psi = generate_testfunctions(Psi);
end

% Estimate variance
V = estimate_variance(Y, dist_params);
if ~isscalar(V)
    V_FFTd = fft2(padarray(V,[ceil(Psi.n/2) ceil(Psi.n/2)]));
end

% Generate t set
tset = generate_tset(Y,Psi.hset,dist_params);

% Compute quantile
qalpha = quantile(max(Psi.pen(Psi.realizations),[],1)',1-alpha);

%% perform test and store detected regions
Y_FFTd = fft2(padarray(Y,[ceil(Psi.n/2) ceil(Psi.n/2)]));

Bnp = cell(numel(Psi.hset),1);

% Transformed quantile
qh = Psi.pen_inv(qalpha);
no_boxes = 0;
for i = 1:numel(Psi.hset)
    % Compute statistic for this scale
    if isscalar(V)
        % estimated variance does not depend on t
        Th = mean_values_fixed_scale(Y_FFTd,Psi.Phi_FFTd{i},Psi.n,Psi.hset{i})/sqrt(V)/Psi.Phi_norms{i};
    else
        % Due to heterogeneity, the estimated variance depends on t
        Th = mean_values_fixed_scale(Y_FFTd,Psi.Phi_FFTd{i},Psi.n,Psi.hset{i})./sqrt(mean_values_fixed_scale(V_FFTd,Psi.Phisq_FFTd{i},Psi.n,Psi.hset{i}));
    end
    % Use statistic only on tset{i}
    Th(~tset{i}) = -Inf;
    [t1,t2] = ind2sub([Psi.n Psi.n]-Psi.hset{i}+1,find(Th>qh(i)));
    Bnp{i} = [t1 t2];
    no_boxes = no_boxes + size(Bnp{i},1);
end
fprintf(' --> %i significant boxes found\n',no_boxes);

end