function test_estimation_single_cluster()
% Tests the full algorithm in a simple case (a single cluster)
% Currently the number of molecules n is consistently understimated

close all;
fprintf('Repeated simulation and analysis of a single cluster.\n');

no_repetitions = 1e2;

% Set parameters
lambda = 0;
n0 = 10;
p0 = 0.02;
no_pulses = 3e2;

data.t_STED = no_pulses;
data.t_CONF = no_pulses;

% Computed by a fitting with h from above
params.a = 2;
params.b = 0.0673;
% Background is given here
params.lambda = lambda;
% Other parameters
params.pixelsize_in_nm = 10;
params.FWHM_STED = 7;
params.FWHM_CONF = 23;
params.md = 4;
% Remaining parameters are used with default values

% Level
alpha = 0.1;

% load test (comment out and use empty Psi for creating a new test)
% Psi = [];
load('miscat/test_single_cluster2.mat');

% Call main routine
results = cell(no_repetitions, 1);
pbar = progressbar('estimate molecular map');
for i = 1 : no_repetitions
    % Generate synthetic data ( using the full_forward_simulation_cluster)
    % Confocal
    h  = gaussian([64,64],params.FWHM_CONF);
    Y_CONF = full_forward_simulation_cluster(n0,p0,lambda,h,params.md,no_pulses);
    % STED
    h  = gaussian([64,64],params.FWHM_STED);
    Y_STED = full_forward_simulation_cluster(n0,p0,lambda,h,params.md,no_pulses);
    
    data.Y_CONF = Y_CONF;
    data.Y_STED = Y_STED;
    
    [segments, n, p, conf, aux] = estimate_molecular_map(data,params,Psi,alpha);
    %     [segments, n, p, conf, aux] = estimate_molecular_number_single_segment(data,params,alpha);
    r.segments = segments;
    r.n = n;
    r.p = p;
    r.conf = conf;
    r.aux = aux;
    results{i} = r;
    pbar.set(i / no_repetitions);
end
pbar.close();

% statistics
number_segments = cellfun(@(x) max(x.segments(:)), results);
% get those where number segments == 1
r = results(number_segments == 1);
n_estimated = cellfun(@(x) x.n{1}, r);
p_estimated = cellfun(@(x) x.p{1}, r);
confidence_hit = cellfun(@(x) n0 >= x.conf{1}(1) && n0 <= x.conf{1}(2), r);

% show results
figure();
subplot(2,2,1);
histogram(n_estimated);
title(sprintf('n = %d', n0));
subplot(2,2,2);
histogram(p_estimated);
title(sprintf('p = %.3f', p0));
subplot(2,2,3);
c = mean(confidence_hit);
bar([c, 1-c]);
title(sprintf('alpha = %.2f', alpha));

end