function test_algorithm_single_pixel_cluster()
% Tests the full algorithm in a simple case (a single cluster consisting of a single pixel), tries to
% make everything as simple as possible. (No STED, no background, no
% segmentation.)

no_repetitions = 100;

% Set parameters
lambda = 0;
n0 = 100;
p0 = 0.02;
no_pulses = 1e4;

data.t = no_pulses;
data.lambda = lambda;

% Computed by a fitting with h from above
params.a = 2;
params.b = 0.0673;
% Other parameters
params.pixelsize_in_nm = 4;
params.FWHM_CONF = 8;
params.md = 8;
L = 32;

% Level
alpha = 0.1;

% Call main routine
results = cell(no_repetitions, 1);
pbar = progressbar('estimate molecular map');
for i = 1 : no_repetitions
    % Generate synthetic data ( using the full_forward_simulation_cluster)
    % Confocal
    h  = 1; %gaussian([L,L],params.FWHM_CONF);
    Y_CONF = full_forward_simulation_cluster(n0,p0,lambda,h,params.md,no_pulses);
    
    data.img = Y_CONF;
    
    %     [segments, n, p, conf, aux] = estimate_molecular_map(data,params,Psi,alpha);
    [n, p, conf, aux] = estimate_molecular_number_single_segment(data,params,alpha);
    %r.segments = segments;
    r.n = n;
    r.p = p;
    r.conf = conf;
    r.aux = aux;
    results{i} = r;
    pbar.set(i / no_repetitions);
end
pbar.close();

% statistics
n_estimated = cellfun(@(x) x.n{1}, results);
p_estimated = cellfun(@(x) x.p{1}, results);
confidence_hit = cellfun(@(x) n0 >= x.conf{1}(1) && n0 <= x.conf{1}(2), results);
confidence_range = cellfun(@(x) diff(x.conf{1}), results);

% show results
fig = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1);
histogram(n_estimated, 0:1:2*n0);
title(sprintf('n_0=%d, median(n est.)=%.1f, mean(n est.)=%.1f, std(n est.)=%.1f', n0, median(n_estimated), mean(n_estimated), std(n_estimated)));
subplot(1,2,2);
histogram(p_estimated, p0-0.02:0.001:p0+0.02);
title(sprintf('p_0=%.3f, median(p est.)=%.3f, mean(p est.)=%.3f', p0, median(p_estimated), mean(p_estimated)));
annotation(fig, 'textbox', [0.2 0.96 0.6 0.03],'String', {sprintf('single cluster, %d pulses, %d repetitions, alpha=%.1f, n0 in confidence interval %.2f, median(confidence range)=%.1f', no_pulses, no_repetitions, alpha, mean(confidence_hit), median(confidence_range))}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
saveas(fig, 'cluster_sim.jpg');

end