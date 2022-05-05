function test_algorithm_cluster2()
% Tests the full algorithm in a simple case (a single cluster), tries to
% make everything as simple as possible. (No STED, no background, no
% segmentation.)

% close all;

no_repetitions = 1000;

% Set parameters
lambda = 0;
n0 = 10;
p0 = 0.02;
no_pulses = 1e3*4;

data.t_CONF = no_pulses;

% Computed by a fitting with h from above
params.a = 2;
params.b = 0.0673;
% Other parameters
params.pixelsize_in_nm = 4;
params.FWHM_CONF = 8/2;
params.md = 4;
L = 32;

% Level
alpha = 0.1;

% Call main routine
results = cell(no_repetitions, 1);
pbar = progressbar('estimate molecular map');
for i = 1 : no_repetitions
    % Generate synthetic data ( using the full_forward_simulation_cluster)
    % Confocal
    h  = gaussian([L,L],params.FWHM_CONF);
    img = full_forward_simulation_cluster(n0,p0,lambda,h,params.md,no_pulses);
    
    data.img = img;
    data.t = no_pulses;
    data.lambda = zeros(size(img{1})) + lambda;
    
    %     [segments, n, p, conf, aux] = estimate_molecular_map(data,params,Psi,alpha);
    [n, p, conf, iterations] = estimate_molecular_number_single_segment(data,params,alpha);
    r.n = n;
    r.p = p;
    r.conf = conf;
    r.aux = iterations;
    results{i} = r;
    pbar.set(i / no_repetitions);
end
pbar.close();

% statistics
n_estimated = cellfun(@(x) x.n{1}, results);
p_estimated = cellfun(@(x) x.p{1}, results);
% confidence_hit = cellfun(@(x) n0 >= x.conf{1}(1) && n0 <= x.conf{1}(2), results);
% confidence_range = cellfun(@(x) diff(x.conf{1}), results);

% show results
fig = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1);
histogram(n_estimated, 0:1:40);
title(sprintf('n_0=%d, median(n est.)=%.1f, mean(n est.)=%.1f, std(n est.)=%.1f', n0, median(n_estimated), mean(n_estimated), std(n_estimated)));
subplot(1,2,2);
histogram(p_estimated, p0-0.02:0.001:p0+0.02);
title(sprintf('p_0=%.3f, median(p est.)=%.3f, mean(p est.)=%.3f, std(p est.)=%.3f', p0, median(p_estimated), mean(p_estimated), std(p_estimated)));
% annotation(fig, 'textbox', [0.2 0.96 0.6 0.03],'String', {sprintf('single cluster, %d pulses, %d repetitions, alpha=%.1f, n0 in confidence interval %.2f, median(confidence range)=%.1f', no_pulses, no_repetitions, alpha, mean(confidence_hit), median(confidence_range))}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
annotation(fig, 'textbox', [0.2 0.96 0.6 0.03],'String', {sprintf('single cluster, %d pulses, %d repetitions', no_pulses, no_repetitions)}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
saveas(fig, 'cluster_sim.jpg');

% fig = figure();
% np_estimated = n_estimated .* p_estimated;
% histogram(np_estimated, 0.1:0.005:0.3);
% xlim([0.1,0.3]);
% xlabel('np');
% title(sprintf('n0*p0=%.2f, mean(n-est*p-est)=%.3f, mean(n-est)=%.1f, mean(p-est)=%.3f', n0*p0, mean(np_estimated), mean(n_estimated), mean(p_estimated)));
end