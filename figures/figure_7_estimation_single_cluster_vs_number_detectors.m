function figure_7_estimation_single_cluster_vs_number_detectors()
% runs single clusters vs. number of detectors and n,p for a high number of
% repetitions to detect possible biases

simulation_run();
display_results();

end

function display_results()

file = 'figure_7_results.mat';

t = tic();
results = load(file);
results = results.results;
toc(t);

N = numel(results);
z = zeros(N, 13);
for i = 1 : N
    r = results{i};
    values = analyse_single_result(r.results, r.params(1));
    z(i, :) = [r.params, values];
    % order of params: n, p, md
end

nmax = max(z(:, 1));
R = [0, nmax];

md = sort(unique(z(:, 3)));
Nd = numel(md);
p = sort(unique(z(:, 2)));

fig = figure();
set(fig,'defaultAxesTickLabelInterpreter', 'latex');
set(fig,'defaultBubblelegendInterpreter', 'latex');
set(fig,'defaultColorbarTickLabelInterpreter', 'latex');
set(fig,'defaultConstantlineInterpreter', 'latex');
set(fig,'defaultGraphplotInterpreter', 'latex');
set(fig,'defaultLegendInterpreter', 'latex');
set(fig,'defaultPolaraxesTickLabelInterpreter', 'latex');
set(fig,'defaultTextInterpreter', 'latex');
set(fig,'defaultTextarrowshapeInterpreter', 'latex');
set(fig,'defaultTextboxshapeInterpreter', 'latex');
fig.Position = [100, 100, 1800, 500];
hold on;
% diagonal line
plot(R, R, 'k-');
k = [0.5, 1, 0.5];
k = k / sum(k);

pi = p;
h = zeros(Nd, 1);
for j = 1 : Nd
    mdj = md(j);
    idx = z(:, 2) == pi & z(:, 3) == mdj;
    zi = z(idx, :);
    n0 = zi(:, 1);
    % median
    nest = zi(:, 5);
    %     mean
    %     nest = zi(:, 4);
    x = [0; n0];
    y = [0; nest];
    ys = sgolayfilt(y, 2, 9);
    y(y > 80) = ys(y>80);
    h(j) = plot(x, ys, '-', 'LineWidth', 2, 'DisplayName', sprintf('$m_d=%d$', mdj));
end

xlabel('N');
ylabel('$\hat N$');
xlim(R);
ylim(R);
xticks(0:50:200);
yticks(0:50:200);
daspect([1,1,1]);
grid on;
box on;
ax = gca;
ax.FontSize = 12;
legend(h, 'location', 'NorthWest', 'FontSize', 16);
exportgraphics(fig, 'figure_7_single_cluster_vs_number_detectors.pdf');
end


function r = analyse_single_result(results, n0)
% output order is: mean n-estimated, quantile n-estimated 0.15, 0.85, mean
% p-estimated, quantiles of p-estimates, mean confidence interval hists (or
% true n), mean relative confidence interval size (extent divided by mean)

n_est = cellfun(@(x) x.n{1}, results);
p_est = cellfun(@(x) x.p{1}, results);
confidence_hit = cellfun(@(x) n0 >= x.conf{1}(1) && n0 <= x.conf{1}(2), results);
rel_conf_int = cellfun(@(x) (x.conf{1}(2) - x.conf{1}(1)) / (x.conf{1}(2) + x.conf{1}(1)) / 2, results);

r = [mean(n_est), median(n_est), std(n_est), quantile(n_est, [0.15, 0.85]), mean(p_est), quantile(p_est, [0.15, 0.85]), mean(confidence_hit), mean(rel_conf_int)];

end



function simulation_run()

file = 'figure_7_results.mat';
warning('off');

n = 5:5:200;
pi = 0.02;
md = [4, 6, 8];
no_pulses = 2e4;

% results = [];
if exist(file, 'file')
    load(file);
else
    results = {};
end
N = numel(n) * numel(md);
c = 0;
for ni = n
    for mdi = md
        c = c + 1;
        fprintf('%d / %d - %s\n', c, N, datetime('now'));
        r.params = [ni, pi, mdi];
        already_done = false;
        for i = 1 : numel(results)
            if results{i}.params == r.params
                already_done = true;
                break;
            end
        end
        if already_done == true
            continue;
        end
        r.results = single_run(ni, pi, mdi, no_pulses);
        results = [results; {r}];
        save(file, 'results');
    end
end

end

function results = single_run(n0, p0, md, no_pulses)

% fundamental parameters
no_repetitions = 1000; % number of repetitions for the statistic

params.md = md; % number of detection channels
params.FWHM_CONF = 4;

% Level
alpha = 0.1;

% generate parameters for the forward simulation
R = 8;
dims = (2*R+1)*[1,1];
lambda = zeros(dims);
object = zeros(dims);
object(R+1,R+1)=n0;
psf = gaussian(dims,params.FWHM_CONF);

% Call main routine
results = cell(no_repetitions, 1);
pbar = progressbar(sprintf('n=%d, p=%.2f, md=%d', n0, p0, md));
for i = 1 : no_repetitions
    
    % perform forward simulation
    img = full_forward_simulation(object, p0, lambda, psf, params.md, no_pulses);
    data.img = img;
    data.t = no_pulses;
    data.lambda = lambda; % must be whole array
    
    % count
    [n, p, conf] = estimate_molecular_number_single_segment(data,params,alpha);
    
    % store results
    r.n = n;
    r.p = p;
    r.conf = conf;
    results{i} = r;
    pbar.set(i / no_repetitions);
end
pbar.close();

end