function test_algorithm_regular_grid()
% tests only the n-p estimation on a regular grid

close all;

% general parameter
params.md = 6;
params.mD = params.md;
params.mQ = params.md;
params.mS = params.md;
params.ms = params.md;
params.mp = params.md;
params.mc = 2;
params.maxit_preprocessing = 10;
params.maxit_counting = 10;
params.tol_preprocessing = 1e-4;
params.tol_counting = 1e-4;
% params.preprocessing_inversion_method = 'fmincon';
params.FWHM_CONF = 5;
params.no_neighboring_px = 0;

lambda = 0;
pc = 0.02;
nc = 100; % number of molecules per cluster
D = 3 * params.FWHM_CONF;
g = round(D+(0:4)*D); % cluster localizations
L = ceil(max(g)+D);
object = zeros(L);
T = 2e4;
% T = 1e4;
alpha = 0.1;
R = 10;

% read results
r = load('test_regular_grid_results_t=4e3.mat');
r = r.r;

N = numel(r);
ch = cellfun(@(x) all(x.ch == 1), r);
fprintf('in %.1f %% of %d cases all n within confidence\n', mean(ch)*100, N);

n = cellfun(@(x) x.n, r, 'UniformOutput', false);
n = cat(1, n{:});
p = cellfun(@(x) x.p, r, 'UniformOutput', false);
p = cat(1, p{:});
c = cellfun(@(x) x.c, r, 'UniformOutput', false);
c = cat(1, c{:});

fig = figure();
fig.Position = [100, 100, 900, 300];

subplot(1,3,1);
h = histogram(n(:), 40);
hold on;
plot(nc*[1,1], [0, max(h.Values(:))], 'r-', 'LineWidth', 1.5);
xlabel('est. n');
ylabel('occ.');
title(sprintf('median %.1f', median(n)));

subplot(1,3,2);
h = histogram(p(:), 40);
hold on;
plot(pc*[1,1], [0, max(h.Values(:))], 'r-', 'LineWidth', 1.5);
xlabel('est. p');
ylabel('occ.');
title(sprintf('median %.2f', median(p)));

subplot(1,3,3);
h = histogram(c(:, 1), 40);
hold on;
h = histogram(c(:, 2), 40);
xlabel('conf. intervals');
ylabel('occ.');
title(sprintf('medians %.1f, %.1f', median(c)));

exportgraphics(fig, 'test_regular_grid_20_simulations_t=4e3.png');

psf = gaussian([13, 13], params.FWHM_CONF);

% create object
% a = [1/2,1,1/2];
% a = a .* a.';
% a = a / sum(a(:)) * nc;
for i = g
    for j = g
        %         object(i+(-1:1), j+(-1:1)) = round(poissrnd(a));
        object(i, j) = nc;
    end
end
fprintf('object with %d molecules\n', sum(object(:)));

% segmentation
Rs = floor(min(diff(g))/2);
segments = zeros(L);
l = 0;
for i = g
    for j = g
        l = l + 1;
        segments(i-Rs:i+Rs, j-Rs:j+Rs) = l;
    end
end
Ns = l;

r = cell(R, 1);
for o = 1 : R
    fprintf('iteration %d/%d\n', o, R);
    
    % perform forward simulation
    Y = full_forward_simulation(object, pc, lambda, psf, params.md, T, [], true);
    
    % preprocess data
    lambda = lambda + zeros(L);
    [S_CONF,Cov] = preprocess_data(Y,T,lambda,params);
    
    % estimate number and brightness
    [n,p,conf] = estimate_np(S_CONF,Cov,segments,alpha/2,params); % Bonferroni adjustment
    ni = cellfun(@(x) x, n);
    pi = cellfun(@(x) x, p);
    c1 = cellfun(@(x) x(1), conf);
    c2 = cellfun(@(x) x(2), conf);
    
    % confidence hits
    ch = zeros(Ns, 1);
    for i = 1 : Ns
        ch(i) = ni(i) >= c1(i) && ni(i) <= c2(i);
    end
    
    r{o} = struct('n', ni, 'p', pi, 'c', [c1, c2], 'ch', ch);
    save('test_regular_grid_results_t=4e3.mat', 'r');
    
    if o == 1
        fig = figure();
        fig.Position = [100, 100, 1600, 600];
        subplot(2,5,1);
        imagesc(object);
        axis image;
        colorbar();
        title('true object');
        for i = 2 : 5
            subplot(2,5,i);
            img = Y{i};
            imagesc(img);
            axis image;
            title(sprintf('%d-photons', i-1));
            colorbar();
        end
        subplot(2,5,6);
        imagesc(segments);
        axis image;
        colorbar();
        title('segmentation');
        for i = 1 : 4
            subplot(2,5,6+i);
            img = S_CONF{i};
            imagesc(img);
            axis image;
            title(sprintf('S_%d', i));
            colorbar();
        end
        colormap(hot);
        exportgraphics(fig, 'test_regular_grid_single_simulation_data_t=4e3.png');
        
        fig = figure();
        
        subplot(2,2,1);
        h = histogram(ni(:), 40);
        hold on;
        plot(nc*[1,1], [0, max(h.Values(:))], 'r-', 'LineWidth', 1.5);
        xlabel('est. n');
        ylabel('occ.');
        title(sprintf('conf. hit prob. %.1f%%', mean(ch)*100));
        
        subplot(2,2,2);
        h = histogram(pi(:), 40);
        hold on;
        plot(pc*[1,1], [0, max(h.Values(:))], 'r-', 'LineWidth', 1.5);
        xlabel('est. p');
        ylabel('occ.');
        
        subplot(2,2,3);
        h = histogram(c1(:), 40);
        xlabel('left conf.');
        ylabel('occ.');
        title(sprintf('median %.1f', median(c1)));
        
        subplot(2,2,4);
        h = histogram(c2(:), 40);
        xlabel('right conf.');
        ylabel('occ.');
        title(sprintf('median %.1f', median(c2)));
        
        exportgraphics(fig, 'test_regular_grid_single_simulation_np_t=4e3.png');
    end
    
end

end