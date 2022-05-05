function test_clustered_cluster()
% One test with clustered cluster

close all;

%% general parameters
p = 0.02;
params.md = 4;
data.t_STED = 5e2;
data.t_CONF = 3e3;

params.FWHM_STED = 4; % 1px = 15nm
params.FWHM_CONF = 12;

% PSFs (the same for all examples here)
h_sted = Gaussian([11, 11], params.FWHM_STED);
h_conf = Gaussian([39, 39], params.FWHM_CONF);

%% cluster simulation
o.type = 'random cluster';
o.border = max(floor(size(h_conf)/2)); % border
B = o.border - 7;
o.size = [140, 140];
o.density = 40/prod(o.size);
o.cluster_mean = 6;
o.cluster_std = 2;
o.cluster_xy_std = 1;
o.cluster_min = 2;

% lambda = 0.01;
lambda = 0;

rng_seed = 1;
object = generate_object(o, rng_seed);
% figure; subplot(1,2,1);imagesc(object); axis image;
object(93:end, (1:35)+4) = object(85:end-8, 1:35);
object(1:60, :) = object(9:68, :);
object(61:68, :) = 0;
% subplot(1,2,2); imagesc(object); axis image;
fprintf('object type %s with %d molecules\n', o.type, sum(object(:)));

% blur object and save
g = exp(-4*log(2)*(-10:10).^2/5^2);
g = g' .* g;
o = conv2(object, g, 'same');
o = -o;
o = (o - min(o(:))) / (max(o(:)) - min(o(:)));
o = uint8(o*255);
imwrite(o, gray, 'map_truth.png');

% do not blur and show
fig = figure();
hold on;
o = [];
for i = 1 : size(object, 1)
    for j = 1 : size(object, 2)
       for k = 1 : object(i, j)
           o = [o; [i, j] + rand(1, 2) - 0.5];
       end
    end
end
plot(o(:, 2), o(:, 1), 'k.', 'MarkerSize', 8);
xlim([1, size(object, 1)]);
ylim([1, size(object, 2)]);
pbaspect([1,1,1]);
box on;
xticks([]);
yticks([]);
exportgraphics(fig, 'clustered_cluster_truth.emf');

save('clustered_cluster.mat','object')

% perform forward simulation
t = tic();
% STED
data.Y_STED = full_forward_simulation(object, p, lambda, h_sted, params.md, data.t_STED, [], true);
% confocal
data.Y_CONF = full_forward_simulation(object, p, lambda, h_conf, params.md, data.t_CONF, [], true);
toc(t);

% store as images
cm = hot(255);
for i = 1:4
    img = data.Y_CONF{i+1};
    img = uint8(img / max(img(:)) * 254);
    img = upsample(img);
    imwrite(img, cm, sprintf('conf_%d.png', i));
        img = data.Y_STED{i+1};
    img = uint8(img / max(img(:)) * 254);
    img = upsample(img);
    imwrite(img, cm, sprintf('sted_%d.png', i));
end

fig = show_results(data.Y_CONF, data.Y_STED, object, B);
exportgraphics(fig, 'clustered_cluster-raw_data.png');
drawnow;

%% estimate molecular map
alpha = 0.1;
% load('test_2021-May-07_10-45-34.mat');
% load('test_2021-Jun-01_08-05-28.mat');
% load('test_2021-Jun-16_12-54-31.mat');
load('test_2021-Jun-17_13-48-34.mat');
% load('test_2021-Jun-18_08-46-39.mat');
% Psi = []; % need to compute a new test before
params.lambda = lambda;
params.a = 2; % matched by hand
params.b = 0.04;
params.segmentation.threshold = 0.05;
params.segmentation.smooth_width = params.FWHM_STED;

t = tic();
[segments, n, p, conf, aux] = estimate_molecular_map(data,params,Psi,alpha);
% "Making the segmentation valid......done!" - can take very long
toc(t);

% plot S_CONF
fig = figure();
fig.Position = [100, 100, 1400, 400];
for i = 1 : 4
    subplot(1,4,i);
    img = aux.S_CONF{i};
    img = img(1+B:end-B,1+B:end-B);
    imagesc(img);
    axis image;
    colormap(hot);
    colorbar();
    title(sprintf('S_%d', i));
end
exportgraphics(fig, 'clustered_cluster-estimated_s.png');

% show segmentation
fig = figure();
subplot(1,2,1);
img = data.Y_STED{2};
img = img(1+B:end-B,1+B:end-B);
imagesc(img);
axis image;
colormap(hot);
title('Y_1^S');
subplot(1,2,2);
s = segments(1+B:end-B,1+B:end-B);
imagesc(s);
axis image;
% hold on;
colormap(jet);
% imagesc(h_sted*numel(conf));
title(sprintf('%d segments', max(segments(:))));
exportgraphics(fig, 'clustered_cluster-segmentation.png');

% fig = figure();
fig = figure();
fig.Position = [100, 100, 1200, 500];
ax = subplot(1,2,1);
ob = object(1+B:end-B,1+B:end-B);
imagesc(ob);
axis image;
title(sprintf('truth n=%d (%d in segments)', sum(object(:)), sum(object(segments>0))));
colormap(ax, jet);
ax = subplot(1,2,2);
imagesc(s);
axis image;
colormap(ax, parula);
hold on;
[xi, yi] = ndgrid(1:size(s, 1), 1:size(s, 2));
hit = zeros(numel(conf), 1);
for i = 1 : numel(conf)
    seg = s == i;
    % center of mass
    xc = sum(xi(:) .* seg(:)) ./ sum(seg(:));
    yc = sum(yi(:) .* seg(:)) ./ sum(seg(:));
    no = sum(ob(seg));
    if no >= conf{i}(1) && no <= conf{i}(2)
        hit(i) = 1;
        text(yc, xc, sprintf('%d (%d,%d)', no, conf{i}), 'HorizontalAlignment', 'center', 'Color', 'w');
    else
        text(yc, xc, sprintf('%d (%d,%d)', no, conf{i}), 'HorizontalAlignment', 'center', 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 12);
    end
%     text(yc, xc, sprintf('%d (%d,%d)', no, conf{i}), 'HorizontalAlignment', 'center', 'Color', 'k');
end
title(sprintf('%d segments, confidence hits=%.2f', numel(conf), mean(hit)));
exportgraphics(fig, 'clustered_cluster-molecular_map.png');

% show density as gray scale plot
Nseg = double(max(segments(:)));
density = zeros(Nseg, 1);
for i = 1 : Nseg
    seg = segments == i;
    density(i) = sum(object(seg)) / sum(seg(:));
end
density_max = max(density);

fig = figure();
hold on;
seg = segments;
cm = gray(1000);
cm = cm(end:-1:1, :);
for i = 1 : Nseg
    B = bwboundaries(seg == i, 'noholes');
    B = B{1};
    c = cm(round(density(i)/density_max*999), :);
    patch(B(:, 2), -B(:, 1), c, 'EdgeColor', 'k')
end
xlim([1, size(seg, 2)]);
ylim([-size(seg, 1), -1]);
pbaspect([1,1,1]);
box on;
xticks([]);
yticks([]);
exportgraphics(fig, 'clustered_cluster_estimated_density.emf');


end

function y = upsample(x)

y = zeros(2*size(x));
y(1:2:end, 1:2:end) = x;
y(2:2:end, 1:2:end) = x;
y(1:2:end, 2:2:end) = x;
y(2:2:end, 2:2:end) = x;

end

function fig = show_results(y_conf, y_sted, object, B)
%

fig = figure();
fig.Position = [100, 100, 1400, 500];
for i = 1 : 4
    img = y_conf{i+1};
    img = img(1+B:end-B,1+B:end-B);
    ax = subplot(2, 4, i);
    imagesc(img);
    axis image;
    colormap(ax, hot);
    colorbar();
    caxis([0, max(img(:))]);
    title(sprintf('Y_%d^C \\Sigma ph. = %d', i, sum(img(:))));
end
for i = 1 : 3
    img = y_sted{i+1};
    img = img(1+B:end-B,1+B:end-B);
    ax = subplot(2, 4, 4+i);
    imagesc(img);
    axis image;
    title(sprintf('Y_%d^S \\Sigma ph. = %d', i, sum(img(:))));
    colormap(ax, hot);
    caxis([0, max(img(:))]);
    colorbar();
end

% object
ax = subplot(2,4,8);
o = object(1+B:end-B,1+B:end-B);
imagesc(o);
axis image;
colormap(ax, jet);
colorbar();
title(sprintf('truth %d mol.', sum(object(:))));

end