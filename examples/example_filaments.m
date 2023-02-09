function example_filaments()
% One test run with segmentation on simulated filaments.

fprintf('statistical analysis on simulated filaments (takes ~10 minutes)\n');

%% general parameters
p = 0.02;
params.md = 4; 
data.t_STED = 5e2;
data.t_CONF = 2e3;

params.FWHM_STED = 7;
params.FWHM_CONF = 23;

% PSFs (the same for all examples here)
h_sted = gaussian([31, 31], params.FWHM_STED);
h_conf = gaussian([63, 63], params.FWHM_CONF);

%% filament simulation
o.type = 'filaments';
o.border = max(floor(size(h_conf)/2));
o.size = [200, 200];
o.length = 3000;
o.alpha = 0.95;
o.filaments_sigma = 0.5;
o.number_molecules = 10000;

% lambda = 0.01;
lambda = 0;

rng_seed = 4;
[object, filaments] = generate_object(o, rng_seed);
object = round(object / 13);
fprintf('object type %s with %d molecules\n', o.type, sum(object(:)));

% save('filament.mat','object')

% perform forward simulation
t = tic();
% STED
data.Y_STED = full_forward_simulation(object, p, lambda, h_sted, params.md, data.t_STED, [], true);
% confocal
data.Y_CONF = full_forward_simulation(object, p, lambda, h_conf, params.md, data.t_CONF, [], true);
toc(t);

show_results(data.Y_CONF, data.Y_STED, object);
drawnow;

%% estimate molecular map
alpha = 0.1;
load('miscat/test_filaments.mat'); %Psi = []; % % need to compute a new test before
params.lambda = lambda;
params.a = 2; % matched by hand
params.b = 0.04;
t = tic();
[segments, n, p, conf, aux] = estimate_molecular_map(data,params,Psi,alpha);
% "Making the segmentation valid......done!" - takes very long
toc(t);

% Generating the segmentation...
%  --> 208529 significant boxes found
% Discarding all boxes which contain smaller ones......done!
%  --> 1995 boxes left
% Make valid --> 34 boxes left

fig = figure();
subplot(1,2,1);
imagesc(data.Y_STED{2});
axis image;
title('STED 1-ph');
subplot(1,2,2);
imagesc(segments);
axis image;
hold on;
imagesc(h_sted*numel(conf));
title('segmentation');


% fig = figure();
fig = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1);
imagesc(object);
axis image;
title(sprintf('truth n=%d (%d in segments)', sum(object(:)), sum(object(segments>0))));
subplot(1,2,2);
imagesc(segments);
axis image;
hold on;
[xi, yi] = ndgrid(1:size(segments, 1), 1:size(segments, 2));
hit = zeros(numel(conf), 1);
for i = 1 : numel(conf)
    seg = segments == i;
    % center of mass
    xc = sum(xi(:) .* seg(:)) ./ sum(seg(:));
    yc = sum(yi(:) .* seg(:)) ./ sum(seg(:));
    no = sum(object(seg));
    text(yc, xc, sprintf('%d (%d,%d)', no, conf{i}), 'HorizontalAlignment', 'center');
    if no >= conf{i}(1) && no <= conf{i}(2)
        hit(i) = 1;
    end
end
title(sprintf('%d segments, confidence hits=%.2f', numel(conf), mean(hit)));


% write to files
imwrite_tiff(flip(filaments, 1), 'fil-example_filaments.tif');
imwrite_tiff(flip(data.Y_STED{2}, 1), 'fil-example_sted1.tif', 'hot');
imwrite_tiff(flip(data.Y_CONF{2}, 1), 'fil-example_conf1.tif', 'hot');
imwrite_tiff(flip(data.Y_CONF{3}, 1), 'fil-example_conf2.tif', 'hot');
f = figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(flip(segments,1));
axis image;
hold on;
[xi, yi] = ndgrid(1:size(segments, 1), 1:size(segments, 2));
hit = zeros(numel(conf), 1);
for i = 1 : numel(conf)
    seg = segments == i;
    % center of mass
    xc = sum(xi(:) .* seg(:)) ./ sum(seg(:));
    yc = sum(yi(:) .* seg(:)) ./ sum(seg(:));
    no = sum(object(seg));
    text(yc, 200-xc, sprintf('%d (%d,%d)', no, conf{i}), 'HorizontalAlignment', 'center');
    if no >= conf{i}(1) && no <= conf{i}(2)
        hit(i) = 1;
    end
end
xticks([]);
yticks([]);
colormap(jet);
saveas(f, 'fil-example_counting.tif');

end

function show_results(y_conf, y_sted, object)
%

fig = figure();
for i = 1 : 4
    img = y_conf{i+1};
    subplot(2, 4, i);
    imagesc(img);
    axis image;
    title(sprintf('conf. %d-events: %d', i, sum(img(:))));
end
for i = 1 : 3
    img = y_sted{i+1};
    subplot(2, 4, 4+i);
    imagesc(img);
    axis image;
    title(sprintf('STED %d-events: %d', i, sum(img(:))));
end

subplot(2,4,8);
imagesc(object);
axis image;
title('true object');

end