function example_forward_simulations()
% Example code showing how to do forward simulations of complex objects
%
% - Generate object and PSFs
% - Perform forward simulation
% - Show results

close all;
fprintf('Simulating photon incidence images for various objects and imaging modes, can take a while!\n\n');

%% general parameters
lambda = 0.01; 
p = 0.02;
no_channels = 4; 
no_pulses = 1e2;

% PSFs (the same for all examples here)
h_conf = gaussian([41, 41], 16);
h_sted = gaussian([13, 13], 4);

% simulation become considerably faster if we set psfs below 1 / 100 to
% zero (~20% faster)
h_conf(h_conf < 0.01 * max(h_conf(:))) = 0;
h_sted(h_sted < 0.01 * max(h_sted(:))) = 0;


%% cluster simulation
o.type = 'random cluster';
o.border = max(floor(size(h_conf)/2));
o.size = [100, 100];
o.density = 0.01;
o.cluster_min = 1;
o.cluster_mean = 10;
o.cluster_std = 3;
o.cluster_xy_std = 1;

object = generate_object(o);
fprintf('object "%s" with %d molecules\n', o.type, sum(object(:)));

% perform forward simulation
t0 = tic();
% STED
y_sted = full_forward_simulation(object, p, lambda, h_sted, no_channels, no_pulses, [], true);
% confocal
y_conf = full_forward_simulation(object, p, lambda, h_conf, no_channels, no_pulses, [], true);
dt = toc(t0);
fprintf(' confocal + STED simulation took %.2fs\n', dt);

show_results(y_conf, y_sted, object, o.type);

%% the same cluster simulation but with non-homogeneous background

[~, yi] = ndgrid(1:size(object, 1), 1:size(object, 2));
lambda = yi / size(object,2) * 0.2;

% perform forward simulation
t0 = tic();
% STED
y_sted = full_forward_simulation(object, p, lambda, h_sted, no_channels, no_pulses, [], true);
% confocal
y_conf = full_forward_simulation(object, p, lambda, h_conf, no_channels, no_pulses, [], true);
dt = toc(t0);
fprintf(' confocal + STED simulation with inhomogeneous background took %.2fs\n', dt);

show_results(y_conf, y_sted, object, o.type);

%% filament simulation
lambda = 0.01;

o.type = 'filaments';
o.border = max(floor(size(h_conf)/2));
o.size = [100, 100];
o.length = 1500;
o.alpha = 0.95;
o.filaments_sigma = 0.5;
o.number_molecules = 300;

object = generate_object(o);
fprintf('object "%s" with %d molecules\n', o.type, sum(object(:)));

% perform forward simulation
t0 = tic();
% STED
y_sted = full_forward_simulation(object, p, lambda, h_sted, no_channels, no_pulses, [], true);
% confocal
y_conf = full_forward_simulation(object, p, lambda, h_conf, no_channels, no_pulses, [], true);
dt = toc(t0);
fprintf(' confocal + STED simulation took %.2fs\n', dt);

show_results(y_conf, y_sted, object, o.type);

%% now simulation with variable brightness

o.type = 'random cluster';
o.border = max(floor(size(h_conf)/2));
o.size = [100, 100];
o.density = 0.01;
o.cluster_mean = 10;
o.cluster_std = 3;

object = generate_object(o);
fprintf('object "%s" with %d molecules\n', o.type, sum(object(:)));

% background increases from left to right
[~, yi] = ndgrid(1:size(object, 1), 1:size(object, 2));
lambda = yi / size(object,2) * 0.2;

% molecular brightness decreases from left to right
pxy = 2 * p * (1 - yi / max(yi(:)));

% perform forward simulation
t0 = tic();
% STED
y_sted = full_forward_simulation(object, pxy, lambda, h_sted, no_channels, no_pulses, [], true);
% confocal
y_conf = full_forward_simulation(object, pxy, lambda, h_conf, no_channels, no_pulses, [], true);
dt = toc(t0);
fprintf(' confocal + STED simulation took %.2fs\n', dt);

show_results(y_conf, y_sted, object, o.type);

end

function show_results(y_conf, y_sted, object, name)
% display photon incidence images

fig = figure();
fig.Position = [100, 100, 1200, 500];
fig.Name = name;
for i = 1 : 4
    img = y_conf{i+1};
    subplot(2, 4, i);
    imagesc(img);
    axis image;
    title(sprintf('Y_%d^C (\\Sigma=%d)', i, sum(img(:))));
end
for i = 1 : 3
    img = y_sted{i+1};
    subplot(2, 4, 4+i);
    imagesc(img);
    axis image;
    title(sprintf('Y_%d^S (\\Sigma=%d)', i, sum(img(:))));
end

subplot(2,4,8);
imagesc(object);
axis image;
title('true object');

end