function example_forward_simulation_single_cluster()
% tests the forward simulations

% Set parameters
lambda = 0.01;
n = 10;
p = 0.02;
no_channels = 4;
no_pulses = 1e3;

%% Generate synthetic data ( using the full_forward_simulation_cluster)

t = tic();
% Confocal
h  = gaussian([128,128],23);
Y_CONF = full_forward_simulation_cluster(n,p,lambda,h,no_channels,no_pulses);
% STED
h  = gaussian([128,128],7);
Y_STED = full_forward_simulation_cluster(n,p,lambda,h,no_channels,no_pulses);
toc(t);

%% Generate synthetic data (using the full_forward_simulation)

object = zeros(127,127);
object(63,63) = n;

t = tic();
% Confocal
h = gaussian([65, 65], 23);
Y_CONF2 = full_forward_simulation(object, p, lambda, h, no_channels,no_pulses);
% STED
h  = gaussian([33,33], 7);
Y_STED2 = full_forward_simulation(object, p, lambda, h, no_channels,no_pulses);
toc(t);

%% Compare results

% forward_cluster
fig = figure();
for i = 1 : 4
    img = Y_CONF{i+1};
    subplot(2, 4, i);
    imagesc(img);
    axis image;
    title(sprintf('Y_%d^C (\\Sigma=%d)', i, sum(img(:))));
end
for i = 1 : 4
    img = Y_STED{i+1};
    subplot(2, 4, 4+i);
    imagesc(img);
    axis image;
    title(sprintf('Y_%d^S (\\Sigma=%d)', i, sum(img(:))));
end

% forward (general)
fig = figure();
for i = 1 : 4
    img = Y_CONF2{i+1};
    subplot(2, 4, i);
    imagesc(img);
    axis image;
    title(sprintf('Y_%d^C (\\Sigma=%d)', i, sum(img(:))));
end
for i = 1 : 4
    img = Y_STED2{i+1};
    subplot(2, 4, 4+i);
    imagesc(img);
    axis image;
    title(sprintf('Y_%d^S (\\Sigma=%d)', i, sum(img(:))));
end

end