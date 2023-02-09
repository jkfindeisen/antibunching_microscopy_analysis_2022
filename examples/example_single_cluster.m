function example_single_cluster()
% analyses synthetic data, simulation of a single cluster

close all;
fprintf('Simulate and analyze single cluster data (takes < 1 minute).\n');

%% simulation parameters
lambda = 0.01;
n_true = 10;
p = 0.02;
no_channels = 4;
no_pulses = 1e3;

%% Generate synthetic data ( using the full_forward_simulation_cluster)

% Confocal
h  = gaussian([128,128],23);
Y_CONF = full_forward_simulation_cluster(n_true,p,lambda,h,no_channels,no_pulses);
% STED
h  = gaussian([128,128],7);
Y_STED = full_forward_simulation_cluster(n_true,p,lambda,h,no_channels,no_pulses);

data.Y_CONF = Y_CONF;
data.Y_STED = Y_STED;
data.t_STED = no_pulses;
data.t_CONF = no_pulses;

%% analysis parameters

% Computed by a fitting with h from above
params.a = 2;
params.b = 0.0673;
% Other parameters
params.pixelsize_in_nm = 10;
params.FWHM_STED = 7;
params.FWHM_CONF = 23;
params.md = 4;
% Remaining parameters are used with default values

% Level
alpha = 0.1;
% load test (comment out and use empty Psi for creating a new test)
load('miscat/test_single_cluster.mat');

%% estimate molecular map
[segments,n_est,p,conf,aux] = estimate_molecular_map(data,params,Psi,alpha);

%% plot molecular map
plot_molecular_map(segments,p,conf,alpha,n_true);

end