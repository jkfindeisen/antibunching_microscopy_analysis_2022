function example_experimental_data_origami()
% Analyses the experimental data for DNA origami

close all;
fprintf('Will analyze and show experiment DNA origami data. (Takes 1-2 minutes)\n');

%% load confocal data
load('data/origami_data_confocal_image_data.mat');
% this data does not contain the zero-count data, so we reconstruct that
% first
Y_CONF = cell(5,1);
Y_CONF(2:5) = img;
t_CONF = parameters.dw;
Y_CONF{1} = t_CONF*ones(size(img{1}))-cellsum(img);

%% load STED data
load('data/origami_data_STED_image_data.mat');
Y_STED = cell(5,1);
Y_STED(2:5) = img;
t_STED = parameters.dw;
Y_STED{1} = t_STED*ones(size(img{1}))-cellsum(img);

data.Y_CONF = Y_CONF;
data.Y_STED = Y_STED;
data.t_STED = t_STED;
data.t_CONF = t_CONF;

%% analysis parameters

% Computed by a fitting with the experimental STED psf:
params.a = 2;
params.b = 0.0149;
% Other parameters
params.pixelsize_in_nm = 10;
params.FWHM_STED = parameters.fwhmsted/params.pixelsize_in_nm;
params.FWHM_CONF = parameters.fwhmconf/params.pixelsize_in_nm;
params.md = 4;
% Remaining parameters are used with default values

% Level
alpha = 0.1;

% load test (comment out and use empty Psi for creating a new test)
load('miscat/test_origami.mat');

%% estimate molecular map
[segments,n,p,conf] = estimate_molecular_map(data,params,Psi,alpha);

%% plot molecular map
plot_molecular_map(segments,p,conf,alpha);

end