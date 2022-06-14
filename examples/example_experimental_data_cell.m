function example_experimental_data_cell()
% analyses the Origami data, shows a small cutout of it

close all;
fprintf('Will analyze and show experimental cell data. May take a while!\n');

%% load confocal data
load('data/cell_data_confocal_image_data.mat');
% This data does not contain the zero-count data, so we reconstruct that
% first
Y_CONF = cell(5,1);
Y_CONF(2:5) = img;
t_CONF = parameters.dw;
Y_CONF{1} = t_CONF*ones(size(img{1}))-cellsum(img);

%% load STED data

load('data/cell_data_STED_image_data.mat');
Y_STED = cell(5,1);
Y_STED(2:5) = img;
t_STED = parameters.dw;
Y_STED{1} = t_STED*ones(size(img{1}))-cellsum(img);

fig = figure();
for i = 1 : 5
    subplot(2,5,i);
    imagesc(Y_STED{i}); axis image;
    subplot(2,5,i+5);
    imagesc(Y_CONF{i}); axis image;
end

window1 = 431:580;
window2 = 331:480;
for i=1:5
    Y_STED{i} = Y_STED{i}(window1,window2);
    Y_CONF{i} = Y_CONF{i}(window1,window2);
end

data.Y_CONF = Y_CONF;
data.Y_STED = Y_STED;
data.t_STED = t_STED;
data.t_CONF = t_CONF;

% Computed by a fitting with the experimental STED psf:
params.a = 2;
params.b = 0.0530;
% Other parameters
params.pixelsize_in_nm = 10;
params.FWHM_STED = 75/params.pixelsize_in_nm; %Not given in the parameters!!
params.FWHM_CONF = 240/params.pixelsize_in_nm; %Taken from parameters
params.md = 4; % Remaining parameters are used with default values
params.preprocessing_inversion_method = 'newton'; %JAN: CHANGE HERE

% Level
alpha = 0.1;

% load test (comment out and use empty Psi for creating a new test)
load('miscat/test_cell.mat'); % Could be modified, but should also work

%% estimate molecular map
[segments,n,p,conf] = estimate_molecular_map(data,params,Psi,alpha);

%% plot molecular map
plot_molecular_map(segments,p,conf,alpha);

end