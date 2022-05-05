function plot_cell_detail()
% analyses the Origami data, shows a small cutout of it

% Confocal data:
load('data/02_T20160127T223150_rep001_exc00uw_sted00mw_cntr0_0_0nm_pix20_20nm_lpc0mw_780nm_msr(0)_mgate_bchimages.mat');
% This data does not contain the zero-count data, so we reconstruct that
% first
Y_CONF = cell(5,1);
Y_CONF(2:5) = img;
t_CONF = parameters.dw;
Y_CONF{1} = t_CONF*ones(size(img{1}))-cellsum(img);

% STED data
load('data/02_T20160127T223454_rep002_exc00uw_sted00mw_cntr0_0_0nm_pix20_20nm_lpc850mw_780nm_msr(0)_mgate_bchimages.mat');
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

% Parameters for test:
% n =150 + 2*13 = 176;
% a = 2;
% b = 0.0530

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
% Test
load('miscat/test_origami.mat'); %Could be modified, but should also work

% Call main routine
[segments,n,p,conf] = estimate_molecular_map(data,params,Psi,alpha);

% Plot

plot_molecular_map(segments,p,conf,alpha);

end