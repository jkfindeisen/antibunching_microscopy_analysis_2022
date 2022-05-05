function MAIN_origami()

% Confocal data:
load('data/01_T20130411T231851_rep001_exc_0005.0uW_lpc0_0_pix_10_10_0nm_scan_bchimages.mat');
% This data does not contain the zero-count data, so we reconstruct that
% first
Y_CONF = cell(5,1);
Y_CONF(2:5) = img;
t_CONF = parameters.dw;
Y_CONF{1} = t_CONF*ones(size(img{1}))-cellsum(img);

% STED data
load('data/01_T20130411T232014_rep002_exc_0005.0uW_lpc402_780_pix_10_10_0nm_scan_bchimages.mat');
Y_STED = cell(5,1);
Y_STED(2:5) = img;
t_STED = parameters.dw;
Y_STED{1} = t_STED*ones(size(img{1}))-cellsum(img);

data.Y_CONF = Y_CONF;
data.Y_STED = Y_STED;
data.t_STED = t_STED;
data.t_CONF = t_CONF;

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
% Test
load('miscat/test_origami.mat');

% Call main routine
[segments,n,p,conf] = estimate_molecular_map(data,params,Psi,alpha);

% Plot

plot_molecular_map(segments,p,conf,alpha);

end