function [segments,n,p,conf,aux] = estimate_molecular_map(data,params,Psi,alpha)
% Main function to estimate the molecular map from data

% (c) Frank Werner, IMS Uni Göttingen, 25.02.2019 or later

% INPUTS: 
% data --> struct containing the fields
%       Y_STED --> antibunching images from STED experiment
%       Y_CONF --> antibunching images from confocal experiment
%       t_STED --> dwell time in STED experiment
%       t_CONF --> dwell time in confocal experiment
% params --> struct containing the fields 
%       FWHM_STED --> full width at half maximum in pixels for STED experiment
%       FWHM_CONF --> full width at half maximum in pixels for confocal experiment
%       a --> shape parameter a for approximate STED psf (precomputed)
%       b --> scale parameter b for approximate STED psf (precomputed)
%       lambda --> OPTIONAL: Estimate of background Poisson level
%       pixelsize_in_nm --> pixelsize in nm
%       md --> number of detectors used for measurements
%       mD --> number of channels used in reconstruction (<= md, DEFAULT: md)
%       mQ --> number of Q's used in approximation (free, DEFAULT: md)
%       mS --> number of S's used in approximation formula for Q (<= ms, >= mp, DEFAULT: md)
%       ms --> number of s's to be determined (>=mS, DEFAULT: md)
%       mp --> number of p's or lambda's used in formula for Q (<=mS, DEFAULT: md)
%       mc --> number of s's used for counting (<= ms, DEFAULT:2)
%       maxit_preprocessing --> Number of maximal Newton iterations for preprocessing (DEFAULT: 10)
%       maxit_counting --> Number of maximal Newton iterations for counting (DEFAULT: 10)
%       tol_preprocessing --> Relative tolerance for Newton's method for preprocessing (DEFAULT: 1e-4)
%       tol_copunting --> Relative tolerance for Newton's method for counting (DEFAULT: 1e-4)
%       no_neighboring_px --> number of neighboring pixels used for counting (DEFAULT: FWHMP_CONF)
% Psi --> filename of test used for multiscale location (or [] if unknown)
% alpha --> final error level (1-alpha is the final confidence level)

% OUTPUTS:
% segments --> regions containing markers with uniform prob >= 1-alpha
% n --> estimated number of markers in the segments
% p --> estimated brightness of markers in the segments
% conf --> confidence interval for number of markers in the segments, coverage >= 1-alpha
% aux --> further information for debugging

%% Set default values
if ~isfield(params,'md')
    fprintf('Number of channels not given in params, taking default value!\n');
    params.md = numel(data.Y_STED)-1;
end
if ~isfield(params,'mD')
    params.mD = params.md;
end
if ~isfield(params,'mQ')
    params.mQ = params.md;
end
if ~isfield(params,'mS')
    params.mS = params.md;
end
if ~isfield(params,'ms')
    params.ms = params.md;
end
if ~isfield(params,'mp')
    params.mp = params.md;
end
if ~isfield(params,'mc')
    params.mc = 2;
end
if ~isfield(params,'maxit_preprocessing')
    params.maxit_preprocessing = 10;
end
if ~isfield(params,'maxit_counting')
    params.maxit_counting = 10;
end
if ~isfield(params,'tol_preprocessing')
    params.tol_preprocessing = 1e-4;
end
if ~isfield(params,'tol_counting')
    params.tol_counting = 1e-4;
end
if ~isfield(params,'no_neighboring_px')
    params.no_neighboring_px = params.FWHM_CONF;
end
if ~isfield(params, 'segmentation')
    params.segmentation.smooth_width = 4;
    params.segmentation.threshold = 0.15;
end

%% Check data consistency
Y_STED = data.Y_STED;
Y_CONF = data.Y_CONF;
t_CONF = data.t_CONF;
t_STED = data.t_STED;

fprintf('Preparations...');
if size(Y_STED,1) ~= size(Y_CONF,1) || size(Y_STED,2) ~=1 || size(Y_CONF,2) ~=1 || size(Y_STED{1},1)~= size(Y_CONF{1},1) || size(Y_STED{1},2) ~= size(Y_CONF{1},2)
    error('Inconsistent data sizes!');
end

% Compute necessary padding
pad_STED = ceil(4*params.FWHM_STED/2/sqrt(2*log(2)));
% pad_STED = ceil(4*params.FWHM_CONF/2/sqrt(2*log(2)));

% Check if padded image is suitable for test in Psi or generate new test
% TODO is checking the size sufficient to decide if the test in Psi is suitable?
if isempty(Psi) || Psi.n ~= size(Y_STED{1},1) + 2*pad_STED || Psi.n ~= size(Y_STED{1},2) + 2*pad_STED
    % In this case, no test has been specified, or the one specified does
    % not apply. We have to generate a new one.
    hset = generate_set_of_scales(ceil(params.FWHM_STED/2),4*params.FWHM_STED,2);
%     hset = generate_set_of_scales(ceil(params.FWHM_STED),3*params.FWHM_STED,2);
    [k,otf] = generate_psf(size(Y_STED{1},1) + 2*pad_STED,params.b,params.a);
    kernel_params.k = k;
    kernel_params.a = params.a;
    kernel_params.b = params.b;
    Psi = setup_test(size(Y_STED{1},1) + 2*pad_STED,kernel_params,otf,[2*params.a,2*params.a],hset,['test_',datestr(clock,'yyyy-mmm-dd_HH-MM-SS'),'.mat']);
    Psi = simulate_gaussian_approximation(Psi.filename,1e3);    
end
fprintf('...done!\n');

%% Detect background?
if ~isfield(params,'lambda')
    fprintf('Detecting background...');
    lambda_STED = detect_background(Y_STED{2}, 2*params.FWHM_STED);
    lambda_CONF = detect_background(Y_CONF{2}, 2*params.FWHM_CONF);
    fprintf('...done!\n');
else
    lambda_STED = t_STED*params.lambda*ones(size(Y_STED{1}));
    lambda_CONF = t_CONF*params.lambda*ones(size(Y_CONF{1}));
end

%% Preprocess data
fprintf('Preprocessing data...');
S_STED = preprocess_data(Y_STED,t_STED,lambda_STED,params);
[S_CONF,Cov] = preprocess_data(Y_CONF,t_CONF,lambda_CONF,params);
fprintf('...done!\n');

%% Locate with multiscale test
fprintf('Generating the segmentation...\n');
dist_params.name = 'Binomial';
dist_params.t = t_STED;
dist_params.bg = 0;
Bnp = multiscale_test(padarray(max(S_STED{1},0),[pad_STED pad_STED]),Psi,dist_params,alpha/2); % Bonferroni adjustment

%% Merge segmentation with estimated locations
segments = boxes_to_segments(max(S_STED{1},0),Bnp,Psi.hset,pad_STED, params.segmentation);
fprintf('...done!\n');

%% Estimate number and brightness on each segment
fprintf('Counting...');
[n,p,conf] = estimate_np(S_CONF,Cov,segments,alpha/2,params); % Bonferroni adjustment
fprintf('...done!\n');

if nargout == 5
    aux.lambda_STED = lambda_STED;
    aux.lambda_CONF = lambda_CONF;
    aux.S_CONF = S_CONF;
end
end