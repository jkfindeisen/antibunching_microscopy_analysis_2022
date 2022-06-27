function [n,p,conf,iterations] = estimate_molecular_number_single_segment(data,params,alpha, segments)
% Similar to estimate_molecular_map with the limitation that it only
% regards the confocal data of a simulation and does no segmentation but
% assumes that there is only one segment spanning the whole given area.
% Also there is no background detection.
%
% Used for testing the algorithm and running without segmentation.

%% check that all parameters are there
if ~isfield(params,'md')
    fprintf('Number of channels not given in params, taking default value!\n');
    params.md = numel(data.img)-1;
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
    params.maxit_counting = 20;
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
if nargin < 4 || isempty(segments)
    segments = ones(size(data.img{1}), 'uint32');
end

%% Preprocess data
[S,Cov] = preprocess_data(data.img,data.t,data.lambda,params);

%% Estimate number and brightness on each segment
[n,p,conf,iterations] = estimate_np(S,Cov,segments,alpha/2,params); % Bonferroni adjustment

end