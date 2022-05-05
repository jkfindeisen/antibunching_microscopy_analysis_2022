function background = detect_background(image, feature_size)
% Given an 2D image and an idea about the interesting feature size
% detects the background.
%
% The estimated background is a weighted smoothing of the original data,
% where weights are smaller the larger the discrepance between the current
% background estimation at a pixel and the data value. Iterative procedure.
%
% That means:
%   - peaks in the data are nicely preserved
%   - homogeneous areas in the data will mistakenly identified as background
%
% Jan Keller-Findeisen (jan.keller@mpibpc.mpg.de)

assert(nargin == 2, 'Not enough arguments');
image = double(image);

% parameter
N = 3; % the number of iterations

% compute Gaussian convolution kernel with certain FWHM
B = 2 * ceil(feature_size);
g = -B:B;
[x, y] = ndgrid(g, g);
kernel = exp(-4*log(2)*(x.^2+y.^2)/feature_size^2);
kernel = kernel / sum(kernel(:));

% initial estimate of the background is the mean of the image
% background = zeros(size(image)) + mean(image(:));
% background = conv2(image, kernel, 'same');
% background = zeros(size(image)) + median(image(:));
background = 0.5 * conv2(image, kernel, 'same');

% perform a certain number of iterations
for i = 1 : N
    
    % ensure background is larger than zero
    background = max(background, 1e-6);
    
    % compute weights
    w = exp(-(image - background).^2 ./ background.^2);
    
    % background weighted smoothing of the image
    background = conv2(w .* image, kernel, 'same') ./ max(conv2(w, kernel, 'same'), 1e-6);
end

end