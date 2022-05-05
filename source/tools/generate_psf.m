function [k, otf] = generate_psf(n, b, a)

assert(nargin == 3, 'Not enough arguments');

% we define our kernel in fourier space
[xi1,xi2] = meshgrid(-n:2:n-1,-n:2:n-1);
otf = (1./(1 + b^2 *( xi1.^2 + xi2.^2))).^a;
k = ifftshift(ifft2(fftshift(otf)));

end