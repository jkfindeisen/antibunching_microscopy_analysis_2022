function Th = mean_values_fixed_scale(Y,Phi,n,h)
% Computes the single-scale statistic
% Th (t) = sum_{i} Y_{i} Phi_h(p_{i}-t)
% for all t in {1,...,n}^d
%
% The argument Phi has already been rotated by 180 degrees, continued
% periodically and fftd. Consequently, the above sum is a convolution, and
% Th = ifftn(Y.*Phi);

assert(nargin == 4, 'Not enough arguments');

tmp = ifftn(Y.*Phi);

m = ceil(n/2);
Th = tmp(m:end-m-h(1), m:end-m-h(2));

end