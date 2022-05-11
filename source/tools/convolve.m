function res = convolve(f, otf)
% Performs periodic convolution with already FFT-ed kernel

res = ifftn(fftn(f) .* fftshift(otf));

end