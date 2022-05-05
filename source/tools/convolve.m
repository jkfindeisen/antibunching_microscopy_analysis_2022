function res = convolve(f,otf)
% performs periodic convolution with already FFTd kernel

res = ifft2(fft2(f) .* fftshift(otf));

end