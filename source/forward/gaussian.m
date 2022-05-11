function k = gaussian(N, fwhmp)
% Returns unnormalized Gaussian PSFs in 2D or 3D
%
%   fwhmp   FWHM of kernel in pixels
%   N       total number of pixels

assert(nargin >=2, 'Not enough arguments!');
assert(all(mod(N, 2) == 1) || all(mod(N, 2) == 0), 'N must all be even or all be odd');

switch numel(N)
    case 2
        if all(mod(N, 2) == 0)
            Np = N/2;
            [X, Y] = ndgrid(-Np(1):Np(1)-1, -Np(2):Np(2)-1);
        else
            Np = (N-1)/2;
            [X, Y] = ndgrid(-Np(1):Np(1), -Np(2):Np(2));
        end
        k = power(2., -(X.^2+Y.^2) / (fwhmp / 2)^2);
    case 3
        if all(mod(N, 2) == 0)
            Np = N/2;
            [X, Y, Z] = ndgrid(-Np(1):Np(1)-1, -Np(2):Np(2)-1, -Np(3):Np(3)-1);
        else
            Np = (N-1)/2;
            [X, Y, Z] = ndgrid(-Np(1):Np(1), -Np(2):Np(2), -Np(3):Np(3));
        end
        k = power(2., -(X.^2+Y.^2+Z.^2) / (fwhmp / 2)^2);
    otherwise
        error('Dimension not implemented!');
end

end