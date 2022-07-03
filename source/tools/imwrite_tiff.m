function imwrite_tiff(data, name, type)
% very simple helper function using the fire colormap

assert(nargin >= 2, 'Not enough arguments!');

if nargin == 2
    type = 'gray';
end

switch type
    case 'gray'
        % as gray scale for ImageJ
        if max(data(:)) > 2^16
            data = data / max(data(:)) * 2^16;
        end
        imwrite(uint16(data), name);
    case 'hot'
        % as hot
        N = 256;
        data = double(data);
        data = data / max(data(:)) * N;
        imwrite(data, hot(N), name);
    otherwise
        error('Unknown type!');
end

end