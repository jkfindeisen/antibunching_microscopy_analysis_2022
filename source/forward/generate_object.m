function [object, additional_output] = generate_object(p, rng_seed)
% Generates an random object according to a type and some options which can
% be used in a forward simulation
%
% Inputs
%   p   parameters struct, should at least contain a field called type and
%       then type specific parameters
%
% Output
%   object  an array with integer numbers indicating how many molecules
%           there were at every position

assert(nargin >= 1, 'Not enough parameters!');

% random number generator seed
if nargin >= 2 && ~isempty(rng_seed)
    rng(rng_seed);
else
    rng('shuffle');
end

additional_output = [];

% switch for type
switch p.type
    case 'random cluster'
        % randomly distributed clusters with a certain mean density of
        % clusters and a certain distribution of molecules per cluster
        %
        % p.density = density of clusters in space
        % p.cluster_mean = mean number of molecules per cluster
        % p.cluster_std = std. of number of molecules per cluster
        
        object = zeros(p.size);
        
        [xc, yc] = ind2sub(p.size, find(rand(p.size) < p.density));
        for i = 1:numel(xc)
            xi = xc(i);
            yi = yc(i);
            ni = max(p.cluster_min, round(randn() * p.cluster_std + p.cluster_mean));
            for j = 1 : ni
                xj = xi + round(randn() * p.cluster_xy_std);
                xj = max(1, min(p.size(1), xj));
                yj = yi + round(randn() * p.cluster_xy_std);
                yj = max(1, min(p.size(1), yj));
                object(xj, yj) = object(xj, yj) + 1;
            end
        end
        
        % clear the border
        object = clear_border(object, p.border);
        
    case 'filaments'
        % filament simulation using a correlated random work which is
        % periodically entering the area of interest
        %
        % p.length = determines length of the filaments (in ~ pixel)
        % p.alpha = how flexible are the filaments (1 = straight) typical
        % values between 0.8 and 1
        % p.filaments_sigma = thickness of the filaments
        % p.number_molecules = number of molecules drawn from the filament
        % distribution
        
        % create filaments
        filaments = zeros(p.size);
        points = zeros(p.length, 2);
        points(2, :) = randn(1, 2);
        % start correlated random walk
        for i = 3 : p.length
            d = points(i-1, :) - points(i-2, :);
            d = d / norm(d);
            points(i, :) = points(i-1,:) + p.alpha * d + (1-p.alpha) * randn(1, 2);
        end
        % points = round(points);
        x = mod(points(:, 1) - 1, p.size(1)-1) + 1;
        y = mod(points(:, 2) - 1, p.size(2)-1) + 1;
        %         idx = sub2ind(p.size, x, y);
        %         filaments(idx) = 1;
        
        % weighted adding
        xp = floor(x);
        yp = floor(y);
        dx = x - xp;
        dy = y - yp;
        idx = sub2ind(p.size, xp, yp);
        filaments(idx) = max(filaments(idx), (1-dx)+(1-dy));
        idx = sub2ind(p.size, xp+1, yp);
        filaments(idx) = max(filaments(idx), dx+(1-dy));
        idx = sub2ind(p.size, xp, yp+1);
        filaments(idx) = max(filaments(idx), (1-dx)+dy);
        idx = sub2ind(p.size, xp+1, yp+1);
        filaments(idx) = max(filaments(idx), dx+dy);
        filaments = imgaussfilt(filaments, p.filaments_sigma);
        
        % clear the border of the filaments
        filaments = clear_border(filaments, p.border);
        
        % draw N random positions from the filaments distribution
        filaments = filaments / max(filaments(:));
        object = zeros(p.size);
        length = 0;
        while length < p.number_molecules
            x = randi(p.size(1));
            y = randi(p.size(2));
            if rand() < filaments(x, y)
                object(x, y) = object(x, y) + 1;
                length = length + 1;
            end
        end
        
        additional_output = filaments;
        
    otherwise
        error('Unknown type');
end

end

function [object] = clear_border(object, D)
% sets all the object within a distance to the border to 0

object(1:D, :) = 0;
object(end-D+1:end, :) = 0;
object(:, 1:D) = 0;
object(:, end-D+1:end) = 0;

end