function hset = generate_set_of_scales(hmin,hmax,hstep)
% Use scales of size hmin:hstep:hmax pixels and all possible cross-products

assert(nargin == 3, 'Not enough arguments');

hset = [];

if (nargin < 3)
    hstep = 1;
end
if (nargin < 2)
    hmax = 15;
end
if (nargin < 1)
    hmin = 8;
end

for i=hmin:hstep:hmax
    for j=hmin:hstep:hmax
        hset = [hset; i j];
    end
end

hset = num2cell(hset,2);

end