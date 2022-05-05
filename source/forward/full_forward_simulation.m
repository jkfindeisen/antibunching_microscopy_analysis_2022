function img = full_forward_simulation(object, p, lambda, psf, d, t, rng_seed, show_progressbar, molecules_close_to_border_okay)
% Performs the full forward simulation for many molecules distributed on a
% grid
%
% Note: Only 2D implemented so far!
%
% Inputs
%   object  grid with integer numbers indicating the number of molecules
%           at each grid position
%   p       Molecular brightness of a molecule at the center of the PSF,
%           either a grid as big as object or a single value (will then
%           used for all)
%   lambda  Mean background probability (can also be a grid)
%   psf     PSF (can be different grid, but we prefer odd numbers of pixels here)
%   d       Number of detection channels
%   t       Number of repetitions / pulses
%   rng_seed    Random number generator seed (optional)
%   show_progressbar    Shows a progress bar (optional)
%
% Output
%   img     Cell array with d+1 entries for the 0..d photon event images

% safety checks
assert(nargin >= 6, 'Not enough parameters!');
assert(numel(size(object)) == 2, 'object must be a 2D array');
assert(all(object(:) >= 0), 'object must be >= 0');
assert(all(object(:) == round(object(:))), 'object must only contain integer values');
assert(numel(p) == 1 || isequal(size(object), size(p)), 'p must either be a single value or an array of the same size as object');
assert(all(p(:) >= 0), 'p must be >= 0');
assert(numel(lambda) == 1 || isequal(size(object), size(lambda)), 'lambda must either be a single value or an array of the same size as object');
assert(all(lambda(:) >= 0), 'lambda must be >= 0');
assert(numel(size(psf)) == 2, 'psf must be a 2D array');
assert(all(psf(:) >= 0), 'psf cannot be negative');
assert(all(mod(size(psf), 2) == 1), 'odd number of pixels in PSF preferred');
assert(numel(d) == 1 && d > 0 && d == round(d), 'd must be positive integer');
assert(numel(t) == 1 && t > 0 && t == round(t), 'no_pulses must be positive integer');

% check that there is no object within psf_radius from the border
if nargin < 9 || isempty(molecules_close_to_border_okay)
    molecules_close_to_border_okay = false;
end
psf_radius = floor(size(psf) / 2);
l = object(1:psf_radius(1), :);
r = object(end-psf_radius(1)+1:end, :);
u = object(:, 1:psf_radius(2));
b = object(:, end-psf_radius(2)+1:end);
assert(molecules_close_to_border_okay || sum(l(:)) + sum(r(:)) + sum(u(:)) + sum(b(:)) == 0, 'there should not be any molecules close to the border (half the psf radius)');

% random number generator seed
if nargin >= 7 && ~isempty(rng_seed)
    rng(rng_seed);
else
    rng('shuffle');
end

% progress bar
if nargin < 8
    show_progressbar = false;
end

% if p constant, expand to size of object
if numel(p) == 1
    p = zeros(size(object)) + p;
end

% if lambda constant, expand to size of object
if numel(lambda) == 1
    lambda = zeros(size(object)) + lambda;
end

% normalize PSF so that it is 1 at maximum
psf = psf / max(psf(:));

% initialize output variable
dims = size(object);
img = cell(d+1,1);
for i = 1:d+1
    img{i} = zeros(dims);
end

% now perform a scan
if show_progressbar
    bar = progressbar('forward simulation');
end
for x = 1 : dims(1)
    for y = 1 : dims(2)
        % (x, y) is the actual scan position
        
        % collect all the relevant n's and p's
        nps = [];
        % loop over all object positions (i, j)
        for i = max(1, x - psf_radius(1)) : min(dims(1), x + psf_radius(1))
            for j = max(1, y - psf_radius(2)) : min(dims(2), y + psf_radius(2))
                % number of molecules at position (i, j)
                nij = object(i, j);
                % brightness of molecules at position (i, j) times psf value
                pij =  p(i, j) * psf(i - x + psf_radius(1) + 1, j - y + psf_radius(2) + 1);
                if pij > 0 && nij > 0
                    nps = [nps; [nij, pij]];
                end
            end
        end
        
        % TODO if speedup is necessary we could replace some nps with very
        % small p into lambda here
        
        % get number of events for each repetition (and sort them into output)
        no_events = draw_events(nps, lambda(x, y), d, t);
        for i = 1:d+1
            img{i}(x, y) = no_events(i);
        end
        
        if show_progressbar
            bar.set(((x-1) * dims(2) + y) / prod(dims));
        end
        
    end
end
if show_progressbar
    bar.close();
end

end

function [no_events, photons] = draw_events(nps, lambda, d, t)
% For a list of numbers and corresponding brightnesses and a mean
% background value and a number of detectors and pulses, draw a random
% sample of number of active detectors per pulse, do it in parallel for
% no_pulses

% first get number of background photons per pulse
bg_photons = poissrnd(lambda, t, 1);

% now get number of photons from molecules
mol_photons = zeros(t, 1);
for i = 1 : size(nps, 1)
    mol_photons = mol_photons + binornd(nps(i, 1), nps(i, 2), t, 1);
end

% total number of photons is the sum of both
photons = bg_photons + mol_photons;

% distribute on detectors
no_events = distribute_on_detectors(photons, d);

end

function no_events = distribute_on_detectors(photons, d)

% distribute on detectors, we only need to treat those above one photon
mi = photons > 1;
pmi = photons(mi);

if ~isempty(pmi)
    
    % get a list of 1..., 2..., 3... where each index is repeated pmi times
    x = zeros(sum(pmi), 1);
    x(cumsum([1; pmi(1:end-1)])) = 1;
    x = cumsum(x);
    
    % draw on which detector they have fallen
    di = randi(d, length(x), 1);
    
    % determine sum of active detectors
    dd = zeros(d, length(pmi));
    idx = sub2ind(size(dd), di, x);
    dd(idx) = 1;
    dd = sum(dd);
    
    % write back to photons per repetition, becomes active detectors
    photons(mi) = dd;
end

% calculate empiric occurences
no_events = zeros(d+1, 1);
for i = 0 : d
    no_events(i+1) = sum(photons == i);
end

end