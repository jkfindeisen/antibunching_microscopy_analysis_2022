function example_cell_dummy2()
% tests the cell example, only the analysis

fprintf('Takes 10-20 minutes!\n')

%% load cell example
pattern = load('data/cell_dummy_object.mat');
pattern = pattern.pattern;
% cut off some border
B = 10;
pattern = pattern(1+B:end-B,1+B:end-B);
L = bwlabel(pattern > 0);
A = L == 1;
B = L == 2;
cell = pattern;
cell(A) = max(0, cell(A) - 0.8);% need to subtract a bit of volume
cell(B) = 0;
cell = cell / max(cell(:));

% density in cell
fig = figure();
imagesc(cell);
pbaspect([1,1,1]);
colormap(flip(gray));

%% draw random position from the cell example
rng(1);
number_molecules = 1000;
dims = size(cell);
object = zeros(dims);
n = 0;
while n < number_molecules
    x = randi(dims(1));
    y = randi(dims(2));
    if rand() < cell(x, y)
        object(x, y) = object(x, y) + 1;
        n = n + 1;
    end
end

% number of molecules in object
fig = figure();
imagesc(object);
pbaspect([1,1,1]);
colormap(flip(gray));

%% general parameters
p = 0.02;
params.md = 4; 
data.t_STED = 1e3;
data.t_CONF = 2e4;

params.FWHM_STED = 8;
params.FWHM_CONF = 32;

% PSFs (the same for all examples here)
h_sted = gaussian([31, 31], params.FWHM_STED);
h_conf = gaussian([63, 63], params.FWHM_CONF);

lambda = 0; %.01;

%% perform forward simulation
t = tic();
% STED
data.Y_STED = full_forward_simulation(object, p, lambda, h_sted, params.md, data.t_STED, [], true);
% confocal
data.Y_CONF = full_forward_simulation(object, p, lambda, h_conf, params.md, data.t_CONF, [], true);
toc(t);

%% estimate molecular map
alpha = 0.1;
load('miscat/test_cell_dummy.mat'); %Psi = []; % % need to compute a new test before
params.lambda = lambda;
params.a = 2; % matched by hand
params.b = 0.031;
t = tic();
[segments, n, p, conf] = estimate_molecular_map(data,params,Psi,alpha);
% "Making the segmentation valid......done!" - takes very long
toc(t);

%% check if in confidence interval
Nseg = double(max(segments(:)));
for i = 1 : Nseg
    seg = segments == i;
    n_true = sum(object(seg));
    ci = conf{i};
    fprintf('seg %d \ttrue n %d \test. n %.1f \test. p %.2f \tconf [%.1f, %.1f] \tin conf %d\n', i, n_true, n{i}, p{i}, ci, n_true >= ci(1) & n_true < ci(2));
end

end