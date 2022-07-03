function example_cell_dummy()
% Shows another simulation with a cell dummy (quite crowded but not enough
% fine structure).

fprintf('Takes very long!\n');

% load cell example
pattern = load('data/cell_dummy_object.mat');
pattern = pattern.pattern;
% cut off some border
B = 10;
pattern = pattern(1+B:end-B,1+B:end-B);
cell = max(0, pattern - 0.6); % need to subtract a bit of volume
cell = cell / max(cell(:));

% save as transparent image
% A = uint8(cell/max(cell(:))*255);
% cm = gray(255);
% imwrite(A, cm(end:-1:1,:), 'cell.png'); % transparentcolor doesn't work with png

% draw random position from the cell example
rng('shuffle');
number_molecules = 2000;
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

figure();
subplot(1,2,1);
imagesc(cell);
axis image;
subplot(1,2,2);
imagesc(object);
axis image;

%% general parameters
p = 0.02;
params.md = 4;
data.t_STED = 1e3;
data.t_CONF = 3e2;

params.FWHM_STED = 7;
params.FWHM_CONF = 23;

% PSFs (the same for all examples here)
h_sted = gaussian([31, 31], params.FWHM_STED);
h_conf = gaussian([63, 63], params.FWHM_CONF);

lambda = 0; %.01;

% perform forward simulation
t = tic();
% STED
data.Y_STED = full_forward_simulation(object, p, lambda, h_sted, params.md, data.t_STED, [], true);
% confocal
data.Y_CONF = full_forward_simulation(object, p, lambda, h_conf, params.md, data.t_CONF, [], true);
toc(t);

show_results(data.Y_CONF, data.Y_STED, object, cell);

% store as images
% cm = hot(255);
% for i = 1:2
%     img = data.Y_CONF{i+1};
%     img = uint8(img / max(img(:)) * 255);
%     imwrite(img, cm, sprintf('conf_%d.png', i));
%         img = data.Y_STED{i+1};
%     img = uint8(img / max(img(:)) * 255);
%     imwrite(img, cm, sprintf('sted_%d.png', i));
% end

%% estimate molecular map
alpha = 0.1;
load('miscat/test_cell_dummy.mat'); %Psi = []; % % need to compute a new test before
% Psi = [];
params.lambda = lambda;
params.a = 2; %matched by hand
params.b = 0.031;
t = tic();
[segments, n, p, conf, aux] = estimate_molecular_map(data,params,Psi,alpha);
% "Making the segmentation valid......done!" - takes very long
toc(t);

fig = figure();
hold on;
axis off;
box on;
seg = d.segmentation;
Nseg = double(max(seg(:)));
cm = jet(Nseg);
for i = 1 : Nseg
    B = bwboundaries(seg == i, 'noholes');
    B = B{1};
    patch(B(:, 2), -B(:, 1), cm(i, :), 'EdgeColor', 'k')
end
xlim([1, size(seg, 2)]);
ylim([-size(seg, 1), -1]);
pbaspect([1,1,1]);
exportgraphics(fig, 'cell_segmentation.emf');

fig = figure();
hold on;
axis off;
box on;
N = size(d.boxes, 1);
% sort boxes by area
a = (diff(d.boxes(:, [1,3]), 1, 2) + 1) .* (diff(d.boxes(:, [2,4]), 1, 2) + 1);
[a, ix] = sort(a, 'desc');
boxes = d.boxes(ix, :);
a = sqrt(a); % we use the sqrt of the area as color indicator
R = [min(a), max(a)];
% cm = jet(1000);
cm = parula(1000);
% cm = img_colormap_dawn(1000);
% cm = bone(1000);
% cm = winter(1000);
for i = 1 : N
    b = boxes(i, :);
    ci = floor((a(i) - R(2))/(R(1)-R(2))*999)+1;
    %     c = floor(cm(ci, :)*255);
    c = cm(ci, :);
    patch(b([2, 2, 4, 4, 2]), -b([1,3,3,1,1]), c, 'EdgeColor', 'k');
end
xlim([1, size(seg, 2)]);
ylim([-size(seg, 1), 1]);
pbaspect([1,1,1]);
exportgraphics(fig, 'cell_boxes.emf');


% A = double(d.watershed);
% A = uint8(A / max(A(:)) * 255);
% cm = jet(255);
% cm(1, :) = [1,1,1];
% imwrite(A, cm, 'cell_watershed.png');
%
% A = double(d.segmentation);
% A = uint8(A / max(A(:)) * 255);
% cm = jet(255);
% cm(1, :) = [1,1,1];
% imwrite(A, cm, 'cell_segmentation.png');
%

N = size(d.boxes, 1);
A = zeros([size(d.segmentation), 3], 'uint8') + 255; % white at the beginning
% sort boxes by area
a = (diff(d.boxes(:, [1,3]), 1, 2) + 1) .* (diff(d.boxes(:, [2,4]), 1, 2) + 1);
[a, ix] = sort(a, 'desc');
boxes = d.boxes(ix, :);
a = sqrt(a); % we use the sqrt of the area as color indicator
R = [min(a), max(a)];
% cm = jet(1000);
cm = parula(1000);
% cm = img_colormap_dawn(1000);
% cm = bone(1000);
% cm = winter(1000);
for i = 1 : N
    b = boxes(i, :);
    ci = floor((a(i) - R(2))/(R(1)-R(2))*999)+1;
    c = floor(cm(ci, :)*255);
    A(b(1):b(3),b(2):b(4),1) = c(1);
    A(b(1):b(3),b(2):b(4),2) = c(2);
    A(b(1):b(3),b(2):b(4),3) = c(3);
end
imshow(A);
% colormap(parula);
% colormap(img_colormap_dawn);
colormap(bone);
colormap(winter);
colorbar();
% exportgraphics(gcf, 'cell_boxes_parula.png');
% exportgraphics(gcf, 'cell_boxes_dawn.png');
% exportgraphics(gcf, 'cell_boxes_bone.png');
% exportgraphics(gcf, 'cell_boxes_winter.png');
imwrite(A, 'cell_boxes.png');

% molecular map (color coded for density)
N = max(d.segmentation(:));
out = zeros(size(d.segmentation));
for i = 1 : N
    seg = d.segmentation == i;
    no = sum(object(seg));
    A = sum(seg(:));
    out(seg) = no / A;
end
out = uint8(out / max(out(:)) * 255);
out(out == 0) = 255;
imwrite(out, gray, 'map_truth.png');

% molecular map
fig = figure();
imagesc(segments);
axis image;
cm = jet(256);
cm(1, :) = [1,1,1];
colormap(cm);
xticks([]);
yticks([]);
hold on;
[xi, yi] = ndgrid(1:size(segments, 1), 1:size(segments, 2));
hit = zeros(numel(conf), 1);
for i = 1 : numel(conf)
    seg = segments == i;
    % center of mass
    xc = sum(xi(:) .* seg(:)) ./ sum(seg(:));
    yc = sum(yi(:) .* seg(:)) ./ sum(seg(:));
    no = sum(object(seg));
    %     text(yc, xc, sprintf('%d (%d,%d)', no, conf{i}), 'HorizontalAlignment', 'center');
    text(yc, xc, sprintf('%d', no), 'HorizontalAlignment', 'center', 'FontSize', 14);
    if no >= conf{i}(1) && no <= conf{i}(2)
        hit(i) = 1;
    end
end
exportgraphics(fig, 'cell_molecular_map.png');

fig = figure();
subplot(1,2,1);
imagesc(data.Y_STED{2});
axis image;
title('STED 1-ph');
subplot(1,2,2);
imagesc(segments);
axis image;
hold on;
imagesc(h_sted*numel(conf));
title('segmentation');

% fig = figure();
fig = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1);
imagesc(object);
axis image;
title(sprintf('truth n=%d (%d in segments)', sum(object(:)), sum(object(segments>0))));
subplot(1,2,2);
imagesc(segments);
axis image;
hold on;
[xi, yi] = ndgrid(1:size(segments, 1), 1:size(segments, 2));
hit = zeros(numel(conf), 1);
for i = 1 : numel(conf)
    seg = segments == i;
    % center of mass
    xc = sum(xi(:) .* seg(:)) ./ sum(seg(:));
    yc = sum(yi(:) .* seg(:)) ./ sum(seg(:));
    no = sum(object(seg));
    text(yc, xc, sprintf('%d (%d,%d)', no, conf{i}), 'HorizontalAlignment', 'center');
    if no >= conf{i}(1) && no <= conf{i}(2)
        hit(i) = 1;
    end
end
title(sprintf('%d segments, confidence hits=%.2f', numel(conf), mean(hit)));

end

function show_results(y_conf, y_sted, object, cell)
%

fig = figure();
for i = 1 : 4
    img = y_conf{i+1};
    subplot(2, 4, i);
    imagesc(img);
    axis image;
    title(sprintf('conf. %d-events: %d', i, sum(img(:))));
end
for i = 1 : 2
    img = y_sted{i+1};
    subplot(2, 4, 4+i);
    imagesc(img);
    axis image;
    title(sprintf('STED %d-events: %d', i, sum(img(:))));
end

subplot(2,4,7);
imagesc(object);
axis image;
title(sprintf('true mol.=%d',sum(object(:))));

subplot(2,4,8);
imagesc(cell);
axis image;
title('true object');

end