function res = plot_molecular_map(segmentation,br,conf,alpha,n_true)
% Plot results

% Plot segments, counted numbers and uncertainty

if nargin < 5
    n_true = [];
end

scrsz = get(0,'ScreenSize');
figure('Position', scrsz);

res = zeros(size(segmentation));
for i=1:max(max(segmentation))
    res = res + br{i}*(segmentation==i);
end

imagesc(res);
title(['Segmentation and molecule numbers with uncertainty at uniform level ',num2str(100*alpha),'%']);
axis square;
colorbar;
% jet colormap with white background
cmap = jet(64);
cmap(1, :) = 1;
colormap(cmap);

for i=1:max(max(segmentation))
    % Compute center of mass of set{i}
    tmp = regionprops(segmentation==i,'Centroid');
    center = tmp.Centroid;
    if ~isempty(n_true)
        if conf{i}(1) <= n_true(i) && n_true(i) <= conf{i}(2)
            text(center(1)-5,center(2)-2, sprintf('\\bf[%d, %d]: %d', conf{i}, n_true(i)));
        else
            text(center(1)-5,center(2)-2, sprintf('\\bf[%d, %d]: %d', conf{i}, n_true(i)),'Color','r');
        end
    else
        text(center(1)-5,center(2)-2, sprintf('\\bf[%d, %d]', conf{i}));
    end
end

