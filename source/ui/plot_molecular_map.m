function res = plot_molecular_map(segmentation,br,conf,alpha)
% Plot results

% (c) Frank Werner, IMS Uni Göttingen, 03.01.2019 or later

%assert(nargin == 10, 'Not enough arguments');

% n = size(Y,1);
% 
% % Plot Data
% 
% figure('Position', scrsz);
% 
% imagesc(Y);
% axis square;
% colorbar;
% title('Data');

% Plot segments, counted numbers and uncertainty

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
    text(center(1)-5,center(2)-2,['\bf[',num2str(conf{i}(1)),',',num2str(conf{i}(2)),']']);
%     if nargin==10
%         if conf{i}(1) <= N_true(i) && N_true(i) <= conf{i}(2)
%             text(center(1)-1,center(2)+2,num2str(N_true(i)));
%         else
%             text(center(1)-1,center(2)+2,num2str(N_true(i)),'Color','r');
%         end
%     end
end

