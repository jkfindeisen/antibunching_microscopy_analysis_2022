function segments = boxes_to_segments(Y,Bnp,hset,pad,options)
% Converts boxes obtained by MISCAT to segments using the Watershed
% segmentation algorithm to obtain more natural segments

% INPUT:
% Y = unpadded data used for MISCAT
% Bnp = cell with rejections of the MISCAT algorithm
% hset = cell of used scales for MISCAT
% pad = padding used for MISCAT to avoind boundary effects

% OUTPUT:
% segments = segmentation of the image

% (c) Frank Werner, IMS Uni Göttingen, 03.01.2019 or later

%% First generate a watershed segmentation of the image
% first smooth the image
image_smoothed = gaussian_smooth(Y, options.smooth_width);
% apply relative threshold to determine bg
T = options.threshold * max(image_smoothed(:)) + (1 - options.threshold) * min(image_smoothed(:));
bg = image_smoothed < T;
image_smoothed(bg) = T;
% apply watershed
segments = watershed(-image_smoothed);
% set on background to 0
segments(bg) = 0;
% expand by one
for j = 1 : 2
    for i = 1 : max(segments(:))
        seg = segments == i;
        seg = imdilate(seg, ones(3));
        seg = seg & segments == 0;
        segments(seg) = i;
    end
end

% TODO: Parellelize what follows!

%% Now reduce complexity of boxes found by MISCAT
Bnp_red = reduce_complexity(Bnp, hset);

%% Generate boxes encoded in Bnp-red, but only those which are not in the padded region

n1 = size(Y,1);
n2 = size(Y,2);
boxes = [];
for i=1:numel(hset)
    for j=1:size(Bnp_red{i},1)
        ulx = Bnp_red{i}(j,1); % x-coordinate of upper left point of box
        uly = Bnp_red{i}(j,2); % y-coordinate of upper left point of box
        lrx = Bnp_red{i}(j,1)+hset{i}(1)-1; % x-coordinate of lower right point of box
        lry = Bnp_red{i}(j,2)+hset{i}(2)-1; % y-coordinate of lower right point of box
        % Candidate for next box, corrected for padding
        aux = [min(max(ulx-pad,1),n1) min(max(uly-pad,1),n2) min(max(lrx-pad,1),n1) min(max(lry-pad,1),n2)];
        if (aux(3)-aux(1)) * (aux(4)-aux(2))>0
            boxes = [boxes; aux];
        end
    end
end

% make segmenation valid
segments = make_segmentation_valid(boxes, segments);

end

function smoothed = gaussian_smooth(image, width)
% width = FWHM of gaussian kernel in pixel

% a grid
L = ceil(3*width);
[x,y] = ndgrid(-L:L,-L:L);

% gaussian kernel
k = exp(-4*log(2) * (x.^2 + y.^2) / width^2);
k = k / sum(k(:));

% convolve
smoothed = conv2(image, k, 'same');

end