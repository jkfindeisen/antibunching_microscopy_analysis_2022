function Bnp = reduce_complexity(Bnp, hset)
% Reduces complexity by discarding all boxes with contain smaller ones

assert(nargin == 2, 'Not enough arguments');

fprintf('Discarding all boxes which contain smaller ones...');

% get box rectangles as list
boxes = [];
for i=1:numel(hset)
    boxes = [boxes; Bnp{i} repmat(hset{i},size(Bnp{i},1),1)];
end

% loop over box size
no_boxes = 0;
for i=1:numel(hset)
    current_size = hset{i}; % box size
    % smaller boxes (not larger and smaller in at least one dimension)
    I = (boxes(:,3) <= current_size(1) & boxes(:,4) <= current_size(2)) & (boxes(:,3) < current_size(1) | boxes(:,4) < current_size(2));
    smaller_boxes = boxes(I, :);
    current_boxes = Bnp{i}; % boxes of this size
    % for all boxes of this size
    for j=1:size(current_boxes,1)
        % get rectangle
        current = [current_boxes(j,1) current_boxes(j,2) current_boxes(j,1)+current_size(1)-1 current_boxes(j,2)+current_size(2)-1];
        % test if within smaller boxes
        if any(smaller_boxes(:,1)>= current(1) & smaller_boxes(:,2)>= current(2) & smaller_boxes(:,1)+smaller_boxes(:,3)-1<= current(3) & smaller_boxes(:,2)+smaller_boxes(:,4)-1<= current(4))
            % Smaller box found, discard this one
            current_boxes(j,:) = [0 0];
        end
    end
    current_boxes = current_boxes(current_boxes~=0);
    current_boxes = reshape(current_boxes,[],2);
    no_boxes = no_boxes + size(current_boxes,1);
    
    Bnp{i} = current_boxes;
end

fprintf('...done!\n --> %i boxes left\n',no_boxes);
end