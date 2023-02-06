function segmentation = make_segmentation_valid(boxes,segmentation)
% Checks for all segments if they contain at least one box and if not,
% joins the nearest box into them if possible.
% Input: boxes = K x 4 matrix where each row is of the form
%                [lower left x, lower left y, upper right x, upper right y]
%        segmentation = integer matrix, where pixel value indicates to
%                       which segment the pixel belongs
% Output: segmentation = integer matrix, where pixel value indicates to
%                        which segment the pixel belongs
%
% Details:
% First we create a list, which box makes which segment valid. Out
% of this list, we first choose those boxes which are contained in a
% segment. These segments are valid.
% Then we remove those segments, which are not connected to any box (those
% are considered as artifacts).
% Now only segments are left which have to be enlarged. For each segment i
% and each box k which could make the segment i valid, we check what would
% happen if we use the box k. If the box k intersects another segment j or
% another box l which makes another segment j valid, this is stored in the
% list problems. If not, the pair i k is stored in solutions, together with
% the amount of enlargement.
% Afterwards, out of the solutions, the ones with smallest enlargement are
% applied.
% Now we check for each remaining segment i if there is a box k which can be
% merged such that any other remaining segment j still has another box to
% make it valid. If so, the corresponding box k is merged to the segment i
% All remaining segments have to be merged somehow to make them valid. All
% possibilities are computed and the final size of the corresponding
% segment is stored. Those are then solved iteratively by always creating
% the smallest possible joint segment.
% Finally, the values in segmentation have to be adjusted

fprintf('Making the segmentation valid...');

%% Number of segments and list of segments, which have to be validated
M = max(max(segmentation));
todo = uint32(1:M);

%% make a list which box could make which segment valid
candidates = uint32.empty;
for i=1:size(boxes,1)
    for j=todo
        % Check all remaining segments how much must be added by the box
        contained = sum(sum(segmentation(boxes(i,1):boxes(i,3),boxes(i,2):boxes(i,4)) ==j));
        total = (boxes(i,4)+1-boxes(i,2))*(boxes(i,3)+1-boxes(i,1));
        if contained > 0
            candidates = [candidates; i uint32(j) contained total];
        end
    end
    % Here we could store those boxes, which do not intersect with any
    % segments. If wanted, we could check at the end for which this is
    % still true and create new segments out of those.
end

%% find those segments which contain at least one box (those are easy)
I = candidates(:,3)==candidates(:,4);
todo = setdiff(todo,unique(candidates(I,2)));
% now all candidates which are contained in such a segment can be neglected
candidates = candidates(~I,:);

%% now remove all segments for which no candidate at all was found
I = setdiff(todo,unique(candidates(:,2)));
for i=I
    segmentation(segmentation==i) = 0;
end
todo = setdiff(todo,I);

if ~isempty(todo)
    %% create lists of problems and solutions
    solutions = double.empty;
    problems = [];
    % for each segment and possible box
    for i=todo
        I = candidates(:,2) == i; % candidates for i
        for k=unique(candidates(I,1))'
            % Now check for any other segment and any other box if box k can be
            % merged to segment i without trouble
            match = 1;
            for j=setdiff(uint32(1:M),i)
                % if segment j overlaps with box k, then box is problematic
                if any(candidates(:,1) == k & candidates(:,2) == j)
                    problems = [problems; i k j 0];
                    match = 0;
                else
                    % if not, then check for all possible boxes
                    J = candidates(:,2) == j & candidates(:,2) ~=i; %candidates for j
                    % note that boxes in I overlap with segment i but not with segment j
                    for l = unique(candidates(J,1))'
                        % check if boxes overlap
                        if boxes(k,1) <= boxes(l,3) && boxes(k,3) >= boxes(l,1) && boxes(k,2) <= boxes(l,4) && boxes(k,4) >= boxes(l,2)
                            % FOUND a match which cause trouble
                            problems = [problems; i k j l];
                            match = 0;
                        end
                    end
                end
            end
            if match == 1
                % box k fits
                contained = sum(sum(segmentation(boxes(k,1):boxes(k,3),boxes(k,2):boxes(k,4)) ==i));
                total = (boxes(k,4)+1-boxes(k,2))*(boxes(k,3)+1-boxes(k,1));
                solutions = [solutions; double(k) double(i) contained/total]; %Here the values are assigned
            end
        end
    end
    % now solutions is a list with all boxes which can be joined to segments
    % without causing trouble (and the corresponding value)
    % problems is a list which contains all matchings segments to boxes which
    % might cause trouble as either the boxes overlap or the box interesects
    % more than one segment
    
    %% Now choose for each segment for which we have a solution the best one
    if ~isempty(solutions)
        for i = unique(solutions(:,2))'
            % find best possible box
            I = solutions(:,2)==i;
            val = max(solutions(I,3));
            k = solutions(find(solutions(:,2)==i & solutions(:,3)==val,1),1);
            segmentation(boxes(k,1):boxes(k,3),boxes(k,2):boxes(k,4)) = i;
            todo = setdiff(todo,i);
            % remove i from problems list
            if ~isempty(problems)
                problems(problems(:,1)==i,:) = -1;
                problems = reshape(problems(problems~=-1),[],4);
            end
        end
    end
    
    %% now handle those segments which can be made valid without merging
    for i=todo
        for k=unique(problems(problems(:,1)==i,2))'
            problems2 = problems;
            % Do we loose segments by joining box k to segment i?
            tmp = problems(problems(:,1)==i & problems(:,2)==k,3:4);
            if all(tmp(:,2)~=0)
                problems2(ismember(problems(:,1),tmp(:,1)) & ismember(problems(:,2),tmp(:,2)),:) = -1;
                problems2 = reshape(problems2(problems2~=-1),[],4);
                if numel(setdiff(todo,unique(problems2(:,1)))) == 0
                    % no problem!
                    segmentation(boxes(k,1):boxes(k,3),boxes(k,2):boxes(k,4)) = i;
                    todo = setdiff(todo,i);
                    % update problems list
                    problems(problems(:,1)==i,:) = -1;
                    problems(problems(:,3)==i & problems(:,4)==k,4) = 0;
                    problems = reshape(problems(problems~=-1),[],4);
                    break;
                end
            end
        end
    end
    clear problems2;
    
    %% Compute all possible solutions for the remaining segments by merging
    solutions = [];
    sets = cell(0,1);
    for i=todo
        % Check with which other segments this one could be merged and compute
        % how large the new segment would have to be
        for j=unique(problems(problems(:,1)==i,3))'
            config = [0 0 Inf];
            for k=unique(problems(problems(:,1)==i & problems(:,3)==j,2))'
                for l=unique(problems(problems(:,1)==i & problems(:,3)==j & problems(:,2)==k,4))'
                    tmp = zeros(size(segmentation));
                    % create a list of segments which is connected to box k
                    list = [];
                    for m=1:M
                        if sum(sum(segmentation(boxes(k,1):boxes(k,3),boxes(k,2):boxes(k,4)) ==m))>0
                            tmp(segmentation==m) = 1;
                            list = [list m];
                        end
                    end
                    tmp(boxes(k,1):boxes(k,3),boxes(k,2):boxes(k,4)) = 1;
                    if l>0
                        for m=setdiff(1:M,list)
                            if sum(sum(segmentation(boxes(l,1):boxes(l,3),boxes(l,2):boxes(l,4)) ==m))>0
                                tmp(segmentation==m) = 1;
                                list = [list m];
                            end
                        end
                        tmp(boxes(l,1):boxes(l,3),boxes(l,2):boxes(l,4)) = 1;
                    end
                    val = sum(sum(tmp));
                    if val<config(3)
                        config =[k l val];
                        current_list = list;
                    end
                end
            end
            solutions = [solutions; i j config];
            sets = [sets current_list];
        end
    end
    % solutions is now a list which contains all possible matchings of segments
    % i with segments j, which boxes are necessary to obtain this matching, and
    % how large the resulting segment would be
    % sets is a cell of the same size as solutions, and each entry is a list of
    % those segments which would have to be merged if this solution is chosen
    
    %% Now solve iteratively, starting with the smallest resulting segment
    while ~isempty(solutions)
        % find cheapest solution
        I = find(solutions(:,5) == min(solutions(:,5)),1);
        % Merge segments as computed
        list = sets{I};
        val = min(list);
        for i=list
            segmentation(segmentation==i) = val;
        end
        % join computed box
        k=solutions(I,3);
        segmentation(boxes(k,1):boxes(k,3),boxes(k,2):boxes(k,4)) = val;
        % join other box if necessary
        l=solutions(I,4);
        if l>0
            segmentation(boxes(l,1):boxes(l,3),boxes(l,2):boxes(l,4)) = val;
        end
        % Remove all entries from solutions, which are for segments which have
        % just been solved
        todo = setdiff(todo,uint32(list));
        sets = sets(ismember(solutions(:,1),todo));
        solutions = solutions(ismember(solutions(:,1),todo),:);
        % Update solutions lists
        for j=1:numel(sets)
            if any(ismember(list,sets{j}))
                % corresponding solution has changed...compute new value
                tmp = zeros(size(segmentation));
                tmp(segmentation==val) = 1;
                for m=sets{j}
                    tmp(segmentation==m) = 1;
                end
                k=solutions(j,3);
                tmp(boxes(k,1):boxes(k,3),boxes(k,2):boxes(k,4)) = 1;
                l=solutions(j,4);
                if l>0
                    tmp(boxes(l,1):boxes(l,3),boxes(l,2):boxes(l,4)) = val;
                end
                solutions(j,:) = [solutions(j,1) val k l sum(sum(tmp))];
                sets{j} = union(setdiff(sets{j},list),val);
            end
        end
        % Now we might have some multiple entries, but it doesn't matter, as we
        % always choose the cheapest
    end
end

%% Readjust segmentation values
nhood = [0 1 0; 1 1 1; 0 1 0];
i = 1;
val = min(min(segmentation(segmentation>0)));
while val <= max(max(segmentation))
    % try to fill holes or gaps in the segment, currently only one pixel
    A = (segmentation==val); %Current segment
    A = A | imerode(imdilate(A,true(3)),nhood);
    % set segment with value val to value i
    segmentation(A) = i;
    % update
    i = i+1;
    val = val+1;
    % find next value
    while val < max(max(segmentation)) && all(all(segmentation~=val))
        val = val+1;
    end
end

fprintf('...done!\n');

end