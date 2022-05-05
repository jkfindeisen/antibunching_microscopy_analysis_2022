% s=cellsum(co,ct)
% sum up two cells or all elements of once cell array
% it works iteratively and, therefore, can sum up cells with same level of
% subcells.
% 
% co: cell with numbers as the lowest level
% ct: cell with numbers as the lowest level
% 
% s:  sum
% 
% for example:
% co = cellfun(@(x) {ones(5,5)},num2cell(1:3));
% ct = cellfun(@(x) {rand(5,5)},num2cell(1:3));
% s1=cellsum(co,ct);
% s2=cellsum(ct); % here only sumup one level
% 
% Haisen Ta 20120611 @ goettingen

function s=cellsum(co,ct)
assert(nargin>0 && nargin<3,'invalid number of inputs!');
if nargin<2
   s=co{1};
   for i=2:numel(co)
      s=cellsum(s,co{i});
   end
   return;
end

if ~strcmp( class(co),class(ct) )
   errorfun;
   return;
end
if ndims(co)~=ndims(ct) || any(size(co)~=size(ct))
   errorfun;
   return;
end
   
if iscell(co)
   if isempty(co) && isempty(ct)
      s=cell(size(co));
      return;
   end
   s=cellfun(@(x,y) {cellsum(x,y)},co,ct);
elseif (isnumeric(co) && isnumeric(ct) ) || (islogical(co) && islogical(ct))
   s=co+ct;
   return;
else
   error('not implemented for such data type!');
end
end

function errorfun()
   error('can not sum up cells with different stucture!');
end
