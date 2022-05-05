function res = wh(h,C)

assert(nargin == 2, 'Not enough arguments');

%res = sqrt(2*log(C/prod(h))) + log(2*log(C/prod(h)))/2/sqrt(2*log(C/prod(h)));

inner = sqrt(2*log(C./prod(h,2)));
res = inner + log(inner)./inner;

end