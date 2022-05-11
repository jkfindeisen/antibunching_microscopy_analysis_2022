function hess = l_hess(x,Y,W,ord)
% Computes the hessian of the likelihood function
%
% INPUT:
% x = 2x1-vector with entries n (number of markers) and p (brightness of markers)
% W = weighting matrix
% Y = data
% ord =  order up to which data in Y is used
%
% OUTPUT:
% hess = hessian
%
% (c) Frank Werner, IMS Uni Göttingen, 25.02.2019 or later

tmp = x(1) * x(2).^(1:ord)';

% res =  (tmp-Y(1:ord))' * W(1:ord,1:ord) * (tmp-Y(1:ord))/2;

% grad = [x(2).^(1:ord); x(1)*(1:ord).*(x(2).^(0:ord-1))] * W(1:ord,1:ord) * (tmp-Y(1:ord));

tmp2 = (x(1) *(1:ord).*(x(2).^(0:ord-1)))';

C = x(2).^(1:ord);
D = (1:ord).*(x(2).^(0:ord-1));
E = (0:ord-1).*(1:ord).*x(2).^(-1:ord-2);

hess  = zeros(2,2);
hess(1,1) = C*W(1:ord,1:ord)*C';
hess(1,2) = D*W(1:ord,1:ord)*(tmp-Y(1:ord))+C*W(1:ord,1:ord)*tmp2;
hess(2,1) = hess(1,2);
hess(2,2) = (x(1)*E*W(1:ord,1:ord)*(tmp-Y(1:ord)))+(tmp2'*W(1:ord,1:ord)*tmp2);

end