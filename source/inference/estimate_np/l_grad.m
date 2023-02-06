function grad = l_grad(x,Y,W,ord)
% Computes the gradient of the likelihood function
%
% INPUT:
% x = 2x1-vector with entries n (number of markers) and p (brightness of markers)
% W = weighting matrix
% Y = data
% ord =  order up to which data in Y is used
%
% OUTPUT:
% grad = gradient

tmp = x(1) * x(2).^(1:ord)';

% res =  (tmp-Y(1:ord))' * W(1:ord,1:ord) * (tmp-Y(1:ord))/2;

grad = [x(2).^(1:ord); x(1)*(1:ord).*(x(2).^(0:ord-1))] * W(1:ord,1:ord) * (tmp-Y(1:ord));

end