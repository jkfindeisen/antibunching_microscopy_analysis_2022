function res = l(x,Y,W,ord)
% Computes the likelihood function
% l(n,p,Y) = 1/2 (F(n,p) - Y)^T W (F(n,p) - Y)
% where F(n,p) = (n*p, n*p^2, ... n^p^k)
%
% INPUT:
% x = 2x1-vector with entries n (number of markers) and p (brightness of markers)
% W = weighting matrix
% Y = data
% ord =  order up to which data in Y is used
%
% OUTPUT:
% res = likelihood value
%
% (c) Frank Werner, IMS Uni Göttingen, 25.02.2019 or later

tmp = x(1) * x(2).^(1:ord)';

res =  (tmp-Y(1:ord))' * W(1:ord,1:ord) * (tmp-Y(1:ord))/2;

end