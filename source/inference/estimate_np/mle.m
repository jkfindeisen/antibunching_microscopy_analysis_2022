function [n,p,Cov,it] = mle(Y,Sigma,W,k,maxit,tol, pmin)
% Computes maximum likelihood estimator according to the likelihood
% l(n,p,Y) = 1/2 (F(n,p) - Y)^T W (F(n,p) - Y)
% where F(n,p) = (n*p, n*p^2, ... n^p^k)
%
% INPUT:
% Y = data
% Sigma = covariance matrix of data
% W = weighting matrix for maximum likelihood estimation (ideally: W = Sigma^{-1})
% k = order up to which data in Y is used
% maxit = maximal number of Newton iterations
% tol = tolerance for Newton iteration
%
% OUTPUT:
% n = estimated number
% p = estimated brightness
% Cov = covariance matrix of estimators
% it = number of Newton iterations (-1 if the Newton method did not converge and the initial guess was used)
%
% (c) Frank Werner, IMS Uni Göttingen, 25.02.2019 or later

if nargin < 7
    pmin = 0.001;
end

% use Newton's method to solve l_grad(x) = 0
p0 = max(Y(2)/Y(1), pmin); % we have a minimal p0
n0 = max(Y(1)/p0, 1);
x = [n0; p0];
it = 0;
res = @(x) norm(l_grad(x,Y,W,k),2);
res0 = res(x);
while res(x)/res0 > tol && it<maxit
    %fprintf('Objective value: %f\t Residual: %f\n',obj(x),res(x));
    x = x-l_hess(x,Y,W,k)\l_grad(x,Y,W,k);
    x = max(x, [1; pmin]); % n >= 1 and p >= pmin, whatever we do
    it = it+1;
end
if it == maxit
    % Newtons method did not converge. Use initial guess
    n = n0;
    p = p0;
    it = -1; % set to -1 to indicate no convergence
else
    n = x(1);
    p = x(2);
end

% Compute covariance matrix
% minimzing l is like solving l_grad = 0, hence, the new covariance matrix is
B = [p.^(1:k); n*(1:k).*(p.^(0:k-1))]; % gradient of F
C = B*W*B';
Cov = C\(B*W*Sigma*W'*B')/C';

% I thought the correct formula should be
%A = l_hess(x,zeros(size(Y)),W,k);
%Cov = A\(B*W*Sigma*W'*B')/A';
% Note that l_hess consists of a part like hess_F*F + grad_F^2 and in the
% above formula we hence just neglect the hess_F part....
% For details, see molecule_inference.odf, pages 3-4.
end