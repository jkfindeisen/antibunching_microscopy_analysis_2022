function Psi = setup_test(n,kernel_params,otf,betamax,hset,filename)
% Provides a struct Psi containing all the test information

assert(nargin == 6, 'Not enough arguments');

% General parameters

Psi.n = n;
Psi.otf = otf;
Psi.kernel = kernel_params;
Psi.hset = hset;
Psi.betamax = betamax;
Psi.filename = filename;

% Probe functionals
Psi.phi = @(x,y) 1/Psi.n^2 * betapdf(x,Psi.betamax(1)+1,Psi.betamax(1)+1).*betapdf(y,Psi.betamax(2)+1,Psi.betamax(2)+1);

% Penalization
Psi.C = 3;
Psi.pen = @(m) repmat(wh(cell2mat(Psi.hset)/Psi.n,Psi.C),1,size(m,2)).*(m-repmat(wh(cell2mat(Psi.hset)/Psi.n,Psi.C),1,size(m,2)));
Psi.pen_inv = @(q) q./wh(cell2mat(Psi.hset)/Psi.n,Psi.C) + wh(cell2mat(Psi.hset)/Psi.n,Psi.C);

save(filename,'Psi');

end