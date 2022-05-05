function test_g()
% test_g checks by finite differences if g and dgds are implemented correctly

params.mQ=10;  %number of Q's used in approximation (<=mS)
params.mS=10;  %number of S's used in approximation formula for Q (<= ms, >= mp)
params.ms=10;  %number of s's to be determined (>=mS)
params.mp=10;  %number of p's or lambda's used in formula for Q (<=mS)


params.max = max([params.mQ params.mS params.ms params.mp]);
params.f = factorial(0:params.max);

s = randn(params.ms,1);
h = randn(params.ms,1);
lambda = abs(randn(1));

tvalues = 10.^(-1:-1:-7);

fprintf('\nTESTING g AND dgds BY FINITE DIFFERENCES:\n');
P= g(s,lambda,params);
der = dgds(s,lambda,params);
fprintf('||dgds[s]|| = %4.3e \n',norm(der));
fprintf('t \t\t||dgsg[s]h - (1/t)(g(s+th)-g(s))||\n');
for t = tvalues
    Y = g(s+t*h,lambda,params);
    diff = norm(der*h-(Y-P)/t);
    fprintf('%2.1e \t %4.3e \n',t,diff);
end

end