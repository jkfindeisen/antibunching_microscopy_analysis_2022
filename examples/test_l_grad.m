function test_l_grad()
% test_l checks by finite differences if l_grad and l_hess are implemented correctly
order = 4;

x = randn(2,1);
h = randn(2,1);

Y = randn(order,1);
tmp = randn(order);
Sigma = tmp'*tmp;

tvalues = 10.^(-1:-1:-7);

fprintf('\nTESTING l_grad AND l_hess BY FINITE DIFFERENCES:\n');
val = l_grad(x,Y,Sigma,order);
der = l_hess(x,Y,Sigma,order);
fprintf('||hess l(x)|| = %4.3e \n',norm(der));
fprintf('t \t\t||hess_l*h - (1/t)(grad_l(x+th)-grad_l(x))||\n');
for t = tvalues
    val2 = l_grad(x+t*h,Y,Sigma,order);
    diff = norm(der'*h-(val2-val)/t);
    fprintf('%2.1e \t %4.3e \n',t,diff);
end

end