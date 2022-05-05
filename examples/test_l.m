function test_l()
% test_l checks by finite differences if l and l_grad are implemented correctly
order = 4;

x = randn(2,1);
h = randn(2,1);

Y = randn(order,1);
tmp = randn(order);
Sigma = tmp'*tmp;

tvalues = 10.^(-1:-1:-7);

fprintf('\nTESTING l AND l grad BY FINITE DIFFERENCES:\n');
val = l(x,Y,Sigma,order);
der = l_grad(x,Y,Sigma,order);
fprintf('||grad l(x)|| = %4.3e \n',norm(der));
fprintf('t \t\t||<grad l,h> - (1/t)(l(x+th)-l(x))||\n');
for t = tvalues
    val2 = l(x+t*h,Y,Sigma,order);
    diff = norm(der'*h-(val2-val)/t);
    fprintf('%2.1e \t %4.3e \n',t,diff);
end

end