function ddd()
% Tests the mle algorithm by checking mean and distribution
close all;

no_trials = 1000;

n= 10;
p = 0.01;

k = 4; %number of s's used for counting
maxit = 10;
tol = 1e-4;

%Truth
s = n*p.^(1:k)';
stdev = 0.1 * s;
Sigma = diag(stdev.^2);
W = diag(1./stdev);

nhat = zeros(no_trials,1);
phat = zeros(no_trials,1);
type1 = zeros(no_trials,1);
Cov = cell(no_trials,1);
it = zeros(no_trials,1);

bar = progressbar('test mle');
for i=1:no_trials
    noise = stdev .* randn(k,1);
    [nhat(i),phat(i),Cov{i},it(i)] = mle(s+noise,Sigma,W,k,maxit,tol);
    bar.set(i / no_trials);
end
bar.close();
clear bar;

% Covariance
W = zeros(2,2);
for i=1:no_trials
    W = W+Cov{i};
end
W = W/no_trials;
s = sqrt(diag(W));

figure();
subplot(1,3,1)
data = nhat;
histogram(data,'Normalization','pdf')
hold on;
grid = min(data):(max(data)-min(data))/200:max(data);
f = normpdf(grid,n,s(1));
plot(grid,f);
legend('Empirical errors','Estimated asymptotic errors')
title('relative errors for n')
hold off;

subplot(1,3,2)
data = phat;
histogram(data,'Normalization','pdf')
hold on;
grid = min(data):(max(data)-min(data))/200:max(data);
f = normpdf(grid,p,s(2));
plot(grid,f);
legend('Empirical errors','Estimated asymptotic errors')
title('relative errors for p')
hold off;

subplot(1,3,3)
histogram(it)
title('number of Newton iterations')

end