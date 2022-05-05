% Tests the estimation algorithm by checking mean and distribution
close all;
clear all;

init;

no_trials = 10000;

n= 10;
p = 0.01;
h = 1; %Gaussian([32 32],4);
segments = 1; % zeros(32,32); segments(8:24,8:24) = 1;

alpha = 0.1; %Confidence level

params.FWHM_CONF = 1;
params.mc = 4; %number of s's used for counting
params.no_neighboring_px = 1;
params.maxit_counting = 10;
params.tol_counting = 1e-4;

%Truth
s = n*p.^(1:params.mc)';
stdev = 0.1 * s;
Sigma{1} = diag(stdev.^2);

nhat = zeros(no_trials,1);
phat = zeros(no_trials,1);
conf = zeros(no_trials,2);
type1 = zeros(no_trials,1);
it = zeros(no_trials,1);

bar = progressbar('test estimate np');
for i=1:no_trials
    noise = stdev .* randn(params.mc,1);
    [nest,pest,confest,it(i)] = estimate_np(num2cell(s+noise),Sigma,segments,alpha,params);
    nhat(i) = nest{1};
    phat(i) = pest{1};
    conf(i,:) = confest{1};
    % Check if type 1 error
    if conf(i,1) > n || n > conf(i,2)
        type1(i) = 1;
    end
    bar.set(i / no_trials);
end
bar.close();
clear bar;

fprintf('Empirical type 1 error is %3.2f%%\n',sum(type1)/no_trials*100);

figure(); 
subplot(1,3,1)
histogram(nhat)
title('errors in estimating n')

subplot(1,3,2)
histogram(phat)
title('errors in estimating p')

subplot(1,3,3)
histogram(it)
title('number of Newton iterations')