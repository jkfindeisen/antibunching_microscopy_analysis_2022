function test_preprocessing2()
% tests especially the case of n=100 identical molecules at one position
% with p=0.02

close all;

t = 1e5;
n = 100;
p = 0.02;
lambda = 0;
R = 0;
rng_seed = [];

% create object and psf
C = R+1;
dims = (C+R)* [1,1];

object = zeros(dims);
object(C, C) = n;

psf = zeros(dims);
psf(C, C) = 1;


ord = 4;
params.md=ord;  %number of detectors used for measurements
params.mD=ord;  %number of channels used in reconstruction (<= md)
params.mQ=ord;  %number of Q's used in approximation (<=mS)
params.mS=ord;  %number of S's used in approximation formula for Q (<= ms, >= mp)
params.ms=ord;  %number of s's to be determined (>=mS)
params.mp=ord;  %number of p's or lambda's used in formula for Q (<=mS)
params.maxit_preprocessing = 10; % number of maximal Newton iterations for preprocessing
params.tol_preprocessing = 1e-4; % Relative tolerance for Newton's method for preprocessing
% params.preprocess_inversion_method = 'lsqnonlin';

% theoretical values
for i=1:params.ms
    truth(i) = n*p^i*sum(sum(psf(:).^i));
end
disp(truth);

Y = full_forward_simulation(object, p, lambda, psf, params.md, t, rng_seed);
[S_me,Cov] = preprocess_data2(Y,t,lambda*ones(32,32)*t,params);
for j=1:params.ms
    s_me(j)= sum(sum(S_me{j}));
end
disp(s_me);

end