function test_simulation()      
% tests the simulation part

close all;

% parameters
R = 10; % radius of simulation area
n = 100; % number of molecules
p = 0.02; % probability
lambda = 0; % background
d = 4; % number of detectors
t = 1e5; % number of pulses
rng_seed = 0;

% calculate with our formulas
s2 = zeros(d+1,1);
for i = 1 : d+1
    s2(i) = calc_di(i-1, d, n, p);
end

% create object and psf
C = R+1;
dims = (C+R)* [1,1];

object = zeros(dims);
object(C, C) = n;

psf = zeros(dims);
psf(C, C) = 1;

% forward simulation
img = full_forward_simulation(object, p, lambda, psf, d, t, rng_seed);

% analysis
s = cellfun(@(x) x(C,C), img);
s = s / t; % normalize by number of pulses

% visualization
figure;
plot(s, 'b:o');
hold on;
plot(s2, 'r:s');
legend('sim', 'bino');

end

function di = calc_di(i, d, N, p)
% Di

di = 0;
for j = i : min(N, 8)
    di = di + calc_qjd(j, N, p) * calc_wijd(i, j, d);
end

end

function wijd = calc_wijd(i, j, d)

wijd = stirling2(j, i) * factorial(d-1) / factorial(d-i) / d^(j-1);

end

function qkd = calc_qjd(k, N, p)

% add up qkd
qkd = 0;
for j = 0 : min(8, N) - k
    qkd = qkd + (-1)^j/factorial(j)*calc_Sk(j+k, N, p);
end
qkd = qkd / factorial(k);

end

function Sk = calc_Sk(k, N, p)

if k == 0
    Sk = 1;
    return;
end

% recursive_calculation
Sk = 0;
for j = 1 : k
    Sk = Sk + (-1)^(j+1)*factorial(k-1)/factorial(k-j)*calc_Sk(k-j, N, p)*calc_s(j, N, p);
end

end

function sk = calc_s(k, N, p)
% N equal molecules with prob. p

sk = N * p^k;

end