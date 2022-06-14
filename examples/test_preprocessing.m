function test_preprocessing()
% Tests the preprocessing algorithm by checking mean and distribution

no_trials = 100;

n= 20;
p = 0.02;
lambda = 0;
L = 32;
h = gaussian([L L],8);

no_pulses = 1e3;
params.md=4;  %number of detectors used for measurements
params.mD=4;  %number of channels used in reconstruction (<= md)
ord = 10;
params.mQ=ord;  %number of Q's used in approximation (<=mS)
params.mS=ord;  %number of S's used in approximation formula for Q (<= ms, >= mp)
params.ms=ord;  %number of s's to be determined (>=mS)
params.mp=ord;  %number of p's or lambda's used in formula for Q (<=mS)
params.maxit_preprocessing = 10; % number of maximal Newton iterations for preprocessing
params.tol_preprocessing = 1e-8; % Relative tolerance for Newton's method for preprocessing
params.preprocessing_inversion_method = 'gauss_newton';

s_me = zeros(no_trials,params.ms);
Y = cell(no_trials,params.md+1);
bar = progressbar('test preprocessing');
for i=1:no_trials
    Y(i,:) = full_forward_simulation_cluster(n,p,lambda,h,params.md,no_pulses);
    [S_me,Cov{i},it{i}] = preprocess_data(Y(i,:),no_pulses,lambda*ones(32,32)*no_pulses,params);
    for j=1:params.ms
        s_me(i,j)= sum(sum(S_me{j}));
    end
    bar.set(i / no_trials);
end
bar.close();
clear bar;

%Theoretical values
for i=1:params.ms
    truth(i) = n*p^i*sum(sum(h(:).^i));
end

% Covariance
W = zeros(params.ms,params.ms);
for i=1:no_trials
    for j=1:L
        for k=1:L
            W = W+Cov{i}{j,k};
        end
    end
end
W = W/no_trials;
s = sqrt(diag(W));

diff = s_me-repmat(truth,no_trials,1);

% display relative errors (S_me - truth) / truth
% figure();
fig = figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1 : 4
    subplot(2,2,i)
    data = diff(:,i)/truth(i);
    d = std(data);
    edges = -3*d-d/4:d/2:3*d;
    center = (edges(1:end-1)+edges(2:end))/2;
    histogram(data, edges, 'Normalization', 'pdf')
    f = normpdf(center,0,s(i)/truth(i));
    hold on;
    plot(center, f, 'ro');
    center = -3*d:d/20:3*d;
    f = normpdf(center,0,s(i)/truth(i));
    plot(center, f, 'r-');
    legend('Empirical errors','Estimated asymptotic errors')
    title(sprintf('%d-th order, <rel.err.>=%.2f, [15%%,85%%]=[%.2f,%.2f]', i, mean(data), quantile(data, [0.15, 0.85])));
    xlabel('rel. error');
end


end