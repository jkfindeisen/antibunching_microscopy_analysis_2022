function test_preprocessing3()
% tests both, the case of n=10 and n=100 and p=0.02 of molecules at one
% position on different ways to make the preprocessing and for different
% number of pulses (no background for the moment)

optimized_test();

end

function optimized_test()
% optimized test (for speed)

pp = 0.02;
repetitions = 1e3;

order = 8;
method = 'fmincon';

% comparison fmincon for scaled, unscaled
for n = [10, 40, 100]
    for t = [1e3, 1e4, 1e5]
        test_it(n, pp, t, order, method, [], repetitions);
    end
end

% order = 4;
% method = 'newton';
% 
% % comparison fmincon for scaled, unscaled
% for n = [10, 40, 100]
%     for t = [1e3, 1e4, 1e5]
%         test_it(n, pp, t, order, method, [], repetitions);
%     end
% end

end


function full_test()

close all;

pp = 0.02;
repetitions = 1e3;

% newton
% order = 4;
% method = 'newton';
% for t = [1e3, 1e4, 1e5]
%     test_it(10, p, t, order, method, [], repetitions);
%     test_it(40, p, t, order, method, [], repetitions);
% end
% test_it(100, p, 1e5, order, method, [], repetitions); % won't work

% fmincon
order = 10;
method = 'lsqnonlin';

% % comparison initial_s for fmincon for scaled but not bounded
% p.scale_s = true; % true or false
% p.solver = 'fmincon'; % fmincon or lsqnonlin
% p.bounded = false; % true or false
% for n = [10, 40, 100]
%     p.initial_s = 'one'; % zeros, one, many
%     test_it(n, pp, 1e3, order, method, p, repetitions);
%     p.initial_s = 'zeros'; % zeros, one, many
%     test_it(n, pp, 1e3, order, method, p, repetitions);
%     p.initial_s = 'many'; % zeros, one, many
%     test_it(n, pp, 1e3, order, method, p, repetitions);
% end

% comparison fmincon for scaled and bounded/not bounded
% p.scale_s = true; % true or false
% p.solver = 'fmincon'; % fmincon or lsqnonlin
% p.initial_s = 'zeros'; % zeros, one, many
% for n = [10, 40, 100]
%     p.bounded = true; % true or false
%     test_it(n, pp, 1e3, order, method, p, repetitions);
%     % unbounded already done above (if it does something than for low
%     % amount of pulses)
%     p.bounded = true; % true or false
%     test_it(n, pp, 1e5, order, method, p, repetitions);
% end

% comparison fmincon for scaled, unscaled
p.solver = 'fmincon'; % fmincon or lsqnonlin
p.initial_s = 'zeros'; % zeros, one, many
p.bounded = false;
p.scale_s = false;
for n = [10, 40, 100]
    test_it(n, pp, 1e4, order, method, p, repetitions);
    % unscaled will run later or so
end

% run fmincon for scaled, unbounded
p.solver = 'fmincon'; % fmincon or lsqnonlin
p.initial_s = 'zeros'; % zeros, one, many
p.bounded = false;
p.scale_s = true;
for n = [10, 40, 100]
    for t = [1e3, 1e4, 1e5]
        test_it(n, pp, 1e4, order, method, p, repetitions);
    end
end

% run lsqnonlin for scaled, unbounded
p.solver = 'lsqnonlin'; % fmincon or lsqnonlin
p.initial_s = 'zeros'; % zeros, one, many
p.bounded = false;
p.scale_s = true;
for n = [10, 40, 100]
    for t = [1e3, 1e4, 1e5]
        test_it(n, pp, 1e4, order, method, p, repetitions);
    end
end

% comparison lsqnonlin for scaled, unscaled
p.solver = 'lsqnonlin'; % fmincon or lsqnonlin
p.initial_s = 'zeros'; % zeros, one, many
p.bounded = false;
p.scale_s = false;
for n = [10, 40, 100]
    test_it(n, pp, 1e4, order, method, p, repetitions);
    % unscaled will run later or so
end

% comparison initial_s for lsqnonlin for scaled but not bounded
p.scale_s = true; % true or false
p.solver = 'lsqnonlin'; % fmincon or lsqnonlin
p.bounded = false; % true or false
for n = [10, 40, 100]
    p.initial_s = 'one'; % zeros, one, many
    test_it(n, pp, 1e3, order, method, p, repetitions);
    p.initial_s = 'zeros'; % zeros, one, many
    test_it(n, pp, 1e3, order, method, p, repetitions);
    p.initial_s = 'many'; % zeros, one, many
    test_it(n, pp, 1e3, order, method, p, repetitions);
end

end

function first_run()

% what we have already
% method = 'newton';
% method_params = [];
% single_sweep(method, method_params);

% new ways
method = 'lsqnonlin';
p.scale_s = true; % true or false
p.initial_s = 'one'; % zeros, one, many
p.solver = 'fmincon'; % fmincon or lsqnonlin
p.bounded = false; % true or false
single_sweep(method, p);

% new ways
method = 'lsqnonlin';
p.scale_s = true; % true or false
p.initial_s = 'zeros'; % zeros, one, many
p.solver = 'fmincon'; % fmincon or lsqnonlin
p.bounded = false; % true or false
single_sweep(method, p);

% new ways
method = 'lsqnonlin';
p.scale_s = true; % true or false
p.initial_s = 'many'; % zeros, one, many
p.solver = 'fmincon'; % fmincon or lsqnonlin
p.bounded = false; % true or false
single_sweep(method, p);

% new ways
method = 'lsqnonlin';
p.scale_s = true; % true or false
p.initial_s = 'one'; % zeros, one, many
p.solver = 'fmincon'; % fmincon or lsqnonlin
p.bounded = true; % true or false
single_sweep(method, p);

% new ways
method = 'lsqnonlin';
p.scale_s = false; % true or false
p.initial_s = 'one'; % zeros, one, many
p.solver = 'fmincon'; % fmincon or lsqnonlin
p.bounded = false; % true or false
single_sweep(method, p);

% new ways
method = 'lsqnonlin';
p.scale_s = true; % true or false
p.initial_s = 'one'; % zeros, one, many
p.solver = 'lsqnonlin'; % fmincon or lsqnonlin
p.bounded = false; % true or false
single_sweep(method, p);

% new ways
method = 'lsqnonlin';
p.scale_s = true; % true or false
p.initial_s = 'zeros'; % zeros, one, many
p.solver = 'lsqnonlin'; % fmincon or lsqnonlin
p.bounded = false; % true or false
single_sweep(method, p);

% new ways
method = 'lsqnonlin';
p.scale_s = true; % true or false
p.initial_s = 'many'; % zeros, one, many
p.solver = 'lsqnonlin'; % fmincon or lsqnonlin
p.bounded = false; % true or false
single_sweep(method, p);

end

function single_sweep(method, method_p)

repetitions = 100;
p = 0.02;
for order = [10, 4]
    for t = [1e3, 1e4, 1e5]
        for n = [10, 40, 100]
            test_it(n, p, t, order, method, method_p, repetitions);
        end
    end
end

end

function test_it(n, p, t, order, method, method_p, repetitions)

lambda = 0;
L = 9;
psf = gaussian([L L],4);
object = zeros(L, L);
object(5,5) = n;

params.md=4;  %number of detectors used for measurements
params.mD=4;  %number of channels used in reconstruction (<= md)
params.mQ=order;  %number of Q's used in approximation (<=mS)
params.mS=order;  %number of S's used in approximation formula for Q (<= ms, >= mp)
params.mp=order;  %number of p's or lambda's used in formula for Q (<=mS)
params.ms=4;  %number of s's to be determined (>=mS)
params.maxit_preprocessing = 10; % number of maximal Newton iterations for preprocessing
params.tol_preprocessing = 1e-4; % Relative tolerance for Newton's method for preprocessing
params.preprocessing_inversion_method = method; % newton or lsqnonlin
params.preprocessing_inversion = method_p; % all the method parameters

% Theoretical values
for i=1:params.ms
    s_true(i) = n*p^i*sum(sum(psf(:).^i));
end

s_estimated = zeros(repetitions,params.ms);
Y = cell(repetitions,params.md+1);
bar = progressbar('test preprocessing');
time = 0;
Cov = cell(repetitions, 1);
for i=1:repetitions
    % Y(i,:) = full_forward_simulation_cluster(n,p,lambda,h,params.md,t);
    Y(i,:) = full_forward_simulation(object, p, lambda, psf, params.md, t);
    % preprocessing with time measurement
    tic();
    %     profile off;
    %     profile on;
    [S_me, Cov{i}] = preprocess_data(Y(i,:),t,lambda*ones(32,32)*t,params);
    %     profile viewer;
    time = time + toc();
    for j=1:params.ms
        s_estimated(i,j)= sum(sum(S_me{j}));
    end
    bar.set(i / repetitions);
end
time = time / repetitions;
bar.close();
clear bar;

% Covariance
W = zeros(params.ms,params.ms);
for i=1:repetitions
    for j=1:L
        for k=1:L
            W = W+Cov{i}{j,k};
        end
    end
end
W = W / repetitions;
std_theo = sqrt(diag(W));

% display absolute errors
% figure();
fig = figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1 : 4
    subplot(2,2,i)
    s = s_estimated(:, i);
    tr = s_true(i); % the true value
    h = histogram(s, 'Normalization', 'pdf');
    center = (h.BinEdges(1:end-1)+h.BinEdges(2:end))/2;
    f = normpdf(center, tr, std_theo(i));
    hold on;
    plot(s_true(i) * [1,1], ylim(), 'm-');
    plot(center, f, 'ro-');
    xlabel(sprintf('s_%d', i));
    title(sprintf('bias %.2f%%, %.2f%%, emp-rel-std %.2f%%, cov-rel-std %.2f%%', (mean(s) - tr) / tr * 100, (median(s) - tr) / tr * 100, std(s) / tr * 100, std_theo(i) / tr * 100));
end
switch method
    case 'newton'
        s = sprintf('n=%d, t=1e%d, ord=%d, Newton', n, log10(t), order);
    case 'lsqnonlin'
        s = sprintf('n=%d, t=1e%d, ord=%d, solver=%s, scaled=%d, bounded=%d, initial=%s', n, log10(t), order, method_p.solver, method_p.scale_s, method_p.bounded, method_p.initial_s);
    case 'fmincon'
        s = sprintf('n=%d, t=1e%d, ord=%d, fmincon', n, log10(t), order);
    otherwise
        error('unknown method');
end
annotation(fig, 'textbox', [0.2 0.96 0.6 0.03],'String', {sprintf('%s, time=%.3fs', s, time)}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
saveas(fig, ['preprocessing_', s, '.jpg']);
close(fig);

end