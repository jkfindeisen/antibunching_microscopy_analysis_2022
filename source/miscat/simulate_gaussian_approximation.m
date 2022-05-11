function Psi = simulate_gaussian_approximation(filename, number_runs)
%

assert(nargin >= 1, 'Not enough arguments');
if nargin < 2
    number_runs = 1e4;
end

fprintf('-----------------------------------------------------------------\n');
fprintf('Performing %i pure noise runs for the test in %s!\n', number_runs, filename);

fftw('planner', 'exhaustive');

if ~exist(filename,'file')
    error('The file does not exist!\n');
end
load(filename);

if ~exist('Psi','var')
    error('The file does not contain a test!\n');
end

% If simulate_quantiles crashed or was interrupted for some reason, then
% there should be a variable 'realizations'
if exist('realizations','var')
    Psi.realizations = realizations;
end


save(Psi.filename,'Psi');

% Now check if there are already runs and if they suffice
if isfield(Psi,'realizations')
    if (size(Psi.realizations,2)>= number_runs)
        fprintf('Nothing to do! \n');
        return
    else
        lb = size(Psi.realizations,2);
    end
else
    Psi.realizations = [];
    lb = 0;
end


%% Compute precomputable quantities

fprintf('Performing preparations...');

tic;

if ~isfield(Psi,'Phi_FFTd')
    Psi = generate_testfunctions(Psi);
end

k = numel(Psi.hset);

fprintf('...done!\n');
toc;

%% evaluate statistic
tic;
tmp = zeros(k,1);
fprintf('Need to do %i runs in total!\n',number_runs-lb);
for j=lb+1:number_runs
    Y_FFTd = fft2(padarray(randn(Psi.n,Psi.n),[ceil(Psi.n/2) ceil(Psi.n/2)]));
    for i = 1:numel(Psi.hset)
        Th = abs(mean_values_fixed_scale(Y_FFTd,Psi.Phi_FFTd{i},Psi.n,Psi.hset{i}))/Psi.Phi_norms{i};
        tmp(i) = max(Th(:));
    end
    Psi.realizations = [Psi.realizations tmp];
    if (mod(j,100) == 0 || j == number_runs)
        fprintf('%i runs done ...\n',j-lb);
        save(Psi.filename,'-struct','Psi','realizations','-append');
        toc;
    end
end
Psi = rmfield(Psi,'Phi_FFTd');
Psi = rmfield(Psi,'Phisq_FFTd');
Psi = rmfield(Psi,'Phi_norms');

save(Psi.filename,'Psi');

fprintf('Finished!\n');
fprintf('-----------------------------------------------------------------\n');

end