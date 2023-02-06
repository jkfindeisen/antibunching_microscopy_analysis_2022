function [n,p,conf, iterations] = estimate_np(s,Sigma,segments,alpha,params)
% Estimates the number of molecules in each of the sets under the
% assumption that the brightness is constant inside each set.
% Furthermore asymptotic confidence intervals at (uniform) level alpha are given
%
% INPUT:
% s = cell with empirical measurements of (np) * k, (np^2) * k^2, ...
% Sigma = cell of covariance matrices of the measurements (pixelwise)
% segments = segmentation on which n should be estimated. The sets are not allowed to overlap
% alpha = level (1-alpha is the confidence level)
% params
%       no_neighboring_px = number of neighboring pixels used for counting
%       maxit_counting = maximal number of Newton iterations used for counting
%       mc = number of s's used for counting
%       FWHM_CONF = full width at half maximum in pixels for confocal experiment
%       tol_counting = tolerance for Newton iteration
%
% OUTPUT:
% n = estimated number of molecules in the sets (as cell)
% p = estimated brightness of molecules in the sets (as cell)
% conf = confidence radius at level alpha for number of molecules in the sets (as cell)
% iterations = number of iterations (-1 if the Newton method did not converge and the initial guess was used)

ord = min(length(s),params.mc);

k = gaussian(size(s{1}),params.FWHM_CONF);

% First we sum over the corresponding measurements for each segment
% Therefore compute disjoint neighborhoods
K = max(segments(:));
mask = segments;
if K>0
    covered = segments>0;
    for i = 1 : params.no_neighboring_px
        [mask,covered] = extend_by_one(mask,covered);
    end
end

% Now sum over the segments
ssum = zeros(ord,K);
for i=1:ord
    for j=1:K
        ssum(i,j) = sum(sum(s{i}(mask==j)));
    end
    ssum(i,:) = ssum(i,:)/sum(k(:).^i);
end
% ssum = max(ssum, 0); % we know that the summed s cannot be negative (estimates of np^k)

% The same for the covariance matrices exploiting spatial independence
aux = eye(ord);
for i=1:ord
    aux(i,i) = 1/sum(k(:).^i);
end
W = cell(K,1);
for j=1:K
    W{j} = zeros(ord,ord);
    for i=1:size(k,1)
        for l=1:size(k,2)
            if mask(i,l)==j
                W{j} = W{j} + Sigma{i,l}(1:ord,1:ord);
            end
        end
    end
    W{j} = aux*W{j}*aux;
end

n = cell(K,1);
p = cell(K,1);
conf = cell(K,1);

calpha = icdf('Normal',1-alpha/2/double(K),0,1);

iterations = zeros(K, 1);
for i=1:K
    % Now estimate n and p by maxmimum likelihood:
    % [n,p] = argmin s^T W^(-1) s with s = (np,...,np^k)-Y (encoded in l)
    % we choose W to be diag(1./sqrt(diag(W{i}))) %% HERE %%
    [nhat,p{i},Cov,iterations(i)] = mle(ssum(:,i),W{i},diag(1./sqrt(diag(W{i}))),ord,params.maxit_counting,params.tol_counting);
    % n is encoded as first variable, hence the confidence interval is...
    sig = sqrt(Cov(1,1));
    lb = nhat - calpha*sig;
    ub = nhat + calpha*sig;
    % store n (and use rounded n for confidence interval)
    n{i} = max(nhat,1);
    conf{i} = [min(max(ceil(lb),1),round(n{i})) max(max(floor(ub),1),round(n{i}))];
end

end

function [mask, covered] = extend_by_one(mask, covered)
for i=1:max(mask(:))
    extension = imdilate(mask==i, true(3));
    mask(extension & ~covered) = i;
    covered = covered | extension;
end
end