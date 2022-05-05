function V = estimate_variance(Y, dist_params)
% Estimate variance from data under the given distribution assumptions

assert(nargin == 2, 'Not enough parameters');

switch dist_params.name
    case {'Gauss','Students t'}
        if isfield(dist_params,'sigma')
            % If variance is known, use this one
            V =  dist_params.sigma^2;
        else
            % Otherwise estimate variance
            V = var(Y(:));
        end
    case 'Poisson'
        assert(all(Y(:) >= 0), 'Data cannot be negative for Poisson distributed noise');
        V = (Y + dist_params.bg) / dist_params.t;
    case 'Binomial'
        assert(all(Y(:) >= 0), 'Data cannot be negative for Binomial distributed noise');
        V = (Y + dist_params.bg) .* (1-(Y + dist_params.bg)) / dist_params.t;
    case 'Mixture'
        if isfield(dist_params,'g')
            % If variance is known, use this one
            V = (dist_params.g+dist_params.bg)/dist_params.t + dist_params.sigma^2;
        else
            % Otherwise estimate variance
            V = max((Y+dist_params.bg)/dist_params.t + dist_params.sigma^2,dist_params.bg/dist_params.t + dist_params.sigma^2);
        end
    otherwise
        error('Not yet implemented!');
end

end