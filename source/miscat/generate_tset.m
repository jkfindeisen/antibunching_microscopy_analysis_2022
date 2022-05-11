function tset = generate_tset(Y,hset,dist_params)
% 

assert(nargin == 3, 'Not enough arguments');

tset = cell(numel(hset),1);

switch dist_params.name
    case {'Gauss','Students t','Mixture'}
        for i=1:numel(hset)
            tset{i} = true(size(Y));
        end
    case {'Poisson','Binomial'}
        n = size(Y,1);
        Y_FFTd = fft2(padarray(Y + dist_params.bg,[ceil(n/2) ceil(n/2)]));
        for i=1:numel(hset)
            tmp = zeros(size(Y));
            tmp(1:hset{i}(1),1:hset{i}(2)) = 1;
            tmp = mean_values_fixed_scale(Y_FFTd,fft2(repmat(rot90(tmp,2),2,2)),size(Y,1),hset{i});
            % use only those t which are in a box in which there was at
            % least 1 count
            tset{i} = tmp>=1/dist_params.t;
        end
    otherwise
        error('Not yet implemented!');
end

end