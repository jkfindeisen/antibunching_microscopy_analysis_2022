function Q = g(s, lambda, params)
% Computes the nonlinear transform
% s_1, ..., s_{params.ms},lambda -> Q
% where Q are the probabilities that exactly k photons (k=1,...,params.mQ) are
% emitted and recorded during a single excitation pulse
% Q{k} is given by exp(-lambda) sum_{l=0}^k lambda^(k-l)/l!/(k-l)! sum_{j=0}^{ord-k} (-1)^j/j! S_{j+l}
%
% params.f(i) = (i-1)! (precomputed to save time)

% (c) Frank Werner, IMS Uni Göttingen, 03.01.2019 or later

% Compute the S's
S = zeros(params.max+1,1);
S(1) = 1;
for k=2:params.mS+1
    for j=1:k-1
        % S(k) = S(k) + (-1)^(j+1)*factorial(k-2)/factorial(k-1-j).*s(j).*S(k-j);
        S(k) = S(k) + (-1)^(j+1)*params.f(k-1)/params.f(k-j) .* s(j) .* S(k-j);
    end
end

% Compute the Q's
Q = zeros(params.max+1,1);
for k=1:params.mQ+1
    for l=0:k-1
        for j=0:max(0,params.mp-(k-1))
            % Q(k) = Q(k) + exp(-lambda)*lambda^(k-1-l)/factorial(l)/factorial(k-1-l)*(-1)^j/factorial(j)*S(j+l+1);
            Q(k) = Q(k) + exp(-lambda)*lambda^(k-1-l)/params.f(l+1)/params.f(k-l)*(-1)^j/params.f(j+1)*S(j+l+1);
        end
    end
end

end