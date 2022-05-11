function dQ = dgds(s, lambda, params)
% Computes the derivative of the nonlinear transform
% s_1, ..., s_{params.ms},lambda -> Q
% w.r.t s
%
% params.f(i) = (i-1)! (precomputed to save time)
%
% (c) Frank Werner, IMS Uni Göttingen, 03.01.2019 or later

% Compute the S's
S = zeros(params.max+1,1);
S(1) = 1;
for k=2:params.mS+1
    for j=1:k-1
        % S(k) = S(k) + (-1)^(j+1)*factorial(k-2)/factorial(k-1-j).*s(j).*S(k-j);
        S(k) = S(k) + (-1)^(j+1)*params.f(k-1)/params.f(k-j).*s(j).*S(k-j);
    end
end

% Compute dSds
dS = zeros(params.max+1,params.max);
for k=2:params.mS+1
    for i=1:params.max
        for j=1:k-1
            if j==i
                % dS(k,i) = dS(k,i) + (-1)^(j+1)*factorial(k-2)/factorial(k-1-j).*S(k-j) + (-1)^(j+1)*factorial(k-2)/factorial(k-1-j).*s(j).*dS(k-j);
                dS(k,i) = dS(k,i) + (-1)^(j+1)*params.f(k-1)/params.f(k-j).*S(k-j) + (-1)^(j+1)*params.f(k-1)/params.f(k-j).*s(j).*dS(k-j,i);
            else
                % dS(k,i) = dS(k,i) + (-1)^(j+1)*factorial(k-2)/factorial(k-1-j).*s(j).*dS(k-j,i);
                dS(k,i) = dS(k,i) + (-1)^(j+1)*params.f(k-1)/params.f(k-j).*s(j).*dS(k-j,i);
            end
        end
    end
end

% Compute dQds
dQ = zeros(params.max+1,params.max);
for k=1:params.mQ+1
    for l=0:k-1
        for j=0:max(0,params.mp-(k-1))
            % dQ(k,:) = dQ(k,:) + exp(-lambda)*lambda^(k-1-l)/factorial(l)/factorial(k-1-l)*(-1)^j/factorial(j)*dS(j+l+1,:);
            dQ(k,:) = dQ(k,:) + exp(-lambda)*lambda^(k-1-l)/params.f(l+1)/params.f(k-l)*(-1)^j/params.f(j+1)*dS(j+l+1,:);
        end
    end
end

end