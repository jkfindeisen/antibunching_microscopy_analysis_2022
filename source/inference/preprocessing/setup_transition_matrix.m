function A = setup_transition_matrix(ord,m)
% Setup the transition matrix, which converts probabilities of multiple
% photon emissions to probabilities of multiple photon detections using m
% detectors
% ord is the order up to which multiple photon emissions are modeled

A = zeros(m+1,ord+1);
A(1,1) = 1;
for k=2:ord+1
    A(2:min(k,m+1),k) = P(m,k-1);
end

end

function res = P(m,k)
% Computes the probability that exactly l out of m detectors are hit by k
% photons

res = zeros(min(m,k),1);

for l=1:min(m,k)
    res(l) = nchoosek(m,l)*(l/m)^k*B(k,l);
end

end

function res = B(k,l)
% Computes probability to hit all l urns with k balls (k>=l)
% Calculated recursively

if k<l
    res = 0;
    return;
end

if l==1
    if k>=1
        res = 1;
    else
        res = 0;
    end
    return;
end

if k==1
    if l>1
        res = 0;
    else
        res = 1;
    end
    return;
end

if k==l
    res = factorial(l)/l^l;
else
    res = 0;
    for j=1:k-l+1
        res = res + nchoosek(k,j)*(l-1)^(k-j)*B(k-j,l-1);
    end
    res = res/l^(k);
end

end