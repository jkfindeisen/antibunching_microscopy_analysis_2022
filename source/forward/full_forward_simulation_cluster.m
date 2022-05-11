function img = full_forward_simulation_cluster(n,p,lambda,psf,d,no_pulses)
% Performs the full forward simulation for a single cluster located at one
% pixel with n molecules of brightness p
%
% (c) Frank Werner, IMS Uni Göttingen, 03.01.2019 or later

img = cell(d+1,1);
for i=1:d+1
    img{i} = zeros(size(psf));
end

for i=1:no_pulses
    res = single_excitation(n,p,psf,lambda,d);
    for k=1:d+1
        img{k} = img{k} + (res==k-1);
    end
end

end

function res = single_excitation(n,p,psf,lambda,d)
% Simulates one single excitation pulse at each pixel

% photons from n molecules and from the background
no_photons = binornd(n,p*psf) + poissrnd(lambda*ones(size(psf)));
m = max(no_photons(:));
if m>1
    detectors = zeros(size(psf));
    for i=1:m
        tmp = m.^(randi(d,size(psf))-1); %transforms 1:d to [1 m m^2 ... m^(d-1)]
        detectors = detectors + tmp.*(no_photons>=i);
        %detectors stores the detectors hit in a m-base expansion to ensure
        %distinghuishability
    end
    res = zeros(size(psf));
    for i=1:d
        % we want to know if the ith detector was hit. Therefore we have to
        % check if the ith coefficient in the m-base expansion of detectors is
        % non-zero. The ith coefficient can be computed as
        % rem(floor(detectors*m.^(i-d))),m)
        res = res + (rem(floor(detectors*m.^(i-d)),m)>0);
    end
else
    res = no_photons;
end

end