function Psi = generate_testfunctions(Psi)
% 

assert(nargin == 1, 'Not enough arguments');

[y,x] = meshgrid(1/(2*Psi.n):1/Psi.n:1-1/(2*Psi.n));

phih = @(h) 1/prod(h) * Psi.phi(x/h(1),y/h(2));

% Phi_h = F^{-1}(F phih ./ F k)  with the Fourier transform F

Phih = @(h) convolve(phih(h), 1./Psi.otf);

Psi.Phisq_FFTd = cell(numel(Psi.hset),1);
Psi.Phi_FFTd = cell(numel(Psi.hset),1);
Psi.Phi_norms = cell(numel(Psi.hset),1);

for i=1:numel(Psi.hset)
    tmp = Phih(Psi.hset{i}/Psi.n);
    Psi.Phi_FFTd{i} = fft2(repmat(rot90(tmp,2),2,2));
    Psi.Phisq_FFTd{i} = fft2(repmat(rot90(tmp.^2,2),2,2));
    Psi.Phi_norms{i} = sqrt(sum(sum(tmp.^2)));
end

end