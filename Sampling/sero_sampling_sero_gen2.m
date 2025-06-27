function [posi, b, sthi] = sero_sampling_sero_gen2(nShot, ar, nz, bv, in_seed)
% function [pos, b, sthi] = sero_sampling_sero_gen2(nShot, ar, nz, bv, in_seed)
%
% Second generation SERO sampling protocol that only produces positions
% (pos), b-values (b) and slice thicknesses.

if nargin < 4
    bv = [1.2 0.8 0.4 0.1];
end

if nargin < 5
    in_seed = 1234;
end
rng(in_seed);


b  = repmat(bv', ceil(nShot/numel(bv)), 1);
b  = b(1:nShot);
b  = b(randperm(numel(b)));

n_history = 10;

xold = ones(1,n_history)*inf;

x0 = repmat(1:(nz-ar+1), ceil(nShot/(nz-ar+1)),1);

x0 = x0(:);

x0l = x0(randperm(numel(x0)));

posi = zeros(nShot, 1);
sthi = zeros(nShot, 1);

for i = 1:nShot

    x0 = x0l(i);

    for j = 1:n_history
        if any( abs(x0-xold)<(ar/2) )
            tmp =  x0l(i:end);
            x0l(i:end) = tmp(randperm(numel(tmp)));
            x0 = x0l(i);
        end
    end

    ind = x0:(x0+ar-1);

    posi(i) = (mean(ind)-1/2)-nz/2;
    sthi(i) = ar;

    xold(1)=[];
    xold(end+1) = x0;
end
