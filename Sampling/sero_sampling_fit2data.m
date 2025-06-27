function S = sero_sampling_fit2data(T, tr, b, do_t1, W, te)
% function S = sero_sampling_fit2data(T, tr, b, do_t1, W, te)

if nargin < 6
    te = [];
end

switch do_t1
    case 0 % No T1 effects
        T1W = double(tr>0);

    case 1 % Simplistic T1 effect
        T1W = 1-exp(-tr./T(:,3)');

    case 2 % Accurate T1 effect (TR is time from 180 pulse).
        sp  = 1-exp(-te./2./T(:,3)');
        T1W = -sp + (1+sp).*(1-exp(-tr./T(:,3)'));

    otherwise
        error('mode not recognized!')
end

S = W .* T1W .* exp( -b*T(:,2)' + b.^2*T(:,4)'/2 ) * T(:,1);
S = abs(S);