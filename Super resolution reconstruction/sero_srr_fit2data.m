function S = sero_srr_fit2data(T, tr, b, do_t1, do_v, W, te)
% function S = sero_srr_fit2data(T, tr, b, do_t1, W, te)

if nargin < 7
    te = [];
end

if ~do_v
    T(:,4) = 0;
end

switch do_t1
    case 0 % No T1 effects
        T1W = double(tr>0);

    case 1 % Simplistic T1 effect
        T1W = 1-exp(-tr./T(:,3)');

    case 2 % Accurate T1 effect (including recovery on both sides of 180 degree pulse).
        te = 0.08; %seconds
        % sp  = 1-exp(-te./2./T(:,3)');
        % T1W = -sp + (1+sp).*(1-exp(-tr./T(:,3)'));
        
        ti = tr - te/2;
        T1W = 1 - 2*exp(-ti./T(:,3)') + exp(-tr./T(:,3)');

    otherwise
        error('mode not recognized!')
end

S = W .* T1W .* exp( -b*T(:,2)' + b.^2*T(:,4)'/2 ) * T(:,1);
S = abs(S);