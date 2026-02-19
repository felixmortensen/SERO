function x0 = sero_srr_setGuess(S, TR, B, W, do_t1, do_v, sm_siz)
% function x0 = srp_setGuess(S, TR, B, W, do_t1, sm_siz)

z = size(TR,2);
TRnan = TR;
TRnan(TRnan==0) = nan;

if nargin < 4 || isempty(W)
    W = TR>0;
end

if nargin < 5 || isempty(do_t1)
    do_t1 = 1;
end

if nargin < 6 || isempty(do_v)
    do_v = 1;
end

if nargin < 7 || isempty(sm_siz)
    sm_siz = z/5;
end

x0 = zeros(z,4);

for i = 1:z
    ind = TR(:,i)>0;
    w   = sum(W(ind,:), 2);
    s   = S(ind)./w;
    tr  = mean(TRnan(ind,:), 2, 'omitnan');
    b   = B(ind);

    try
        x0(i,:) = sero_srr_data2fit_1d_reg_v2(s, tr, b, 1, do_t1, do_v, 0, 0);
    catch
    end
end


if sm_siz
    for i = 1:size(x0,2)
        x0(:,i) = smooth(x0(:,i), sm_siz);
    end
end