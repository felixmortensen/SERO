function WGHT = sero_srr_mask2wght(M, TR, mode)

% This function needs some love since it behaves in unpredictable ways

if nargin < 3
    mode = 2;
end

WGHT = zeros(size(M,1), size(M,2), size(TR,1));

for i = 1:size(M,1)
    for j = 1:size(M,2)
        m = squeeze(M(i,j,:));

        if ~any(m)
            continue
        end

        switch mode
            case 0
                % do nothing

            case 1
                m = smooth(m,4);

            case 2
                m = m*0+1;
                WGHT(i,j,:) = m;
                return

        end

        WGHT(i,j,:) = (TR>0)*m./sum(TR>0,2)*0.9 + 0.1;
    end
end