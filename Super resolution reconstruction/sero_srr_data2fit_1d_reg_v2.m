function m = sero_srr_data2fit_1d_reg_v2(S, tr, b, w, do_t1, do_v, reg_str, mode, x0, opt, wght, te)
    % function m = sero_srr_data2fit_1d(S, tr, b, w, do_t1, opt, x0, wght, te)
    
    if nargin < 7 || isempty(reg_str)
        % Define the regularization parameter
        reg_str=0;
    end

    if nargin < 8 || isempty(mode)
        %choose regularisation mode
        mode = 1;
    end

    if nargin < 10 || isempty(opt)
        opt = sero_srr_data2fit_opt();
    end

    if nargin < 9 || isempty(x0)
        do_guess = 1;
    else
        do_guess = 0;
    end

    if nargin < 11 || isempty(wght)
        wght = sero_srr_weight_f(tr);
        wghtInput = wght;
        wghtInput(end + 1) = 1;
    end

    if nargin < 12 || isempty(te)
        te = [];
    end

    z = size(tr, 2);
    
    ar = sum(ceil(w(1,:)));
    reg_w = opt.lambda.brain;

    lb = ones(z, 1) * opt.lb;
    ub = ones(z, 1) * opt.ub;
    lg = ones(z, 1) * opt.lg;
    ug = ones(z, 1) * opt.ug;

    if do_guess
        x0 = rand(z, 4) .* (ug - lg) + lg;
    end

    nrm = max(S);
    if mode == 1
        input = [S/nrm; zeros(4 * (z - 1), 1)];
        wghtInput = ones(size(input));
        wghtInput(1:size(wght, 1)) = wght;
    else
        input = [S/nrm; 0];
    end

    if ~do_t1
        x0(:, 3) = 0;
        lb(:, 3) = 0;
        ub(:, 3) = 0;
    end

    if ~do_v
        x0(:, 4) = 0;
        lb(:, 4) = 0;
        ub(:, 4) = 0;
    end

    m = lsqcurvefit(@par2signal, x0, [], input.*wghtInput, lb, ub, opt.lsq);

    m(:, 1) = m(:, 1) * nrm;

    if ~do_t1
        m(:, 3) = nan;
    end

    if ~do_v
        m(:, 4) = nan;
    end

 
    function S = par2signal(p, varargin)
        % Compute the base model signal
        S = sero_srr_fit2data(p, tr, b, do_t1, do_v, w, te) .* wght;

        switch mode
            case 0
                regularization_term = 0; %do nothing
    
            case 1
                %Smoothing regularisation basd on neighbouring differences
                param_diff = diff(p, 1, 1);
                param_diff = sero_srr_f_params(param_diff, ar); %default is linear ramp filter
                weighted_diff = param_diff .* repmat(reg_w, size(param_diff, 1), 1);
                regularization_term = reg_str*weighted_diff(:);
    
            otherwise
                error('Unknown regularization mode.');
        end
    
        S = [S; regularization_term];
    end

end
