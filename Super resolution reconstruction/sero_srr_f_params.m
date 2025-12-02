function param_diff_w = sero_srr_f_params(param_diff, num_edge, lower_val, rampType)
%param_diff_w = sero_srr_f_params(param_diff, num_edge, lower_val, rampType)

    if nargin < 3 || isempty(lower_val)
        lower_val = 1e-3;
    end
    if nargin < 4 || isempty(rampType)
        rampType = 'linear';
    end

    N = size(param_diff, 1);
    
    weights = ones(N, 1);
    
    switch lower(rampType)
        case 'linear'
            %Linearly interpolate from lower_val to 1 for the top edge.
            weights(1:num_edge) = linspace(lower_val, 1, num_edge);
            weights(end-num_edge+1:end) = linspace(1, lower_val, num_edge);

        case 'exponential'
            %Exponentially ramp from lower_val to 1.
            ramp = lower_val * ((1/ lower_val) .^ ((0:num_edge-1) / (num_edge-1)));
            weights(1:num_edge) = ramp;
            weights(end-num_edge+1:end) = fliplr(ramp);
        otherwise
            error('Unknown rampType. Use "linear" or "exponential".');
    end
    
    param_diff_w = param_diff .* weights;
end