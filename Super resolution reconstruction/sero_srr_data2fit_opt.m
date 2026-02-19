function opt = sero_srr_data2fit_opt(opt)

if nargin < 1
    opt.initialized = 1;
end

opt.lsq = optimoptions('lsqcurvefit');
opt.lsq.Display = 'off';

% Fitting bounds
opt.lb = [0    0.3    .3    0 ];
opt.ub = [4    4    5     2 ];

% Guess bounds
opt.lg = [0.8  0.7    1   0.0 ];
opt.ug = [1.2  2    3   0.1 ];

opt.lambda.prost = [1 3 3 4];
opt.lambda.brain = [1 1 2 4];