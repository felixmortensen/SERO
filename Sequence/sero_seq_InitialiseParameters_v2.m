function params = sero_seq_InitialiseParameters_v2(varargin)
    % pulseq_InitialiseParameters: Initialize and pack pulse sequence parameters
    % with default values that can be overridden by user input.
    % Example usage: pulseq_InitialiseParameters('FOV', 50e-3)
    
    % Set default parameter values
    defaults = struct(...
        'ts', 0.15, ...          % Time per shot [s]
        'nShot', 1000, ...          % Number of shots
        'nPos', 50, ...             % Number of positions/slices
        'ar', 4, ...                % Aspect ratio
        'dz', 1.5e-3, ...           % Step length in z [m]
        'FOV', 180e-3, ...          % In-plane FOV [m]
        'Nx', 120, ...              % Image matrix size in x
        'Ny', 120, ...              % Image matrix size in y (same as Nx by default)
        'SliceWidth', 6e-3, ...     % Slice thickness [m]
        'sliceGap', 0, ...          % Slice gap factor
        'Max_b_value', 2000, ...    % b-value [s/mm^2]
        'TE', 80e-3, ...            % Echo time [s]
        'SpoilerFactor', 20, ...    % Spoiler strength
        'pe_enable', 1, ...         % Phase encoding enable flag
        'ro_os', 1, ...             % Frequency oversampling factor
        'ReadoutTime', 6.0e-4, ...  % Readout time per line [s]
        'PFFactor', 0.5, ...        % Partial Fourier factor
        'rf_faFactor', 1, ...       % RF flip angle factor
        'do_fatsat', 0, ...         % Fat saturation flag
        'do_rfSpoil', 0 ...         % RF spoiling flag
    );
    
    % Parse input arguments
    p = inputParser;
    paramNames = fieldnames(defaults);
    for i = 1:length(paramNames)
        addParameter(p, paramNames{i}, defaults.(paramNames{i}));
    end
    parse(p, varargin{:});
    
    % Pack parameters into a structure
    params = struct();
    for i = 1:length(paramNames)
        params.(paramNames{i}) = p.Results.(paramNames{i});
    end

    % Ensure consistency between related parameters
    params.Ny = params.Nx;  % Default Ny = Nx unless overridden

    % Return the packed parameters structure
end
