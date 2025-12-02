function sero_srr_simulation_1D_v2(snrRange, nObj, nSims, ar, reg_str, mode, checkpointFileInput)
    % sero_srr_simulation_1D_v2 runs simulations over a range of SNR values, objects, and inner simulations.
    %
    % Inputs:
    %   snrRange    - Array of SNR values to simulate.
    %   nObj        - Number of simulated objects.
    %   nSims       - Number of noise realisations per object.
    %   ar          - Aspect ratio parameter for the sampling scheme.
    %   reg_str     - Regularization strength for the fitting process.
    %   mode        - Determines the cost function: modes 1 L2 norm,
    %                 otherwise 'None'.
    %   checkpointFileInput - (Optional) Filename of a previously saved checkpoint.
    %                         If provided, the simulation will resume from that file.
    %                         If not provided, any existing default checkpoint file is deleted
    %                         and a new simulation starts from scratch.
    
    %% Determine checkpoint file and resume indices
    if nargin < 7 || isempty(checkpointFileInput)
        checkpointFile = 'temp_simulation_checkpoint.mat';
        % If starting from scratch, delete any existing default checkpoint file.
        if exist(checkpointFile, 'file')
            delete(checkpointFile);
        end
        resumeSNR = 1;
        resumeObject = 1;
    else
        checkpointFile = checkpointFileInput;
        if exist(checkpointFile, 'file')
            load(checkpointFile, 'checkpointData');
            % Retrieve resume indices from the loaded checkpoint file, if available.
            if isfield(checkpointData, 'lastSNRIndex') && isfield(checkpointData, 'lastObject')
                resumeSNR = checkpointData.lastSNRIndex;
                resumeObject = checkpointData.lastObject + 1;
                if resumeObject > nObj
                    resumeSNR = resumeSNR + 1;
                    resumeObject = 1;
                end
            else
                resumeSNR = 1;
                resumeObject = 1;
            end
        else
            resumeSNR = 1;
            resumeObject = 1;
        end
    end

    %% Select cost function based on mode using switch-case
    switch mode
        case 1
            costFunctionType = 'L2norm';
        otherwise
            costFunctionType = 'None';
            mode = 0;
    end
    
    %% Setup problem parameters and generate object structure
    numPositions = 50;   % Number of positions
    numSamples   = 1000; % Total number of samples
    timeStep     = 0.15; % Time step
    
    [S, D, T1, V] = sero_srr_generate_structure(numPositions, nObj);
    S = S ./ mean(S, 1);  % Normalize object intensity
    
    %% Create sampling schemes for each method
    [samplingDirect.t, samplingDirect.tr, samplingDirect.b, samplingDirect.w] = sero_srr_sampling_direct(timeStep, numSamples, numPositions);
    [samplingSlice.t, samplingSlice.tr, samplingSlice.b, samplingSlice.w] = sero_srr_sampling_conventional_slice_shift(timeStep, numSamples, ar, numPositions);
    [samplingSero.t, samplingSero.tr, samplingSero.b, samplingSero.w] = sero_srr_sampling_v3(timeStep, numSamples, ar, numPositions);
    
    TR_direct  = nanmedian(samplingDirect.tr(:));
    TR_slice   = nanmedian(samplingSlice.tr(:));
    
    samplingDirect.tr(isnan(samplingDirect.tr)) = 0;
    samplingSlice.tr(isnan(samplingSlice.tr))   = 0;
    samplingSero.tr(isnan(samplingSero.tr))     = 0;
    
    %% Preallocate or load previous Bias and Precision results
    if exist('checkpointData','var')
        biasDirect      = checkpointData.biasDirect;
        biasSlice       = checkpointData.biasSlice;
        biasSero        = checkpointData.biasSero;
        precisionDirect = checkpointData.precisionDirect;
        precisionSlice  = checkpointData.precisionSlice;
        precisionSero   = checkpointData.precisionSero;
    else
        biasDirect      = zeros(4, numel(snrRange), nObj);
        biasSlice       = zeros(4, numel(snrRange), nObj);
        biasSero        = zeros(4, numel(snrRange), nObj);
        precisionDirect = zeros(4, numel(snrRange), nObj);
        precisionSlice  = zeros(4, numel(snrRange), nObj);
        precisionSero   = zeros(4, numel(snrRange), nObj);
    end
    
    %% Setup initial checkpointing structure (if not loaded)
    if ~exist('checkpointData','var')
        checkpointData = struct();
        checkpointData.biasDirect      = biasDirect;
        checkpointData.biasSlice       = biasSlice;
        checkpointData.biasSero        = biasSero;
        checkpointData.precisionDirect = precisionDirect;
        checkpointData.precisionSlice  = precisionSlice;
        checkpointData.precisionSero   = precisionSero;
        checkpointData.snr             = snrRange;
        checkpointData.reg_str         = reg_str;
        checkpointData.ar              = ar;
        checkpointData.mode            = mode;
    end
    
    %% Simulation Loop
    tic;
    for snrIndex = resumeSNR : numel(snrRange)
        currentSNR = snrRange(snrIndex);
        % For the first resumed SNR level, start at resumeObject; otherwise start at object 1.
        if snrIndex == resumeSNR
            startObject = resumeObject;
        else
            startObject = 1;
        end
        
        for objIndex = startObject : nObj
            try
                % Preallocate temporary storage for inner-loop results for the current object.
                residualsDirect   = zeros(numPositions, 4, nSims);
                predictionsDirect = zeros(numPositions, 4, nSims);
                residualsSlice    = zeros(numPositions, 4, nSims);
                predictionsSlice  = zeros(numPositions, 4, nSims);
                residualsSero     = zeros(numPositions, 4, nSims);
                predictionsSero   = zeros(numPositions, 4, nSims);
                
                % Inner simulation loop.
                for simIndex = 1:nSims
                    fprintf('SNR: %f, Object: %d/%d, Simulation: %d/%d\n', ...
                        currentSNR, objIndex, nObj, simIndex, nSims);
                    
                    % Assemble true parameter vector for the current object:
                    % [intensity, diffusivity, T1, diffusional variance]
                    trueValues = [S(:, objIndex), D(:, objIndex), T1(:, objIndex), V(:, objIndex)];
                    
                    % Generate simulated signals using the three sampling schemes.
                    dataDirect = sero_srr_fit2data(trueValues, samplingDirect.tr, samplingDirect.b, 1, samplingDirect.w);
                    dataSlice  = sero_srr_fit2data(trueValues, samplingSlice.tr, samplingSlice.b, 1, samplingSlice.w);
                    dataSero   = sero_srr_fit2data(trueValues, samplingSero.tr, samplingSero.b, 1, samplingSero.w);
                    
                    % Add noise based on the current SNR value.
                    noisyDataDirect = sero_srr_add_noise(double(dataDirect), currentSNR);
                    noisyDataSlice  = sero_srr_add_noise(double(dataSlice), currentSNR);
                    noisyDataSero   = sero_srr_add_noise(double(dataSero), currentSNR);
                    
                    % Perform fitting using the provided regularization strength and regularization mode.
                    tic;
                    fitDirect = sero_srr_data2fit_1d_reg_v2(noisyDataDirect, samplingDirect.tr, samplingDirect.b, samplingDirect.w, 0, reg_str, mode);
                    fitSlice  = sero_srr_data2fit_1d_reg_v2(noisyDataSlice, samplingSlice.tr, samplingSlice.b, samplingSlice.w, 0, reg_str, mode);
                    fitSero   = sero_srr_data2fit_1d_reg_v2(noisyDataSero, samplingSero.tr, samplingSero.b, samplingSero.w, 1, reg_str, mode);
                    toc;
                    
                    % Correct for T1 error in the first parameter.
                    fitDirect(:, 1) = fitDirect(:, 1) ./ (1 - exp(-TR_direct ./ trueValues(:, 3)));
                    fitSlice(:, 1)  = fitSlice(:, 1)  ./ (1 - exp(-TR_slice  ./ trueValues(:, 3)));
                    
                    % Calculate residuals (difference between true values and the fit).
                    residualsDirect(:, :, simIndex)  = trueValues - fitDirect;
                    predictionsDirect(:, :, simIndex) = fitDirect;
                    residualsSlice(:, :, simIndex)   = trueValues - fitSlice;
                    predictionsSlice(:, :, simIndex)  = fitSlice;
                    residualsSero(:, :, simIndex)    = trueValues - fitSero;
                    predictionsSero(:, :, simIndex)   = fitSero;
                end
                
                % Compute Bias (mean squared error) and Precision (variance) over inner simulations.
                biasDirect(:, snrIndex, objIndex)  = mean(mean(residualsDirect .^ 2, 3), 1)';
                biasSlice(:, snrIndex, objIndex)   = mean(mean(residualsSlice .^ 2, 3), 1)';
                biasSero(:, snrIndex, objIndex)    = mean(mean(residualsSero .^ 2, 3), 1)';
                precisionDirect(:, snrIndex, objIndex) = mean(var(predictionsDirect, [], 3), 1);
                precisionSlice(:, snrIndex, objIndex)  = mean(var(predictionsSlice, [], 3), 1);
                precisionSero(:, snrIndex, objIndex)   = mean(var(predictionsSero, [], 3), 1);
                
                % Update checkpoint structure with results and current indices.
                checkpointData.biasDirect(:, snrIndex, objIndex)      = biasDirect(:, snrIndex, objIndex);
                checkpointData.biasSlice(:, snrIndex, objIndex)       = biasSlice(:, snrIndex, objIndex);
                checkpointData.biasSero(:, snrIndex, objIndex)        = biasSero(:, snrIndex, objIndex);
                checkpointData.precisionDirect(:, snrIndex, objIndex) = precisionDirect(:, snrIndex, objIndex);
                checkpointData.precisionSlice(:, snrIndex, objIndex)  = precisionSlice(:, snrIndex, objIndex);
                checkpointData.precisionSero(:, snrIndex, objIndex)   = precisionSero(:, snrIndex, objIndex);
                checkpointData.lastSNRIndex = snrIndex;
                checkpointData.lastObject   = objIndex;
                
                % Save a checkpoint after finishing the current object's simulation at this SNR level.
                save(checkpointFile, 'checkpointData', '-v7.3');
                
            catch ME
                % In case of an error, save progress, display the error, then rethrow.
                checkpointData.lastSNRIndex = snrIndex;
                checkpointData.lastObject   = objIndex;
                save(checkpointFile, 'checkpointData', '-v7.3');
                fprintf('Error at SNR = %f, Object = %d: %s\n', currentSNR, objIndex, ME.message);
                rethrow(ME);
            end
        end
    end
    toc;
    
    %% Structure final results for output
    data.direct.bias      = biasDirect;
    data.direct.precision = precisionDirect;
    data.slice.bias       = biasSlice;
    data.slice.precision  = precisionSlice;
    data.sero.bias        = biasSero;
    data.sero.precision   = precisionSero;
    data.snr              = snrRange;
    
    % Convert regularization strength to scientific notation for the filename.
    regStr = sprintf('%.0e', reg_str);
    
    %% Save final results to output directory
    outputDir = '';
    cd(outputDir);
    finalFileName = sprintf('Sero_Slice_Dir_acc_prec_nObj%d_nSims%d_ar%d_regstr%s_RegType%s.mat', ...
                            nObj, nSims, ar, regStr, costFunctionType);
    save(finalFileName, 'data');
    disp(['Simulation file saved in ' outputDir ' as ' finalFileName]);
end


