function sero_seq_gui_recon
% Create the main figure
fig = uifigure('Position', [100, 100, 400, 300], 'Name', 'Recon of TWIX data');

% Persistent variable to store last used folder
persistent lastUsedFolder;
if isempty(lastUsedFolder)
    lastUsedFolder = pwd; % Default to the current working directory
end

% Data file components
lbl1 = uilabel(fig, 'Position', [30, 250, 80, 22], 'Text', 'Data:');
txt1 = uieditfield(fig, 'text', 'Position', [100, 250, 200, 22]);
btn1 = uibutton(fig, 'push', 'Text', 'Browse', 'Position', [310, 250, 70, 22],...
    'ButtonPushedFcn', @(btn, event) selectFile(txt1, '*.dat'));

% Sequence file components
lbl2 = uilabel(fig, 'Position', [30, 220, 80, 22], 'Text', 'Sequence:');
txt2 = uieditfield(fig, 'text', 'Position', [100, 220, 200, 22]);
btn2 = uibutton(fig, 'push', 'Text', 'Browse', 'Position', [310, 220, 70, 22],...
    'ButtonPushedFcn', @(btn, event) selectFile(txt2, '*.seq'));

% XPS file components
lbl3 = uilabel(fig, 'Position', [30, 190, 80, 22], 'Text', 'XPS (opt):');
txt3 = uieditfield(fig, 'text', 'Position', [100, 190, 200, 22], 'HorizontalAlignment', 'right');
btn3 = uibutton(fig, 'push', 'Text', 'Browse', 'Position', [310, 190, 70, 22],...
    'ButtonPushedFcn', @(btn, event) selectFile(txt3, '*xps.mat'));

% Radio button group
bg  = uibuttongroup(fig, 'Position', [100, 80, 200, 100], 'Title', 'Recon type');
rb1 = uiradiobutton(bg,  'Position', [10, 50, 120, 22],   'Text', 'AdaptiveCombine');
rb2 = uiradiobutton(bg,  'Position', [10, 30, 120, 22],   'Text', 'SumOfSquares');
rb3 = uiradiobutton(bg,  'Position', [10, 10, 120, 22],   'Text', 'NoCombine');

% Run Recon button
btnRun = uibutton(fig, 'push', 'Text', 'Run Recon', 'Position', [150, 30, 100, 30],...
    'ButtonPushedFcn', @(btn, event) runRecon(txt1.Value, txt2.Value, bg, txt3.Value));

% Function to select files
    function selectFile(editField, filter)
        % Open file selection dialog
        [file, path] = uigetfile({filter, 'All Files'}, 'Select a File', lastUsedFolder);
        if file ~= 0  % If a file was selected
            editField.Value = fullfile(path, file);
            lastUsedFolder = path;
        end
    end

% Function to run "recon"
    function runRecon(fn_data, fn_seq, buttonGroup, fn_xps)
        % Get selected radio button
        selectedOption = buttonGroup.SelectedObject.Text;

        fn_nii = pulseq_fn_data2nii(fn_data);

        if ~isempty(fn_xps)
            copyfile(fn_xps, mdm_fn_nii2xps(fn_nii));
        end

        % Validate inputs before calling recon function
        if isempty(fn_data) || isempty(fn_seq)
            uialert(fig, 'Please select both files before running.', 'Error');
        else
            switch selectedOption
                case 'AdaptiveCombine'
                    pulseq_data2nii_adaptiveCombine(fn_data, fn_seq, fn_nii, 1);

                case 'SumOfSquares'
                    pulseq_data2nii_sumOfSquares(fn_data, fn_seq, fn_nii, 1);

                case 'NoCombine'
                    pulseq_data2nii_multicoil(fn_data, fn_seq, fn_nii, 1);
            end
        end
    end
end
