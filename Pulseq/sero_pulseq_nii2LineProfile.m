function sero_pulseq_nii2LineProfile(niftiFile, sliceNum)

% Function allowing the user to choose a column from a nii image and
% plot the line profile of that column.

try
    %Read the NIfTI image with double precision
    nii = mdm_nii_read(niftiFile);
    
    %Display the nii image
    if nargin<2
        %get center slice if no slice is manually selected
        sliceNum = round(size(nii, 3)/2);
        imshow(squeeze(nii(:, :, sliceNum)), []);
    else
        imshow(squeeze(nii(:, :, sliceNum)), []);
    end

    %Allow the user to select a column
    disp('Click on a column to plot its line profile.');
    [x, ~] = ginput(1);

    %Round the x-coordinate to the nearest integer
    x = round(x);

    %Extract the selected column from the selected x-coordinate and slice
    selectedColumn = squeeze(nii(:, x, sliceNum));

    %Plot the line profile
    figure;
    plot(abs(selectedColumn));
    title(['Line Profile of Column ' num2str(x) ' of Slice ' num2str(sliceNum)]);
    xlabel(['Pixel Number of Column ' num2str(x)]);
    ylabel('Intensity');
    grid on;

catch
    %My nifti files never work :))))))
    disp('Error reading or processing the NIfTI file.');
end

end