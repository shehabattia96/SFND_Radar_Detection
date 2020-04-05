% SFND Radar Detection - Shehab Attia
% Example input: % offset the threshold by SNR value in dB
% thresholdOffset = 5; %dB
% trainingCellSize = [5, 2]; % [doppler, range], axis flipped because it makes sense for moving window
% guardCellSize = [2, 1]; % [doppler, range]
function [signal_cfar] = TwoDimensionalCFAR(RDM, trainingCellSize, guardCellSize, thresholdOffset, isDisplaySuftPlot)
    [numRangeIndecies, numDopplerIndecies] = size(RDM);
    threshold_cfar = zeros(size(RDM));
    signal_cfar = zeros(size(RDM));

    trainingAndguardCellSizeize = trainingCellSize + guardCellSize;

    startingDopplerIndex = trainingAndguardCellSizeize(1) + 1;
    endingDopplerIndex = numDopplerIndecies - trainingAndguardCellSizeize(1) -1;
    startingRangeIndex = trainingAndguardCellSizeize(2) + 1;
    endingRangeIndex = numRangeIndecies - trainingAndguardCellSizeize(2) - 1;

    %design a loop such that it slides the CUT across range doppler map by
    %giving margins at the edges for Training and Guard Cells.
    for CUTDopplerIndex = startingDopplerIndex:endingDopplerIndex
        for CUTRangeIndex = startingRangeIndex:endingRangeIndex
            gridTrainingAndGuardDopplerIndecies = CUTDopplerIndex - trainingAndguardCellSizeize(1): CUTDopplerIndex + trainingAndguardCellSizeize(1);
            gridTrainingAndGuardRangeIndecies = CUTRangeIndex - trainingAndguardCellSizeize(2): CUTRangeIndex + trainingAndguardCellSizeize(2);

            gridTrainingAndGuard = RDM(gridTrainingAndGuardRangeIndecies, gridTrainingAndGuardDopplerIndecies);
            gridTraining = gridTrainingAndGuard;


            gridGuardDopplerIndecies = startingDopplerIndex - guardCellSize(1): startingDopplerIndex + guardCellSize(1);
            gridGuardRangeIndecies = startingRangeIndex - guardCellSize(2): startingRangeIndex + guardCellSize(2);

            gridTraining(gridGuardRangeIndecies, gridGuardDopplerIndecies) = nan;

            trainingAverage = nanmean(gridTraining, 'all');
            trainingAverageWithThreshold = trainingAverage + thresholdOffset;
            threshold_cfar(CUTRangeIndex, CUTDopplerIndex) = trainingAverageWithThreshold;

            CUTValue = RDM(CUTRangeIndex, CUTDopplerIndex);

            if (CUTValue >= trainingAverageWithThreshold)
                signal_cfar(CUTRangeIndex, CUTDopplerIndex) = CUTValue;
            end
        end
    end

    if ~isDisplaySuftPlot
        return
    end
    
    threshold_cfar_without_zeroPadding = threshold_cfar(startingRangeIndex:endingRangeIndex, startingDopplerIndex:endingDopplerIndex);

    doppler_axis = linspace(-100,100,numDopplerIndecies);
    range_axis = linspace(-200,200,numRangeIndecies)*((numRangeIndecies)/400);

    figure,surf(doppler_axis,range_axis,signal_cfar); hold on
    surf(doppler_axis(startingDopplerIndex:endingDopplerIndex),range_axis(startingRangeIndex:endingRangeIndex),threshold_cfar_without_zeroPadding, 'EdgeColor', 'none');
    colorbar;

    title('Signal after 2D-CFAR filtering')
    xlabel("Doppler (m/s)")
    ylabel("Range (m)")
    zlabel("Intensity (dBW)")
end

