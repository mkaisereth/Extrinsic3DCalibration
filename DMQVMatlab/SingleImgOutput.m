classdef SingleImgOutput
    % Class for Single image evaluation output (results)
    
    properties
        % Evaluation mode, e.g. cboard, cylinder, ...
        Mode
        % The version
        Version
        % Mean of absolute binary data depth image
        AbsoluteMean
        % Std of absolute binary data depth image
        AbsoluteStd
        % Mean of target absolute binary data depth image
        AbsoluteTargetMean
        % Std of target absolute binary data depth image
        AbsoluteTargetStd
        % Mean of relative difference (actual - target) of binary data
        % depth image
        RelativeMean
        % Std of relative difference (actual - target) of binary data depth
        % image
        RelativeStd
        % Mean of best fit object difference (actual - target)
        BestFitDiffMean
        % Std of best fit object difference (actual - target)
        BestFitDiffStd
        % Mean of target object difference (actual - target)
        TargetDiffMean
        % Std of target object difference (actual - target)
        TargetDiffStd
        % The unit for depth image
        Unit
        % Root Mean square Error for best fit object difference for point cloud
        BestFitPCDiffRMSE
        % Root Mean square Error for registration for point cloud
        RegistrationPCDiffRMSE
        % Point cloud density (region) in pts/mm2
        PtDensityMM
        % Mean difference (region) after registration in mm
        ZMeanMM
        % Std (region) after registration in mm
        ZStdMM
        % some description
        Description
        % whether we used transform from metadata or cpd/icp
        usedMetaDataTransform
        % camera pair combination (optional)
        Combination
    end
end

