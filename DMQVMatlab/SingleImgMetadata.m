classdef SingleImgMetadata
    %SINGLEIMGMETADATA Metadata for a single image evaluation
    %   Metadata read from JSON File
    %   At least Mode and SingleDepthImgPath are required
    
    properties
        % Evaluation mode (e.g. cboard, sphere, cylinder, 3dstaticback, ...)
        Mode
        % Metadata version
        Version
        % The path to the binary depth image to evaluate
        SingleDepthImgPath
        % The path to the target (e.g. ground truth) depth image
        SingleDepthTargetImgPath
        % The path to the point cloud file (e.g. pcd or ply)
        SinglePointCloudImgPath
        % The path to the target point cloud
        SinglePointCloudTargetImgPath
        % A geometric representation of a target point cloud
        SinglePointCloudTargetObject
        % The height in pixels for the depth image
        SingleDepthImgHeight
        % The width in pixels for the depth image
        SingleDepthImgWidth
        % The unit (e.g. m, mm, cm, um) for the depth image
        SingleDepthImgUnit
        % The minimum distance in SingleDepthImgUnit for foreground
        % segmentation
        SingleDepthImgForegroundMin
        % The maximum distance in SingleDepthImgUnit for foreground
        % segmentation
        SingleDepthImgForegroundMax
        % Whether to do M360 evaluation
        M360Eval
        % Optional target object to compare to
        SingleTargetObject
        % some description
        Description
        % path to 3d checkerboard calibration
        Calib3dPath
        % camera combination pair (optional)
        Combination
    end
    
    methods(Static)
        
        function obj = FromJson(jsonContent)
            % Read the metadata from a jsonContent
            
            % create a new SingleImgMetadata
            obj = SingleImgMetadata;
            % read all fields (if available), else set some defaults
            if(isfield(jsonContent, "Mode"))
                obj.Mode = jsonContent.Mode;
            else
                obj.Mode = "None";
            end
            if(isfield(jsonContent, "Version"))
                obj.Version = jsonContent.Version;
            else
                obj.Version = 0;
            end
            if(isfield(jsonContent, "SingleDepthImgPath"))
                obj.SingleDepthImgPath = jsonContent.SingleDepthImgPath;
            else
                obj.SingleDepthImgPath = "";
            end
            if(isfield(jsonContent, "SingleDepthTargetImgPath"))
                obj.SingleDepthTargetImgPath = jsonContent.SingleDepthTargetImgPath;
            else
                obj.SingleDepthTargetImgPath = "";
            end
            if(isfield(jsonContent, "SinglePointCloudImgPath"))
                obj.SinglePointCloudImgPath = jsonContent.SinglePointCloudImgPath;
            else
                obj.SinglePointCloudImgPath = "";
            end
            if(isfield(jsonContent, "SinglePointCloudTargetImgPath"))
                obj.SinglePointCloudTargetImgPath = jsonContent.SinglePointCloudTargetImgPath;
            else
                obj.SinglePointCloudTargetImgPath = "";
            end
            if(isfield(jsonContent, "SinglePointCloudTargetObject"))
                obj.SinglePointCloudTargetObject = jsonContent.SinglePointCloudTargetObject;
            else
                obj.SinglePointCloudTargetObject = [];
            end
            if(isfield(jsonContent, "SingleDepthImgHeight"))
                obj.SingleDepthImgHeight = jsonContent.SingleDepthImgHeight;
            else
                obj.SingleDepthImgHeight = 480;
            end
            if(isfield(jsonContent, "SingleDepthImgWidth"))
                obj.SingleDepthImgWidth = jsonContent.SingleDepthImgWidth;
            else
                obj.SingleDepthImgWidth = 640;
            end
            if(isfield(jsonContent, "SingleDepthImgUnit"))
                obj.SingleDepthImgUnit = jsonContent.SingleDepthImgUnit;
            else
                obj.SingleDepthImgUnit = "mm";
            end
            if(isfield(jsonContent, "SingleDepthImgForegroundMin"))
                obj.SingleDepthImgForegroundMin = jsonContent.SingleDepthImgForegroundMin;
            else
                obj.SingleDepthImgForegroundMin = 500;
            end
            if(isfield(jsonContent, "SingleDepthImgForegroundMax"))
                obj.SingleDepthImgForegroundMax = jsonContent.SingleDepthImgForegroundMax;
            else
                obj.SingleDepthImgForegroundMax = 2000;
            end
            if(isfield(jsonContent, "M360Eval"))
                obj.M360Eval = jsonContent.M360Eval;
            else
                obj.M360Eval = "";
            end
            if(isfield(jsonContent, "SingleTargetObject"))
                obj.SingleTargetObject = jsonContent.SingleTargetObject;
            else
                obj.SingleTargetObject = [];
            end
            if(isfield(jsonContent, "Description"))
                obj.Description = jsonContent.Description;
            else
                obj.Description = [];
            end
            if(isfield(jsonContent, "Calib3dPath"))
                obj.Calib3dPath = jsonContent.Calib3dPath;
            else
                obj.Calib3dPath = [];
            end
            if(isfield(jsonContent, "Combination"))
                obj.Combination = jsonContent.Combination;
            else
                obj.Combination = [];
            end
        end
    end
end

