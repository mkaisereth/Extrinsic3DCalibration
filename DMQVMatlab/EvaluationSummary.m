classdef EvaluationSummary
   
    methods (Static)
        
        function h = DoJson(jsonBaseFolderPath, jsonFilePath)
            % only evaluate results (output)
            if(contains(upper(jsonFilePath), "\OUTPUT\"))
                
                jsonFilePath = strrep(jsonFilePath, "\\", "\");
                jsonFilePath = strrep(jsonFilePath, "\", "\\");
                fprintf(jsonFilePath);
                addpath('..\Utils');
                if(~exist('ReadJson'))
                    msgbox('You need to clone utils: https://gitlab.ti.bfh.ch/4d-spine/utils on relative path ..\Utils');
                end
                jsonContent = ReadJson.Read(jsonFilePath);
                
                obj = [];
                % Mode
                if(isfield(jsonContent, "Mode"))
                    obj.Mode = jsonContent.Mode;
                else
                    obj.Mode = "None";
                end
                % object
                pathFolders = split(jsonFilePath, "\\");
                folderInd = -1;
                for i=1:length(pathFolders)
                   if(contains(upper(pathFolders(i)), "SHUTTERS"))
                       folderInd = i;
                       break;
                   end
                end
                % if not found, search for folder with mode
                for i=1:length(pathFolders)
                   if(contains(upper(pathFolders(i)), upper(jsonContent.Mode)))
                       folderInd = i;
                       break;
                   end
                end
                underlineInd = strfind(pathFolders(folderInd), "_");
                objectName = extractBetween(pathFolders(folderInd), underlineInd(end)+1, strlength(pathFolders(folderInd)));
                obj.ObjectName = objectName;
                
                % PC registration RMSE
                if(isfield(jsonContent, "RegistrationPCDiffRMSE"))
                    obj.RegistrationPCDiffRMSE = jsonContent.RegistrationPCDiffRMSE;
                else
                    obj.RegistrationPCDiffRMSE = [];
                end
                % Shutters
                if(contains(jsonFilePath, "OpenShutters"))
                    obj.Shutters = "Open";
                elseif(contains(jsonFilePath, "ClosedShutters"))
                    obj.Shutters = "Closed";
                else
                    obj.Shutters = [];
                end
                
                % distance
                cmInd = strfind(jsonFilePath, "cm_");
                cmUnderlineInd = -1;
                if ~isempty(cmInd)
                    cmInd = cmInd(1);
                    underlineInd = strfind(jsonFilePath, "_");
                    for i=1:length(underlineInd)
                        if(underlineInd(i)<cmInd)
                            cmUnderlineInd = underlineInd(i);
                        end
                    end
                end
                if cmUnderlineInd > -1
                    distance = extractBetween(jsonFilePath, cmUnderlineInd+1, cmInd+1);
                    obj.Distance = distance;
                else
                    obj.Distance = [];
                end
                % resolution
                if(contains(jsonFilePath, "maxRes"))
                    obj.Resolution = "maxRes";
                elseif(contains(jsonFilePath, "640x480"))
                    obj.Resolution = "640x480";
                else
                    underlineInd = strfind(pathFolders(folderInd), "_");
                    resolution = extractBetween(pathFolders(folderInd), underlineInd(2)+1, underlineInd(3)-1);
                    if contains(resolution, "x")
                        obj.Resolution = resolution;
                    else
                        obj.Resolution = [];
                    end
                end
                
                % get used combination
                pathFolders = split(jsonFilePath, "\\");
                
%                 if(pathFolders(7) == "Output")
%                     obj.Combination = "Single";
%                 elseif(pathFolders(7) == "Single")
%                     obj.Combination = "Single";
%                 else
%                     obj.Combination = pathFolders(7);
%                 end
                combinationFolderName = pathFolders(length(pathFolders)-2);
                cFNunderlineInd = strfind(combinationFolderName, "_");
                if(length(cFNunderlineInd) > 0)
                    combinationFolderName = extractBetween(combinationFolderName, 1, cFNunderlineInd(2)-1); %% TODO think about this
                end

                obj.Combination = combinationFolderName;
                
                % get camera name
                [~,fileName,~] = fileparts(jsonFilePath);
                underlineInd = strfind(fileName, "_");
                cameraName = extractBetween(fileName, underlineInd(1)+1, underlineInd(3)-1); % TODO
                cameraInd = extractBetween(fileName, underlineInd(3), underlineInd(3)+1);
                obj.CameraName = cameraName + cameraInd;
                
                % get point density
                if(isfield(jsonContent, "PtDensityMM"))
                    obj.PtDensityMM = jsonContent.PtDensityMM;
                else
                    obj.PtDensityMM = [];
                end
                % get mean difference in z (actual - target) after
                % registration
                if(isfield(jsonContent, "ZMeanMM"))
                    obj.ZMeanMM = jsonContent.ZMeanMM;
                else
                    obj.ZMeanMM = [];
                end
                % get std difference in z (actual - target) after
                % registration
                if(isfield(jsonContent, "ZStdMM"))
                    obj.ZStdMM = jsonContent.ZStdMM;
                else
                    obj.ZStdMM = [];
                end
                
                if(isfield(jsonContent, "Description"))
                    obj.Description = jsonContent.Description;
                else
                    obj.Description = [];
                end
                
                if(isfield(jsonContent, "usedMetaDataTransform"))
                    obj.usedMetaDataTransform = jsonContent.usedMetaDataTransform;
                else
                    obj.usedMetaDataTransform = [];
                end
                
                disp(' ')
                disp(obj.Mode);
                disp(obj.RegistrationPCDiffRMSE);
                disp(obj.Shutters);
                disp(obj.Distance);
                disp(obj.Description)
                
                outputFileName = jsonBaseFolderPath + "\\EvaluationSummary.csv";
                
                fid = fopen(outputFileName, 'a+');
                fprintf(fid, "%s,%s,%0.2f,%s,%s,%s,%s,%s,%0.2f,%0.2f,%0.2f,%s,%d\n",...
                    obj.Mode, obj.ObjectName, obj.RegistrationPCDiffRMSE, obj.Shutters, obj.Distance, obj.Resolution, obj.Combination, obj.CameraName, obj.PtDensityMM, obj.ZMeanMM, obj.ZStdMM, obj.Description, obj.usedMetaDataTransform);
                fclose(fid);
                
                h = 1;
            else
                h = 0;
            end
        end

        function h = DoJsonMinimal(jsonBaseFolderPath, jsonFilePath)
            % only evaluate results (output)
            if(contains(upper(jsonFilePath), "\OUTPUT\"))
                
                jsonFilePath = strrep(jsonFilePath, "\\", "\");
                jsonFilePath = strrep(jsonFilePath, "\", "\\");
                fprintf(jsonFilePath);
                addpath('..\Utils');
                if(~exist('ReadJson'))
                    msgbox('You need to clone utils: https://gitlab.ti.bfh.ch/4d-spine/utils on relative path ..\Utils');
                end
                jsonContent = ReadJson.Read(jsonFilePath);
                
                obj = [];
                % Mode
                if(isfield(jsonContent, "Mode"))
                    obj.Mode = jsonContent.Mode;
                else
                    obj.Mode = "None";
                end
                % object
                pathFolders = split(jsonFilePath, "\\");
                
                objectName = pathFolders(end-2);
                obj.ObjectName = objectName;
                
                % PC registration RMSE
                if(isfield(jsonContent, "RegistrationPCDiffRMSE"))
                    obj.RegistrationPCDiffRMSE = jsonContent.RegistrationPCDiffRMSE;
                else
                    obj.RegistrationPCDiffRMSE = [];
                end

                % TODO
                obj.Shutters = [];
                obj.Distance = [];
                obj.Resolution = [];
                obj.Combination = jsonContent.Combination;
                obj.CameraName = [];
                
                % get point density
                if(isfield(jsonContent, "PtDensityMM"))
                    obj.PtDensityMM = jsonContent.PtDensityMM;
                else
                    obj.PtDensityMM = [];
                end
                % get mean difference in z (actual - target) after
                % registration
                if(isfield(jsonContent, "ZMeanMM"))
                    obj.ZMeanMM = jsonContent.ZMeanMM;
                else
                    obj.ZMeanMM = [];
                end
                % get std difference in z (actual - target) after
                % registration
                if(isfield(jsonContent, "ZStdMM"))
                    obj.ZStdMM = jsonContent.ZStdMM;
                else
                    obj.ZStdMM = [];
                end
                
                if(isfield(jsonContent, "Description"))
                    obj.Description = jsonContent.Description;
                else
                    obj.Description = [];
                end
                
                if(isfield(jsonContent, "usedMetaDataTransform"))
                    obj.usedMetaDataTransform = jsonContent.usedMetaDataTransform;
                else
                    obj.usedMetaDataTransform = [];
                end
                
                disp(' ')
                disp(obj.Mode);
                disp(obj.RegistrationPCDiffRMSE);
                disp(obj.Shutters);
                disp(obj.Distance);
                disp(obj.Description)
                
                outputFileName = jsonBaseFolderPath + "\\EvaluationSummary.csv";
                
                fid = fopen(outputFileName, 'a+');
                fprintf(fid, "%s,%s,%0.2f,%s,%s,%s,%s,%s,%0.2f,%0.2f,%0.2f,%s,%d\n",...
                    obj.Mode, obj.ObjectName, obj.RegistrationPCDiffRMSE, obj.Shutters, obj.Distance, obj.Resolution, obj.Combination, obj.CameraName, obj.PtDensityMM, obj.ZMeanMM, obj.ZStdMM, obj.Description, obj.usedMetaDataTransform);
                fclose(fid);
                
                h = 1;
            else
                h = 0;
            end
        end
    end
end