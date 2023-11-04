function DoDMQV(pc1Path, pc1FileName, pc2_refPath, tFormPath, comboName, basePath)

addpath("DMQVMatlab");
addpath('Utils');
json1FileName = replace(pc1FileName, ".ply", ".json");
%json1FileName = append("*", json1FileName);
jsonFilePath = replace(pc1Path, pc1FileName, json1FileName);
jsonFilePath2 = dir(jsonFilePath);
if isempty(jsonFilePath2)
    % create default json if missing
    jsonObject.SingleDepthImgHeight = 0;
    jsonObject.SingleDepthImgWidth = 0;
    jsonObject.SingleDepthImgPath = "";
    jsonObject.SingleDepthTargetImgPath = "";
    jsonObject.Version = 1;
    jsonObject.Description = "";
else
    jsonFilePath = append(jsonFilePath2(1).folder, "\", jsonFilePath2(1).name);
    jsonString = fileread(jsonFilePath);
    jsonObject = jsondecode(jsonString);
end

% adjust some parameters and save back
jsonObject.Mode = "3dStaticBack";
jsonObject.SingleDepthImgForegroundMin = [-1,-1,0.1];
jsonObject.SingleDepthImgForegroundMax = [1,1,2];
jsonObject.SingleDepthImgUnit = "m";
jsonObject.SinglePointCloudImgPath = replace(pc1Path, "\", "/");
jsonObject.SinglePointCloudTargetImgPath = replace(pc2_refPath, "\", "/");
jsonObject.Calib3dPath = replace(tFormPath, "\", "/");
comboNameString = convertStringsToChars(comboName);
%%% TODO this belongs to the other json from DMQV
jsonObject.Combination = comboNameString(1:end-9);
jsonString = jsonencode(jsonObject, "PrettyPrint", true);

fid = fopen(jsonFilePath, "wt");
fprintf(fid, jsonString);
fclose(fid);

showPcRoi = 0;
basePathTemp = basePath;
Run3DStaticBack;
basePath = basePathTemp;

end