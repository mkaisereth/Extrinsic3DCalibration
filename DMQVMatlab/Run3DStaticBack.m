% Script to run 3dStaticBack evaluation
% Requires JSON metadata given by 'jsonFilePath'

% get timestamp
ts = RunObject.GetTimestamp();

% initialization
RunObject.Init();

% Read metadata from json file
if(~exist('jsonFilePath', 'var') || isempty(jsonFilePath))
    jsonFilePath = 'Test_3dstaticback.json';
end

% log file path
[jsonFolder, jsonFile, ~] = fileparts(jsonFilePath);

% diary
if(isempty(jsonFolder))
    outputFolder = "Output";
else
    outputFolder = sprintf("%s\\Output", jsonFolder);
end
diaryFileName = sprintf("%s\\%s_log_%s.txt", outputFolder, jsonFile, ts);
if(~exist(outputFolder, 'dir'))
    mkdir(outputFolder);
end
diary(diaryFileName);
diary on;
disp('Run3DStaticBack ...');

metaData = RunObject.ReadJson(jsonFilePath);
% clear the path (will not be cleared during init)
clear jsonFilePath;

if ~isempty(metaData.SingleDepthImgPath)
% read the absolute data (depth image)
absBinData = RunObject.ReadAbsBinaryData(metaData);
[msData.absFgMean, msData.absFgStd] = RunObject.CalcAbsFgMeanAndStd(absBinData, metaData);

% show figure for absolute data
[f1, absImg] = RunObject.ShowFigure(absBinData, '3DStaticBack absolute depth image', metaData);
end
% TODO this is just to generate some testing target data
%binDataArr = ReadBinary.ReadBinaryInt32(metaData.SingleDepthImgPath);
%binDataArr = SegmentForeground.SegmentMinMax(binDataArr, metaData.SingleDepthImgForegroundMin, metaData.SingleDepthImgForegroundMax);
%binDataArr(binDataArr > 0) = 900; %TODO
%ReadBinary.WriteBinaryInt32(binDataArr,'TestOut_3dStaticBack.bin');

if ~isempty(metaData.SingleDepthTargetImgPath)
% read the relative data (depth image data - target depth image data)
[absTargetBinData, relDiffBinData] = RunObject.ReadAbsTargetBindaryData(metaData, absBinData);
[msData.absTargetFgMean, msData.absTargetFgStd] = RunObject.CalcAbsTargetFgMeanAndStd(absTargetBinData, metaData);
[msData.relFgMean, msData.relFgStd] = RunObject.CalcRelFgMeanAndStd(relDiffBinData, metaData);

% show figure for relative data
[f2, relImg] = RunObject.ShowFigure(relDiffBinData, '3DStaticBack relative depth image (actual - target)', metaData);
end

if ~isempty(metaData.SingleDepthImgPath)
% best fitting plane
[f3, bestFitDiffData, bestFitDiffImg] = RunObject.BestFitObject(absBinData, 'poly11', metaData, '3DStaticBack depth image with best fit plane');
[msData.bestFitDiffFgMean,msData.bestFitDiffFgStd] = RunObject.CalcBestFitDiffFgMeanAndStd(bestFitDiffData, metaData);
end

% best fitting sphere for point cloud
pcfitplane2 = @(ptCloud) pcfitplane(ptCloud, 0.02);
[ptCloud, mFactor] = RunObject.ReadAbsoluteFgPointCloud(metaData.SinglePointCloudImgPath, metaData, showPcRoi);

% try to align point cloud to axes for best fitting
ptCloudAligned = ptCloud; %RunObject.AlignPCWithYCoordSym(ptCloud, mFactor);

%[f5, ptCloudSmooth, planeModel, rmse] = RunObject.PCBestFitObject(ptCloudAligned, pcfitplane2, 10);
f5 = [];
rmse = [];
msData.bestFitPCDiffRMSE = rmse;

%try to estimate width and height of point cloud
RunObject.GetPcSize(ptCloud, mFactor);

% pointcloud registration with target point cloud
tPtCloud = RunObject.ReadTargetFgPointCloud(metaData.SinglePointCloudTargetImgPath, metaData, showPcRoi);
% registration
[ptCloudTrans, tForm, rmse, f6, f7, usedMetaDataTransform] = RunObject.RegisterPointCloud(ptCloud, tPtCloud, [], metaData);
msData.registrationPCDiffRMSE = rmse;
msData.usedMetaDataTransform = usedMetaDataTransform;

% do another evaluation of smaller region (to remove borderline issues)
[msData.ptDensityMM, msData.zMeanMM, msData.zStdMM, f8] = RunObject.EvaluateRegion(ptCloud, tForm, tPtCloud);

% save image & json output
basePath = outputFolder+"\\"+ jsonFile;

if ~isempty(metaData.SingleDepthImgPath)
    saveas(absImg, basePath + "_abs_"+ts+".png");
end
if ~isempty(metaData.SingleDepthTargetImgPath)
    saveas(relImg, basePath + "_rel_"+ts+".png");
end
if ~isempty(metaData.SingleDepthImgPath)
    if ~isempty(f3)
        saveas(f3, basePath + "_bestFit_"+ts+".png");
        savefig(f3, basePath + "_bestFit_"+ts+".fig");
    end
    saveas(bestFitDiffImg, basePath + "_bestFitDiff_"+ts+".png");
end

if ~isempty(f5)
    saveas(f5, basePath + "_pcBestFit_"+ts+".png");
    savefig(f5, basePath + "_pcBestFit_"+ts+".fig");
end
if ~isempty(f6)
    savefig(f6, basePath + "_pcDiff_"+ts+".fig");
end
if ~isempty(f7)
    savefig(f7, basePath + "_pcRegistration_"+ts+".fig");
end
if ~isempty(f8)
    savefig(f8, basePath + "_pcRegionEval_"+ts+".fig");
end

siOutput = RunObject.PrepareSiOutput(metaData, msData, mFactor);
RunObject.SaveSiOutput(siOutput, basePath, ts);

diary off;