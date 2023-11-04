close all; clear all; clc;

% script to merge sagittal bending captures from two photoneos
global printLevel;
printLevel = 1;

doUseIcp = true;
doUseUserInput = false;
letUserSelectMarkers = false;
directVersion = false;
withTexture = false;

% basefolder where photoneo_*.ply are found
baseFolder = "I:\Temp\3DcbData\Study2Data\Dynamic\08";
cb3dBaseFolder = "I:\Temp\3DcbData\Study2Data\Dynamic\50";
genericcb3dFolderPaths = [];

patientNums = 1:18;

%%
[rmseVals1s, stdVals1s] = Apply3DcbDynamic_func(baseFolder, cb3dBaseFolder, genericcb3dFolderPaths, doUseIcp, doUseUserInput, letUserSelectMarkers, patientNums, directVersion, withTexture, printLevel);

%% some checks on rmse and std
rmseValses = [];
stdValses = [];
for patNum = patientNums
    if patNum>length(rmseVals1s) || isempty(rmseVals1s(patNum))
        disp("Warning: data for patient missing: " + string(patNum));
        continue;
    end
    disp("Max of patient " + string(patNum) + " is: " + string(1000*max(rmseVals1s{patNum})) + " +- " + string(1000*max(stdVals1s{patNum})))
    rmseValses = [rmseValses; rmseVals1s{patNum}];
    stdValses = [stdValses; stdVals1s{patNum}];
end

%% plot eval with Boxplots
if ~exist("rmseVals1s")
load('C:\Users\ETH\OneDrive - ETH Zurich\Papers\BFH\3dCheckerboard\Resources\Balgrist_18samples (2)\rmseValses_Phneo2-Phneo.mat')
load('C:\Users\ETH\OneDrive - ETH Zurich\Papers\BFH\3dCheckerboard\Resources\Balgrist_18samples (2)\stdValses_Phneo2-Phneo.mat')
if doUseIcp
load('C:\Users\ETH\OneDrive - ETH Zurich\Papers\BFH\3dCheckerboard\Resources\Balgrist_18samples (2)\rmseValses_Phneo2-Phneo_ICP.mat')
load('C:\Users\ETH\OneDrive - ETH Zurich\Papers\BFH\3dCheckerboard\Resources\Balgrist_18samples (2)\stdValses_Phneo2-Phneo_ICP.mat')
end
end
outputPath = append(baseFolder, "\Output");
if ~exist(outputPath, "dir")
    mkdir(outputPath)
end
combination = "Photoneo2-Photoneo";
if doUseIcp
    combination = "Photoneo2-Photoneo_ICP";
end

f2 = figure;
boxplot(rmseValses*1000, 'Labels', combination);
title("RMSE [mm]")
ylabel("RMSE [mm]")
saveas(f2, append(outputPath, "\ZMeanDiff_", combination, ".fig"));
ax = gca;
exportgraphics(ax,append(outputPath, "\ZMeanDiff_", combination, ".png"),"Resolution",600)
f3 = figure;
boxplot(stdValses*1000, 'Labels', combination);
title("SD [mm]")
ylabel("SD [mm]")
saveas(f3, append(outputPath, "\ZStd_", combination, ".fig"));
ax = gca;
exportgraphics(ax,append(outputPath, "\ZStd_", combination, ".png"),"Resolution",600)
% plot together
f4 = figure;
boxplot([rmseValses*1000, stdValses*1000], 'Labels', ["RMSE", "SD"]);
ylim([0,19])
title("MotionCam 2 (above) - MotionCam 1 (behind)")
ylabel("[mm]")
saveas(f4, append(outputPath, "\ZMeanStd_", combination, ".fig"));
ax = gca;
exportgraphics(ax,append(outputPath, "\ZMeanStd_", combination, ".png"),"Resolution",600)