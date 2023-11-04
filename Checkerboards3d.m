clear all; close all; clc;

%dataSetName = "Larsen";
dataSetName = "Balgrist";
dynamicData = true;
directVersion = false;

% only do the first two
numOfCheckerboards = inf;
numOfCheckerboards = 2;

printLevel = 1;

% Larsen
if strcmp(dataSetName, "Larsen")
    basePath = "I:\Temp\3DcbData\Study1Data";
    baseFoldercb3d = [];
end
% Balgrist
if strcmp(dataSetName, "Balgrist")
    basePath = "I:\Temp\3DcbData\Study2Data\Static";
    baseFoldercb3d = [];
    basePath = "I:\Temp\3DcbData\Study2Data\Dynamic\08";
    baseFoldercb3d = "I:\Temp\3DcbData\Study2Data\Dynamic\50";
end

% Larsen study has patiens from 5 to 75
if strcmp(dataSetName, "Larsen")
    patientNumscb3d = 5:75;
end
if strcmp(dataSetName, "Balgrist")
    patientNumscb3d = 1:22;
end

% systemCombos to evaluate
if ~dynamicData
    systemCombos1 = ["Photoneo", "Astida" "Sltida", "Astra_5", "D415_1", "Astra_6", "D415_2", "Astra_6", "D415_2"];
    systemCombos2 = ["Photoneo", "Photoneo", "Photoneo", "Photoneo", "Photoneo", "Photoneo", "Photoneo", "Astra_5", "D415_1"];
    systemNames = ["Photoneo", "Astida", "Sltida", "Astra", "Astra", "D415", "D415"];
    systemNames2 = ["", "", "", "_5", "_6", "_1", "_2"];
else
    %dynamic
    systemCombos1 = ["Photoneo2"];
    systemCombos2 = ["Photoneo"];
    systemNames = ["Photoneo", "Photoneo2"];
    systemNames2 = ["", ""];
end

Checkerboards3d_func(basePath, baseFoldercb3d, patientNumscb3d, systemCombos1, systemCombos2, systemNames, systemNames2, directVersion, printLevel, numOfCheckerboards)