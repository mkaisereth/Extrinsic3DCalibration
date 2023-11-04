% Do the DMQV summary for Balgrist and Larsen study
clear all; close all; clc;

dataSetNames = ["Larsen","Balgrist"];
combinedOutputPath = "I:\Temp\";

boxplotTablePtDensityMMs = [];
boxplotTableZMeanMMs = [];
boxplotTableZStdMMs = [];
uniqueCombinationss = [];
for dataSetName=dataSetNames
    if strcmp(dataSetName, "Larsen")
        basePath = "I:\Temp\3DcbData\Study1Data";
    end
    if strcmp(dataSetName, "Balgrist")
        basePath = "I:\Temp\3DcbData\Study2Data\Static";
    end
    
    % Larsen study has patiens from 5 to 75
    if strcmp(dataSetName, "Larsen")
        patientNumscb3d = 5:75;
    end
    if strcmp(dataSetName, "Balgrist")
        patientNumscb3d = 1:22;
    end
    
    DoEvalForCombination = ["Photoneo-Photoneo"; "Sltida-Photoneo"; "Astida-Photoneo"; "Astra-Photoneo"; "D415-Photoneo"];

    global printLevel;
    printLevel = 0;
    
    % calculate all transformations for all systems
    systemNames = ["Photoneo", "Astida", "Sltida", "Astra", "Astra", "D415", "D415"];
    systemNames2 = ["", "", "", "_5", "_6", "_1", "_2"];
    
    systemCombos1 = ["Photoneo", "Sltida", "Astida", "Astra_5", "D415_1", "Astra_6", "D415_2"];
    systemCombos2 = ["Photoneo", "Photoneo", "Photoneo", "Photoneo", "Photoneo", "Astra_5", "D415_1"];
    
    %
    addpath("DMQVMatlab")
    addpath("Utils")
    outputPath = append(basePath, "\Output");
    outputFileName = outputPath + "\\EvaluationSummary.csv";
    % show plots
    evalSummaryCsvPath = outputFileName;
    [boxplotTablePtDensityMM, boxplotTableZMeanMM, boxplotTableZStdMM, uniqueCombinations] = BoxPlotEvaluationSummary(evalSummaryCsvPath, DoEvalForCombination);
    if isempty(uniqueCombinationss)
        uniqueCombinationss = uniqueCombinations;
    else
        for uci=1:length(uniqueCombinationss)
            if strcmp(uniqueCombinationss(uci), uniqueCombinations(uci)) ~= 1
                disp("Warning: unique combinations MUST match!")
                keyboard;
            end
        end
    end
    boxplotTablePtDensityMMs = [boxplotTablePtDensityMMs; boxplotTablePtDensityMM];
    boxplotTableZMeanMMs = [boxplotTableZMeanMMs; boxplotTableZMeanMM];
    boxplotTableZStdMMs = [boxplotTableZStdMMs; boxplotTableZStdMM];
end

%% plot the combined results
uniqueCombinations = uniqueCombinationss;
filepath = combinedOutputPath;

f1 = figure;
boxplot(boxplotTablePtDensityMMs, 'Labels', uniqueCombinations);
title("Point density 1/[mm^2]")
ylabel("1/[mm^2]")
saveas(f1, append(filepath, "\PointDensity_", join(uniqueCombinations, "_"), ".fig"));
ax = gca;
exportgraphics(ax,append(filepath, "\PointDensity_", join(uniqueCombinations, "_"), ".png"),"Resolution",600)

f2 = figure;
boxplot(boxplotTableZMeanMMs, 'Labels', uniqueCombinations);
title("Z Mean difference [mm]")
ylabel("[mm]")
saveas(f2, append(filepath, "\ZMeanDiff_", join(uniqueCombinations, "_"), ".fig"));
ax = gca;
exportgraphics(ax,append(filepath, "\ZMeanDiff_", join(uniqueCombinations, "_"), ".png"),"Resolution",600)
% do it again with DataLim
f2 = figure;
boxplot(boxplotTableZMeanMMs, 'Labels', ["MC1", "TIDA", "BoofCV", "Astra", "D415"], 'DataLim',[-30,30]);
title("RMSE")
ylabel("[mm]")
saveas(f2, append(filepath, "\ZMeanDiff_", join(uniqueCombinations, "_"), "_DataLim.fig"));
ax = gca;
exportgraphics(ax,append(filepath, "\ZMeanDiff_", join(uniqueCombinations, "_"), "_DataLim.png"),"Resolution",600)

f3 = figure;
boxplot(boxplotTableZStdMMs, 'Labels', uniqueCombinations);
title("Z Std [mm]")
ylabel("[mm]")
saveas(f3, append(filepath, "\ZStd_", join(uniqueCombinations, "_"), ".fig"));
ax = gca;
exportgraphics(ax,append(filepath, "\ZStd_", join(uniqueCombinations, "_"), ".png"),"Resolution",600)
% do it again with DataLim
f3 = figure;
boxplot(boxplotTableZStdMMs, 'Labels', ["MC1", "TIDA", "BoofCV", "Astra", "D415"], 'DataLim',[-Inf,30]);
title("SD")
ylabel("[mm]")
saveas(f3, append(filepath, "\ZStd_", join(uniqueCombinations, "_"), "_DataLim.fig"));
ax = gca;
exportgraphics(ax,append(filepath, "\ZStd_", join(uniqueCombinations, "_"), "_DataLim.png"),"Resolution",600)

%% do it for UseCase 1: left/right camera pairs

clearvars boxplotVals;
fnss = [];
for dataSetName=dataSetNames
    if strcmp(dataSetName, "Larsen")
        basePath = "I:\Temp\LarsenStudyData_02_20230424_1\cut";
    end
    if strcmp(dataSetName, "Balgrist")
        basePath = "I:\Temp\BalgristStudyData_01_20230424_1\cut";
    end
    outputPath = append(basePath, "\Output");
    % and summary for overlapping regions
    evalOverClust = load(append(outputPath, "\EvalOverClust.mat"));
    evalOverClust = evalOverClust.evalOverClust;
    fns = fieldnames(evalOverClust);
    if isempty(fnss)
        fnss = fns;
    else
        for uci=1:length(fnss)
            if strcmp(fnss(uci), fns(uci)) ~= 1
                disp("Warning: field names MUST match!")
                keyboard;
            end
        end
    end
    for fni=fns'
        systemCombo12 = string(fni);
        if ~exist('boxplotVals', 'var') || ~isfield(boxplotVals, systemCombo12)
            boxplotVals.(systemCombo12) = [];
        end
        boxplotVals.(systemCombo12) = [boxplotVals.(systemCombo12);[evalOverClust.(systemCombo12).rmseVals*1000, evalOverClust.(systemCombo12).stdVals*1000]];
    end
end
for fni=fnss'
    systemCombo12 = string(fni);
    f286 = figure;
    boxplot(boxplotVals.(systemCombo12), 'Labels', ["RMSE", "SD"]);
    ylabel("[mm]")
    title(replace(replace(replace(replace(evalOverClust.(systemCombo12).Name, "_6-", " (right) - "), "_5", " (left)"), "_2-", " (right) - "), "_1", " (left)"));
    saveas(f286, append(combinedOutputPath, "\Eval_", evalOverClust.(systemCombo12).Name, ".fig"));
    ax = gca;
    exportgraphics(ax,append(combinedOutputPath, "\Eval_", evalOverClust.(systemCombo12).Name, ".png"),"Resolution",600)
    % do a version with cut axis
    f134 = figure;
    boxplot(boxplotVals.(systemCombo12), 'Labels', ["RMSE", "SD"], 'DataLim',[-Inf, 8]);
    ylabel("[mm]")
    title(replace(replace(replace(replace(evalOverClust.(systemCombo12).Name, "_6-", " (right) - "), "_5", " (left)"), "_2-", " (right) - "), "_1", " (left)"));
    saveas(f134, append(combinedOutputPath, "\Eval_", evalOverClust.(systemCombo12).Name, "_DataLim.fig"));
    ax = gca;
    exportgraphics(ax,append(combinedOutputPath, "\Eval_", evalOverClust.(systemCombo12).Name, "_DataLim.png"),"Resolution",600)
end