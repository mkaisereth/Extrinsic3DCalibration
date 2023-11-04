function [boxplotTablePtDensityMM, boxplotTableZMeanMM, boxplotTableZStdMM, uniqueCombinations] = BoxPlotEvaluationSummary(evalSummaryCsvPath, DoEvalForCombination)

if ~exist('evalSummaryCsvPath', 'var') || isempty(evalSummaryCsvPath)
    close all; clear all; clc;
    evalSummaryCsvPath = "G:\Temp\LarsenStudyData_02_20230424_1\cut\Output\EvaluationSummary.csv";
end

csvTable = readtable(evalSummaryCsvPath);

% make and show a boxplot for each combination
combinations = string(csvTable.Combination);
uniqueCombinations = unique(combinations);

boxplotTablePtDensityMM = [];
boxplotTableZMeanMM = [];
boxplotTableZStdMM = [];
if ~isempty(DoEvalForCombination)
    uniqueCombinations = DoEvalForCombination;
end
for uci=1:size(uniqueCombinations,1)
    comboInds = startsWith(combinations, uniqueCombinations(uci, :));
    ptDensityData = csvTable{comboInds,'PtDensityMM'};
    zMeanData = csvTable{comboInds,'ZMeanMM'};
    zStdData = csvTable{comboInds,'ZStdMM'};
    if ~isempty(boxplotTablePtDensityMM)
        matSize = size(boxplotTablePtDensityMM,1);
        dataSize = size(ptDensityData, 1);
        if dataSize<matSize
            % append nans to data
            tempd = nan(matSize, 1);
            tempd(1:dataSize, 1) = ptDensityData;
            ptDensityData = tempd;
            tempd = nan(matSize, 1);
            tempd(1:dataSize, 1) = zMeanData;
            zMeanData = tempd;
            tempd = nan(matSize, 1);
            tempd(1:dataSize, 1) = zStdData;
            zStdData = tempd;
        elseif dataSize>matSize
            % append nans to data
            tempd = nan(dataSize, size(boxplotTablePtDensityMM, 2));
            tempd(1:matSize, 1:size(boxplotTablePtDensityMM, 2)) = boxplotTablePtDensityMM;
            boxplotTablePtDensityMM = tempd;
            tempd = nan(dataSize, size(boxplotTableZMeanMM, 2));
            tempd(1:matSize, 1:size(boxplotTableZMeanMM, 2)) = boxplotTableZMeanMM;
            boxplotTableZMeanMM = tempd;
            tempd = nan(dataSize, size(boxplotTableZStdMM, 2));
            tempd(1:matSize, 1:size(boxplotTableZStdMM, 2)) = boxplotTableZStdMM;
            boxplotTableZStdMM = tempd;
        end
    end
    boxplotTablePtDensityMM = [boxplotTablePtDensityMM, ptDensityData];
    boxplotTableZMeanMM = [boxplotTableZMeanMM, zMeanData];
    boxplotTableZStdMM = [boxplotTableZStdMM, zStdData];
end

[filepath, ~, ~] = fileparts(evalSummaryCsvPath);

f1 = figure;
boxplot(boxplotTablePtDensityMM, 'Labels', uniqueCombinations);
title("Point density 1/[mm^2]")
ylabel("1/[mm^2]")
saveas(f1, append(filepath, "\PointDensity_", join(uniqueCombinations, "_"), ".fig"));
saveas(f1, append(filepath, "\PointDensity_", join(uniqueCombinations, "_"), ".png"));

f2 = figure;
boxplot(boxplotTableZMeanMM, 'Labels', uniqueCombinations);
title("Z Mean difference [mm]")
ylabel("[mm]")
saveas(f2, append(filepath, "\ZMeanDiff_", join(uniqueCombinations, "_"), ".fig"));
saveas(f2, append(filepath, "\ZMeanDiff_", join(uniqueCombinations, "_"), ".png"));

f3 = figure;
boxplot(boxplotTableZStdMM, 'Labels', uniqueCombinations);
title("Z Std [mm]")
ylabel("[mm]")
saveas(f3, append(filepath, "\ZStd_", join(uniqueCombinations, "_"), ".fig"));
saveas(f3, append(filepath, "\ZStd_", join(uniqueCombinations, "_"), ".png"));

end