function DoDMQVSummary(basePath, patientNumscb3d, DoEvalForCombination)

outputPath = append(basePath, "\Output");
if ~exist(outputPath, "dir")
    mkdir(outputPath)
end
outputFileName = outputPath + "\\EvaluationSummary.csv";
if(isfile(outputFileName))
    delete(outputFileName);
end
% write header
fid = fopen(outputFileName, 'a+');
fprintf(fid, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s, %s\n", "Mode", "ObjectName", "RegistrationPCDiffRMSE", "Shutters", "Distance", "Resolution", "Combination", "CameraName", "PtDensityMM", "ZMeanMM", "ZStdMM", "Description", "usedMetaDataTransform");
fclose(fid);

for pi=patientNumscb3d
    subjectFolderPath = append(basePath, "\", string(pi));
    if ~exist(subjectFolderPath, 'dir')
        continue;
    end
    % get all eval json
    jsonFiles = dir(append(subjectFolderPath, "\Output\*.json"));
    for jfi=1:length(jsonFiles)
        % check whether eval required
        doEval = false;
        for eci=1:length(DoEvalForCombination)
            sysName = split(DoEvalForCombination(eci), "-");
            sysName = sysName(1);
            if startsWith(jsonFiles(jfi).name, sysName)
                doEval = true;
                break;
            end
        end
        if doEval
            jsonFilePath = append(jsonFiles(jfi).folder, "\", jsonFiles(jfi).name);
            EvaluationSummary.DoJsonMinimal(outputPath, jsonFilePath);
        end
    end
end

% show plots
evalSummaryCsvPath = outputFileName;
BoxPlotEvaluationSummary(evalSummaryCsvPath, DoEvalForCombination);

end

