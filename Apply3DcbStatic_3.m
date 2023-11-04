clear all; close all; clc;

dataSetName = "Larsen";
dataSetName = "Balgrist";

if strcmp(dataSetName, "Larsen")
    cb3dPath = "I:\Temp\3DcbData\Study1Data\cb3d\Output";
    basePath = "I:\Temp\3DcbData\Study1Data";
end
if strcmp(dataSetName, "Balgrist")
    cb3dPath =  "I:\Temp\3DcbData\Study2Data\Static\cb3d\Output";
    basePath = "I:\Temp\3DcbData\Study2Data\Static";
end

% Larsen study has patiens from 5 to 75
if strcmp(dataSetName, "Larsen")
    patientNumscb3d = 5:75;
end
if strcmp(dataSetName, "Balgrist")
    patientNumscb3d = 1:22;
end

doDMQV = true;
global printLevel;
printLevel = 0;
showRef = true;

doUseIcp = true;
doUseUserInput = false;

% first look for corresponding transform (Larsen same day, Balgrist same
% patient), if not found whether to use the nearest transform
doSearchForNearestTransform = true;

% calculate all transformations for all systems

systemNames = ["Photoneo", "Sltida", "Astra", "Astra", "D415", "D415"];
systemNames2 = ["", "", "_5", "_6", "_1", "_2"];

systemCombos1 = ["Astra_5", "D415_1"];
systemCombos2 = ["Astra_6", "D415_2"];
systemCombosBase = ["Photoneo", "Photoneo"];

for sni=1:length(systemNames)
    fullSystemNames(sni) = append(systemNames(sni), systemNames2(sni));
    plyDict.(fullSystemNames(sni)) = [];
end

% loop all patients
for pi=patientNumscb3d
    subjectFolderPath = append(basePath, "\", string(pi));
    if ~exist(subjectFolderPath, 'dir')
        continue;
    end

    % get all files and sort
    plyList = dir(append(subjectFolderPath, "\*.ply"));
    for i=1:length(plyList)
        for sni=1:length(systemNames)
            if startsWith(plyList(i).name, systemNames(sni)) && endsWith(plyList(i).name, append(systemNames2(sni), ".ply"))
                % add id (day)
                plyList(i).id = extractAfter(plyList(i).name, strlength(systemNames(sni))+1);
                plyList(i).id = string(plyList(i).id(1:8));
                plyDict.(fullSystemNames(sni)) = [plyDict.(fullSystemNames(sni)), plyList(i)];
            end
        end
    end
end

% loop all system combos and apply transformation
if printLevel > 0
    f1 = figure;
    f2 = figure;
end
for sci=1:length(systemCombos1)
    systemCombo12 = append(systemCombos2(sci), systemCombos1(sci));
    evalOverClust.(systemCombo12) = [];
    ids1 = [plyDict.(systemCombos1(sci)).id];
    ids2 = [plyDict.(systemCombos2(sci)).id];
    idsBase = [plyDict.(systemCombosBase(sci)).id];
    if length(ids1) ~= length(idsBase) || length(ids2) ~= length(idsBase)
        disp('Warning: number of ids different')
        keyboard
    end
    for ii=1:length(ids1)
        disp(string(ii) + " / " + string(length(ids1)));
        close all;
        if strcmp(ids1(ii), idsBase(ii)) == 0 || strcmp(ids2(ii), idsBase(ii)) == 0
            disp('Warning: order of ids different')
            keyboard
        end
        % now do the calibration for each combo
        plyDictCombo1 = plyDict.(systemCombos1(sci));
        plyDictCombo2 = plyDict.(systemCombos2(sci));
        plyDictComboBase = plyDict.(systemCombosBase(sci));
        
        subjectFolderPath = plyDictCombo1(ii).folder;

        pc1FileName = plyDictCombo1(ii).name;
        pc1Path = append(plyDictCombo1(ii).folder, "\", pc1FileName);
        pc1 = pcread(pc1Path);
        pc2FileName = plyDictCombo2(ii).name;
        pc2Path = append(plyDictCombo2(ii).folder, "\", pc2FileName);
        pc2 = pcread(pc2Path);
        
        [pc1, pc2] = PcsSegDist(pc1, pc2);
        
        % second is reference coordinate system
        pcBase_refFileName = plyDictComboBase(ii).name;
        pcBase_refPath = append(plyDictComboBase(ii).folder, "\", pcBase_refFileName);
        pcBase_ref = pcread(pcBase_refPath);        
        comboName1 = append(systemCombos1(sci), "-", systemCombosBase(sci), "-", ids1(ii));
        comboName2 = append(systemCombos2(sci), "-", systemCombosBase(sci), "-", ids2(ii));
        camName = convertStringsToChars(systemCombos1(sci));
        camName = camName(1:end-2);
        comboName = append(camName, "-", systemCombosBase(sci), "-", ids1(ii));
        % read the transformation
        tFormPath1 = append(cb3dPath, "\", comboName1, "-tForm2s.mat");
        if ~exist(tFormPath1, 'file')
            if doSearchForNearestTransform
                % try to find a close transform (e.g. from another time or
                % date)
                for diffCounter=1:max(abs(ii-1),abs(length(ids1)-ii))
                    diffUp = ii+diffCounter;
                    if diffUp<=length(ids1)
                        comboName1 = append(systemCombos1(sci), "-", systemCombosBase(sci), "-", ids1(diffUp));
                        tFormPath1 = append(cb3dPath, "\", comboName1, "-tForm2s.mat");
                        if exist(tFormPath1, 'file')
                            % OK, found one, continue
                            break;
                        end
                    end
                    diffDown = ii-diffCounter;
                    if diffDown>=1
                        comboName1 = append(systemCombos1(sci), "-", systemCombosBase(sci), "-", ids1(diffDown));
                        tFormPath1 = append(cb3dPath, "\", comboName1, "-tForm2s.mat");
                        if exist(tFormPath1, 'file')
                            % OK, found one, continue
                            break;
                        end
                    end
                end
                % check again whether found
                if ~exist(tFormPath1, 'file')
                    disp('Warning: no suitable tForm found!!!')
                    %keyboard;
                    continue;
                end
            else
                disp('Warning: no suitable tForm found!!!')
                %keyboard;
                continue;
            end
        end
        tForm2s1 = load(tFormPath1).tForm2s;
        % apply it and show the result
        pc1t = pctransform(pc1, tForm2s1);
        % read the transformation
        tFormPath2 = append(cb3dPath, "\", comboName2, "-tForm2s.mat");
        if ~exist(tFormPath2, 'file')
            if doSearchForNearestTransform
                % try to find a close transform (e.g. from another time or
                % date)
                for diffCounter=1:max(abs(ii-1),abs(length(ids1)-ii))
                    diffUp = ii+diffCounter;
                    if diffUp<=length(ids1)
                        comboName2 = append(systemCombos2(sci), "-", systemCombosBase(sci), "-", ids1(diffUp));
                        tFormPath2 = append(cb3dPath, "\", comboName2, "-tForm2s.mat");
                        if exist(tFormPath2, 'file')
                            % OK, found one, continue
                            break;
                        end
                    end
                    diffDown = ii-diffCounter;
                    if diffDown>=1
                        comboName2 = append(systemCombos2(sci), "-", systemCombosBase(sci), "-", ids1(diffDown));
                        tFormPath2 = append(cb3dPath, "\", comboName2, "-tForm2s.mat");
                        if exist(tFormPath2, 'file')
                            % OK, found one, continue
                            break;
                        end
                    end
                end
                % check again whether found
                if ~exist(tFormPath2, 'file')
                    disp('Warning: no suitable tForm found!!!')
                    %keyboard;
                    continue;
                end
            else
                disp('Warning: no suitable tForm found!!!')
                %keyboard;
                continue;
            end
        end
        tForm2s2 = load(tFormPath2).tForm2s;
        % apply it and show the result
        pc2t = pctransform(pc2, tForm2s2);
        % transform into pc1 coordinate system (for DMQV)
        pc2tt = pctransform(pc2t, invert(tForm2s1));
        
        if printLevel > 0
            figure(1);
            hold off;
            pcshow(pc1t);
            hold on;
            pcshow(pc2t);
            if showRef
                pcshow(pcBase_ref);
            end
            view(0,-89);
            figure(2);
            hold off;
            pcshow(pc1t.Location, 'y');
            hold on;
            pcshow(pc2t.Location, 'g');
            if showRef
                hold on;
                pcshow(pcBase_ref.Location, 'r');
            end
            view(0,-89);
            disp('')
            pause(0.1);
        end
        
        % save the new (pc1 and pc2 merged) pointcloud as pc3
        pc3FileName = plyDictCombo1(ii).name;
        pc3FileName = convertStringsToChars(pc3FileName);
        pc3FileName = append(pc3FileName(1:end-6), ".ply");
        pc3Path = append(plyDictCombo1(ii).folder, "\", pc3FileName);
        pc3 = pointCloud([pc1.Location;pc2tt.Location]);
        pcwrite(pc3, pc3Path, 'Encoding', 'binary');
        % now do DMQV
        if doDMQV
            json1FileName = replace(pc1FileName, ".ply", ".json");
            json1FilePath = replace(pc1Path, pc1FileName, json1FileName);
            json3FileName = replace(pc3FileName, ".ply", ".json");
            json3FilePath = replace(pc1Path, pc1FileName, json3FileName);
            if ~exist(json3FilePath, 'file')
                copyfile(json1FilePath, json3FilePath);
            end
            DoDMQV(pc3Path, pc3FileName, pcBase_refPath, tFormPath1, comboName, basePath);
            % find overlapping clusters
            systemCombo12Translation = append(systemCombos2(sci), "-", systemCombos1(sci));
            distRejection = 0.006; % distance to surface normal must no be larger than 6mm, otherwise not overlapping
            [rmseVals1,stdVals1] = EvaluateOverlappingClusters(pc1, pc2tt, subjectFolderPath, systemCombo12Translation, false, true, distRejection, doUseIcp, doUseUserInput);
            
            evalOverClust.(systemCombo12).Name = systemCombo12Translation;
            if ~isfield(evalOverClust.(systemCombo12), 'rmseVals')
                evalOverClust.(systemCombo12).rmseVals = [];
            end
            if ~isfield(evalOverClust.(systemCombo12), 'stdVals')
                evalOverClust.(systemCombo12).stdVals = [];
            end
            evalOverClust.(systemCombo12).rmseVals = [evalOverClust.(systemCombo12).rmseVals; rmseVals1];
            evalOverClust.(systemCombo12).stdVals = [evalOverClust.(systemCombo12).stdVals; stdVals1];
        end
    end
end

%% do DMQV summary
if doDMQV
    outputPath = append(basePath, "\Output");
    if ~exist(outputPath, "dir")
        mkdir(outputPath)
    end
    % and summary for overlapping regions
    fns = fieldnames(evalOverClust);
    for fni=fns'
        systemCombo12 = string(fni);
        f286 = figure;
        boxplot([evalOverClust.(systemCombo12).rmseVals*1000, evalOverClust.(systemCombo12).stdVals*1000], 'Labels', ["ZMean", "ZStd"]);
        ylabel("[mm]")
        title(evalOverClust.(systemCombo12).Name);
        saveas(f286, append(outputPath, "\Eval_", evalOverClust.(systemCombo12).Name, ".fig"));
    end
    save(append(outputPath, "\EvalOverClust.mat"), 'evalOverClust');
end

%%
function [systemName, systemName2] = GetSystemNames(systemCombos)

splitStr = split(systemCombos, "_");
if length(splitStr)<=1
    systemName = systemCombos;
    systemName2 = "";
    return;
end
systemName = splitStr(1);
systemName2 = append("_", splitStr(2));
end