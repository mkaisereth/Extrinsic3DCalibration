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

DoEvalForCombination = ["Photoneo-Photoneo"; "Sltida-Photoneo"; "Astida-Photoneo"; "Astra-Photoneo"; "D415-Photoneo"];

doDMQV = true;
global printLevel;
printLevel = 0;

% first look for corresponding transform (Larsen same day, Balgrist same
% patient), if not found whether to use the nearest transform
doSearchForNearestTransform = true;

% calculate all transformations for all systems
systemNames = ["Photoneo", "Astida", "Sltida", "Astra", "Astra", "D415", "D415"];
systemNames2 = ["", "", "", "_5", "_6", "_1", "_2"];

systemCombos1 = ["Photoneo", "Sltida", "Astida", "Astra_5", "D415_1", "Astra_6", "D415_2"];
systemCombos2 = ["Photoneo", "Photoneo", "Photoneo", "Photoneo", "Photoneo", "Astra_5", "D415_1"];

%%
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

totalNum = 0;
for sci=1:length(systemCombos1)
    ids1 = [plyDict.(systemCombos1(sci)).id];
    for ii=1:length(ids1)
        totalNum = totalNum+1;
    end
end
inc = 0;
startTime = now;

% loop all system combos and apply transformation
if printLevel > 0
f1 = figure;
f2 = figure;
end
for sci=1:length(systemCombos1)
    ids1 = [plyDict.(systemCombos1(sci)).id];
    ids2 = [plyDict.(systemCombos2(sci)).id];
    if length(ids1) ~= length(ids2)
        disp('Warning: number of ids different')
        keyboard
    end
    for ii=1:length(ids1)

        inc = inc+1;
        if inc > 1
            eta = ((now-startTime)/(inc-1)*(totalNum-inc))*24*60;
        else
            eta = 0;
        end
        disp(inc + "/" + totalNum + " eta: " + round(eta,1) + " min");

        if strcmp(ids1(ii), ids2(ii)) == 0
            disp('Warning: order of ids different')
            keyboard
        end
        % now do the calibration for each combo
        plyDictCombo1 = plyDict.(systemCombos1(sci));
        plyDictCombo2 = plyDict.(systemCombos2(sci));

        pc1FileName = plyDictCombo1(ii).name;
        pc1Path = append(plyDictCombo1(ii).folder, "\", pc1FileName);
        pc1 = pcread(pc1Path);
        % second is reference coordinate system
        pc2_refFileName = plyDictCombo2(ii).name;
        pc2_refPath = append(plyDictCombo2(ii).folder, "\", pc2_refFileName);
        pc2_ref = pcread(pc2_refPath);
        
        [pc1, pc2_ref] = PcsSegDist(pc1, pc2_ref);
        
        comboName = append(systemCombos1(sci), "-", systemCombos2(sci), "-", ids1(ii));
        % read the transformation
        tFormPath = append(cb3dPath, "\", comboName, "-tForm2s.mat");
        if ~exist(tFormPath, 'file')
            if doSearchForNearestTransform
                % try to find a close transform (e.g. from another time or
                % date)
                for diffCounter=1:max(abs(ii-1),abs(length(ids1)-ii))
                    diffUp = ii+diffCounter;
                    if diffUp<=length(ids1)
                        comboName = append(systemCombos1(sci), "-", systemCombos2(sci), "-", ids1(diffUp));
                        tFormPath = append(cb3dPath, "\", comboName, "-tForm2s.mat");
                        if exist(tFormPath, 'file')
                            % OK, found one, continue
                            break;
                        end
                    end
                    diffDown = ii-diffCounter;
                    if diffDown>=1
                        comboName = append(systemCombos1(sci), "-", systemCombos2(sci), "-", ids1(diffDown));
                        tFormPath = append(cb3dPath, "\", comboName, "-tForm2s.mat");
                        if exist(tFormPath, 'file')
                            % OK, found one, continue
                            break;
                        end
                    end
                end
                % check again whether found
                if ~exist(tFormPath, 'file')
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
        tForm2s = load(tFormPath).tForm2s;
        % apply it and show the result
        pc1t = pctransform(pc1, tForm2s);
        
        if printLevel > 0
            figure(1);
            hold off;
            pcshow(pc2_ref);
            hold on;

            % shift a bit for better visualization
            pc1tLoc = pc1t.Location;
            pc1tLoc(:,3) = pc1tLoc(:,3)-0.02;

            pcshow(pc1tLoc, 'MarkerSize', 64);
            view(0,-89);
            figure(2);
            hold off;
            pcshow(pc2_ref.Location, 'r');
            hold on;
            pcshow(pc1t.Location, 'g');
            view(0,-89);
            disp('')
            pause(0.1);
        end

        % now do DMQV
        if doDMQV
            DoDMQV(pc1Path, pc1FileName, pc2_refPath, tFormPath, comboName, basePath);
        end
    end
end

%% do DMQV summary
if doDMQV
    addpath("dmqvmatlab")
    addpath("Utils")
    DoDMQVSummary(basePath, patientNumscb3d, DoEvalForCombination);
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