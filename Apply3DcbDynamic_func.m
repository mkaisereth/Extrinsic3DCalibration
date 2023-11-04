function [rmseVals1s, stdVals1s] = Apply3DcbDynamic_func(baseFolder, cb3dBaseFolder, genericcb3dFolderPaths, doUseIcp, doUseUserInput, letUserSelectMarkers, patientNums, directVersion, withTexture, lPrintLevel)

% script to merge sagittal bending captures from two photoneos
global printLevel;
printLevel = lPrintLevel;

mergedPath = append(baseFolder, "\Merged");
if ~exist(mergedPath, "dir")
    mkdir(mergedPath);
end

numTotal = length(patientNums);
startTime = now;

% do this for all subfolders
clearvars rmseVals1s;
clearvars stdVals1s;
rmseVals1s = cell(numTotal,1);
stdVals1s = cell(numTotal,1);
% parallel can only be used without printLevel > 0
if lPrintLevel <= 0
    numWorkers = parcluster;
    numWorkers = numWorkers.NumWorkers;
else
    numWorkers = 0;
end
parfor (patInd=patientNums, numWorkers)
    disp("Patient " + string(patInd))
    num = patInd;
    if num > 1
        eta = ((now-startTime)/(num-1)*(numTotal-num))*24*60;
    else
        eta = 0;
    end
    disp(num + "/" + numTotal + " eta: " + round(eta,1) + " min");

    if lPrintLevel>0
    close all;
    end
    % get the patient folder
    if directVersion
        patientFolderPath = baseFolder;
        patientcb3dFolderPath = cb3dBaseFolder;
    else
        patientFolderPath = append(baseFolder, "\", string(patInd));
        patientcb3dFolderPath = append(cb3dBaseFolder, "\", string(patInd));
    end

    % get id to check for corresponding transformation file
    firstPlyFileName = dir(append(patientFolderPath, "\Photoneo_*.ply"));
    if isempty(firstPlyFileName)
        keyboard;
    end
    firstPlyFileName = firstPlyFileName(1).name;
    comboId = replace(firstPlyFileName, "Photoneo_", "");
    comboId = comboId(1:8);
    comboName = append("Photoneo2-Photoneo-", comboId);
    % read the transformation

    % search for 3d checkerboard transform in patient specific tranform folder
    tFormPath = append(patientcb3dFolderPath, "\Output\", comboName, "-tForm2s.mat");
    % backup, search in general folder
    if ~exist(tFormPath, "file")
        for genericcb3dFolderPath=genericcb3dFolderPaths
            tFormPath = append(genericcb3dFolderPath, "\", comboName, "-tForm2s.mat");
            if exist(tFormPath, "file")
                break;
            end
        end
    end
    
    % filter for first and second photoneo captures
    phneoFilter = "Photoneo_*.ply";
    phneo2Filter = "Photoneo2_*.ply";
        
    phneoFiles = dir(append(patientFolderPath, "\", phneoFilter));
    phneo2Files = dir(append(patientFolderPath, "\", phneo2Filter));
    
    numOfPhneoFiles = length(phneoFiles);
    numOfPhneo2Files = length(phneo2Files);
    % number of captures should be equals (HW triggered)
    if numOfPhneoFiles ~= numOfPhneo2Files
        disp("Superbad, now make matching!!!");
    end
    numOfFiles = min(numOfPhneoFiles, numOfPhneo2Files);
    % read all point clouds
    pcs = pointCloud.empty(0, numOfFiles);
    pcs2 = pointCloud.empty(0, numOfFiles);
    for i=1:numOfFiles
        pcs(i) = pcread(append(phneoFiles(i).folder, "\", phneoFiles(i).name));
        pcs2(i) = pcread(append(phneo2Files(i).folder, "\", phneo2Files(i).name));
    end
    
    %% transform second photoneo point clouds
    
    % sync offset from first photoneo (w.r.t. second)
    syncOffset = 0;
    % check whether json file has syncOffset set
    phneoJsonFiles = dir(append(patientFolderPath, "\*_", replace(phneoFilter, ".ply", ".json")));
    if length(phneoJsonFiles)>0
        fid = fopen(append(phneoJsonFiles(1).folder, "\", phneoJsonFiles(1).name)); 
        raw = fread(fid,inf); 
        str = char(raw'); 
        fclose(fid); 
        jsonData2 = jsondecode(str);
        if isfield(jsonData2, 'PhotoneoSyncDiff')
            if isempty(jsonData2.PhotoneoSyncDiff)
                disp("Warning: this patient has been declared as failed: " + string(patInd))
                continue;
            end
            syncOffset = jsonData2.PhotoneoSyncDiff;
        end
    end
    if syncOffset<0
        syncOffset2 = 1;
        syncOffset = 0;
    else
        syncOffset2 = 0;
    end
    % load tForm (from checkerboard calibration)
    tForm = load(tFormPath);
    tForm = tForm.tForm2s;
    
    % transform all point clouds from 2nd photoneo to first
    % -1 and +1 because captures don't seem to be synchronized in time (maybe
    % first photoneo captured one capture more in beginning)
    pcst = pointCloud.empty(0,numOfFiles);
    pcs2t = pointCloud.empty(0,numOfFiles);
    for i=1:numOfFiles
        syncInd = i+syncOffset;
        syncInd2 = i+syncOffset2;
        if syncInd > numOfFiles || syncInd2 > numOfFiles
            continue;
        end
        % cut stuff away
        pcCols = pcs(syncInd).Color; % +1 sync
        pcLoc = pcs(syncInd).Location; % +1 sync
        inds = pcLoc(:,2) > 0.7 | pcLoc(:,2) < -0.7 | pcLoc(:,1) > 0.5 | pcLoc(:,1) < -0.5 ...
            | pcLoc(:,3) < 0.1; % cut away floor, headrest, ...
        pcCols(inds,:) = [];
        pcLoc(inds,:) = [];
        pcst(i) = pointCloud(pcLoc, 'Color', pcCols);
        
        pcs2t(i) = pctransform(pcs2(syncInd2),tForm); % transform P2 to P1
        % cut stuff away
        pcCols = pcs2t(i).Color;
        pcLoc = pcs2t(i).Location;
        %inds = pcLoc(:,2) > 0.7 | pcLoc(:,2) < -0.7 | pcLoc(:,1) > 0.5 | pcLoc(:,1) < -0.5 ...
        %    | pcLoc(:,3) < 0.1;
        %pcCols(inds,:) = [];
        %pcLoc(inds,:) = [];
        pcs2t(i) = pointCloud(pcLoc, 'Color', pcCols);
    end
    %% find overlapping clusters
    distRejection = 0.006; % distance to surface normal must no be larger than 6mm, otherwise not overlapping
    [rmseVals1,stdVals1, icptForms, closeRatios] = EvaluateOverlappingClusters(pcst,pcs2t, patientFolderPath, "Photoneo2-Photoneo", true, false, distRejection, doUseIcp, doUseUserInput);    
%     if patInd<=length(rmseVals1s) && ~isempty(rmseVals1s{patInd})
%         disp("Warning: this code is now different!!!")
%     end
    rmseVals1s{patInd} = rmseVals1;
    stdVals1s{patInd} = stdVals1;
    for i=1:length(pcst)
        if doUseIcp
            pcs2t(i) = pctransform(pcs2t(i), icptForms(i));
        end
    end

    %% let user select some markers
    if letUserSelectMarkers
        f139 = figure;
        f146 = figure;
        for i=1:length(pcst)
            if closeRatios(i)>0.5
                figure(f139.Number);
                pcshow(pcst(i));
                view(0,-30)
                figure(f146.Number);
                pcshow(pcs2t(i));
                view(0,-30)
                circlePts1 = [];
                circlePts2 = [];
                circlePtInds1 = [];
                circlePtInds2 = [];
                for mit=1:3
                    figure(f139.Number);
                    [cpInd1, circlePt1] = SelectPointFromPc(pcst(i).Location', true);
                    circlePts1 = [circlePts1;circlePt1];
                    circlePtInds1 = [circlePtInds1;cpInd1];
                    hold on;
                    h1 = plot3(circlePt1(:,1), circlePt1(:,2), circlePt1(:,3), 'r.', 'MarkerSize', 12);
                    figure(f146.Number);
                    [cpInd2, circlePt2] = SelectPointFromPc(pcs2t(i).Location', true);
                    circlePts2 = [circlePts2;circlePt2];
                    circlePtInds2 = [circlePtInds2;cpInd2];
                    hold on;
                    h2 = plot3(circlePt2(:,1), circlePt2(:,2), circlePt2(:,3), 'r.', 'MarkerSize', 12);
                    pDiff = norm(circlePt1-circlePt2)
                end
                % estimate another transformation
                pc2Shift = mean(circlePts1-circlePts2);
                pc2RigidTForm = rigid3d(eye(3), pc2Shift);
                pcs2t(i) = pctransform(pcs2t(i), pc2RigidTForm);
                for mit=1:3
                    circlePt1 = pcst(i).Location(circlePtInds1(mit), :);
                    circlePt2 = pcs2t(i).Location(circlePtInds2(mit), :);
                    pDiff = norm(circlePt1-circlePt2);
                end
            end
        end
        distRejection = 0.006; % distance to surface normal must no be larger than 6mm, otherwise not overlapping
        [rmseVals1,stdVals1, icptForms, closeRatios] = EvaluateOverlappingClusters(pcst,pcs2t, patientFolderPath, "Photoneo2-Photoneo", true, false, distRejection, doUseIcp, doUseUserInput);
        for i=1:length(pcst)
            if doUseIcp
                pcs2t(i) = pctransform(pcs2t(i), icptForms(i));
            end
        end
    end

    %% save merged pcs
    for i=1:length(pcst)
        pcLocs = [pcst(i).Location;pcs2t(i).Location];
        pcCols = [pcst(i).Color;pcs2t(i).Color];
        mergedPc = pointCloud(pcLocs, 'Color', pcCols);
        pcwrite(mergedPc, append(mergedPath, "\Photoneo12_", string(i)), "Encoding", "binary");
    end

    %% visualize as video
    if lPrintLevel>0
    rotateCamera = false;
    rotateCamAngle = 0;
    
    xMin=inf;
    xMax=-inf;
    yMin=inf;
    yMax=-inf;
    zMin=inf;
    zMax=-inf;
    for i=1:length(pcst)
        [xMin, xMax] = GetMinMax(pcst(i), "XLimits", xMin,xMax);
        [yMin, yMax] = GetMinMax(pcst(i), "YLimits", yMin,yMax);
        [zMin, zMax] = GetMinMax(pcst(i), "ZLimits", zMin,zMax);
    end
    for i=1:length(pcs2t)
        [xMin, xMax] = GetMinMax(pcs2t(i), "XLimits", xMin,xMax);
        [yMin, yMax] = GetMinMax(pcs2t(i), "YLimits", yMin,yMax);
        [zMin, zMax] = GetMinMax(pcs2t(i), "ZLimits", zMin,zMax);
    end
    f = figure;
    end
    for i=1:length(pcst)
        if lPrintLevel>0
            if withTexture
                pcshow(pcst(i), 'MarkerSize', 12);
                hold on;
                pcshow(pcs2t(i), 'MarkerSize', 12);
            else
                pcshow(pcst(i).Location, 'r', 'MarkerSize', 12);
                hold on;
                pcshow(pcs2t(i).Location, 'g', 'MarkerSize', 12);
            end
            
            xlim([xMin xMax]);
            ylim([yMin yMax]);
            zlim([zMin zMax]);
            view(0,-89.9);
            title("");
            set(gcf,'color','w');
            axis off;
            camorbit(-89 + rotateCamAngle,0,'camera');
            hold off;
            drawnow;
            pause(0.1);
            if ~ishandle(f)
                %break;
            end
            if rotateCamera
                rotateCamAngle = rotateCamAngle + 10;
            end
        end
        if lPrintLevel>0
            if ishandle(f)
                title("pause");
            else
                %break;
            end
            pause(1);
        end
    end
end
end

function [minX,maxX] = GetMinMax(pc, DimLimits, minX, maxX)

tempMinX = pc.(DimLimits)(1);
if tempMinX<minX
    minX = tempMinX;
end
tempMaxX = pc.(DimLimits)(2);
if tempMaxX>maxX
    maxX = tempMaxX;
end

end