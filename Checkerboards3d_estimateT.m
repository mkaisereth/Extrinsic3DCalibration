function Checkerboards3d_estimateT(pc1Path, pc2Path, pc1_ref, pc2, cb3dParams, printLevel, outputFolderPath, outputFileName)

% check whether all holes are visible (no config file), else use
% information from config file
jsonFilePath1 = replace(pc1Path, ".ply", ".json");
jsonFilePath2 = replace(pc2Path, ".ply", ".json");
cb3dParams1 = cb3dParams;
cb3dParams2 = cb3dParams;
rotMatZ = [];
if exist(jsonFilePath1, 'file') && exist(jsonFilePath2, 'file')
% if one has less, the other has to have less as well
% TODO
elseif exist(jsonFilePath1, 'file')
    fid = fopen(jsonFilePath1); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    jsonData1 = jsondecode(str);
    if isfield(jsonData1, "LeftColsVisible")
        cb3dParams1.LeftColsVisible = jsonData1.LeftColsVisible;
        cb3dParams1.RightColsVisible = jsonData1.RightColsVisible;
        % adjust holes, and columns
        cb3dParams1.NumOfCols = cb3dParams1.LeftColsVisible+cb3dParams1.RightColsVisible;
        cb3dParams1.NumOfHoles = 2*cb3dParams1.NumOfCols*cb3dParams1.NumOfRows;
        % if one has less, the other has to have less as well
        cb3dParams2.LeftColsVisible = cb3dParams1.LeftColsVisible;
        cb3dParams2.RightColsVisible = cb3dParams1.RightColsVisible;
        cb3dParams2.NumOfCols = cb3dParams1.NumOfCols;
    end
    if isfield(jsonData1, "TopRowsVisible")
        cb3dParams1.TopRowsVisible = jsonData1.TopRowsVisible;
        cb3dParams1.BottomRowsVisible = jsonData1.BottomRowsVisible;
        % adjust holes, and columns
        cb3dParams1.NumOfRows = cb3dParams1.TopRowsVisible+cb3dParams1.BottomRowsVisible;
        cb3dParams1.NumOfHoles = 2*cb3dParams1.NumOfCols*cb3dParams1.NumOfRows;
        % if one has less, the other has to have less as well
        cb3dParams2.TopRowsVisible = cb3dParams1.TopRowsVisible;
        cb3dParams2.BottomRowsVisible = cb3dParams1.BottomRowsVisible;
        cb3dParams2.NumOfRows = cb3dParams1.NumOfRows;
    end
elseif exist(jsonFilePath2, 'file')
    fid = fopen(jsonFilePath2); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    jsonData2 = jsondecode(str);
    if isfield(jsonData2, "LeftColsVisible")
        cb3dParams2.LeftColsVisible = jsonData2.LeftColsVisible;
        cb3dParams2.RightColsVisible = jsonData2.RightColsVisible;
        % adjust holes, and columns
        cb3dParams2.NumOfCols = cb3dParams2.LeftColsVisible+cb3dParams2.RightColsVisible;
        cb3dParams2.NumOfHoles = 2*cb3dParams2.NumOfCols*cb3dParams2.NumOfRows;
        rotZAngle = jsonData2.RotationAroundZ;
        rotMatZ = [cosd(rotZAngle) -sind(rotZAngle) 0; sind(rotZAngle) cosd(rotZAngle) 0; 0 0 1];
        pcLoc = pc2.Location;
        pcLoc = (rotMatZ*pcLoc')';
        pc2_orig = pc2;
        pc2 = pointCloud(pcLoc);
        % if one has less, the other has to have less as well
        cb3dParams1.LeftColsVisible = cb3dParams2.LeftColsVisible;
        cb3dParams1.RightColsVisible = cb3dParams2.RightColsVisible;
        cb3dParams1.NumOfCols = cb3dParams2.NumOfCols;
    end
    if isfield(jsonData2, "TopRowsVisible")
        cb3dParams2.TopRowsVisible = jsonData2.TopRowsVisible;
        cb3dParams2.BottomRowsVisible = jsonData2.BottomRowsVisible;
        % adjust holes, and columns
        cb3dParams2.NumOfRows = cb3dParams2.TopRowsVisible+cb3dParams2.BottomRowsVisible;
        cb3dParams2.NumOfHoles = 2*cb3dParams2.NumOfCols*cb3dParams2.NumOfRows;
        % if one has less, the other has to have less as well
        cb3dParams1.TopRowsVisible = cb3dParams2.TopRowsVisible;
        cb3dParams1.BottomRowsVisible = cb3dParams2.BottomRowsVisible;
        cb3dParams1.NumOfRows = cb3dParams2.NumOfRows;
    end
end

[reconstructed2, ~, holeMedians2] = Checkerboard3d(pc2, cb3dParams2, printLevel, 2.5);
[reconstructed1, ~, holeMedians1] = Checkerboard3d(pc1_ref, cb3dParams1, printLevel);

%% registration
if printLevel > 1
    figure;
    pcshow(pc1_ref.Location, 'c');
    hold on;
    pcshow(pc2.Location, 'm');
    pcshow(reconstructed1, 'r', 'MarkerSize', 24);
    pcshow(reconstructed2, 'g', 'MarkerSize', 24);
    title("checkerboard back in 3d space");
end

%% another approach, try to find correspondence points (with sparse)

[largestPointInd, smallestPointInd, middle1PointInd] = FindTopRightPoint(holeMedians1);
[largestPoint2Ind, smallestPoint2Ind, middle1Point2Ind] = FindTopRightPoint(holeMedians2);

% visualize
if printLevel > 1
    figure;
    pcshow(holeMedians1, 'r', 'MarkerSize', 24);
    hold on;
    pcshow(holeMedians1(largestPointInd, :), 'b', 'MarkerSize', 24);
    pcshow(holeMedians1(smallestPointInd, :), 'y', 'MarkerSize', 24);
    pcshow(holeMedians1(middle1PointInd, :), 'm', 'MarkerSize', 24);
    pcshow(holeMedians2, 'g', 'MarkerSize', 24);
    pcshow(holeMedians2(largestPoint2Ind, :), 'b', 'MarkerSize', 24);
    pcshow(holeMedians2(smallestPoint2Ind, :), 'y', 'MarkerSize', 24);
    pcshow(holeMedians2(middle1Point2Ind, :), 'm', 'MarkerSize', 24);
end

% now estimate transformation
tPts1 = [holeMedians1(largestPointInd, :); holeMedians1(smallestPointInd, :); holeMedians1(middle1PointInd, :)];
tPts2 = [holeMedians2(largestPoint2Ind, :); holeMedians2(smallestPoint2Ind, :); holeMedians2(middle1Point2Ind, :)];

tFormEst = estimateGeometricTransform3D(double(tPts2),double(tPts1),'rigid');

% visualize
if printLevel > 1
    figure;
    tPts2t = pctransform(pointCloud(holeMedians2), tFormEst);
    tPts2t = tPts2t.Location;
    pcshow(holeMedians1, 'r', 'MarkerSize', 24);
    hold on;
    pcshow(holeMedians1(largestPointInd, :), 'b', 'MarkerSize', 24);
    pcshow(holeMedians1(smallestPointInd, :), 'y', 'MarkerSize', 24);
    pcshow(tPts2t, 'g', 'MarkerSize', 24);
    pcshow(tPts2t(largestPoint2Ind, :), 'b', 'MarkerSize', 24);
    pcshow(tPts2t(smallestPoint2Ind, :), 'y', 'MarkerSize', 24);
end

% visualize the real checkerboard
ptCloudTrans2Est = pctransform(pc2, tFormEst);
if printLevel > 1
    figure;
    pcshow(pc1_ref.Location, 'r', 'MarkerSize', 24);
    hold on;
    pcshow(ptCloudTrans2Est.Location, 'g', 'MarkerSize', 24);
    title("checkerboards calibrated (estimate)");
end

% do with sparse data
tForm2s = tFormEst;
% save it
if ~isempty(outputFolderPath)
    % if there was a transformation, add the inverse again
    if ~isempty(rotMatZ)
        rotMatZT = rigid3d(rotMatZ', [0 0 0]);
        overallT = rotMatZT.T*tForm2s.T;
        tForm2s = rigid3d(overallT);
        % revert to original now
        pc2 = pc2_orig;
    end
    save(convertStringsToChars(append(outputFolderPath, "\", outputFileName, "-tForm2s.mat")), 'tForm2s');
end
% transform the original point cloud to see the result
ptCloudTrans2s = pctransform(pc2, tForm2s);
if printLevel > 0
    fig157 = figure;
    pcshow(pc1_ref.Location, 'r', 'MarkerSize', 24);
    hold on;
    pcshow(ptCloudTrans2s.Location, 'g', 'MarkerSize', 24);
    title("checkerboards calibrated: " + outputFileName);
    view(0,-89.9);
    saveas(fig157, convertStringsToChars(append(outputFolderPath, "\", outputFileName, "-tForm2s.fig")));
end
pause(0.5);

end