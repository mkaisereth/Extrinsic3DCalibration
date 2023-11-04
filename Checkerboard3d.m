function [reconstructed, sparseRecon, holeMedians] = Checkerboard3d(pc, cb3dParams, printLevel, holecondition)

if ~exist('holecondition')
    holecondition = 0.85;
end

pcLocation = pc.Location;
pc = pointCloud(pcLocation);

% fit plane
[model,inlierInds,~] = pcfitplane(pc,0.02); %0.01% 1cm variation allowed

plane1 = select(pc,inlierInds);
[labels,numClusters] = pcsegdist(plane1,0.01);

if printLevel > 1
    figure;
    pcshow(pc);
    title("original pointcloud")

    figure;
    pcshow(pc);
    hold on;
    plot(model);
    title("original pointcloud with fitted plane")

    figure;
    pcshow(plane1);
    title("fitted plane only (inliers)");
    
    pcshow(plane1.Location,labels)
    colormap(hsv(numClusters))
    title('Point Cloud Clusters')
end

values = ones(length(labels),1);
labelsTable = table(labels, values);
groupsum = groupsummary(labelsTable,'labels');
[~,maxgroupInd] = max(groupsum{:,'GroupCount'});
plane1Label = groupsum{maxgroupInd,'labels'};

% get plane only
p2Location = plane1.Location;
p2Location(labels~=plane1Label, :) = [];
plane2 = pointCloud(p2Location);

if printLevel > 1
    figure;
    pcshow(plane2);
    title("Plane only according to clusters")
end

% get center
xMean = 0.5*(plane2.XLimits(1)+plane2.XLimits(2));
yMean = 0.5*(plane2.YLimits(1)+plane2.YLimits(2));
zMean = 0.5*(plane2.ZLimits(1)+plane2.ZLimits(2));

% cut edges
if printLevel > 1
    uv=null(model.Normal(:).');
    figure;
    pcshow(plane2);
    hold on;
    tempN = 0.03*model.Normal;
    tempUv = 0.03*uv;
    quiver3(xMean,yMean,zMean, tempN(1),tempN(2),tempN(3));
    quiver3(xMean,yMean,zMean, tempUv(1,1),tempUv(2,1),tempUv(3,1));
    quiver3(xMean,yMean,zMean, tempUv(1,2),tempUv(2,2),tempUv(3,2));
    title("plane with normal and nullspace")
end
% pca
[coeff,score,~,~,~,mu] = pca(plane2.Location);

u = coeff(:,1);
v = coeff(:,2);
if printLevel > 1
    figure;
    pcshow(plane2);
    hold on;
    tempN = 0.03*model.Normal;
    tempU = 0.03*u;
    tempV = 0.03*v;
    quiver3(xMean,yMean,zMean, tempN(1),tempN(2),tempN(3));
    quiver3(xMean,yMean,zMean, tempU(1),tempU(2),tempU(3));
    quiver3(xMean,yMean,zMean, tempV(1),tempV(2),tempV(3));
    title("plane with normal and basis of pca");
end
% map to 2d
%uvPoints = plane2.Location*[u,v];
uvPoints = score(:,1:2);
plane2uv = zeros(size(uvPoints,1),3);
plane2uv(:,1:2) = uvPoints;
plane2uv = pointCloud(plane2uv);
if printLevel > 1
    figure;
    pcshow(plane2uv);
    zlim([-0.01,0.01]);
    title("plane projected into 2d space")
end

% regular grid
xGrid = plane2uv.XLimits(1):cb3dParams.Resolution:plane2uv.XLimits(2);
yGrid = plane2uv.YLimits(1):cb3dParams.Resolution:plane2uv.YLimits(2);
[xGrid,yGrid] = meshgrid(xGrid,yGrid);
zGrid = nan(size(xGrid));
for yi=1:size(xGrid,1)
   for xi=1:size(xGrid,2)
       [~,dists] = findNearestNeighbors(plane2uv, [xGrid(yi,xi),yGrid(yi,xi),0],1);
       if dists<cb3dParams.Resolution
           zGrid(yi,xi)=1;
       end
   end
end

if printLevel > 1
    figure;
    surf(xGrid,yGrid,zGrid);
    title("rasterized plane (grid)");
end

% align properly
zGridInv = zGrid;
zGridInv(isnan(zGrid)) = 1;
zGridInv(~isnan(zGrid)) = 0;
xInv = xGrid(:);
yInv = yGrid(:);

CC = bwconncomp(zGridInv);
ccCols = zeros(length(zGridInv(:)),1);
numOfObj = zeros(CC.NumObjects,1);
for i=1:CC.NumObjects
    numOfObj(i) = length(CC.PixelIdxList{i});
end

[~,si] = sort(numOfObj, 'descend');
holeXMedian = zeros(cb3dParams.NumOfHoles,1);
holeYMedian = zeros(cb3dParams.NumOfHoles,1);
ri=1;
for i=1:CC.NumObjects
    if ri == cb3dParams.NumOfHoles+1
        break;
    end
    tempXMedian = median(xInv(CC.PixelIdxList{si(i)}));
    tempYMedian = median(yInv(CC.PixelIdxList{si(i)}));
    tempXMin = min(xInv(CC.PixelIdxList{si(i)}));
    tempYMin = min(yInv(CC.PixelIdxList{si(i)}));
    tempXMax = max(xInv(CC.PixelIdxList{si(i)}));
    tempYMax = max(yInv(CC.PixelIdxList{si(i)}));
    
    % ignore border (either close to border or not square)
    if abs(tempXMedian-xGrid(1))<0.01 || abs(tempXMedian-xGrid(end))<0.01 ...
            || abs(tempYMedian-yGrid(1))<0.01 || abs(tempYMedian-yGrid(end))<0.01
        continue;
    end
    if abs((tempYMax-tempYMin)/(tempXMax-tempXMin)-1)>holecondition
        continue;
    end
    
    ccCols(CC.PixelIdxList{si(i)}) = ri;
    holeXMedian(ri) = tempXMedian;
    holeYMedian(ri) = tempYMedian;
    ri=ri+1;
end

% if colums are uneven, use less holes
if ~isempty(cb3dParams.LeftColsVisible) && cb3dParams.LeftColsVisible > 0
    leftColsRemainder = cb3dParams.LeftColsVisible-floor(cb3dParams.LeftColsVisible);
    if leftColsRemainder>0.1
        % remove last column
        [~,six] = sort(holeXMedian);
        six(1:cb3dParams.NumOfRows) = []; % TODO this is wrong, should be holes per Column
        holeXMedian = holeXMedian(six);
        holeYMedian = holeYMedian(six);
        cb3dParams.NumOfCols = floor(cb3dParams.LeftColsVisible);
        cb3dParams.NumOfHoles = cb3dParams.NumOfHoles-cb3dParams.NumOfRows; % TODO this is wrong, should be holes per Column
    end
    while cb3dParams.NumOfHoles > 2*cb3dParams.NumOfRows*cb3dParams.NumOfCols
        % remove last column
        [~,six] = sort(holeXMedian);
        six(1:cb3dParams.NumOfRows) = []; % TODO this is wrong, should be holes per Column
        holeXMedian = holeXMedian(six);
        holeYMedian = holeYMedian(six);
        cb3dParams.NumOfCols = floor(cb3dParams.LeftColsVisible);
        cb3dParams.NumOfHoles = cb3dParams.NumOfHoles-cb3dParams.NumOfRows; % TODO this is wrong, should be holes per Column
    end
end
% if colums are uneven, use less holes
if ~isempty(cb3dParams.RightColsVisible) && cb3dParams.RightColsVisible > 0
    rightColsRemainder = cb3dParams.RightColsVisible-floor(cb3dParams.RightColsVisible);
    if rightColsRemainder>0.1
        % remove last column
        [~,six] = sort(holeXMedian, 'descend');
        six(1:cb3dParams.NumOfRows) = []; % TODO this is wrong, should be holes per Column
        holeXMedian = holeXMedian(six);
        holeYMedian = holeYMedian(six);
        cb3dParams.NumOfCols = floor(cb3dParams.RightColsVisible);
        cb3dParams.NumOfHoles = cb3dParams.NumOfHoles-cb3dParams.NumOfRows; % TODO this is wrong, should be holes per Column
    end
    while cb3dParams.NumOfHoles > 2*cb3dParams.NumOfRows*cb3dParams.NumOfCols
        % remove last column
        [~,six] = sort(holeXMedian, 'descend');
        six(1:cb3dParams.NumOfRows) = []; % TODO this is wrong, should be holes per Column
        holeXMedian = holeXMedian(six);
        holeYMedian = holeYMedian(six);
        cb3dParams.NumOfCols = floor(cb3dParams.RightColsVisible);
        cb3dParams.NumOfHoles = cb3dParams.NumOfHoles-cb3dParams.NumOfRows; % TODO this is wrong, should be holes per Column
    end
end
% if rows are uneven, use less holes
if ~isempty(cb3dParams.TopRowsVisible) && cb3dParams.TopRowsVisible > 0
    topRowsRemainder = cb3dParams.TopRowsVisible-floor(cb3dParams.TopRowsVisible);
    if topRowsRemainder>0.1
        % remove last row
        [~,six] = sort(holeYMedian);
        six(1:cb3dParams.NumOfCols) = []; % TODO this is wrong, should be holes per Row
        holeXMedian = holeXMedian(six);
        holeYMedian = holeYMedian(six);
        cb3dParams.NumOfRows = floor(cb3dParams.TopRowsVisible);
        cb3dParams.NumOfHoles = cb3dParams.NumOfHoles-cb3dParams.NumOfCols; % TODO this is wrong, should be holes per Row
    end
    while cb3dParams.NumOfHoles > 2*cb3dParams.NumOfRows*cb3dParams.NumOfCols
        % remove last column
        [~,six] = sort(holeYMedian);
        six(1:cb3dParams.NumOfCols) = []; % TODO this is wrong, should be holes per Row
        holeXMedian = holeXMedian(six);
        holeYMedian = holeYMedian(six);
        cb3dParams.NumOfRows = floor(cb3dParams.TopRowsVisible);
        cb3dParams.NumOfHoles = cb3dParams.NumOfHoles-cb3dParams.NumOfCols; % TODO this is wrong, should be holes per Row
    end
end
% if rows are uneven, use less holes
if ~isempty(cb3dParams.BottomRowsVisible) && cb3dParams.BottomRowsVisible > 0
    bottomRowsRemainder = cb3dParams.BottomRowsVisible-floor(cb3dParams.BottomRowsVisible);
    if bottomRowsRemainder>0.1
        % remove last row
        [~,six] = sort(holeYMedian, 'descend');
        six(1:cb3dParams.NumOfCols) = []; % TODO this is wrong, should be holes per Row
        holeXMedian = holeXMedian(six);
        holeYMedian = holeYMedian(six);
        cb3dParams.NumOfRows = floor(cb3dParams.BottomRowsVisible);
        cb3dParams.NumOfHoles = cb3dParams.NumOfHoles-cb3dParams.NumOfCols; % TODO this is wrong, should be holes per Row
    end
    while cb3dParams.NumOfHoles > 2*cb3dParams.NumOfRows*cb3dParams.NumOfCols
        % remove last column
        [~,six] = sort(holeYMedian, 'descend');
        six(1:cb3dParams.NumOfCols) = []; % TODO this is wrong, should be holes per Row
        holeXMedian = holeXMedian(six);
        holeYMedian = holeYMedian(six);
        cb3dParams.NumOfRows = floor(cb3dParams.BottomRowsVisible);
        cb3dParams.NumOfHoles = cb3dParams.NumOfHoles-cb3dParams.NumOfCols; % TODO this is wrong, should be holes per Row
    end
end

if printLevel > 1
    figure;
    pcshow([xGrid(:),yGrid(:),0.1*ccCols/CC.NumObjects], ccCols);
    title("grid with detected holes");
    view(0,60)
end
pause(0.5);

if printLevel > 0
    figure;
    pcshow([xGrid(:),yGrid(:),zGrid(:)]);
    hold on;
    pcshow([holeXMedian, holeYMedian, ones(length(holeXMedian),1)], 'r', 'MarkerSize', 24);
    zlim([1-0.01,1+0.01]);
    title("grid with median of holes")
end

%% fit lines through center of holes
[~,six] = sort(holeYMedian);
sortedHoleYMedian1 = holeYMedian(six);
sortedHoleXMedian1 = holeXMedian(six);
ci=1;
j=1;
for i=1:cb3dParams.NumOfCols:cb3dParams.NumOfHoles %1:5:30
    if j+(cb3dParams.NumOfCols-1)>cb3dParams.NumOfHoles
        break;
    end
    addH = 0;
    if cb3dParams.SkipEven == 1 && mod(ci,2)==0
       addH = 1;
    end
    if cb3dParams.SkipOdd == 1 && mod(ci,2)==1
       addH = 1;
    end
    psx{ci} = polyfit(sortedHoleXMedian1(j:j+(cb3dParams.NumOfCols-1+addH))',...
        sortedHoleYMedian1(j:j+(cb3dParams.NumOfCols-1+addH))',1);
    ci=ci+1;
    if addH==1
        j=j+1;
    end
    j = j+cb3dParams.NumOfCols;
end

%% plot the lines
if printLevel > 1
    figure;
    pcshow([xGrid(:),yGrid(:),zGrid(:)]*1000, 'MarkerSize', 24);
    hold on;
    pcshow([holeXMedian, holeYMedian, ones(length(holeXMedian),1)]*1000, 'r', 'MarkerSize', 36);
    zlim([1000-0.01,1000+0.01]);
    ci=1;
    j=1;
    for i=1:cb3dParams.NumOfCols:cb3dParams.NumOfHoles %1:5:30
        if j+(cb3dParams.NumOfCols-1)>cb3dParams.NumOfHoles
            break;
        end
        addH = 0;
        if cb3dParams.SkipEven == 1 && mod(ci,2)==0
           addH = 1;
        end
        if cb3dParams.SkipOdd == 1 && mod(ci,2)==1
           addH = 1;
        end
        % arrow from smallest X value to largest X value

        smallestInd=j;
        smallestXVal = sortedHoleXMedian1(j);
        largestInd=j;
        largestXVal = sortedHoleXMedian1(j+cb3dParams.NumOfCols-1+addH);
        for smi=j:j+cb3dParams.NumOfCols-1
            if sortedHoleXMedian1(smi)<smallestXVal
                smallestInd=smi;
                smallestXVal = sortedHoleXMedian1(smi);
            end
            if sortedHoleXMedian1(smi)>largestXVal
                largestInd=smi;
                largestXVal = sortedHoleXMedian1(smi);
            end
        end
        arrow = 1.1*diff([sortedHoleXMedian1([smallestInd,largestInd]),...
            polyval(psx{ci},sortedHoleXMedian1([smallestInd,largestInd]))]);
        quiver3(sortedHoleXMedian1(smallestInd)*1000,sortedHoleYMedian1(smallestInd)*1000,1000,arrow(1)*1000,arrow(2)*1000,0);
        ci=ci+1;
        if addH==1
            j=j+1;
        end
        j = j+cb3dParams.NumOfCols;
    end
    view(0,-90)
    ylabel("y [mm]")
    xlabel("x [mm]")
    title("grid with horizontal lines through rows")
end

%% find the middle lines between the holes (edges)
if printLevel > 2
    figure;
    pcshow([xGrid(:),yGrid(:),zGrid(:)]);
    hold on;
    pcshow([holeXMedian, holeYMedian, ones(length(holeXMedian),1)], 'r', 'MarkerSize', 24);
    zlim([1-0.01,1+0.01]);
end
ci=1;
j=1;
hAvgsx1 = zeros(length(1:cb3dParams.NumOfCols:(cb3dParams.NumOfHoles-cb3dParams.NumOfCols)),1);
hAvgsx2 = zeros(length(1:cb3dParams.NumOfCols:(cb3dParams.NumOfHoles-cb3dParams.NumOfCols)),1);
hAvgsy1 = zeros(length(1:cb3dParams.NumOfCols:(cb3dParams.NumOfHoles-cb3dParams.NumOfCols)),1);
hAvgsy2 = zeros(length(1:cb3dParams.NumOfCols:(cb3dParams.NumOfHoles-cb3dParams.NumOfCols)),1);
for i=1:cb3dParams.NumOfCols:(cb3dParams.NumOfHoles-cb3dParams.NumOfCols) %1:5:(30-5)
    if j+(2*cb3dParams.NumOfCols-1)>cb3dParams.NumOfHoles
        break;
    end
    addH = 0;
    if cb3dParams.SkipEven == 1 && mod(ci,2)==0
       addH = 1;
    end
    if cb3dParams.SkipOdd == 1 && mod(ci,2)==1
       addH = 1;
    end
    xMin1 = min(sortedHoleXMedian1(j:j+(cb3dParams.NumOfCols-1+addH)));
    xMax1 = max(sortedHoleXMedian1(j:j+(cb3dParams.NumOfCols-1+addH)));
    xMin2 = min(sortedHoleXMedian1(j+cb3dParams.NumOfCols+addH:j+(2*cb3dParams.NumOfCols-1+addH)));
    xMax2 = max(sortedHoleXMedian1(j+cb3dParams.NumOfCols+addH:j+(2*cb3dParams.NumOfCols-1+addH)));
    
    xAvg1 = (xMin1+ xMin2)/2;
    yAvg1 = 0.5*(polyval(psx{ci},xAvg1)+polyval(psx{ci+1},xAvg1));    
    xAvg2 = (xMax1+xMax2)/2;
    yAvg2 = 0.5*(polyval(psx{ci},xAvg2)+polyval(psx{ci+1},xAvg2));
    hAvgsx1(ci) = xAvg1;
    hAvgsx2(ci) = xAvg2;
    hAvgsy1(ci) = yAvg1;
    hAvgsy2(ci) = yAvg2;
    if printLevel > 1
        line([xAvg1,xAvg2],[yAvg1,yAvg2],[1,1],'Color','r');
    end
    % make another line model
    psdx{ci} = polyfit([xAvg1,xAvg2],[yAvg1,yAvg2],1);
    ci=ci+1;
    if addH==1
       j=j+1; 
    end
    j = j+cb3dParams.NumOfCols;
end
if printLevel > 1
    title("grid with horizontal lines for checkerboard");
end

%% fit lines through center of holes
[~,six] = sort(holeXMedian);
sortedHoleXMedian2 = holeXMedian(six);
sortedHoleYMedian2 = holeYMedian(six);
ci=1;
j=1;
for i=1:cb3dParams.NumOfRows:cb3dParams.NumOfHoles %1:3:30
    if j+(cb3dParams.NumOfRows-1)>cb3dParams.NumOfHoles
        break;
    end
    psy{ci} = polyfit(sortedHoleYMedian2(j:j+(cb3dParams.NumOfRows-1))',sortedHoleXMedian2(j:j+(cb3dParams.NumOfRows-1))',1);
    ci=ci+1;
    j = j+cb3dParams.NumOfRows;
end

%% plot the lines
if printLevel > 1
    figure;
    pcshow([xGrid(:),yGrid(:),zGrid(:)]);
    hold on;
    pcshow([holeXMedian, holeYMedian, ones(length(holeXMedian),1)], 'r', 'MarkerSize', 24);
    zlim([1-0.01,1+0.01]);
    ci=1;
    for i=1:cb3dParams.NumOfRows:cb3dParams.NumOfHoles %1:3:30
        if i+(cb3dParams.NumOfRows-1)>cb3dParams.NumOfHoles
            break;
        end
        arrow = 1.1*diff([polyval(psy{ci},sortedHoleYMedian2([i,i+(cb3dParams.NumOfRows-1)])), sortedHoleYMedian2([i,i+(cb3dParams.NumOfRows-1)])]);
        quiver3(sortedHoleXMedian2(i),sortedHoleYMedian2(i),1,arrow(1),arrow(2),0);
        ci=ci+1;
    end
    title("grid with vertical lines through columns");
end

%% find the middle lines between the holes (edges)
if printLevel > 2
    figure;
    pcshow([xGrid(:),yGrid(:),zGrid(:)]);
    hold on;
    pcshow([holeXMedian, holeYMedian, ones(length(holeXMedian),1)], 'r', 'MarkerSize', 24);
    zlim([1-0.01,1+0.01]);
end
ci=1;
vAvgsx1 = zeros(length(1:cb3dParams.NumOfRows:(cb3dParams.NumOfHoles-cb3dParams.NumOfRows)),1);
vAvgsx2 = zeros(length(1:cb3dParams.NumOfRows:(cb3dParams.NumOfHoles-cb3dParams.NumOfRows)),1);
vAvgsy1 = zeros(length(1:cb3dParams.NumOfRows:(cb3dParams.NumOfHoles-cb3dParams.NumOfRows)),1);
vAvgsy2 = zeros(length(1:cb3dParams.NumOfRows:(cb3dParams.NumOfHoles-cb3dParams.NumOfRows)),1);
for i=1:cb3dParams.NumOfRows:(cb3dParams.NumOfHoles-cb3dParams.NumOfRows) %1:3:(30-3)
    xMin1 = min(sortedHoleYMedian2(i:i+(cb3dParams.NumOfRows-1)));
    xMax1 = max(sortedHoleYMedian2(i:i+(cb3dParams.NumOfRows-1)));
    xMin2 = min(sortedHoleYMedian2(i+cb3dParams.NumOfRows:i+(2*cb3dParams.NumOfRows-1)));
    xMax2 = max(sortedHoleYMedian2(i+cb3dParams.NumOfRows:i+(2*cb3dParams.NumOfRows-1)));
    
    xAvg1 = (xMin1+ xMin2)/2;
    yAvg1 = 0.5*(polyval(psy{ci},xAvg1)+polyval(psy{ci+1},xAvg1));    
    xAvg2 = (xMax1+xMax2)/2;
    yAvg2 = 0.5*(polyval(psy{ci},xAvg2)+polyval(psy{ci+1},xAvg2));
    vAvgsx1(ci) = yAvg1;
    vAvgsx2(ci) = yAvg2;
    vAvgsy1(ci) = xAvg1;
    vAvgsy2(ci) = xAvg2;
    if printLevel > 2
        line([yAvg1,yAvg2],[xAvg1,xAvg2],[1,1],'Color','r');
    end
    % make another line model
    psdy{ci} = polyfit([xAvg1,xAvg2],[yAvg1,yAvg2],1);
    ci=ci+1;
end
if printLevel > 2
    line([hAvgsx1,hAvgsx2],[hAvgsy1,hAvgsy2],[1,1],'Color','r');
    title("grid with horizontal and vertical lines for checkerboard");
end

%% expand borders

% get average slope of linear models and average distance
for ci=1:length(psdx)
    xSlopes(ci) = psdx{ci}(1);
end
xSlope = mean(xSlopes);
for ci=1:length(psdx)-1
    xDiffs(ci) = psdx{ci+1}(2)-psdx{ci}(2);
end
xDiff = mean(xDiffs);
for ci=1:length(psdy)
    ySlopes(ci) = psdy{ci}(1);
end
ySlope = mean(ySlopes);
for ci=1:length(psdy)-1
    yDiffs(ci) = psdy{ci+1}(2)-psdy{ci}(2);
end
yDiff = mean(yDiffs);

% expand one step in all directions
for ci=1:length(psdx)
    psdxt{ci+1} = psdx{ci};
end
for ci=1:length(psdy)
    psdyt{ci+1} = psdy{ci};
end
psdyt{length(psdy)+2} = [];

psdxt{1} = [xSlope,psdxt{2}(2)-xDiff];
psdxt{length(psdx)+2}=[xSlope,psdxt{length(psdx)+1}(2)+xDiff];
psdyt{1} = [ySlope,psdyt{2}(2)-yDiff];
psdyt{length(psdy)+2}=[ySlope,psdyt{length(psdy)+1}(2)+yDiff];

psdx = psdxt;
psdy = psdyt;

%% fill with colors, white and black
% sort median first y, then x

holesLimitPts = [];
holesPts = [];

[~,six] = sort(holeYMedian);
sortedHoleXMedian1 = holeXMedian(six);
sortedHoleYMedian1 = holeYMedian(six);
sortedHoleXMedian2 = zeros(cb3dParams.NumOfHoles,1);
sortedHoleYMedian2 = zeros(cb3dParams.NumOfHoles,1);
yj=1;
ci=1;
for yi=1:cb3dParams.NumOfCols:cb3dParams.NumOfHoles % 1:5:30
    if yj+(cb3dParams.NumOfCols-1)>cb3dParams.NumOfHoles
        break;
    end
    addH = 0;
    if cb3dParams.SkipEven == 1 && mod(ci,2)==0
       addH = 1;
    end
    if cb3dParams.SkipOdd == 1 && mod(ci,2)==1
       addH = 1;
    end
    sortedHoleXMedianTemp = sortedHoleXMedian1(yj:yj+(cb3dParams.NumOfCols-1+addH));
    sortedHoleYMedianTemp = sortedHoleYMedian1(yj:yj+(cb3dParams.NumOfCols-1+addH));
    [~,six] = sort(sortedHoleXMedianTemp);
    sortedHoleXMedian2(yj:yj+(cb3dParams.NumOfCols-1+addH)) = sortedHoleXMedianTemp(six);
    sortedHoleYMedian2(yj:yj+(cb3dParams.NumOfCols-1+addH)) = sortedHoleYMedianTemp(six);
    if addH==1
        yj=yj+1;
    end
    yj=yj+cb3dParams.NumOfCols;
    ci=ci+1;
end

ci=1;
for yi=1:(2*cb3dParams.NumOfRows) % 3
    addH = 0;
    if cb3dParams.SkipEven == 1 && mod(yi,2)==0
       addH = 1;
    end
    if cb3dParams.SkipOdd == 1 && mod(yi,2)==1
       addH = 1;
    end
    for xi=1:cb3dParams.NumOfCols+addH % 5
        if printLevel > 2
            hold on;
            pcshow([sortedHoleXMedian2(ci),sortedHoleYMedian2(ci),1],'g','MarkerSize',24);
            zlim([1-0.01,1+0.01]);
        end
        % get borders
        foundX = 0;
        foundY = 0;
        x1=0;
        x2=0;
        numOfCol2 =(2*cb3dParams.NumOfCols); % 5
        if cb3dParams.SkipEven ==1 || cb3dParams.SkipOdd==1
            numOfCol2 = numOfCol2+1;
        end
        
        for xii=1:numOfCol2
            x1 = polyval(psdy{xii},sortedHoleYMedian2(ci));
            x2 = polyval(psdy{xii+1},sortedHoleYMedian2(ci));
            if sortedHoleXMedian2(ci) > x1 && sortedHoleXMedian2(ci) < x2
                if printLevel > 2
                    hold on;
                    pcshow([x1,sortedHoleYMedian2(ci),1],'b','MarkerSize',24);
                    pcshow([x2,sortedHoleYMedian2(ci),1],'b','MarkerSize',24);
                    zlim([1-0.01,1+0.01]);
                end
                foundX=1;
                break;
            end
        end
        y1=0;
        y2=0;
        for yii=1:(2*cb3dParams.NumOfRows) % 3
            y1 = polyval(psdx{yii},sortedHoleXMedian2(ci));
            y2 = polyval(psdx{yii+1},sortedHoleXMedian2(ci));
            if sortedHoleYMedian2(ci) > y1 && sortedHoleYMedian2(ci) < y2
                if printLevel > 2
                    hold on;
                    pcshow([sortedHoleXMedian2(ci),y1,1],'b','MarkerSize',24);
                    pcshow([sortedHoleXMedian2(ci),y2,1],'b','MarkerSize',24);
                    zlim([1-0.01,1+0.01]);
                end
                foundY=1;
                break;
            end
        end
        if foundX && foundY
            [pX,pY] = meshgrid(x1:cb3dParams.Resolution:x2,y1:cb3dParams.Resolution:y2);
            % append for later
            holesLimitPts = [holesLimitPts; x1,sortedHoleYMedian2(ci); x2,sortedHoleYMedian2(ci); [sortedHoleXMedian2(ci),y1]; [sortedHoleXMedian2(ci),y2];];
            holesPts = [holesPts; [pX(:),pY(:)]];
            pZ = ones(size(pX));
            if printLevel > 2
                hold on;
                pcshow([pX(:),pY(:),pZ(:)], 'r', 'MarkerSize', 24);
                zlim([1-0.01,1+0.01]);
            end
        end
        ci=ci+1;
    end
end

%% bring back to original space and place
% first sparse
holeMedians = [sortedHoleXMedian2 sortedHoleYMedian2];
holeMedians = [holeMedians(:,1), holeMedians(:,2), zeros(size(holeMedians(:,1)))];
holeMedRec1 = holeMedians * coeff';
holeMedians = holeMedRec1 + repmat(mu, length(holeMedRec1), 1);

sparseRecon = holesLimitPts;
sparseRecon = [sparseRecon(:,1), sparseRecon(:,2), zeros(size(sparseRecon(:,1)))];
sparseRec1 = sparseRecon * coeff';
sparseRecon = sparseRec1 + repmat(mu, length(sparseRec1), 1);

% then full
scores2d = [holesPts(:,1), holesPts(:,2), zeros(size(holesPts(:,1)))];
rec1 = scores2d * coeff';
reconstructed = rec1 + repmat(mu, length(rec1), 1);
if printLevel > 1
    figure;
    pcshow(plane2);
    hold on;
    pcshow(reconstructed, 'w', 'MarkerSize', 24);
    title("checkerboard back in 3d space");
end

end