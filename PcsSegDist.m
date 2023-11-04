function [pc1, pc2] = PcsSegDist(pc1, pc2)

% segment point clouds (in case of outliers)
[labels1, numClust1] = pcsegdist(pc1, 0.05); % 5cm
[labels2, numClust2] = pcsegdist(pc2, 0.05); % 5cm
% only keep the biggest one
maxNumClust1 = 0;
maxNumClustInd1 = 0;
for nci=1:numClust1
    if sum(labels1==nci)>maxNumClust1
        maxNumClust1 = sum(labels1==nci);
        maxNumClustInd1 = nci;
    end
end
pcLoc = pc1.Location(labels1==maxNumClustInd1, :);
if isempty(pc1.Color)
    pc1 = pointCloud(pcLoc);
else
    pcCol = pc1.Color(labels1==maxNumClustInd1, :);
    pc1 = pointCloud(pcLoc, 'Color', pcCol);
end
maxNumClust2 = 0;
maxNumClustInd2 = 0;
for nci=1:numClust2
    if sum(labels2==nci)>maxNumClust2
        maxNumClust2 = sum(labels2==nci);
        maxNumClustInd2 = nci;
    end
end
pcLoc = pc2.Location(labels2==maxNumClustInd2, :);
if isempty(pc2.Color)
    pc2 = pointCloud(pcLoc);
else
    pcCol = pc2.Color(labels2==maxNumClustInd2, :);
    pc2 = pointCloud(pcLoc, 'Color', pcCol);
end

end