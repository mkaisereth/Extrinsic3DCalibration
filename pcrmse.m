function [closeCount, rmseOut, stdOut, dirValue] = pcrmse(pc1,pc2, thresh, downsampleLimit, gridStep, dirVector)

if ~isempty(downsampleLimit)
    pc1Size = size(pc1.Location,1);
    if pc1Size>downsampleLimit
        pc1 = pcdownsample(pc1, 'gridAverage', gridStep);
        %disp("pc1 down: " + string(100*size(pc1.Location,1)/pc1Size))
    end
    pc2Size = size(pc2.Location,1);
    if pc2Size>downsampleLimit
        pc2 = pcdownsample(pc2, 'gridAverage', gridStep);
        %disp("pc2 down: " + string(100*size(pc2.Location,1)/pc2Size))
    end
end

distSum = 0;
dists = [];
countSum = 0;
dirValues = [];
for pc1Ind=1:size(pc1.Location,1)
   [ind, dist] = findNearestNeighbors(pc2, pc1.Location(pc1Ind, :), 1);
   dists = [dists; dist];
   distSum = distSum+dist;
   if dist<thresh
       countSum = countSum+1;
       % get the distance w.r.t. dirVector
       if ~isempty(dirVector)
           dotProduct = dot(pc2.Location(ind, :)-pc1.Location(pc1Ind, :), dirVector);
           dirValues = [dirValues, dotProduct];
       end
   end
end
dirValue = mean(dirValues);
for pc2Ind=1:size(pc2.Location,1)
   [~, dist] = findNearestNeighbors(pc1, pc2.Location(pc2Ind, :), 1);
   dists = [dists; dist];
   distSum = distSum+dist;
   if dist<thresh
       countSum = countSum+1;
   end
end
stdOut = std(dists);
rmseOut = distSum/(size(pc1.Location,1)+size(pc2.Location,1));
closeCount = countSum/(size(pc1.Location,1)+size(pc2.Location,1));