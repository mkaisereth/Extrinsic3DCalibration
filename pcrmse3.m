function [rmseOut, stdOut, tForm] = pcrmse3(pc1,pc2, deltadistRejection, doUseIcp)

global printLevel;

rejectionCounter = 0;
%tic
normals1 = pcnormals(pc1, 21);
[idx1,~] = knnsearch(pc2.Location, pc1.Location);

normals2 = pcnormals(pc2, 21);
[idx2,~] = knnsearch(pc1.Location, pc2.Location);

dists = nan(size(pc1.Location,1)+size(pc2.Location,1), 1);
remainingPc1Pts = nan(size(pc1.Location,1),3);
keepCounter = 1;
for pc1Ind=1:size(pc1.Location,1)
   % it only counts if more or less orthogonal to surface (meaning nearest
   % neighbor of nearest neighbor must be very close
   P1 = normals1(pc1Ind,:);
   P2 = pc2.Location(idx1(pc1Ind), :)-pc1.Location(pc1Ind,:);
   sineAngleP2 = norm(cross(P1,P2))/norm(P1);
   distFromNormal = sineAngleP2;
   if abs(distFromNormal)<deltadistRejection
       dist = norm(pc1.Location(pc1Ind, :)-pc2.Location(idx1(pc1Ind), :));
       dists(keepCounter) =dist;
       remainingPc1Pts(keepCounter, :) = pc1.Location(pc1Ind, :);
       keepCounter=keepCounter+1;
   else
       rejectionCounter = rejectionCounter+1;
   end
end

keepPts1 = keepCounter-1;
remainingPc2Pts = nan(size(pc2.Location,1),3);
for pc2Ind=1:size(pc2.Location,1)
   % it only counts if more or less orthogonal to surface (meaning nearest
   % neighbor of nearest neighbor must be very close
   P1 = normals2(pc2Ind,:);
   P2 = pc1.Location(idx2(pc2Ind), :)-pc2.Location(pc2Ind,:);
   sineAngleP2 = norm(cross(P1,P2))/norm(P1);
   distFromNormal = sineAngleP2;
   if abs(distFromNormal)<deltadistRejection
       dist = norm(pc2.Location(pc2Ind, :)-pc1.Location(idx2(pc2Ind), :));
       dists(keepCounter) =dist;
       remainingPc2Pts(keepCounter-keepPts1, :) = pc2.Location(pc2Ind, :);
       keepCounter=keepCounter+1;
   else
       rejectionCounter = rejectionCounter+1;
   end
end
%toc

if printLevel > 1
    figure(2);
    hold off;
    pcshow(pc1);
    hold on;
    pcshow(pc2);
    pcshow(remainingPc1Pts, 'r');
    pcshow(remainingPc2Pts, 'g');
end

stdOut = std(dists, 'omitnan');
rmseOut = sqrt(mean(dists.^2, 'omitnan'));

if doUseIcp
    % do an iterative closest point
    [tForm, ~, rmse] = pcregistericp(pointCloud(remainingPc2Pts), pointCloud(remainingPc1Pts));
    % do it again for remaining points
    remainingPc2Pts = pctransform(pointCloud(remainingPc2Pts), tForm).Location;
    pc2_Z = remainingPc2Pts;
    %pc2_Z(:,3) = 0;
    pc2_Z = pointCloud(pc2_Z);
    pc1_Z = remainingPc1Pts;
    %pc1_Z(:,3) = 0;
    pc1_Z = pointCloud(pc1_Z);
    
    %tic
    [~,dists1] = knnsearch(pc2_Z.Location, pc1_Z.Location);
    [~,dists2] = knnsearch(pc1_Z.Location, pc2_Z.Location);
    dists = [dists1;dists2];
    %toc
else
    tForm = [];
end

if printLevel > 1
    figure(3);
    hold off;
    pcshow(pc1);
    hold on;
    pcshow(pc2);
    pcshow(remainingPc1Pts, 'r');
    pcshow(remainingPc2Pts, 'g');
end

stdOut = std(dists, 'omitnan');
rmseOut = sqrt(mean(dists.^2, 'omitnan'));
%disp("rejectionrate: " + string(rejectionCounter/(size(pc1.Location,1)+size(pc2.Location,1))));