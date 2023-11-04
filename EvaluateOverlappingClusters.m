function [rmseVals,stdVals,icptForms, closeRatios] = EvaluateOverlappingClusters(tmpPcst,tmpPcs2t, patientFolderPath, combination, saveFigures, saveJson, distRejection_t, doUseIcp, doUseUserInput)


global printLevel;
global pcst;
pcst = tmpPcst;
global pcs2t;
pcs2t = tmpPcs2t;
global pcClustLoc1;
global pcClustLoc2;
global pcInd;
global pcInd2;
global distRejection;
distRejection = distRejection_t;

rmseVals = [];
stdVals = [];
closeRatios = [];
for i=1:length(pcst)
    [labels1, labels2, maxClusterNum1, maxClusterNum2, maxLabelInd1, maxLabelInd2] = DoCluster(pcst(i), pcs2t(i));
    % min count of cluster
    %disp(append(string(maxClusterNum1), " ", string(maxClusterNum2)))
    if maxClusterNum2<10000 || maxClusterNum1<10000
        % TODO
        disp("TODO")
    else
        [pcClust1, pcClust2, pcClustLoc1, pcClustLoc2] = GetClusters(pcst(i), pcs2t(i), labels1, maxLabelInd1, labels2, maxLabelInd2);
        [closeRatio, ~] = pcrmse(pcClust1, pcClust2, 0.05, 10000, 0.005, []); % 5cm threshold, 5mm downsampling
        closeRatios(i) = closeRatio;
        % color them
        if printLevel > 0
            f1 = figure(1);
            hold off;
            pcInd = i;
            pcInd2 = i;
            pcshow(pcst(pcInd), 'MarkerSize', 6);
            hold on;
            pcshow(pcs2t(pcInd2), 'MarkerSize', 6);
            pcshow(pcClustLoc1, 'r', 'MarkerSize', 12);
            pcshow(pcClustLoc2, 'g', 'MarkerSize', 12);
            view(0,-89.9);
            camorbit(-89,0,'camera');
            drawnow();
        end
        
        % evaluate only if ... are closer than threshold
        if closeRatio > 0.5 % 30%
            % do evaluation for overlapping region
            [rmseVal2, stdVal2, tForm] = EvalOverlappingRegionFunc(pcClust1, pcClust2, distRejection, doUseIcp);

            if doUseIcp
                icptForms(i) = tForm;
            end
            
            if printLevel > 0
                if doUseUserInput
                    global diff1;
                    diff1 = 0;
                    global diff2;
                    diff2 = 0;
                    disp(append("rmse [mm]: ", string(rmseVal2*1000), " std[mm]: ", string(stdVal2*1000)))
                    addpath("../XrayRegistration/");
                    rotate3d off;
                    set(f1,'KeyPressFcn',@figKeyPress);
                    UpdateClusterErrors(pcst(pcInd),pcs2t(pcInd2));
                    title(append("P1: ", string(diff1), " P2: ", string(diff2)))
                    WaitForSpace();
                end
            end

            rmseVals = [rmseVals; rmseVal2];
            stdVals = [stdVals; stdVal2];
        else
            icptForms(i) = rigid3d;
        end
    end
end

% plot eval with Boxplots
outputPath = append(patientFolderPath, "\Output");
if ~exist(outputPath, "dir")
    mkdir(outputPath)
end
if printLevel > 1
    f2 = figure(4);
    boxplot(rmseVals*1000, 'Labels', combination);
    title("Z Mean difference [mm]")
    ylabel("[mm]")
    if saveFigures
        saveas(f2, append(outputPath, "\ZMeanDiff_", combination, ".fig"));
        saveas(f2, append(outputPath, "\ZMeanDiff_", combination, ".png"));
    end
    f3 = figure(5);
    boxplot(stdVals*1000, 'Labels', combination);
    title("Z Std [mm]")
    ylabel("[mm]")
    if saveFigures
        saveas(f3, append(outputPath, "\ZStd_", combination, ".fig"));
        saveas(f3, append(outputPath, "\ZStd_", combination, ".png"));
    end
end
if saveJson
    clear siOutput;
    siOutput.ZMeanMM = rmseVals*1000;
    siOutput.ZStdMM = stdVals*1000;
    
    jsonOut = jsonencode(siOutput);
    filepath = append(outputPath, "\EvalOverlappingClusters_", combination, ".json");
    fileID = fopen(filepath,'w');
    fprintf(fileID, jsonOut);
    fclose(fileID);
end

tmpPcs2t = pcs2t;
end

function figKeyPress(src, event)
disp(event.Key)
global diff1;
global diff2;
if(strcmp(event.Key,"leftarrow"))
    global pcst;
    global pcs2t;
    global pcClustLoc1;
    global pcClustLoc2;
    global pcInd;
    global pcInd2;

    pcInd = pcInd-1;
    diff1 = diff1-1;
    if pcInd<1
        pcInd=1;
        diff1=diff1+1;
    end

    UpdateClusterErrors(pcst(pcInd),pcs2t(pcInd2));

    hold off;
    pcshow(pcst(pcInd).Location, 'y', 'MarkerSize', 6);
    hold on;
    pcshow(pcs2t(pcInd2).Location, 'c', 'MarkerSize', 6);
    pcshow(pcClustLoc1, 'r', 'MarkerSize', 12);
    pcshow(pcClustLoc2, 'g', 'MarkerSize', 12);
    view(0,-89.9);
    camorbit(-89,0,'camera');
    rotate3d off;
end
if(strcmp(event.Key,"uparrow"))
    global pcst;
    global pcs2t;
    global pcClustLoc1;
    global pcClustLoc2;
    global pcInd;
    global pcInd2;

    diff2=diff2+1;
    pcInd2 = pcInd2+1;
    if pcInd2>length(pcs2t)
        pcInd2=length(pcs2t);
        diff2=diff2-1;
    end

    UpdateClusterErrors(pcst(pcInd),pcs2t(pcInd2));

    hold off;
    pcshow(pcst(pcInd).Location, 'y', 'MarkerSize', 6);
    hold on;
    pcshow(pcs2t(pcInd2).Location, 'c','MarkerSize', 6);
    pcshow(pcClustLoc1, 'r', 'MarkerSize', 12);
    pcshow(pcClustLoc2, 'g', 'MarkerSize', 12);
    rotate3d off;
    view(0,-89.9);
    camorbit(-89,0,'camera');
    rotate3d off;
end
if(strcmp(event.Key,"downarrow"))
    global pcst;
    global pcs2t;
    global pcClustLoc1;
    global pcClustLoc2;
    global pcInd;
    global pcInd2;

    diff2=diff2-1;
    pcInd2 = pcInd2-1;
    if pcInd2<1
        pcInd2=1;
        diff2=diff2+1;
    end

    UpdateClusterErrors(pcst(pcInd),pcs2t(pcInd2));

    hold off;
    pcshow(pcst(pcInd).Location, 'y', 'MarkerSize', 6);
    hold on;
    pcshow(pcs2t(pcInd2).Location, 'c','MarkerSize', 6);
    pcshow(pcClustLoc1, 'r', 'MarkerSize', 12);
    pcshow(pcClustLoc2, 'g', 'MarkerSize', 12);
    rotate3d off;
    view(0,-89.9);
    camorbit(-89,0,'camera');
    rotate3d off;
end
if(strcmp(event.Key,"rightarrow"))
    global pcst;
    global pcs2t;
    global pcClustLoc1;
    global pcClustLoc2;
    global pcInd;
    global pcInd2;

    diff1=diff1+1;
    pcInd = pcInd+1;
    if pcInd>length(pcst)
        pcInd=length(pcst);
        diff1=diff1-1;
    end

    UpdateClusterErrors(pcst(pcInd),pcs2t(pcInd2));

    hold off;
    pcshow(pcst(pcInd).Location, 'y', 'MarkerSize', 6);
    hold on;
    pcshow(pcs2t(pcInd2).Location, 'c','MarkerSize', 6);
    pcshow(pcClustLoc1, 'r', 'MarkerSize', 12);
    pcshow(pcClustLoc2, 'g', 'MarkerSize', 12);
    rotate3d off;
    view(0,-89.9);
    camorbit(-89,0,'camera');
    rotate3d off;
end
title(append("P1: ", string(diff1), " P2: ", string(diff2)))
end

function UpdateClusterErrors(pc1,pc2)

global distRejection;
[labels1, labels2, ~, ~, maxLabelInd1, maxLabelInd2] = DoCluster(pc1,pc2);
[pcClust1, pcClust2] = GetClusters(pc1,pc2, labels1, maxLabelInd1, labels2, maxLabelInd2);
[rmseVal2, stdVal2, ~] = EvalOverlappingRegionFunc(pcClust1, pcClust2, distRejection, false);
disp(append("rmse [mm]: ", string(rmseVal2*1000), " std[mm]: ", string(stdVal2*1000)))

end

function [pcClust1, pcClust2, pcClustLoc1, pcClustLoc2] = GetClusters(pc1,pc2, labels1, maxLabelInd1, labels2, maxLabelInd2)

pcClustLoc1 = pc1.Location(labels1==maxLabelInd1, :);
pcClustLoc2 = pc2.Location(labels2==maxLabelInd2, :);
% check how many points are closer to each other than threshold
pcClust1 = pointCloud(pcClustLoc1);
pcClust2 = pointCloud(pcClustLoc2);

end

function [labels1, labels2, maxClusterNum1, maxClusterNum2, maxLabelInd1, maxLabelInd2] = DoCluster(pc1,pc2)

[labels1, numClusters1] = pcsegdist(pc1, 0.01); %1cm
[labels2, numClusters2] = pcsegdist(pc2, 0.01); %1cm
% keep the largest clusters only
maxLabelInd1=1;
maxClusterNum1 = sum(labels1==1);
for mli=2:numClusters1
    if sum(labels1==mli)>maxClusterNum1
        maxLabelInd1 = mli;
        maxClusterNum1 = sum(labels1==mli);
    end
end
maxLabelInd2=1;
maxClusterNum2 = sum(labels2==1);
for mli=2:numClusters2
    if sum(labels2==mli)>maxClusterNum2
        maxLabelInd2 = mli;
        maxClusterNum2 = sum(labels2==mli);
    end
end

end

function [rmseVal2, stdVal2, tForm] = EvalOverlappingRegionFunc(pcClust1, pcClust2, distRejection, doUseIcp)

% do evaluation for overlapping region
xLimits = [max([pcClust1.XLimits(1), pcClust2.XLimits(1)]),min([pcClust1.XLimits(2), pcClust2.XLimits(2)])];
yLimits = [max([pcClust1.YLimits(1), pcClust2.YLimits(1)]),min([pcClust1.YLimits(2), pcClust2.YLimits(2)])];
zLimits = [max([pcClust1.ZLimits(1), pcClust2.ZLimits(1)]),min([pcClust1.ZLimits(2), pcClust2.ZLimits(2)])];
% cut point clouds and do pcrmse again
pc1Loc = pcClust1.Location;
inds1 = pc1Loc(:,1)>xLimits(2) | pc1Loc(:,1)<xLimits(1) | pc1Loc(:,2)>yLimits(2) | pc1Loc(:,2)<yLimits(1)...
     | pc1Loc(:,3)>zLimits(2) | pc1Loc(:,3)<zLimits(1);
pc1Loc(inds1, :) = [];
pc2Loc = pcClust2.Location;
inds2 = pc2Loc(:,1)>xLimits(2) | pc2Loc(:,1)<xLimits(1) | pc2Loc(:,2)>yLimits(2) | pc2Loc(:,2)<yLimits(1)...
     | pc2Loc(:,3)>zLimits(2) | pc2Loc(:,3)<zLimits(1);
pc2Loc(inds2, :) = [];

[rmseVal2, stdVal2, tForm] = pcrmse3(pointCloud(pc1Loc), pointCloud(pc2Loc), distRejection, doUseIcp); % 5mm NN tolerance for X-Y dist
            
end