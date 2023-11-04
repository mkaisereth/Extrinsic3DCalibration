% Base class for all evaluation runs for e.g. checkerboard, cylinder,
% sphere, 3D static back, ...
classdef RunObject
    
    %static methods
    methods(Static)
       
       function Init()
        % Initalize the Run script
        global printLevel;
        % clear and close all (expect the jsonFilePath)
          close all;
          clear variables -except jsonFilePath showPcRoi printLevel;
       end
       
       
       function jsonFilePath = CheckJsonFilePath(fallbackPath)
           % Check Json File path, use default if none provided
           % fallbackPath - path to be used if jsonFilePath not provided
           % returns the path to the json file
           
          if(~exist('jsonFilePath', 'var') || isempty(jsonFilePath)) %#ok<NODEF> - use default if it does not exist
              jsonFilePath = fallbackPath;
          end
       end
       
       function metaData = ReadJson(filePath)
           % Read the json file given by filePath
           % filePath - path to the json file
           % returns the metaData read from the json file
           
           % read the file
           disp(filePath);
           addpath('..\Utils');
           if(~exist('ReadJson'))
               msgbox('You need to clone utils: https://gitlab.ti.bfh.ch/4d-spine/utils on relative path ..\Utils');
           end
           jsonContent = ReadJson.Read(filePath);
           disp(jsonContent);
           % convert to metadata
           metaData = SingleImgMetadata.FromJson(jsonContent);
           disp(metaData);
       end
       
       function binData = ReadBinaryData(metaData, imgPath)
           % Read the binary data for a depth image (binary Int32 depth
           % data as array (row concat: row1 row2 ...))
           % metaData - the metadata for foreground segmentation
           % imgPath - the path to the binary data depth image
           % returns the binary data as image matrix (height x width)
           
           % read binary file
           binData = ReadBinary.ReadBinaryInt32(imgPath);
           % segment foreground
           binData = SegmentForeground.SegmentMinMax(binData, metaData.SingleDepthImgForegroundMin, metaData.SingleDepthImgForegroundMax);
           % reshape array to image matrix
           binData = reshape(binData,metaData.SingleDepthImgWidth,metaData.SingleDepthImgHeight);
           % transpose to get (height x width)
           binData = binData';
       end
       
       function absBinData = ReadAbsBinaryData(metaData)
           % Read the binary data for a depth image
           % metadata - the metadata with image path
           % returns the binary data as image matrix (height x width)
           
           absBinData = RunObject.ReadBinaryData(metaData, metaData.SingleDepthImgPath);
       end
       
       function [absTargetBinData, relDiffBinData] = ReadAbsTargetBindaryData(metaData, absBinData)
           % Read the binary data for the target depth image
           % metadata - the metadata with target image path
           % absBinData - the binary data for the depth image
           % returns the binary data for the target image matrix and the
           % relative difference image matrix (height x width)
           
           absTargetBinData = RunObject.ReadBinaryData(metaData, metaData.SingleDepthTargetImgPath);
           
           % calculate relative difference
           relDiffBinData = absBinData-absTargetBinData;
       end
       
       function [fig,img] = ShowFigure(binData, imgTitle, metaData)
           % Show a figure for the binary data (depth image)
           % binData - the binary data to show
           % returns the figure and the image of the binary data
           global printLevel;
           if printLevel>0
               fig = figure;
               img = imagesc(binData);
               set(img,'AlphaData',~isnan(binData));
               xlabel('pixel');
               ylabel('pixel');
               title(imgTitle);
               h = colorbar;
               set(get(h,'label'),'string',sprintf('Depth [%s]', metaData.SingleDepthImgUnit));
           else
               fig = [];
               img = [];
           end
       end
       
       function [fgMean, fgStd] = CalcFgMeanAndStd(binData, title, unit)
           fgData = binData(~isnan(binData));
           fgMean = mean2(fgData);
           fgStd = std(fgData);
           fprintf('%s:\n', title);
           fprintf('- Mean: %.4f %s\n', fgMean, unit);
           fprintf('- Std: %.4f %s\n', fgStd, unit);
       end
       
       function [absFgMean, absFgStd] = CalcAbsFgMeanAndStd(absBinData, metaData)
           [absFgMean, absFgStd] = RunObject.CalcFgMeanAndStd(absBinData, 'Absolute depth image', metaData.SingleDepthImgUnit);
       end
       
       function [absTargetFgMean, absTargetFgStd] = CalcAbsTargetFgMeanAndStd(absTargetBinData, metaData)
           [absTargetFgMean, absTargetFgStd] = RunObject.CalcFgMeanAndStd(absTargetBinData, 'Absolute target depth image', metaData.SingleDepthImgUnit);
       end
       
       function [relFgMean, relFgStd] = CalcRelFgMeanAndStd(relDiffBinData, metaData)
           [relFgMean, relFgStd] = RunObject.CalcFgMeanAndStd(relDiffBinData, 'Relative depth image (actual - target)', metaData.SingleDepthImgUnit);
       end
       
       function CheckCreateOutputFolder()
           if ~exist('Output', 'dir')
               mkdir('Output')
           end
       end
       
       function ts = GetTimestamp()
           % Get the current timestamp
           % returns timestamp as string in yyyymmddHHMMSSfff
           
           ts = datestr(now, 'yyyymmddHHMMSSfff'); 
       end
       
       function siOutput = PrepareSiOutput(metaData, msData, mFactor)
           % Prepare the single image output (evaluation results)
           % metaData - metaData from jason file
           % msData - the results
           % returns the single image output
           
           % save all required results
           siOutput = SingleImgOutput;
           siOutput.Mode = metaData.Mode;
           siOutput.Version = metaData.Version;
           siOutput.Description = metaData.Description;
           siOutput.Combination = metaData.Combination;
           if (isfield(msData, 'absFgMean') && isfield(msData, 'absFgStd'))
               siOutput.AbsoluteMean = msData.absFgMean;
               siOutput.AbsoluteStd = msData.absFgStd;
           end
           if (isfield(msData, 'absTargetFgMean') && isfield(msData, 'absTargetFgStd'))
               siOutput.AbsoluteTargetMean = msData.absTargetFgMean;
               siOutput.AbsoluteTargetStd = msData.absTargetFgStd;
           end
           if (isfield(msData, 'relFgMean') && isfield(msData, 'relFgStd'))
               siOutput.RelativeMean = msData.relFgMean;
               siOutput.RelativeStd = msData.relFgStd;
           end
           siOutput.Unit = metaData.SingleDepthImgUnit; 
           % Difference of best fit for depth image is optional
           if(isfield(msData, 'bestFitDiffFgMean') && isfield(msData, 'bestFitDiffFgStd'))
               siOutput.BestFitDiffMean = msData.bestFitDiffFgMean;
               siOutput.BestFitDiffStd = msData.bestFitDiffFgStd;
           end
           % Difference of target object for depth image is optional
           if(isfield(msData, 'targetDiffFgMean') && isfield(msData, 'targetDiffFgStd'))
               siOutput.TargetDiffMean = msData.targetDiffFgMean;
               siOutput.TargetDiffStd = msData.targetDiffFgStd;
           end
           % Mean error of best fit for point cloud is optional
           if(isfield(msData, 'bestFitPCDiffRMSE'))
               siOutput.BestFitPCDiffRMSE = msData.bestFitPCDiffRMSE/mFactor;
           end
           % Mean error of registration for point cloud is optional
           if(isfield(msData, 'registrationPCDiffRMSE'))
               siOutput.RegistrationPCDiffRMSE = msData.registrationPCDiffRMSE/mFactor;
           end
           % region evaluation for point density (pts/mm2), mean difference
           % to target point cloud (mm) and std (mm)
           if(isfield(msData, 'ptDensityMM'))
               siOutput.PtDensityMM = msData.ptDensityMM;
           end
           % Mean difference (region) after registration in mm
           if(isfield(msData, 'zMeanMM'))
               siOutput.ZMeanMM = msData.zMeanMM;
           end
           % Std (region) after registration in mm
           if(isfield(msData, 'zStdMM'))
               siOutput.ZStdMM = msData.zStdMM;
           end
           % whether transform from metadata has been used
           if(isfield(msData, 'usedMetaDataTransform'))
               siOutput.usedMetaDataTransform = msData.usedMetaDataTransform;
           end
           
       end
       
       function SaveSiOutput(siOutput, basePath, ts)
           jsonOut = jsonencode(siOutput);
           filepath = basePath+"_"+ts+".json";
           fileID = fopen(filepath,'w');
           fprintf(fileID, jsonOut);
           fclose(fileID);
       end
       
       function fact = GetUnitFactorToMeter(metaData)
           if(metaData.SingleDepthImgUnit == "um")
               fact = 0.000001;
           elseif(metaData.SingleDepthImgUnit == "mm")
               fact = 0.001;
           elseif(metaData.SingleDepthImgUnit == "cm")
               fact = 0.01;
           else
               fact = 1;
           end   
       end
       
       function [f3, bestFitDiffData, bestFitDiffImg, targetDiffData] = BestFitObject(absBinData, polyStr, metaData, imgTitle)
          % best fitting polynomial plane (defined by polyStr)
          % absBinData - the binary depth data (height x width)
          % polyStr - the polynomial degree, e.g. poly11 for a plane
          % returns the figure, the best fit difference binary data and the
          % best fit difference image
          global printLevel;
          if printLevel>0
              f3 = figure;
              % plot the surface of the binary data
              surf(absBinData);
          end
          
          
          % best fit polynomial plane (regression)
          [m,n] = size(absBinData); 
          [X,Y] = meshgrid(1:n,1:m); 
          Xi = X(:);
          Yi = Y(:);
          Zi = absBinData(:);
          mask = ~isnan(Zi);
          Xim = Xi(mask);
          Yim = Yi(mask);
          Zim = Zi(mask);
          
          %sanity check
          if length(Zim) < 3
              bestFitDiffData = [];
              bestFitDiffImg = imagesc();
              targetDiffData = [];
              return;
          end
          
          P = fit([Xim, Yim], Zim,polyStr);
          % plot the best fit polynomial plane
          if printLevel>0
              bestFitImg = plot(P,[Xim,Yim], Zim);
              title(imgTitle);
              xlabel('x - pixel');
              ylabel('y - pixel');
              zlabel(sprintf('Depth [%s]', metaData.SingleDepthImgUnit));
          end
          % get best fit matrix
          fit2 = P(Xim,Yim);
          absBinData2 = absBinData;
          absBinData2(~isnan(absBinData)) = fit2;
          % get mean and std
          RunObject.CalcFgMeanAndStd(absBinData2, 'Best fit depth image', metaData.SingleDepthImgUnit);

          % additional target object
          if(~isempty(metaData.SingleTargetObject))
              if(metaData.SingleTargetObject.type == "plane")
                  point = metaData.SingleTargetObject.point;
                  normal = metaData.SingleTargetObject.normal;
                  %# a plane is a*x+b*y+c*z+d=0
                  %# [a,b,c] is the normal. Thus, we have to calculate
                  %# d and we're set
                  d = -point'*normal; %'# dot product for less typing

                  [m,n] = size(absBinData);
                  [X,Y] = meshgrid(1:n,1:m);

                  %# calculate corresponding z
                  Z = (-normal(1)*X - normal(2)*Y - d)/normal(3);

                  %# plot the surface
                  if printLevel>0
                      hold on;
                      surf(X,Y,Z)
                      hold off;
                  end

                  % difference between original data and target object
                  targetDiffData = absBinData-Z;
              end
              if(metaData.SingleTargetObject.type == "cylinder")
                  point = metaData.SingleTargetObject.point;
                  normal = metaData.SingleTargetObject.normal;
                  radius = metaData.SingleTargetObject.radius;
                  
                  % get height
                  h = norm(normal);
                  % get default cylinder (in z direction)
                  [X,Y,Z] = cylinder(radius);
                  % scale the height
                  Z = Z*h;
                  % calc rotation
                  aRot = vrrotvec([0,0,1],normal);
                  mRot = axang2rotm(aRot);
                  
                  % rotate
                  Xi = X(:);
                  Yi = Y(:);
                  Zi = Z(:);
                  cyDataVect = [Xi Yi Zi]';
                  rotcyDataVect = mRot*cyDataVect;
                  Xi = rotcyDataVect(1,:);
                  Yi = rotcyDataVect(2,:);
                  Zi = rotcyDataVect(3,:);

                  % translate                  
                  Xi= Xi+point(1);
                  Yi= Yi+point(2);
                  Zi= Zi+point(3);
                  
                  X = reshape(Xi,2,21);
                  Y = reshape(Yi,2,21);
                  Z = reshape(Zi,2,21);
                  
                  %# plot the surface
                  if printLevel>0
                      hold on;
                      surf(X,Y,Z)
                      hold off;
                  end

                  % difference between original data and target object
                  targetDiffData = [];% TODO absBinData-Z;
              end
              if(metaData.SingleTargetObject.type == "sphere")
                  point = metaData.SingleTargetObject.point;
                  radius = metaData.SingleTargetObject.radius;
                  
                  % get default cylinder (in z direction)
                  [X,Y,Z] = sphere;
                  % scale
                  X = X*radius;
                  Y = Y*radius;
                  Z = Z*radius;

                  % translate                  
                  X= X+point(1);
                  Y= Y+point(2);
                  Z= Z+point(3);
                  
                  %# plot the surface
                  if printLevel>0
                      hold on;
                      surf(X,Y,Z)
                      hold off;
                  end

                  % difference between original data and target object
                  targetDiffData = [];% TODO absBinData-Z;
              end
          else
              targetDiffData = [];
          end
          
          % difference between original data and best fit poly plane
          bestFitDiffData = absBinData-absBinData2;
          
          [~,bestFitDiffImg] = RunObject.ShowFigure(bestFitDiffData, sprintf('%s (difference)', imgTitle), metaData);
       end
       
       function [bestFitDiffFgMean, bestFitDiffFgStd] = CalcBestFitDiffFgMeanAndStd(bestFitDiffBinData, metaData)
           [bestFitDiffFgMean, bestFitDiffFgStd] = RunObject.CalcFgMeanAndStd(bestFitDiffBinData, 'Best fit difference depth image (actual - best fit)', metaData.SingleDepthImgUnit);
       end
       
       function distsMean = GetAverageDistanceToKNearestNeighbors(ptCloud, k)
           % Get the average of the distance to the k nearest neighbors for
           % a point cloud
           % ptCloud - the point cloud
           % k - the number o nearest neighbors, e.g. 3
           % returns the average distance to the k nearest neighbors
           
           % get the number of points in the point cloud
           numOfPts = size(ptCloud.Location, 1);
           % array for storing distances
           dists = NaN(1,numOfPts);
           for i = 1:numOfPts
               % find the k nearest neighbors of this point
               [~,dists2] = findNearestNeighbors(ptCloud,ptCloud.Location(i,1:3),k);
               % get the mean and add to array
               distMean = mean(dists2);
               dists(1,i) = distMean;
           end
           % get the mean over all distances
           distsMean = mean(dists);
       end
       
       function [ptCloud, mFactor] = ReadFgPointCloud(fileName, metaData, title, showPcRoi)
           % Read point cloud from file and segment foreground
           % fileName - path to pointcloud file (pcd or ply)
           % metaData - metaData for evaluation
           % returns the point cloud
           
           % read point cloud and get points
           ptCloud = pcread(fileName);
           xyzPts = ptCloud.Location;
           % convert ordered pointcloud to random
           if length(size(xyzPts)) == 3
               xyzPts = reshape(xyzPts, [], 3);
           end
           % check whether in m or mm
           if(max(xyzPts(:,3))>10) % more than 10 must be mm
               xyzPts = xyzPts/1000;
           end
           % get the factor to convert depth image unit to meter (pcd is in
           % meter)
           mFactor = RunObject.GetUnitFactorToMeter(metaData);
           % segment foreground
           if size(metaData.SingleDepthImgForegroundMin, 1) == 1 && size(metaData.SingleDepthImgForegroundMin, 2) == 1
               xyzPts(xyzPts(:,3)>mFactor*metaData.SingleDepthImgForegroundMax, :) = [];
               xyzPts(xyzPts(:,3)<mFactor*metaData.SingleDepthImgForegroundMin, :) = [];
           elseif size(metaData.SingleDepthImgForegroundMin, 1) > 2 || size(metaData.SingleDepthImgForegroundMin, 2) > 2
               xyzPts(xyzPts(:,1)>mFactor*metaData.SingleDepthImgForegroundMax(1), :) = [];
               xyzPts(xyzPts(:,1)<mFactor*metaData.SingleDepthImgForegroundMin(1), :) = [];
               xyzPts(xyzPts(:,2)>mFactor*metaData.SingleDepthImgForegroundMax(2), :) = [];
               xyzPts(xyzPts(:,2)<mFactor*metaData.SingleDepthImgForegroundMin(2), :) = [];
               xyzPts(xyzPts(:,3)>mFactor*metaData.SingleDepthImgForegroundMax(3), :) = [];
               xyzPts(xyzPts(:,3)<mFactor*metaData.SingleDepthImgForegroundMin(3), :) = [];
           end
           % create segmented pointcloud
           ptCloud = pointCloud(xyzPts);
           
           % if pointcloud is huge (e.g. AS), downsample
           if ptCloud.Count>250000
               ptCloud = RunObject.DownSamplePC(ptCloud,250000);
           end
           
           %TODO
           addpath('..\Utils');
           if(~exist('pcRoi'))
               msgbox('You need to clone utils: https://gitlab.ti.bfh.ch/4d-spine/utils on relative path ..\Utils');
           end
           if(exist('showPcRoi') && showPcRoi == 1)
               ptCloud = pcRoi(fileName, ptCloud);
           end
           xyzPts = ptCloud.Location;
           zPts = xyzPts(:,3);
           % get mean and std
           RunObject.CalcFgMeanAndStd(zPts(~isnan(zPts)), title, 'm');
       end
       
       function [ptCloud, mFactor] = ReadAbsoluteFgPointCloud(fileName, metaData, showPcRoi)
           % Read point cloud from file and segment foreground
           % fileName - path to pointcloud file (pcd or ply)
           % metaData - metaData for evaluation
           % returns the point cloud
           
           [ptCloud, mFactor] = RunObject.ReadFgPointCloud(fileName, metaData, 'Absolute point cloud (depth)', showPcRoi);
       end
       
       function [ptCloud, mFactor, gPtCloud] = ReadTargetFgPointCloud(fileName, metaData, showPcRoi)
           % Read point cloud from file and segment foreground
           % fileName - path to pointcloud file (pcd or ply)
           % metaData - metaData for evaluation
           % returns the point cloud
           
           [ptCloud, mFactor] = RunObject.ReadFgPointCloud(fileName, metaData, 'Target point cloud (depth)', showPcRoi);
           
           % read geometric target point cloud (if any)
           if(~isempty(metaData.SinglePointCloudTargetObject))
               if(metaData.SinglePointCloudTargetObject.type == "plane")
                  point = metaData.SinglePointCloudTargetObject.point;
                  normal = metaData.SinglePointCloudTargetObject.normal;
                  sizexy = metaData.SinglePointCloudTargetObject.size;
                  deltaxy = metaData.SinglePointCloudTargetObject.delta;
                  %# a plane is a*x+b*y+c*z+d=0
                  %# [a,b,c] is the normal. Thus, we have to calculate
                  %# d and we're set
                  d = -point'*normal; %'# dot product for less typing

                  n = sizexy(1);
                  m = sizexy(2);
                  [X,Y] = meshgrid(-n/2:deltaxy(1):n/2+point(1),-m/2:deltaxy(2):m/2+point(2));

                  %# calculate corresponding z
                  Z = (-normal(1)*X - normal(2)*Y - d)/normal(3);
                  % add some noise TODO
                  Z = Z + 0.0001*(rand(size(Z,1),size(Z,2))-0.5);
                  
                  gXyzPts(:,:,1) = X;
                  gXyzPts(:,:,2) = Y;
                  gXyzPts(:,:,3) = Z;
                  gPtCloud = pointCloud(gXyzPts);
               end
               if(metaData.SinglePointCloudTargetObject.type == "cylinder")
                  point = metaData.SinglePointCloudTargetObject.point;
                  normal = metaData.SinglePointCloudTargetObject.normal;
                  radius = metaData.SinglePointCloudTargetObject.radius;
                  degree = metaData.SinglePointCloudTargetObject.degree;
                  
                  % get height
                  h = norm(normal);
                  numOfFullPts = 360;
                  
                  numOfMinPts = floor((180-degree)/2);
                  numOfMaxPts = numOfMinPts+degree;
                  numOfPts = degree+1;

                  numOfLinePts = 401;
                  delta = 0.0025;
                  
                  % get default cylinder (in z direction)
                  [X,Y,Z] = cylinder(radius, numOfFullPts);
                  
                  % take only half
                  X2(1,:) = X(1,numOfMinPts:numOfMaxPts);
                  Y2(1,:) = Y(1,numOfMinPts:numOfMaxPts);
                  X = X2;
                  Y = Y2;
                  
                  X = repmat(X(1,:),numOfLinePts,1);
                  Y = repmat(Y(1,:),numOfLinePts,1);
                  Z2(1:numOfLinePts,:) = repmat([0:delta:1]', 1, numOfPts);
                  Z = Z2;
                  
                  % scale the height
                  Z = Z*h;
                  % calc rotation
                  aRot = vrrotvec([0,0,1],normal);
                  mRot = axang2rotm(aRot);
                  
                  % rotate
                  Xi = X(:);
                  Yi = Y(:);
                  Zi = Z(:);
                  cyDataVect = [Xi Yi Zi]';
                  rotcyDataVect = mRot*cyDataVect;
                  Xi = rotcyDataVect(1,:);
                  Yi = rotcyDataVect(2,:);
                  Zi = rotcyDataVect(3,:);

                  % translate                  
                  Xi= Xi+point(1);
                  Yi= Yi+point(2);
                  Zi= Zi+point(3);
                  
                  X = reshape(Xi,numOfLinePts,numOfPts);
                  Y = reshape(Yi,numOfLinePts,numOfPts);
                  Z = reshape(Zi,numOfLinePts,numOfPts);
                  
                  gXyzPts(:,:,1) = X;
                  gXyzPts(:,:,2) = Y;
                  gXyzPts(:,:,3) = Z;
                  gPtCloud = pointCloud(gXyzPts);
               end
               if(metaData.SinglePointCloudTargetObject.type == "sphere")
                  point = metaData.SinglePointCloudTargetObject.point;
                  radius = metaData.SinglePointCloudTargetObject.radius;
                  degree = metaData.SinglePointCloudTargetObject.degree;

                  % get default cylinder (in z direction)
                  [X,Y,Z] = sphere(360);
                  
                  % keep only partial sphere
                  numOfPts = degree+1;
                  X2 = X(1:numOfPts,:);
                  Y2 = Y(1:numOfPts,:);
                  Z2 = Z(1:numOfPts,:);
                  X = X2;
                  Y = Y2;
                  Z = Z2;
                  
                  % scale
                  X = X*radius;
                  Y = Y*radius;
                  Z = Z*radius;

                  % translate                  
                  X= X+point(1);
                  Y= Y+point(2);
                  Z= Z+point(3);
                   
                  gXyzPts(:,:,1) = X;
                  gXyzPts(:,:,2) = Y;
                  gXyzPts(:,:,3) = Z;
                  gPtCloud = pointCloud(gXyzPts);
              end
           else
               gPtCloud = [];
           end
       end
       
       function ptCloud = DownSamplePC(ptCloud, numOfPoints)
           % ptCloud = pcdownsample(ptCloud, 'gridAverage', 0.01);
           
           if ~exist('numOfPoints')
               numOfPoints = 10000;
           end
           
           % adaptive downsampling
           ds = 0.1; % sample down to 10%
           if ptCloud.Count > 0
               ds = 1/(ptCloud.Count/numOfPoints); % sample down to 10'000 pts
               if ds > 1
                   ds = 1;
               end
           end
           ptCloud = pcdownsample(ptCloud, 'random', ds);
       end
       
       function [ptCloudTrans, tForm, rmse, f6, f7, usedMetaDataTransform] = RegisterPointCloud(ptCloud, tPtCloud, gTPtCloud, metaData)
           % Register two point clouds (ptCloud to tPtCloud)
           
           % downsample first
           ptCloudDs = RunObject.DownSamplePC(ptCloud);
           tPtCloudDs = RunObject.DownSamplePC(tPtCloud);
           
           % does not deliver satisfying results
           %tForm = pcregistericp(ptCloud, tPtCloud2);
           
           % did not work at all
           %tForm = pcregistercorr(ptCloud, tPtCloud2, 0.3, 0.01);
           
           % coherent point drift (CPD)
           % pcregistercpd can match the rotation and transformation perfectly ('Rigid')
           %tForm = pcregistercpd(ptCloud, tPtCloud, 'Transform', 'Rigid');
           disp('Going to register point clouds ...');
           
           
           if ~isempty(metaData.Calib3dPath)
               % use transform from file
               
               load(metaData.Calib3dPath);
               if exist("tFormMean","var")
                   tForm = tFormMean.Transform;
               elseif exist("tForm2s", "var")
                   tForm = tForm2s;
               else
                   tForm = [];
               end
               rmse = [];
               
               usedMetaDataTransform = 1;
               fprintf('- Root mean square error (RMSE) of register pcd (3dcb): %.4f m\n', rmse);
           else
               % do it 10 times and pick the best
               tForms = rigid3d.empty(1,0);
               rmses = double.empty(1,0);
               for l=1:1
                   [tForm,~,rmse] = pcregistercpd(ptCloudDs, tPtCloudDs, 'Transform', 'Rigid', 'MaxIterations', 40, 'Tolerance', 1e-10);
                   tForms(l) = tForm;
                   rmses(l) = rmse;
               end
              for l=1:1
                  if(rmse>rmses(l))
                      rmse = rmses(l);
                      tForm = tForms(l);
                  end
              end
               usedMetaDataTransform = 0;
               fprintf('- Root mean square error (RMSE) of register pcd: %.4f m\n', rmse);
           end
           % does not deliver satisfying results
           %tForm = pcregisterndt(ptCloud, tPtCloud2, 0.1);
           
           ptCloudTrans = pctransform(ptCloud, tForm);
           % do again
       if usedMetaDataTransform == 0
           ptCloudDs2 = RunObject.DownSamplePC(ptCloudTrans, 50000);
           tPtCloudDs2 = RunObject.DownSamplePC(tPtCloud, 50000);

           tForm2 = pcregistericp(ptCloudDs2, tPtCloudDs2);
           ptCloudTrans = pctransform(ptCloudTrans, tForm2);
       end
       global printLevel;
       f6 = [];
       f7 = [];
       if printLevel>0
           f6 = figure;
           pcshow(tPtCloud.Location, 'red');
           hold on;
           pcshow(ptCloud.Location, 'green');
           title('point cloud (green) with target point cloud (red)');
           xlabel('x - position [m]');
           ylabel('y - position [m]');
           zlabel('z - depth [m]');
           
           f7 = figure;
           pcshow(tPtCloud.Location, 'red');
           hold on;
           pcshow(ptCloudTrans.Location, 'green');
           title('registered point cloud (green) with target point cloud (red)');
           xlabel('x - position [m]');
           ylabel('y - position [m]');
           zlabel('z - depth [m]');
       end

           % geometric target point cloud (if any)
           if(~isempty(gTPtCloud))
               gTPtCloudDs = RunObject.DownSamplePC(gTPtCloud);

               % coherent point drift (CPD)
               disp('Going to register geometric point clouds ...');
               [gTForm,~,gRmse] = pcregistercpd(ptCloudDs, gTPtCloudDs, 'Transform', 'Rigid');
               fprintf('- Root mean square error (RMSE) of register geometric pc: %.4f m\n', gRmse);

               gPtCloudTrans = pctransform(ptCloud, gTForm);

               if printLevel>0
                   f6_1 = figure;
                   pcshow(gTPtCloud.Location, 'blue');
                   hold on;
                   pcshow(ptCloud.Location, 'green');
                   title('point cloud (green) with geometric target point cloud (blue)');
                   xlabel('x - position [m]');
                   ylabel('y - position [m]');
                   zlabel('z - depth [m]');
    
                   f7_1 = figure;
                   pcshow(gTPtCloud.Location, 'blue');
                   hold on;
                   pcshow(gPtCloudTrans.Location, 'green');
                   title('registered point cloud (green) with geometric target point cloud (blue)');
                   xlabel('x - position [m]');
                   ylabel('y - position [m]');
                   zlabel('z - depth [m]');
               end
           end
       end
       
       function [f6, f7, f8] = RegisterPointCloud3(ptCloud1, ptCloud2, tPtCloud)
           % Register two point clouds (ptCloud to tPtCloud)
           
           % downsample first
           ptCloud1Ds = RunObject.DownSamplePC(ptCloud1);
           ptCloud2Ds = RunObject.DownSamplePC(ptCloud2);
           tPtCloudDs = RunObject.DownSamplePC(tPtCloud);
           
           % coherent point drift (CPD)
           disp('Going to register point clouds ...');
           [tForm1,~,rmse1] = pcregistercpd(ptCloud1Ds, tPtCloudDs, 'Transform', 'Rigid');
           [tForm2,~,rmse2] = pcregistercpd(ptCloud2Ds, tPtCloudDs, 'Transform', 'Rigid');
           fprintf('- Root mean square error (RMSE) of register pcd: %.4f m\n', rmse1);
           fprintf('- Root mean square error (RMSE) of register pcd: %.4f m\n', rmse2);

           ptCloud1Trans = pctransform(ptCloud1, tForm1);
           ptCloud2Trans = pctransform(ptCloud2, tForm2);
           
           f6 = figure;
           pcshow(tPtCloud.Location, 'red');
           hold on;
           pcshow(ptCloud1.Location, 'green');
           pcshow(ptCloud2.Location, 'blue');
           title('point clouds (green, blue) with target point cloud (red)');
           xlabel('x - position [m]');
           ylabel('y - position [m]');
           zlabel('z - depth [m]');
           
           f7 = figure;
           pcshow(tPtCloud.Location, 'red');
           hold on;
           pcshow(ptCloud1Trans.Location, 'green');
           pcshow(ptCloud2Trans.Location, 'blue');
           title('registered point cloud (green, blue) with target point cloud (red)');
           xlabel('x - position [m]');
           ylabel('y - position [m]');
           zlabel('z - depth [m]');
           
           xyzPts = [ptCloud1Trans.Location', ptCloud2Trans.Location', tPtCloud.Location']';
           totalPtCloud = pointCloud(xyzPts);
           totalPtCloud = pcmedian(totalPtCloud, 'Radius', 0.01); % TODO use a value according to average distance to nearest neighbors?
           
           f8 = figure;
           pcshow(totalPtCloud.Location, 'red');
           title('median filtered registered point clouds');
           xlabel('x - position [m]');
           ylabel('y - position [m]');
           zlabel('z - depth [m]');
       end
       
       function [f5, ptCloud2, planeModel, rmse] = PCBestFitObject(ptCloud, pcfitobject, numOfFits)
          % best fitting object (defined by function pcfitobject) for point cloud
          % fileName - Path to point cloud file (pcd or ply)
          % metaData - Metadata from JSON file
          % pcfitobject - a method to fit an object (e.g. pcfitplane, ...)
          % numOfFits - the number of fits (take the best o n tries)
          
          
          % smoothing
%           figure;
%           pcshow(ptCloud3);
%           title('original PC');
          
          % useless (this is only useful if you have noisy outliers)
          %ptCloud2 = pcdenoise(ptCloud3,'NumNeighbors', 4, 'Threshold', 2);
          
          % pcmedian seems useful - median filtering of 3d point cloud data
          % with median over neighbors in Radius
          ptCloud2 = pcmedian(ptCloud, 'Radius', 0.01); % TODO use a value according to average distance to nearest neighbors?
          
          global printLevel;
          f5 = [];
          if printLevel>0
              f5 = figure;
              pcshow(ptCloud2);
              title('median filtered point cloud with best fitting object');
              xlabel('x - position [m]');
              ylabel('y - position [m]');
              zlabel('z - depth [m]');
          end

          %%% wtf? https://ch.mathworks.com/matlabcentral/answers/359659-pcfitcylinder-giving-different-answers
          %rng(0);
          %%%
          % pcfit methods from matlab are quite unstable (e.g.
          % pcfitcylinder), thus use multiple runs and use best one
          planeModel = [];
          meanError = [];
          fprintf('Try to fit best fitting object for %d times:\n', numOfFits);
          for i=1:numOfFits
              try
                  [planeModel2, ~, outInds, meanError2] = pcfitobject(ptCloud2);
                  fprintf('- Mean error: %.4f\n', meanError2);
                  if (isempty(planeModel) || meanError2 < meanError)
                      planeModel = planeModel2;
                      meanError = meanError2;
                  end
                  if(size(outInds,1) ~= 0)
                      disp(' There are some outliers ...')
                  end
              catch e
                  fprintf(' Error: pcfitobject: %s\n',e.message);
              end
          end

          if printLevel>0
          hold on;
              h = plot(planeModel);
              set(h, 'FaceAlpha', 0.2)
              hold off;
          end

          % plot meanError
          rmse = sqrt(meanError);
          fprintf('- Root mean square error (RMSE) of best fit after %d fits: %.4f m\n', numOfFits, rmse);
       end
       
       function [pcWidth, pcHeight, pcDepth] = GetPcSize(ptCloud, mFactor)
           % Get the size of the point cloud
           % ptCloud - point cloud to measure size
           % mFactor - Conversion factor from depthImageUnits to m
           % returns width, height and depth of point cloud
           
           [pcWidth, pcHeight, pcDepth] = RunObject.GetPcSizeK(ptCloud, mFactor, 100);
           
           fprintf('Object size is:\n');
           fprintf('- Width (X): %.2f m\n', pcWidth);
           fprintf('- Height (Y): %.2f m\n', pcHeight);
           fprintf('- Depth (Z): %.2f m\n', pcDepth);
       end
       
       function [pcWidth, pcHeight, pcDepth, pcMiddleX, pcMiddleY, pcMiddleZ] = GetPcSizeK(ptCloud, mFactor, k)
           % Get the size of the point cloud
           % ptCloud - the point cloud to measure the size
           % mFactor - factor for conversion from depthImageUnits to m
           % k - the number of points to get max and min in all dimensions
           % returns the size and middle point of the point cloud in all
           % dimensions
           
           % get max and min in all dimensions
           pcXMax = maxk(ptCloud.Location(:,1), k);
           pcYMax = maxk(ptCloud.Location(:,2), k);
           pcZMax = maxk(ptCloud.Location(:,3), k);
           pcXMin = mink(ptCloud.Location(:,1), k);
           pcYMin = mink(ptCloud.Location(:,2), k);
           pcZMin = mink(ptCloud.Location(:,3), k);
           % difference ist size
           pcWidth = mean2(pcXMax)-mean2(pcXMin);
           pcHeight = mean2(pcYMax)-mean2(pcYMin);
           pcDepth = mean2(pcZMax)-mean2(pcZMin);
           % get the middle
           pcMiddleX = (mean2(pcXMax)+mean2(pcXMin))/2;
           pcMiddleY = (mean2(pcYMax)+mean2(pcYMin))/2;
           pcMiddleZ= (mean2(pcZMax)+mean2(pcZMin))/2;
           
%            fprintf('PcWidth: %.2f mm\n', pcWidth/mFactor);
%            fprintf('PcHeight: %.2f mm\n', pcHeight/mFactor);
%            fprintf('PcDepth: %.2f mm\n', pcDepth/mFactor); 
       end
       
       function [pcXDiff, pcHeight, pcZDiff] = GetPcYOrientation(ptCloud, k)
           [pcXMax, pcXMaxI] = maxk(ptCloud.Location(:,1), k);
           [pcYMax, pcYMaxI] = maxk(ptCloud.Location(:,2), k);
           [pcZMax, pcZMaxI] = maxk(ptCloud.Location(:,3), k);
           [pcXMin, pcXMinI] = mink(ptCloud.Location(:,1), k);
           [pcYMin, pcYMinI] = mink(ptCloud.Location(:,2), k);
           [pcZMin, pcZMinI] = mink(ptCloud.Location(:,3), k);
           
           pcHeight = mean2(pcYMax)-mean2(pcYMin);
           
           % check whether min and max are aligned in different axes
           maxXes = ptCloud.Location(pcYMaxI, 1);
           maxXAvg = mean(maxXes);
           minXes = ptCloud.Location(pcYMinI, 1);
           minXAvg = mean(minXes);
           pcXDiff = maxXAvg-minXAvg;
           
           maxZes = ptCloud.Location(pcYMaxI, 3);
           maxZAvg = mean(maxZes);
           minZes = ptCloud.Location(pcYMinI, 3);
           minZAvg = mean(minZes);
           pcZDiff = maxZAvg-minZAvg;
       end
       
       function ptCloudOut = AlignPCWithYCoord(ptCloud, mFactor)
           
           ptCloudOut = pointCloud(ptCloud.Location);
           % move to origin
           [~, ~, ~, pcMiddleX, pcMiddleY, pcMiddleZ] = RunObject.GetPcSizeK(ptCloudOut, mFactor, 10000);
           ptCloudOut = TransformPC.AffineTransformPC(ptCloudOut, [0 0 0],[-pcMiddleX -pcMiddleY -pcMiddleZ], [1 1 1]);
           
           % Transform x Angle such that Cylinder is aligned with y coordinate (no Rotation)
           [~, pcHeight, pcZDiff] = RunObject.GetPcYOrientation(ptCloudOut, 10000);
           xAngle = atan2(pcZDiff,pcHeight);
           ptCloudOut = TransformPC.AffineTransformPC(ptCloudOut, [xAngle 0 0],[0 0 0], [1 1 1]);
           
           % Transform z Angle such that Cylinder is aligned with y coordinate (no Rotation)
           [pcXDiff, pcHeight, ~] = RunObject.GetPcYOrientation(ptCloudOut, 10000);
           zAngle = atan2(pcXDiff, pcHeight);
           ptCloudOut = TransformPC.AffineTransformPC(ptCloudOut, [0 0 -zAngle],[0 0 0], [1 1 1]); 
       end
       
       function ptCloudOut = AlignPCWithYCoordSym(ptCloud, mFactor)
           
           ptCloudOut = pointCloud(ptCloud.Location);
           numOfPts = floor(0.25 * max(size(ptCloudOut.Location,1),size(ptCloudOut.Location,2)));
           % move to origin
           [~, ~, ~, pcMiddleX, pcMiddleY, pcMiddleZ] = RunObject.GetPcSizeK(ptCloudOut, mFactor, numOfPts);
           ptCloudOut = TransformPC.AffineTransformPC(ptCloudOut, [0 0 0],[-pcMiddleX -pcMiddleY -pcMiddleZ], [1 1 1]);
           
           % Transform x Angle such that Cylinder is aligned with y coordinate (no Rotation)
           [~, pcHeight, pcZDiff] = RunObject.GetPcYOrientation(ptCloudOut, numOfPts);
           xAngle = atan2(pcZDiff,pcHeight);
           ptCloudOut = TransformPC.AffineTransformPC(ptCloudOut, [xAngle 0 0],[0 0 0], [1 1 1]);
           
           % Transform z Angle such that Cylinder is aligned with y coordinate (no Rotation)
           
           minDiff = intmax;
           minPC = ptCloudOut;
           minAngle = 0;
           for i=-15:15
               radAngle = i/180*pi;
               ptCloudTemp = TransformPC.AffineTransformPC(ptCloudOut, [0 0 radAngle],[0 0 0], [1 1 1]);
               [pcXDiff, pcHeight, ~] = RunObject.GetPcYOrientation(ptCloudTemp, numOfPts);
               tempDiff = abs(pcXDiff/pcHeight);
               if(tempDiff < minDiff)
                   minDiff = tempDiff;
                   minPC = ptCloudTemp;
                   minAngle = i;
               end
           end
           ptCloudOut = minPC;
       end
       
       function [ptDensityMM, zMeanMM, zStdMM, f1038] = EvaluateRegion(ptCloud, tForm, tPtCloud)
           disp('Going to evaluate region in point cloud ...');
           regionEvalX = 0.6;
           regionEvalY = 0.8;

           ptCloudTrans = pctransform(pointCloud(ptCloud.Location), tForm);

            % get middle point
%             pcMin1 = min(ptCloudTrans.Location,[],1);
%             pcMax1 = max(ptCloudTrans.Location,[],1);
%             pcMin2 = min(tPtCloud.Location,[],1);
%             pcMax2 = max(tPtCloud.Location,[],1);
%             pcMin = max(pcMin1, pcMin2);
%             pcMax = min(pcMax1, pcMax2);
            xyzrPoints1 = ptCloudTrans.Location;
            xyztPoints = tPtCloud.Location;
%             xyzrPoints1(xyzrPoints1(:,1)>pcMax(1), :) = [];
%             xyzrPoints1(xyzrPoints1(:,1)<pcMin(1), :) = [];
%             xyzrPoints1(xyzrPoints1(:,2)>pcMax(2), :) = [];
%             xyzrPoints1(xyzrPoints1(:,2)<pcMin(2), :) = [];
            P = xyztPoints(:,1:2);
            k = boundary(double(P), 1);
            P_boundary = P(k,:);
            in = inpolygon(xyzrPoints1(:,1), xyzrPoints1(:,2), P_boundary(:,1), P_boundary(:,2));
%             figure;
%             plot(xyzrPoints1(in, 1),xyzrPoints1(in, 2),'r+') % points inside
%             hold on;
%             plot(xyzrPoints1(~in, 1),xyzrPoints1(~in, 2),'bo') % points outside
            xyzrPoints1 = xyzrPoints1(in,:);
%             xyzrPoints2 = ptCloudTrans.Location;
%             xyzrPoints2(xyzrPoints2(:,1)>pcMax(1), :) = [];
%             xyzrPoints2(xyzrPoints2(:,1)<pcMin(1), :) = [];
%             xyzrPoints2(xyzrPoints2(:,2)>pcMax(2), :) = [];
%             xyzrPoints2(xyzrPoints2(:,2)<pcMin(2), :) = [];
            
%             qxmin = quantile([xyzrPoints1(:,1);xyzrPoints2(:,1)], 0.1, 1);
%             qxmax = quantile([xyzrPoints1(:,1);xyzrPoints2(:,1)], 0.9, 1);
%             qymin = quantile([xyzrPoints1(:,2);xyzrPoints2(:,2)], 0.1, 1);
%             qymax = quantile([xyzrPoints1(:,2);xyzrPoints2(:,2)], 0.9, 1);
%             pcMedian = median([xyzrPoints1;xyzrPoints2],1);
%             regionUpperBorderX = qxmax;%pcMedian + regionEval*0.5*(pcMax-pcMin);
%             regionLowerBorderX = qxmin;%pcMedian - regionEval*0.5*(pcMax-pcMin);
%             regionUpperBorderY = qymax;%pcMedian + regionEval*0.5*(pcMax-pcMin);
%             regionLowerBorderY = qymin;%pcMedian - regionEval*0.5*(pcMax-pcMin);
%             regionUpperBorder = [regionUpperBorderX, regionUpperBorderY];
%             regionLowerBorder = [regionLowerBorderX, regionLowerBorderY];
            % get the points in that region
%             xyzrPoints = ptCloudTrans.Location;
%             xyzrPoints(xyzrPoints(:,1)>regionUpperBorder(1), :) = [];
%             xyzrPoints(xyzrPoints(:,1)<regionLowerBorder(1), :) = [];
%             xyzrPoints(xyzrPoints(:,2)>regionUpperBorder(2), :) = [];
%             xyzrPoints(xyzrPoints(:,2)<regionLowerBorder(2), :) = [];

            % how about another approach: use 2D silhouette
            P = xyzrPoints1(:,1:2);
            k = boundary(double(P), 1);
            global printLevel;
            if printLevel>0
                figure;
                plot(P(:,1),P(:,2),'.','MarkerSize',10)
                hold on;
                plot(P(k,1),P(k,2))
            end
            % grab all points with enough distance to the boundary
            [~, dists] = knnsearch(P(k,:), P);
            P_smaller = P;
            indsKeep = dists>0.01;
            P_smaller(~indsKeep, :) = [];
            if printLevel>0
                plot(P_smaller(:,1),P_smaller(:,2),'.','MarkerSize',8, 'Color','r')
            end
            [~,Area] = boundary(double(P_smaller), 1);
            xyzrPoints = xyzrPoints1(indsKeep, :);
            % calculate the point density in that region (pts/m^2)
            %area = (regionUpperBorder(1)-regionLowerBorder(1))*(regionUpperBorder(2)-regionLowerBorder(2));
            %fprintf('- Point cloud area: %.2f mm x %.2f mm\n', (regionUpperBorder(1)-regionLowerBorder(1))*1000, (regionUpperBorder(2)-regionLowerBorder(2))*1000);
            ptDensity = length(xyzrPoints)/Area;
            % pts /mm^2
            ptDensityMM = ptDensity/1000000;
            fprintf('- Point cloud density: %.2f pts/mm^2\n', ptDensityMM);
            % get average difference (after registration transformation)
            rPtCloudTrans = pointCloud(xyzrPoints);
            zMean = 0;
            zDiffs = zeros(length(rPtCloudTrans.Location),1);
            maxNearestNeighbor = 0;
            for i=1:length(rPtCloudTrans.Location)
                [indices,dists] = findNearestNeighbors(tPtCloud,rPtCloudTrans.Location(i,:),1);
                % safety check for x and y (I am only interested in z difference)
                if(dists > maxNearestNeighbor)
                    maxNearestNeighbor = dists;
                end
                % get z distance
                zDiff = rPtCloudTrans.Location(i,3)-tPtCloud.Location(indices,3);
                zDiffs(i) = zDiff;
                zMean = zMean+zDiff;
            end
            fprintf('- Nearest neighbor was %.2f mm away\n', maxNearestNeighbor*1000);
            
            if printLevel > 1
                f1038 = figure;
                pcshow(rPtCloudTrans.Location, 'g');
                hold on;
                pcshow(tPtCloud.Location, 'r');
                pause(1)
                view(0,-60)
            else
                f1038 = [];
            end
            
            % calc the mean difference
            zMean = zMean/length(rPtCloudTrans.Location);
            zMeanMM = zMean*1000;
            fprintf('- Mean difference (actual - target): %.2f mm\n', zMeanMM);
            zStdMM = std(zDiffs*1000);
            fprintf('- Std difference (actual - target): %.2f mm\n', zStdMM);
       end
   end
end