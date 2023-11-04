function LarsenStudyCheckerboards3d_func(basePath, baseFoldercb3d, patientNumscb3d, systemCombos1, systemCombos2, systemNames, systemNames2, directVersion, printLevel, numOfCheckerboards)

cb3dFolders = [append(basePath, "\cb3d"), append(basePath, "\3dcb")];
if directVersion
    cb3dFolders(end+1) = baseFoldercb3d;
else
    if ~isempty(baseFoldercb3d)
        cb3dFolders = [cb3dFolders, append(baseFoldercb3d, "\cb3d"), append(baseFoldercb3d, "\3dcb")];
        for i=patientNumscb3d
            cb3dFolders(end+1) = append(baseFoldercb3d, "\", string(i));
        end
    end
end

for basePath=cb3dFolders

    outputPath = append(basePath, "\Output");
    if ~exist(basePath, "dir")
        continue;
    end
    if ~exist(outputPath, "dir")
        mkdir(outputPath);
    end
    
    cb3dParams = Cboard3dParams;
    cb3dParams = cb3dParams.SetDefaults();
    cb3dParams.NumOfHoles = 18;
    cb3dParams.NumOfRows = 3;
    cb3dParams.NumOfCols = 3;
    
    for sni=1:length(systemNames)
        fullSystemNames(sni) = append(systemNames(sni), systemNames2(sni));
        plyDict.(fullSystemNames(sni)) = [];
    end
    
    % get all files and sort
    plyList = dir(append(basePath, "\*.ply"));
    for i=1:length(plyList)
        for sni=1:length(systemNames)
            if startsWith(plyList(i).name, append(systemNames(sni), "_")) && endsWith(plyList(i).name, append(systemNames2(sni), ".ply"))
                % add id (day)
                plyList(i).id = extractAfter(plyList(i).name, strlength(systemNames(sni))+1);
                plyList(i).id = string(plyList(i).id(1:8));
                plyDict.(fullSystemNames(sni)) = [plyDict.(fullSystemNames(sni)), plyList(i)];
            end
        end
    end
    
    % check uniqueness of id
    for sni=1:length(systemNames)
        idarr = [plyDict.(fullSystemNames(sni)).id];
        uniqueIdArr = unique(idarr);
        if length(idarr) ~= length(uniqueIdArr)
            disp('Warning: ids not unique')
            %keyboard
        end
    end
    
    for sci=1:length(systemCombos1)
        ids1 = [plyDict.(systemCombos1(sci)).id];
        ids2 = [plyDict.(systemCombos2(sci)).id];
        if length(ids1) ~= length(ids2)
            disp('Warning: number of ids different')
            keyboard
        end
        for ii=1:min(numOfCheckerboards, length(ids1))
            close all;
            if strcmp(ids1(ii), ids2(ii)) == 0
                disp('Warning: order of ids different')
                keyboard
            end
            % now do the calibration for each combo
            plyDictCombo1 = plyDict.(systemCombos1(sci));
            plyDictCombo2 = plyDict.(systemCombos2(sci));
            
            pc1Path = append(plyDictCombo1(ii).folder, "\", plyDictCombo1(ii).name);
            pc1 = pcread(pc1Path);
            % second is reference coordinate system
            pc2Path = append(plyDictCombo2(ii).folder, "\", plyDictCombo2(ii).name);
            pc2_ref = pcread(pc2Path);
            comboName = append(systemCombos1(sci), "-", systemCombos2(sci), "-", ids1(ii));
            % first argument is reference system (thus swapped)
            Checkerboards3d_estimateT(pc2Path, pc1Path, pc2_ref, pc1, cb3dParams, printLevel, outputPath, comboName);
            pause(0.1)
        end
    end
end

end