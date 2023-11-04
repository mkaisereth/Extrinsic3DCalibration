function [largestPointInd, smallestPointInd, middle1PointInd] = FindTopRightPoint(holeMedians1)

% fid starting point top right (+x -y)
largestPointInd = 0;
largestPointVal = -Inf;
smallestPointInd = 0;
smallestPointVal = Inf;
middle1PointInd = 0;
middle1PointVal = -Inf;
for i=1:length(holeMedians1)
    currPtVal = holeMedians1(i,1)-holeMedians1(i,2);
    if(currPtVal > largestPointVal)
        largestPointInd = i;
        largestPointVal = currPtVal;
    end
    currPtVal = holeMedians1(i,1)-holeMedians1(i,2);
    if(currPtVal < smallestPointVal)
        smallestPointInd = i;
        smallestPointVal = currPtVal;
    end
    % +x +y => bottom right
    currPtVal = holeMedians1(i,1)+holeMedians1(i,2);
    if(currPtVal > middle1PointVal)
        middle1PointInd = i;
        middle1PointVal = currPtVal;
    end
end

end