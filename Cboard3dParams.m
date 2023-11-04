classdef Cboard3dParams
   properties
       NumOfHoles
       NumOfRows
       NumOfCols
       SkipEven
       SkipOdd
       Resolution
       LeftColsVisible
       RightColsVisible
       TopRowsVisible
       BottomRowsVisible
   end
   methods
       function obj = SetDefaults(obj)
           obj.NumOfHoles = 30;
           obj.NumOfRows = 3; %6/2
           obj.NumOfCols = 5; %10/2
           obj.SkipEven = 0;
           obj.SkipOdd = 0;
           obj.Resolution = 0.003; %2mm
       end
   end
end