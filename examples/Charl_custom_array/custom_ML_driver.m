%clear;
Const = sunem_initialise('custom',false);
Const.FEKOmatfilename          = 'custom.mat'; 
Const.FEKOstrfilename          = 'custom.str';
Const.FEKOrhsfilename          = 'custom.rhs'; % ?
Const.FEKOoutfilename          = 'custom.out'; % 
Const.FEKOefefilename          = 'custom.efe'; % ?
Const.FEKOffefilename          = 'custom.ffe'; % ?


[Const, zMatrices, yVectors, xVectors] = extractFEKOMoMmatrixEq(Const);
[Const, Solver_setup] = parseFEKOoutfile(Const, yVectors);
Const.QUAD_PTS = 12;
%Const.runMLMoMsolver = true;
Const.MLMoMClusterSizeScale = 1;
Const.MLMoMMinPercentImprov = 0;
Const.MLMoMIncludeRealCalc = 1;
Const.MLMoMConstMeshSize = 0;
Const.runMLMoMAddTrianglessolver = true;
Const.SUNEMmlmomstrfilename = 'custom.str';
Const.SUNEMmlmomaddtrianglesstrfilename = 'custom.str';

%Const.indicesStruct = mlmom1.indicesStruct;
%Const.clusterStruct = mlmom1.clusterStruct;
[Solution] = runEMsolvers(Const, Solver_setup, zMatrices, yVectors, xVectors);

mlmomAddTriangles = Solution.mlmomAddTriangles;
%reducedMLMoM = reduceMLMoM(mlmom);





    