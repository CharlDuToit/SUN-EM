%clear;
Const = sunem_initialise('two_triangles',false);
Const.FEKOmatfilename          = 'two_triangles.mat'; 
Const.FEKOstrfilename          = 'two_triangles.str';
Const.FEKOrhsfilename          = 'two_triangles.rhs'; % ?
Const.FEKOoutfilename          = 'two_triangles.out'; % 
Const.FEKOefefilename          = 'two_triangles.efe'; % ?
Const.FEKOffefilename          = 'two_triangles.ffe'; % ?


[Const, zMatrices, yVectors, xVectors] = extractFEKOMoMmatrixEq(Const);
[Const, Solver_setup] = parseFEKOoutfile(Const, yVectors);
Const.QUAD_PTS = 12; %12


Const.MLMoMClusterSizeScale = 1;
Const.MLMoMMinPercentImprov = 0;
Const.MLMoMIncludeRealCalc = 1;
Const.MLMoMConstMeshSize = 1;
Const.SUNEMmlmomstrfilename = 'two_triangles.str';
%Const.SUNEMmomstrfilename = 'two_triangles.str';

%Const.runMLMoMAddTrianglessolver = false;
%Const.runMLMoMsolver = true;
%Const.runMoMsolver = false;

%Const.indicesStruct = mlmom_twoTri.indicesStruct;
%Const.clusterStruct = mlmom_twoTri.clusterStruct;
Const.SUNEMmlmomaddtrianglesstrfilename = 'two_triangles.str';
Const.runMLMoMAddTrianglessolver = true;
[Solution] = runEMsolvers(Const, Solver_setup, zMatrices, yVectors, xVectors);

%mom = Solution.mom;
mlmomAddTriangles = Solution.mlmomAddTriangles;
%mlmomAddTriangles = Solution.mlmomAddTriangles;
%refStruct = []; 
%refStruct.xVectors = xVectors;
%refStruct.yVectors = yVectors;
%refStruct.zMatrices = zMatrices;
%[predictedSetup] = predictSolverSetup(Const,Solver_setup, mlmom,refStruct, 1);


% 
% tic
% [oldAllTerms, oldSingIndices] = fillZmnTermsByEdge(Const,Solver_setup) ;
% oldTermsTime = toc;
% %properties = calcProperties(Solver_setup);
