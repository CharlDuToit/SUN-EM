
Const = sunem_initialise('vivaldi_array',false);

% --------------------------------------------------------------------------------------------------
% Program flow settings
% --------------------------------------------------------------------------------------------------

% Choose the solvers that will be executed
%Const.runMoMsolver       = true;
%C%onst.runCBFMsolver      = true;
Const.runMLMoMsolver      = true;
%Const.runJacobisolver    = false;
%Const.runIFBMoMsolver    = true;

% --------------------------------------------------------------------------------------------------
% Define input files for extracting FEKO data
% --------------------------------------------------------------------------------------------------
Const.FEKOmatfilename          = 'vivaldi_array.mat';
Const.FEKOstrfilename          = 'vivaldi_array.str';
Const.FEKOrhsfilename          = 'vivaldi_array.rhs';
Const.FEKOoutfilename          = 'vivaldi_array.out';
%Const.FEKOefefilename          = 'vivaldi_array.efe';
%Const.FEKOffefilename          = 'vivaldi_array.ffe';

[Const, zMatrices, yVectors, xVectors] = extractFEKOMoMmatrixEq(Const);
[Const, Solver_setup] = parseFEKOoutfile(Const, yVectors);
Const.QUAD_PTS = 12;
Const.MLMoMClusterSizeScale = 1;
Const.MLMoMMinPercentImprov = 0;
Const.MLMoMIncludeRealCalc = 1;
Const.MLMoMConstMeshSize = 0;
Const.SUNEMmlmomstrfilename = 'vivaldi_array.str';
%[Solution] = runEMsolvers(Const, Solver_setup, zMatrices, yVectors, xVectors);
%mlmom = Solution.mlmom;

% first extract .mat files for ref Z and solver setup

 refStruct = [];
 refStruct.xVectors = xVectors;
 refStruct.yVectors = yVectors;
 refStruct.zMatrices = zMatrices;
[predictedSetup] = predictSolverSetup(Const, Solver_setup, reducedMLMoM, refStruct, 1);
%writeSolToFile(Const, predictedSetup);
