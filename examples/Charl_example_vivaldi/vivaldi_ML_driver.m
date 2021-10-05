
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
Const.QUAD_PTS = 1;
Const.MLMoMClusterSizeScale = 1;
Const.MLMoMMinPercentImprov = 0;
Const.MLMoMIncludeRealCalc = 0;
[Solution] = runEMsolvers(Const, Solver_setup, zMatrices, yVectors, xVectors);
mlmom = Solution.mlmom;

% first extract .mat files for ref Z and solver setup

%refZmn_plate = zMatrices_plate.values;
%[predZmn_plate, unityZmn_plate, singInd_plate] = predictSolverSetup(Const,Solver_setup_plate, mlmom,1);
%[comp_real] = compareZmn(refZmn_plate, predZmn_plate,unityZmn_plate,singInd_plate, 1) ;
%[comp_imag] = compareZmn(refZmn_plate, predZmn_plate, unityZmn_plate, singInd_plate, 0);

