%clear;
Const = sunem_initialise('square_plate',false);
Const.FEKOmatfilename          = 'square_plate.mat'; 
Const.FEKOstrfilename          = 'square_plate.str';
Const.FEKOrhsfilename          = 'square_plate.rhs'; % ?
Const.FEKOoutfilename          = 'square_plate.out'; % 
Const.FEKOefefilename          = 'square_plate.efe'; % ?
Const.FEKOffefilename          = 'square_plate.ffe'; % ?


[Const, zMatrices, yVectors, xVectors] = extractFEKOMoMmatrixEq(Const);
[Const, Solver_setup] = parseFEKOoutfile(Const, yVectors);

Const.QUAD_PTS = 1;
Const.runMLMoMsolver = true;
Const.MLMoMClusterSizeScale = 1;
Const.MLMoMMinPercentImprov = 0; %6 for 1 quad pt
Const.MLMoMIncludeRealCalc = 1;
Const.MLMoMConstMeshSize = 1;
Const.SUNEMmlmomstrfilename = 'square_plate.str';

%yVectors.values = ones(343, 1);
%xVectros.Isol = zMatrices.values\yVectors.values;

%xVectros.Isol =ones(343, 1);
%yVectors.values =zMatrices.values * xVectros.Isol;
[Solution] = runEMsolvers(Const, Solver_setup, zMatrices, yVectors, xVectors);
mlmom = Solution.mlmom;

%%========COMPARE TIMING WITH 12 QUAD PTS=========
    %tic
    %Const.QUAD_PTS = 3;
    %[refZMatrices] = FillZMatrixByEdge(Const,Solver_setup) ;
    %refCalcTime = toc;
%  quadX = refZMatrices.values\yVectors.values;
%  [~,quadError] = calcError(xVectors.Isol,quadX);

% [comp_imag] = compareZmn(zMatrices.values, refZMatrices.values, refZMatrices.values, mlmom.singInd, 1, 0);
%%======== COMPARE SAME SETUP AGAIN=========
%Reuse training solver setup
%[comp_real] = compareZmn(zMatrices.values, mlmom.predZmn, mlmom.unityZmn,mlmom.singInd, 1, 1) ;
%[comp_imag] = compareZmn(zMatrices.values, mlmom.predZmn, mlmom.unityZmn, mlmom.singInd, 1, 0);


%%========PREDICT DIFFERENT SETUP=========
%first extract from .mat files
% refStruct = [];
% refStruct.xVectors = xVectors_alt;
% refStruct.yVectors = yVectors_alt;
% refStruct.zMatrices = zMatrices_alt;
% Solver_setup_alt = Solver_setup;
% refStruct.xVectors = xVectors;
% refStruct.yVectors = yVectors;
% refStruct.zMatrices = zMatrices;
%[predictedSetup_plate] = predictSolverSetup(Const,Solver_setup_alt, mlmom,refStruct, 1);
%same setup
%[predictedSetup] = predictSolverSetup(Const,Solver_setup, mlmom,zMatrices, 1); 
% same seteup reduced MlMom
%[predictedSetup] = predictSolverSetup(Const,Solver_setup, reducedMLMoM,zMatrices, 1); 

%%========COMPARE PREDICTIEED SETUP TIMING WITH 12 QUAD PTS=========
%  tic
%  Const.QUAD_PTS = 12;
%  [threeQuadZMatrices] = FillZMatrixByEdge(Const,Solver_setup) ;
%  threeQuadCalcTime = toc;
%  [comp_real] = compareZmn(real(zMatrices_small_plate.values), real(threeQuadZMatrices.values), real(threeQuadZMatrices.values), predictedSetup_plate.singInd);
%  [comp_imag] = compareZmn(imag(zMatrices_small_plate.values), imag(threeQuadZMatrices.values), imag(threeQuadZMatrices.values), predictedSetup_plate.singInd);
%  [comp_complex] = compareZmn(zMatrices_small_plate.values, threeQuadZMatrices.values, threeQuadZMatrices.values, predictedSetup_plate.singInd);
