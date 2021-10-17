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
Const.MLMoMMinPercentImprov = 6;
Const.MLMoMIncludeRealCalc = 0;
[Solution] = runEMsolvers(Const, Solver_setup, zMatrices, yVectors, xVectors);

mlmom = Solution.mlmom;

%%========COMPARE TIMING WITH 12 QUAD PTS=========
% tic
% Const.QUAD_PTS = 12;
% [refZMatrices] = FillZMatrixByEdge(Const,Solver_setup) ;
% refCalcTime = toc;
% [comp_imag] = compareZmn(zMatrices.values, refZMatrices.values, refZMatrices.values, mlmom.singInd, 1, 0);
%%======== COMPARE SAME SETUP AGAIN=========
%Reuse training solver setup
%[comp_real] = compareZmn(zMatrices.values, mlmom.predZmn, mlmom.unityZmn,mlmom.singInd, 1, 1) ;
%[comp_imag] = compareZmn(zMatrices.values, mlmom.predZmn, mlmom.unityZmn, mlmom.singInd, 1, 0);


%%========PREDICT DIFFERENT SETUP=========
%first extract from .mat files

%refZmn_plate = zMatrices_small_plate.values;
%[predictedSetup_plate] = predictSolverSetup(Const,Solver_setup_plate, mlmom,zMatrices_small_plate, 1, 1);
 

%%========COMPARE PREDICTIEED SETUP TIMING WITH 12 QUAD PTS=========
%  tic
%  Const.QUAD_PTS = 12;
%  [threeQuadZMatrices] = FillZMatrixByEdge(Const,Solver_setup_plate) ;
%  threeQuadCalcTime = toc;
%  [comp_real] = compareZmn(real(zMatrices_small_plate.values), real(threeQuadZMatrices.values), real(threeQuadZMatrices.values), predictedSetup_plate.singInd);
%  [comp_imag] = compareZmn(imag(zMatrices_small_plate.values), imag(threeQuadZMatrices.values), imag(threeQuadZMatrices.values), predictedSetup_plate.singInd);
%  [comp_complex] = compareZmn(zMatrices_small_plate.values, threeQuadZMatrices.values, threeQuadZMatrices.values, predictedSetup_plate.singInd);
