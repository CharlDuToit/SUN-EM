%clear;
%Const = sunem_initialise('square_plate',false);
Const.FEKOmatfilename          = 'square_plate.mat'; 
Const.FEKOstrfilename          = 'square_plate.str';
Const.FEKOrhsfilename          = 'square_plate.rhs'; % ?
Const.FEKOoutfilename          = 'square_plate.out'; % 
Const.FEKOefefilename          = 'square_plate.efe'; % ?
Const.FEKOffefilename          = 'square_plate.ffe'; % ?


%[Const, zMatrices, yVectors, xVectors] = extractFEKOMoMmatrixEq(Const);
%[Const, Solver_setup] = parseFEKOoutfile(Const, yVectors);
Const.QUAD_PTS = 12;
%Const.runMLMoMsolver              = true;
Const.runMLMoMAddTrianglessolver = true;
Const.MLMoMClusterSizeScale = 1;
Const.MLMoMMinPercentImprov = 2;
Const.MLMoMIncludeRealCalc = 0;
%[Solution] = runEMsolvers(Const, Solver_setup, zMatrices, yVectors, xVectors);
% 
%mlmomAddTriangles = Solution.mlmomAddTriangles;

%plot(Solver_setup.nodes_xyz(:,1), Solver_setup.nodes_xyz(:,2), '.', 'markerSize', 20);

% %======== Predict different setup
%first extract from .mat files
%predictedSetup_plate = predictSolverSetupAddTriangles(Const,Solver_setup_plate, mlmomAddTriangles,1);

% %======== mlmom addtriangles 3qaud
% tic
% Const.QUAD_PTS = 3;
% [threeQuadZMatrices] = FillZMatrixByEdge(Const,mlmomAddTriangles.new_solver_setup) ;
% threeQuadZMatricesCalcTime = toc;
% %unity = threeQuadZMatrices.values
% %pred = mlmomAddTriangles.predZmn
% [comp_real] = compareZmn(real(mlmomAddTriangles.refZmn), real(mlmomAddTriangles.predZmn), real(threeQuadZMatrices.values), mlmomAddTriangles.newSingInd);
% [comp_imag] = compareZmn(imag(mlmomAddTriangles.refZmn), imag(mlmomAddTriangles.predZmn), imag(threeQuadZMatrices.values), mlmomAddTriangles.newSingInd);
% [comp_complex] = compareZmn(mlmomAddTriangles.refZmn, mlmomAddTriangles.predZmn, threeQuadZMatrices.values, mlmomAddTriangles.newSingInd);
% %===========

%======== predicted plate 3 quad
%tic
%Const.QUAD_PTS = 3;
%[threeQuadZMatrices] = FillZMatrixByEdge(Const,predictedSetup_plate.new_solver_setup) ;
%threeQuadZMatricesCalcTime = toc;
%unity = threeQuadZMatrices.values
%pred = mlmomAddTriangles.predZmn
[comp_real] = compareZmn(real(predictedSetup_plate.refZmn), real(predictedSetup_plate.predZmn), real(threeQuadZMatrices.values), predictedSetup_plate.newSingInd);
[comp_imag] = compareZmn(imag(predictedSetup_plate.refZmn), imag(predictedSetup_plate.predZmn), imag(threeQuadZMatrices.values), predictedSetup_plate.newSingInd);
[comp_complex] = compareZmn(predictedSetup_plate.refZmn, predictedSetup_plate.predZmn, threeQuadZMatrices.values, predictedSetup_plate.newSingInd);
%===========




