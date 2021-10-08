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
Const.MLMoMMinPercentImprov = 0;
Const.MLMoMIncludeRealCalc = 0;
[Solution] = runEMsolvers(Const, Solver_setup, zMatrices, yVectors, xVectors);

mlmom = Solution.mlmom;

%plot(Solver_setup.nodes_xyz(:,1), Solver_setup.nodes_xyz(:,2), '.', 'markerSize', 20);
%[new_solver_setup] = addTriangles(Solver_setup);
%plot(new_solver_setup.nodes_xyz(:,1), new_solver_setup.nodes_xyz(:,2), '.', 'markerSize', 20);

%first extract from .mat files

%refZmn_plate = zMatrices_small_plate.values;
%[predZmn_plate, unityZmn_plate, singInd_plate] = predictSolverSetup(Const,Solver_setup_plate, mlmom,1);
%[comp_real] = compareZmn(refZmn_plate, predZmn_plate,unityZmn_plate,singInd_plate, 1, 1) ;
%[comp_imag] = compareZmn(refZmn_plate, predZmn_plate, unityZmn_plate, singInd_plate, 1, 0);

%Reuse training solver setup
%[comp_real] = compareZmn(zMatrices.values, mlmom.predZmn, mlmom.unityZmn,mlmom.singInd, 1, 1) ;
%[comp_imag] = compareZmn(zMatrices.values, mlmom.predZmn, mlmom.unityZmn, mlmom.singInd, 1, 0);