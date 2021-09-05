clear;
Const = sunem_initialise('square_plate',false);
Const.FEKOmatfilename          = 'square_plate.mat'; 
Const.FEKOstrfilename          = 'square_plate.str';
Const.FEKOrhsfilename          = 'square_plate.rhs'; % ?
Const.FEKOoutfilename          = 'square_plate.out'; % 
Const.FEKOefefilename          = 'square_plate.efe'; % ?
Const.FEKOffefilename          = 'square_plate.ffe'; % ?
Const.runMLMoMsolver              = true;
Const.QUAD_PTS = 3;
[Const, zMatrices, yVectors, xVectors] = extractFEKOMoMmatrixEq(Const);
[Const, Solver_setup] = parseFEKOoutfile(Const, yVectors);
[Solution] = runEMsolvers(Const, Solver_setup, zMatrices, yVectors, xVectors);

mlmom = Solution.mlmom;
predError = mlmom.predError;
unityWeightError = mlmom.unityWeightError;
predZ = mlmom.predNonSingZmn;
unityWeightZ = mlmom.nonSingZmnUnityWeight;
refZ = mlmom.refNonSingZmn;
predRelVal = predZ ./ refZ;
unityWeightRelVal = unityWeightZ ./ refZ;
unityWeightSignErrorCount = numel(find(unityWeightRelVal < 0));
predSignErrorCount = numel(find(predRelVal < 0));

