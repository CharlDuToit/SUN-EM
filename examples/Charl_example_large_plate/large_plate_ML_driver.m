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
predError = mlmom.predMeanError;
unityWeightError = mlmom.unityWeightMeanError;
predZ = mlmom.predNonSingZmn;
unityWeightZ = mlmom.nonSingUnityWeightZmn;
refZ = mlmom.refNonSingZmn;
predRelVal = predZ ./ refZ;
unityWeightRelVal = unityWeightZ ./ refZ;
unityWeightSignErrorCount = numel(find(unityWeightRelVal < 0));
predSignErrorCount = numel(find(predRelVal < 0));
numUnclassified = mlmom.numUnclassified;
varPred = var(predZ);
varUnityWeight = var(unityWeightZ);
varRelPred = var(predRelVal);
varRelUnityWeight = var(unityWeightRelVal);
avgPredRelVal = sum(predRelVal(:))/numel(unityWeightRelVal);
avgUnityWeightRelVal = sum(unityWeightRelVal(:))/numel(unityWeightRelVal);

