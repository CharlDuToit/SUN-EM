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
z = zMatrices.values;

mlmom = Solution.mlmom;
%predError = mlmom.predMeanError;
%unityWeightError = mlmom.unityWeightMeanError;
 predNonSingZ = mlmom.weightModels{1,2}.predNonSingZmn;
 unityNonSingZ = mlmom.weightModels{1,2}.unityNonSingZmn;
 refNonSingZ = mlmom.weightModels{1,2}.refNonSingZmn;
 predNonSingRelVal = predNonSingZ ./ refNonSingZ;
 unityNonSingRelVal = unityNonSingZ ./ refNonSingZ;
% unityWeightSignErrorCount = numel(find(unityWeightRelVal < 0));
% predSignErrorCount = numel(find(predRelVal < 0));
% %numUnclassified = mlmom.numUnclassified;
% varPred = var(predZ);
% varUnityWeight = var(unityWeightZ);
 varRelPred = var(predNonSingRelVal);
 varRelUnity = var(unityNonSingRelVal);
 avgPredRelVal = sum(predNonSingRelVal(:))/numel(predNonSingRelVal);
 avgUnityWeightRelVal = sum(unityNonSingRelVal(:))/numel(unityNonSingRelVal);

%gridSize = 500;
%plotTitle = 'Unity weight error';
%unityWeightDiff = mlmom.unityWeightDiff;
%unityWeightDiff = log(abs(mlmom.nonSingUnityWeightZmn ./ mlmom.refNonSingZmn));
%prop = mlmom.nonSingZmnProp;
%plotError(prop, unityWeightDiff, gridSize, plotTitle);

