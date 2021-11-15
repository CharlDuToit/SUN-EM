%clear;
Const = sunem_initialise('square_plate',false);
Const.FEKOmatfilename          = 'square_plate.mat'; 
Const.FEKOstrfilename          = 'square_plate.str';
Const.FEKOrhsfilename          = 'square_plate.rhs'; % ?
Const.FEKOoutfilename          = 'square_plate.out'; % 
Const.FEKOefefilename          = 'square_plate.efe'; % ?
Const.FEKOffefilename          = 'square_plate.ffe'; % ?
Const.SUNEMmlmomstrfilename    = 'square_plate.str';
Const.runMLMoMsolver              = true;

[Const, zMatrices, yVectors, xVectors] = extractFEKOMoMmatrixEq(Const);
[Const, Solver_setup] = parseFEKOoutfile(Const, yVectors);

Const.MLMoMClusterSizeScale = 1;
Const.QUAD_PTS = 12;
Const.MLMoMMinPercentImprov = 0;
Const.MLMoMIncludeRealCalc = 1;
Const.MLMoMConstMeshSize = 1;
Const.SUNEMmlmomstrfilename = 'square_plate.str';
%[Solution] = runEMsolvers(Const, Solver_setup, zMatrices, yVectors, xVectors);

%mlmom = Solution.mlmom;
 %refStruct = [];
 %refStruct.xVectors = xVectors;
 %refStruct.yVectors = yVectors;
 %refStruct.zMatrices = zMatrices;
%[predictedSetup] = predictSolverSetup(Const, Solver_setup, reducedMLMoM, refStruct, 1);
%[predictedSetup] = predictSolverSetupAddTriangles(Const, Solver_setup, mlmomAddTriangles, 1);
%============= POST PROCESSING ==============

%CURRENT DISTRIBUTION
numFreq = predictedSetup.numSols;
freqSamples = predictedSetup.freqSamples;

% A-A cut x = 0.5, y varies
indNodes = find(Solver_setup.rwg_basis_functions_shared_edge_centre(:,1) == 0.5);
AAnodes = Solver_setup.rwg_basis_functions_shared_edge_centre(indNodes, :);
[~, sortedInd] = sort(AAnodes(:,2) );
indNodes = indNodes(sortedInd);
sortedAANodes = Solver_setup.rwg_basis_functions_shared_edge_centre(indNodes, :);
%numAAsamples = numel(indNodes);
yCoord = sortedAANodes(:, 2);
%numFreq = 1;
% for f = 1:numFreq
%     FEKOsamples = 120*pi*abs(xVectors.Isol(indNodes,f));
%     MLMoMSamples  = 120*pi*abs(predictedSetup.Isol(indNodes,f));
%     RWGsamples = 120*pi*abs(predictedSetup.unityX(indNodes,f));
%     
%     figure;
%     hold on;
%     grid on;
%     plot(yCoord ,MLMoMSamples,'-o', 'LineWidth',2);
%     plot(yCoord ,FEKOsamples,'-*','LineWidth',2);
%     plot(yCoord ,RWGsamples,'-s','LineWidth',2);
%     title(strcat('AA cut Current f=',num2str(freqSamples(f), 8),'Hz'));
%     legend('MLMoM','FEKO', 'RWG');
%     set(get(gca, 'XLabel'), 'String', ('y [m]'));
%     set(get(gca, 'YLabel'), 'String', '|Jx|/|Hinc|');
%     hold off;
% end
% 
% % B-B cut x = y = 0.5, x varies
% indNodes = find(Solver_setup.rwg_basis_functions_shared_edge_centre(:,2) == 0.5);
% BBnodes = Solver_setup.rwg_basis_functions_shared_edge_centre(indNodes, :);
% [~, sortedInd] = sort(BBnodes(:,1) );
% indNodes = indNodes(sortedInd);
% sortedBBNodes = Solver_setup.rwg_basis_functions_shared_edge_centre(indNodes, :);
% xCoord = sortedBBNodes(:, 1);
% %numFreq = 1;
% for f = 1:numFreq
%     FEKOsamples = 120*pi*abs(xVectors.Isol(indNodes,f));
%     MLMoMSamples  = 120*pi*abs(predictedSetup.Isol(indNodes,f));
%     RWGsamples = 120*pi*abs(predictedSetup.unityX(indNodes,f));
%     
%     figure;
%     hold on;
%     grid on;
%     plot(xCoord ,MLMoMSamples,'-o', 'LineWidth',2);
%     plot(xCoord ,FEKOsamples,'-*','LineWidth',2);
%     plot(xCoord ,RWGsamples,'-s','LineWidth',2);
%     title(strcat('BB-Current f=',num2str(freqSamples(f), 8),'Hz'));
%     legend('MLMoM','FEKO', 'RWG');
%     set(get(gca, 'XLabel'), 'String', ('x [m]'));
%     set(get(gca, 'YLabel'), 'String', '|Jy|/|Hinc|');
%     hold off;
%     
% end

