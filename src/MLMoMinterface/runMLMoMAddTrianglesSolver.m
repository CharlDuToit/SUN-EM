function [mlmomAddTriangles] = runMLMoMAddTrianglesSolver(Const, Solver_setup)
    %runMLMoMsolver_addTriangles
    %   Usage:
    %       [mlmom] = runMLMoMsolver_addTriangles(Const, Solver_setup, zMatrices)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       Solver_setup
    %           Solver specific struct, e.g. frequency range, basis function details, geometry details    
    %       zMatrices
    %           The reference Z-matrix data. This MUST be from FEKO (extracted from the *.mat file)
    %           or some other accurate reference simulation data.
    %
    %   Output Arguments:
    %       mlmom
    %           Structs containing trained ML-MoM 
    %
    %   Description:
    %
    %       Runs the ML-MoM solution based on the Z that was read / parsed.
    %       Geometric assumptions:
    %       1. Planar body
    %       2. Consistent mesh size
    %       Algorithm:
    %       For each frequency, split into real and imaginary MLR models
    %       1. Performs quadrature integration and stores 8 terms
    %       2. Group data by 2 , 3  and 4 unique triangles
    %       3. For 4 triangles, performs clustering
    %       4. Use MLR to find optimal term weights of each group and cluster
    %
    %       Trained model can be used to predict a different Solver_setup,
    %       but with the same frequency (and mesh size? to be tested)
    %
    %   =======================
    %   Written by Charl du Toit on August 20, 2021.
    %   Stellenbosch University
    %   Email: 21708886@sun.ac.za

    %narginchk(5,5);
    message_fc(Const,' ');
    message_fc(Const,'------------------------------------------------------------------------------------');
    message_fc(Const,sprintf('Running ML-MoM-add-triangles solver'));
    
    % ======================= VARIABLE CONVENTION =======================

    % ======================= CONSTANTS =======================
    
    %clusterSizeScale = Const.MLMoMClusterSizeScale;
    %minPercentImprov = Const.MLMoMMinPercentImprov;
    %includeRealCalc = Const.MLMoMIncludeRealCalc;
    
    %numFreq = zMatrices.numFreq;
    numFreq = Solver_setup.frequencies.freq_num;
    
    %edgeLengths = Solver_setup.rwg_basis_functions_length_m;
    %numOldEdges = Solver_setup.num_mom_basis_functions;
    
    addBias = 0;
    minClusterSize = 8;
    singDataThresh = 1e-44; %1e-44 %13-60
    maxIter = 6;
    minTol = 0.02;
    
    %useProjectedEdges = 0;

    % ======================= EXTRACT OLD ======================= 
    % Extract old terms and properties
    %message_fc(Const, sprintf('  Calculating old solver setup matrices'));
    tic;
    oldCentreDistanceProperties = createCentreDistanceProperties(Solver_setup );
    [oldTerms, oldSingInd] = fillZmnTermsByEdge(Const,Solver_setup, oldCentreDistanceProperties);
    %[oldTerms, oldSingInd] = fillZmnTermsByEdge(Const,Solver_setup);
    oldEdgeCentreProperties = createEdgeCentreProperties(Solver_setup);
    oldEdgeLengths = Solver_setup.rwg_basis_functions_length_m;
    oldRhoProperties = createRhoProperties(Solver_setup);
    [~ , ~, numTerms] = size(oldTerms);
    oldZmnTime = toc;
    
    % ======================= ADD TRIANGLES ======================= 
    % Add triangles to solver setup
    message_fc(Const, sprintf('  Creating new solver setup '));
    tic;
    [new_solver_setup, newEdgeLinkOldEdge, newEdgeParallelExternalEdgeLinkOldInternalEdge] = addTriangles(Solver_setup);
    newEdgeLengths = new_solver_setup.rwg_basis_functions_length_m;
    numNewEdges = new_solver_setup.num_mom_basis_functions;
    addTrianglesTime = toc;
    
     % ======================= PROJECT ======================= 
     % Project old terms to new terms
     message_fc(Const, sprintf('  Projecting old solver setup matrices '));
     tic;
     newCentreDistancesProperties = createCentreDistanceProperties(new_solver_setup );
     [newTerms, newSingInd, newLinkOld, groupIndices] = ...
         projectOldSolverSetup(new_solver_setup, Solver_setup,newEdgeLinkOldEdge, newEdgeParallelExternalEdgeLinkOldInternalEdge, oldTerms, oldSingInd, newCentreDistancesProperties, oldCentreDistanceProperties );
     newEdgeCentreProperties = createEdgeCentreProperties(new_solver_setup);
     newRhoProperties = createRhoProperties(new_solver_setup);
     projectTime = toc;
    
     % ======================= INIT CLUSTERS =======================
     % Initialise cluster means from new and old properties
     message_fc(Const, sprintf('  Initializing clusters '));
     tic;
     groupMeans = initGroupMeans(groupIndices, newLinkOld, newEdgeCentreProperties, oldEdgeCentreProperties, newEdgeLengths, oldEdgeLengths ,newRhoProperties,oldRhoProperties, minClusterSize);
     initClustersTime = toc;
     
     % ======================= CALC CLUSTERS =======================
     % Find nearest cluster for each new matrix entry
     message_fc(Const, sprintf('  Calculating clusters '));
     tic;
     [groupMeans, groupIndices] =calcGroupMeans(groupMeans, groupIndices, maxIter, minTol);
     calcClustersTime = toc;

     % ======================= CALC REFERENCE ZMN =======================
     %message_fc(Const, sprintf('  Calculating reference matrices '));
     tic;
     quadPts = Const.QUAD_PTS;
     Const.QUAD_PTS = 12;
     [refZMatrices] = FillZMatrixByEdge(Const,new_solver_setup) ;
     refZmn = refZMatrices.values;
     Const.QUAD_PTS = quadPts;
     refZmnTime = toc;

     % ======================= CALC WEIGHTS, X, Y =======================
     % Calculate weights for real and imag at each frequency
     message_fc(Const, sprintf('  Calculating regression models '));
     tic;
     %vectors
     y = zeros(numNewEdges, numFreq);
     refX = zeros(numNewEdges, numFreq);
     predX = zeros(numNewEdges, numFreq);
     %erros
     predXError = zeros(numFreq, 3);
     predFrobNorms = zeros(numFreq, 3);
     projFrobNorms = zeros(numFreq, 3);
     refFrobNorms = zeros(numFreq, 3);
     predRelNormPercentError = zeros(numFreq, 3);
     errorSummary = cell(numFreq, 3); 
     predRowRelNormPercentError = zeros(numNewEdges, numFreq);
     %weights
     weightModels = cell(numFreq, 2); % real and imaginary models
     predZmn = zeros(numNewEdges, numNewEdges,numFreq );
     projZmn = zeros(numNewEdges, numNewEdges,numFreq );
     EMag = 1;
     theta_0 = 0;
     phi_0 = 0;
             
     %unityX(:,f) = unityZmn(:,:,f)\yVectors.values(:,f);
        
     for f = 1:numFreq

         %weights
         [weightModels{f, 1},projZmnReal, predZmnReal ] = calcGroupWeights(groupIndices, real(newTerms(:,:,:,f)), real(refZmn(:,:,f)), singDataThresh);
         [weightModels{f, 2},projZmnImag ,predZmnImag ] = calcGroupWeights(groupIndices, imag(newTerms(:,:,:,f)), imag(refZmn(:,:,f)), singDataThresh);
         projZmn(:, :, f) = projZmnReal + 1i * projZmnImag;
         predZmn(:, :, f) = predZmnReal + 1i * predZmnImag;
         
         %vectors
         y(:,f) = FillVVector(Const, new_solver_setup, EMag,theta_0,phi_0);
         predX(:,f) = predZmn(:,:,f)\y(:,f);
         refX(:,f) = refZmn(:,:,f)\y(:,f);
         [~,predXError(f,1)] = calcError(real(refX(:,f)),real(predX(:,f)));
         [~,predXError(f,2)] = calcError(imag(refX(:,f)),imag(predX(:,f)));
         [~,predXError(f,3)] = calcError(refX(:,f),predX(:,f));
         
         %errors
         [compReal] = compareZmn(real(refZmn(:, :, f)), real(predZmn(:, :, f)), real(projZmn(:, :, f)), newSingInd);
         [compImag] = compareZmn(imag(refZmn(:, :, f)), imag(predZmn(:, :, f)), imag(projZmn(:, :, f)), newSingInd);
         [compComplex] = compareZmn(refZmn(:, :, f), predZmn(:, :, f), projZmn(:, :, f), newSingInd);
         errorSummary{f,1} = compReal;
         errorSummary{f,2} = compImag;
         errorSummary{f,3} = compComplex;
         
         predRelNormPercentError(f,1) = compReal.predRelNormPercentError;
         predRelNormPercentError(f,2)= compImag.predRelNormPercentError;
         [~, predRelNormPercentError(f,3)] = calcError(refZmn(:, :, f), predZmn(:, :, f));

         for mm = 1:numNewEdges
             [~, predRowRelNormPercentError(mm)] = calcError(refZmn(mm, :, f), predZmn(mm, :, f));
         end
         
         predFrobNorms(f,1) = compReal.predFrobNorm;
         predFrobNorms(f,2) = compImag.predFrobNorm;
         predFrobNorms(f,3) = calcFrobNorm(predZmn(:,:,f));
         
         projFrobNorms(f,1) = compReal.unityFrobNorm;
         projFrobNorms(f,2) = compImag.unityFrobNorm;
         projFrobNorms(f,3) = calcFrobNorm(projZmn(:,:,f));
         
         refFrobNorms(f,1) = compReal.refFrobNorm;
         refFrobNorms(f,2) = compImag.refFrobNorm;
         refFrobNorms(f,3) = calcFrobNorm(refZmn(:,:,f));
         

     end
     regressionTime = toc;
     message_fc(Const, sprintf('  Finished ML-MOM training '));
     
     
     %======================= PREDICTION TIMING  =======================
     message_fc(Const, sprintf('  Obtaining prediction timing '));
     %==================================================================
     
     % ======================= ASSIGN  ======================= 
     tic
     groupIndices = assignGroupClusters(groupMeans, groupIndices,newLinkOld, newEdgeCentreProperties, oldEdgeCentreProperties,newEdgeLengths, oldEdgeLengths, newRhoProperties, oldRhoProperties);
     assignTime = toc;
     
     %======================= MULTIPLY WEIGHTS  ======================= 
     tic;
     predZmn = zeros(numNewEdges, numNewEdges,numFreq );
     projZmn = zeros(numNewEdges, numNewEdges,numFreq );
     for f = 1:numFreq
         [~, projZmnReal, predZmnReal] = applyGroupWeights(weightModels{f,1}, groupIndices, real(newTerms(:,:,:,f)), real(refZmn(:,:,f)), 0);
         [~, projZmnImag, predZmnImag] = applyGroupWeights(weightModels{f,2}, groupIndices, imag(newTerms(:,:,:,f)), imag(refZmn(:,:,f)), 0);
         projZmn(:, :, f) = projZmnReal + 1i * projZmnImag;
         predZmn(:, :, f) = predZmnReal + 1i * predZmnImag;
     end
     multiplyTime = toc;
     
    % ======================= DISPLAY TIMING DATA =======================
    trainingTime = oldZmnTime + addTrianglesTime + projectTime + initClustersTime + calcClustersTime + refZmnTime + regressionTime;
    predictTime = oldZmnTime + addTrianglesTime + projectTime + assignTime + multiplyTime;

    message_fc(Const, sprintf('  Input z-Matrices calculation time: %.3f s',oldZmnTime));
    message_fc(Const, sprintf('  Add triangles time: %.3f s',addTrianglesTime));
    message_fc(Const, sprintf('  Projection time: %.3f s',projectTime));
    message_fc(Const, sprintf('  Cluster initialization time: %.3f s',initClustersTime));
    message_fc(Const, sprintf('  Cluster calculation  time: %.3f s',calcClustersTime));
    message_fc(Const, sprintf('  Reference z-Matrices time: %.3f s',refZmnTime));
    message_fc(Const, sprintf('  Regression time: %.3f s',regressionTime));
    message_fc(Const, sprintf('  Training time: %.3f s',trainingTime));
    message_fc(Const, sprintf('  Prediction time: %.3f s',predictTime));
    
     
%     % ======================= UPDATE STRUCT  ======================= 
      mlmomAddTriangles = [];
      mlmomAddTriangles.name = 'MLMOM-ADDTRIANGLES';
      mlmomAddTriangles.numSols = numFreq;
      mlmomAddTriangles.freqSamples = Solver_setup.frequencies.samples;
      
%     %cells
      mlmomAddTriangles.weightModels =weightModels;
      mlmomAddTriangles.errorSummary = errorSummary;
      
      %structures
      mlmomAddTriangles.new_solver_setup = new_solver_setup;
      mlmomAddTriangles.groupIndices = groupIndices;
      mlmomAddTriangles.groupMeans = groupMeans;
      
%     %matrices
      %mlmomAddTriangles.frequencies = 
      mlmomAddTriangles.refZmn = refZmn;
      mlmomAddTriangles.projZmn = projZmn;
      mlmomAddTriangles.predZmn = predZmn;
      
      mlmomAddTriangles.y = y;
      mlmomAddTriangles.refX = refX;
      mlmomAddTriangles.Isol = predX;
      mlmomAddTriangles.predXError = predXError;
      
      mlmomAddTriangles.predRelNormPercentError = predRelNormPercentError;
      mlmomAddTriangles.refFrobNorms = refFrobNorms;
      mlmomAddTriangles.projFrobNorms = projFrobNorms;
      mlmomAddTriangles.predFrobNorms = predFrobNorms;
      mlmomAddTriangles.predRowRelNormPercentError = predRowRelNormPercentError;
      
      mlmomAddTriangles.newLinkOld = newLinkOld;
      
      mlmomAddTriangles.newProperties = newEdgeCentreProperties;
      mlmomAddTriangles.newTerms = newTerms;
      mlmomAddTriangles.newSingInd = newSingInd;

      mlmomAddTriangles.oldProperties = oldEdgeCentreProperties;
      mlmomAddTriangles.oldTerms = oldTerms;
      mlmomAddTriangles.oldSingInd = oldSingInd;
      
      %scalars
      
      %constants
      mlmomAddTriangles.numFreq = numFreq;
      mlmomAddTriangles.quadPts = quadPts;
      mlmomAddTriangles.singDataThresh = singDataThresh;
      %timing
      mlmomAddTriangles.oldZmnTime = oldZmnTime;
      mlmomAddTriangles.addTrianglesTime = addTrianglesTime;
      mlmomAddTriangles.projectTime = projectTime;
      mlmomAddTriangles.initClustersTime = initClustersTime;
      mlmomAddTriangles.calcClustersTime = calcClustersTime;
      mlmomAddTriangles.refZmnTime = refZmnTime;
      mlmomAddTriangles.regressionTime = regressionTime;
      mlmomAddTriangles.assignTime = assignTime;
      mlmomAddTriangles.multiplyTime = multiplyTime;
      mlmomAddTriangles.trainingTime = trainingTime;
      mlmomAddTriangles.predictTime = predictTime;
      

%     mlmom.edgeLengths = edgeLengths;
%     mlmom.nonSingZmnProp = nonSingZmnProp;
%     mlmom.nonSingEdgeLabels = nonSingEdgeLabels;
%     mlmom.maxClusterError = maxClusterError;
%     mlmom.clusterMaxEdgeLength = clusterMaxEdgeLength;
%     mlmom.clusterCounts = clusterCounts;
% 
%     %scalars
%     mlmom.includeRealCalc = includeRealCalc;
%     mlmom.clusterSizeScale = clusterSizeScale;

%     mlmom.totalClusterError = totalClusterError;
%     mlmom.clusterNumIter = clusterNumIter;
%     mlmom.clusterActualTol = clusterActualTol;
%     mlmom.numClusters = numClusters;
%     mlmom.maxDist = maxDist;
%     mlmom.threshDist = threshDist;

%     mmInd = mlmomAddTriangles.groupIndices.twoUnique_pos_indices(:,1);
%     nnInd =  mlmomAddTriangles.groupIndices.twoUnique_pos_indices(:,2);
%     ref = zeros(numel(mmInd),1);
%     proj = zeros(numel(mmInd),1);
%     for k = 1:numel(mmInd)
%         ref(k) = mlmomAddTriangles.refZmn(mmInd(k), nnInd(k), 1);
%         proj(k) = mlmomAddTriangles.projZmn(mmInd(k), nnInd(k), 1);
%     end
       
   % if (~isempty(Const.SUNEMmlmomstrfilename))
        %writeSolToFile(Const, mlmomAddTriangles);
   % end%if
    
    message_fc(Const,...
        '------------------------------------------------------------------------------------');

end