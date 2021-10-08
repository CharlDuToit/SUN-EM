function [mlmomAddTriangles] = runMLMoMAddTrianglesSolver(Const, Solver_setup, zMatrices)
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
    message_fc(Const,sprintf('Running ML-MoM solver'));
    
    % ======================= VARIABLE CONVENTION =======================

    % ======================= CONSTANTS =======================
    
    %clusterSizeScale = Const.MLMoMClusterSizeScale;
    %minPercentImprov = Const.MLMoMMinPercentImprov;
    %includeRealCalc = Const.MLMoMIncludeRealCalc;
    numFreq = zMatrices.numFreq;
    %edgeLengths = Solver_setup.rwg_basis_functions_length_m;
    %numOldEdges = Solver_setup.num_mom_basis_functions;
    
    addBias = 0;
    minClusterSize = 8;
    singDataThresh = 1e-63; %1e-44 
    maxIter = 6;
    minTol = 0.02;
    
    %useProjectedEdges = 0;

    % ======================= EXTRACT OLD ======================= 
    % Extract old terms and properties
    %message_fc(Const, sprintf('  Calculating old solver setup matrices'));
    tic;
    [oldTerms, oldSingInd] = fillZmnTermsByEdge(Const,Solver_setup);
    oldCentreDistances = calcCentreDistance(Solver_setup );
    oldProperties = calcProperties(Solver_setup);
    oldEdgeLengths = Solver_setup.rwg_basis_functions_length_m;
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
     newCentreDistances = calcCentreDistance(new_solver_setup );
     [newTerms, newSingInd, newLinkOld, groupIndices] = ...
         projectOldSolverSetup(new_solver_setup, Solver_setup,newEdgeLinkOldEdge, newEdgeParallelExternalEdgeLinkOldInternalEdge, oldTerms, oldSingInd, newCentreDistances, oldCentreDistances );
     newProperties = calcProperties(new_solver_setup);
     projectTime = toc;
    
     % ======================= INIT CLUSTERS =======================
     % Initialise cluster means from new and old properties
     message_fc(Const, sprintf('  Initializing clusters '));
     tic;
     groupMeans = initGroupMeans(newLinkOld, newProperties, oldProperties, newEdgeLengths, oldEdgeLengths, groupIndices, minClusterSize);
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
     [refZMatrices] = FillZMatrixByEdge(Const,new_solver_setup) ;
     refZmn = refZMatrices.values;
     refZmnTime = toc;

     % ======================= CALC WEIGHTS =======================
     % Calculate weights for real and imag at each frequency
     message_fc(Const, sprintf('  Calculating regression models '));
     tic;
     
     weightModels = cell(numFreq, 2); % real and imaginary models
     predZmn = zeros(numNewEdges, numNewEdges,numFreq );
     projZmn = zeros(numNewEdges, numNewEdges,numFreq );
     for f = 1:numFreq
         [weightModels{f, 1},projZmnReal, predZmnReal ] = calcGroupWeights(groupIndices, real(newTerms(:,:,:,f)), real(refZmn(:,:,f)), singDataThresh);
         [weightModels{f, 2},projZmnImag ,predZmnImag ] = calcGroupWeights(groupIndices, imag(newTerms(:,:,:,f)), imag(refZmn(:,:,f)), singDataThresh);
         projZmn(:, :, f) = projZmnReal + 1i * projZmnImag;
         predZmn(:, :, f) = predZmnReal + 1i * predZmnImag;
     end
     regressionTime = toc;
     message_fc(Const, sprintf('  Finished ML-MOM training '));
     
     % ======================= FROBENIUS NORMS  =======================
     %real imag complex
     predFrobNorms = zeros(numFreq, 3);
     projFrobNorms = zeros(numFreq, 3);
     refFrobNorms = zeros(numFreq, 3);
     for i = 1:numFreq
         predFrobNorms(numFreq,1) = calcFrobNorm(real(predZmn(:,:,numFreq)));
         predFrobNorms(numFreq,2) = calcFrobNorm(imag(predZmn(:,:,numFreq)));
         predFrobNorms(numFreq,3) = calcFrobNorm(predZmn(:,:,numFreq));
         
         projFrobNorms(numFreq,1) = calcFrobNorm(real(projZmn(:,:,numFreq)));
         projFrobNorms(numFreq,2) = calcFrobNorm(imag(projZmn(:,:,numFreq)));
         projFrobNorms(numFreq,3) = calcFrobNorm(projZmn(:,:,numFreq));
         
         refFrobNorms(numFreq,1) = calcFrobNorm(real(refZmn(:,:,numFreq)));
         refFrobNorms(numFreq,2) = calcFrobNorm(imag(refZmn(:,:,numFreq)));
         refFrobNorms(numFreq,3) = calcFrobNorm(refZmn(:,:,numFreq));
     end
     
%     % ======================= UPDATE STRUCT  ======================= 
      mlmomAddTriangles = [];
      
%     %cells
      mlmomAddTriangles.weightModels =weightModels;

      %structures
      mlmomAddTriangles.new_solver_setup = new_solver_setup;
      mlmomAddTriangles.groupIndices = groupIndices;
      mlmomAddTriangles.groupMeans = groupMeans;
      
%     %matrices
      mlmomAddTriangles.refZmn = refZmn;
      mlmomAddTriangles.projZmn = projZmn;
      mlmomAddTriangles.predZmn = predZmn;
      
      mlmomAddTriangles.refFrobNorms = refFrobNorms;
      mlmomAddTriangles.projFrobNorms = projFrobNorms;
      mlmomAddTriangles.predFrobNorms = predFrobNorms;
      
      mlmomAddTriangles.newLinkOld = newLinkOld;
      
      mlmomAddTriangles.newProperties = newProperties;
      mlmomAddTriangles.newTerms = newTerms;
      mlmomAddTriangles.newSingInd = newSingInd;

      mlmomAddTriangles.oldProperties = oldProperties;
      mlmomAddTriangles.oldTerms = oldTerms;
      mlmomAddTriangles.oldSingInd = oldSingInd;
      
      %scalars
      
      %constants
      mlmomAddTriangles.numFreq = numFreq;
      mlmomAddTriangles.quadPts = Const.QUAD_PTS;
      mlmomAddTriangles.singDataThresh = singDataThresh;
      %timing
      mlmomAddTriangles.oldZmnTime = oldZmnTime;
      mlmomAddTriangles.addTrianglesTime = addTrianglesTime;
      mlmomAddTriangles.projectTime = projectTime;
      mlmomAddTriangles.initClustersTime = initClustersTime;
      mlmomAddTriangles.calcClustersTime = calcClustersTime;
      mlmomAddTriangles.refZmnTime = refZmnTime;
      mlmomAddTriangles.regressionTime = regressionTime;
      

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

%     
%     % ------------ CALCULATE Z MATRICES  ------------
%     message_fc(Const, sprintf('  Predicting z-Matrices'));
%     tic;
%     %includeRealCalc = 1;
%     returnUnity = 1;
%     [predZmn, unityZmn] = predictExtractedTerms( mlmom, selfZmnTerms, triZmnTerms,...
%         nonSingZmnTerms, nonSingZmnProp, oldSingInd,edgeLengths,returnUnity);
%     mlmom.predZmn = predZmn;
%     mlmom.unityZmn = unityZmn;
%     predictCalcTime = toc;


%     % ======================= DISPLAY TIMING DATA =======================

%     message_fc(Const, sprintf('  Input z-Matrices calculation time: %.3f s',oldUnityZmnCalcTime));
%     message_fc(Const, sprintf('  Clustering calculation time: %.3f s',clusterCalcTime));
%     message_fc(Const, sprintf('  Regression calculation time: %.3f s',regressionCalcTime));
%     message_fc(Const, sprintf('  Prediction calculation time: %.3f s',predictCalcTime));
%     message_fc(Const, sprintf('  Input z-Matrices + Prediction calculation time: %.3f s',oldUnityZmnCalcTime + predictCalcTime));
%        
%     message_fc(Const,...
%         '------------------------------------------------------------------------------------');

end