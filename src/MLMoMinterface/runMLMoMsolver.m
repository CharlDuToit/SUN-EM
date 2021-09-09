function [mlmom] = runMLMoMsolver(Const, Solver_setup, zMatrices)
    %runMLMoMsolver
    %   Usage:
    %       [mlmom] = runMoMsolver(Const, Solver_setup, zMatrices)
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
    
    % ------------ INIT MLMOM ------------
    
    % Variable convention:
    % dataset + singularism + ?
    %
    % dataset : ref, pred, unity
    % singularism : nonSing, sing ( 2 or 3 triangles), self, tri
    % ? : Context should be clear, first 2 parts could be missing as well
    %
    % cluster variables always involves non singular triangles
   
     
    mlmom = [];
    addBias = 1;
    singDataThresh = 1e-46;
    numFreq = zMatrices.numFreq;
    weightModels = cell(numFreq, 2); % real and imaginary models
    
    % ------------ EXTRACT  ------------ 
    
    % Extract training data
    
    [selfZmnTerms, triZmnTerms, nonSingZmnTerms, nonSingZmnProp, singInd] = extractZmnInfo(Const, Solver_setup);
    [refSelfZmn, refTriZmn, refNonSingZmn] = extractRefZmn(zMatrices.values, singInd );
    [~ , numTerms] = size(nonSingZmnTerms); 
    
    %message_fc(Const, sprintf('  Finished Z-matrix training data calculation '));
    
    % Extract geometric properties
    
    edgeLengths = Solver_setup.rwg_basis_functions_length_m;   
    maxDist = max(nonSingZmnProp(:,1));
    avgEdgeLength = sum(edgeLengths)/numel(edgeLengths);
    threshDist = 0.03584296*maxDist + 1.769794721407624*avgEdgeLength;
    
    % Configure clustering
    
    minExtraInitialClusterPoints = 0;
    minInitialClusterPoints = numTerms + minExtraInitialClusterPoints + addBias;
    minInitialClusterPointsForDynamicRegion = minInitialClusterPoints - 0;
    clusterMinTol = 0.02;
    clusterMaxIter = 6;
    clusterPropScale = [3 1 1];
    
    % Cluster geometric properties
    
    message_fc(Const, sprintf('  Calculating cluster means '));
    
    [clusterMeans, numInitClusterPoints, numInitUnclusteredPoints] = initClusterMeans(nonSingZmnProp,threshDist, minInitialClusterPoints, minInitialClusterPointsForDynamicRegion);
    [clusterMeans, clusterInd, clusterError, clusterNumIter, clusterActualTol] =...
        calcClusterMeans(nonSingZmnProp,clusterMeans, clusterPropScale, clusterMaxIter, clusterMinTol);
    [numClusters , ~] = size(clusterMeans);
    
    % ------------ REAL AND IMAG MODEL FOR EACH FREQ  ------------ 
    message_fc(Const, sprintf('  Calculating regression models '));
    
    for i = 1:numFreq
        message_fc(Const, sprintf('    Processing frequency %d of %d  ',i,numFreq))
        
        weightModels{i, 1} = calcWeightsModel(refSelfZmn(:, i), refTriZmn(:, i), refNonSingZmn(:, i),...
    selfZmnTerms(:,:, i), triZmnTerms(:,:, i), nonSingZmnTerms(:,:, i), nonSingZmnProp,clusterInd,singDataThresh, addBias, 1);

        weightModels{i, 2} = calcWeightsModel(refSelfZmn(:, i), refTriZmn(:, i), refNonSingZmn(:, i),...
    selfZmnTerms(:,:, i), triZmnTerms(:,:, i), nonSingZmnTerms(:,:, i), nonSingZmnProp,clusterInd,singDataThresh, addBias, 0);

    end
    

    
    % ------------ UPDATE STRUCT  ------------ 
    
    %cells
    mlmom.weightModels =weightModels;
    
    %matrices
    mlmom.clusterMeans = clusterMeans;
    mlmom.clusterInd = clusterInd;
    mlmom.clusterPropScale = clusterPropScale;
    
    
    %scalars
    mlmom.numInitClusterPoints = numInitClusterPoints;
    mlmom.numInitUnclusteredPoints = numInitUnclusteredPoints;
    mlmom.numFreq = numFreq;
    mlmom.clusterError = clusterError;
    mlmom.clusterNumIter = clusterNumIter;
    mlmom.clusterActualTol = clusterActualTol;
    mlmom.numClusters = numClusters;
    mlmom.maxDist = maxDist;
    mlmom.threshDist = threshDist;
    
    % ------------ CALCULATE Z MATRICES  ------------
    
    includeRealCalc = 0;
    [zMatrices] = predictExtractedTerms( mlmom, selfZmnTerms, triZmnTerms, nonSingZmnTerms, nonSingZmnProp, singInd, includeRealCalc);
    mlmom.zMatrices = zMatrices;
    
    message_fc(Const, sprintf('  Finished ML-MOM training '));
    message_fc(Const,...
        '------------------------------------------------------------------------------------');
    
end

% Threshold regression
%                 in = [0.1 0.005; 0.1 0.01; 0.1 0.015 ; 0.2 0.01 ; 0.2 0.02; 0.2 0.03;...
%                     0.4 0.02; 0.4 0.04; 0.4 0.06; 0.8 0.04; 0.8 0.08; 0.8 0.12; 1.6 0.08; 1.6 0.16; 1.6 0.24];
%                 out = [0.015; 0.025; 0.03; 0.04; 0.058; 0.068; 0.06; 0.1;0.131;0.12;0.19;0.23;0.2;0.38;0.5];
%                 w = ( (in' *in) \ in')* out ;
%                 w = [0.043777696318019  1.769794721407624 0.006263888888889];
%                 OR w = [0.049472140762463 1.769794721407624]
%                 Divide w(1) by sqrt(2)
    

    