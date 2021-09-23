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
    addBias = 0;
    singDataThresh = 1e-44;
    useProjectedEdges = 0;
    minPercentImprov = 0;
    numFreq = zMatrices.numFreq;
    weightModels = cell(numFreq, 2); % real and imaginary models
    
    % ------------ EXTRACT  ------------ 
    
    % Extract training data
    
    tic;
    [selfZmnTerms, triZmnTerms, nonSingZmnTerms, nonSingZmnProp,nonSingEdgeLabels, singInd] = extractZmnInfo(Const, Solver_setup);
    unityZmnCalcTime = toc;
    [refSelfZmn, refTriZmn, refNonSingZmn] = extractRefZmn(zMatrices.values, singInd );
    [~ , numTerms] = size(nonSingZmnTerms);
    
    % Extract geometric properties
    
    edgeLengths = Solver_setup.rwg_basis_functions_length_m;
    numEdges = Solver_setup.num_mom_basis_functions;
    maxDist = max(nonSingZmnProp(:,1));
    avgEdgeLength = sum(edgeLengths)/numel(edgeLengths);
    threshDist = 0.03584296*maxDist + 1.769794721407624*avgEdgeLength;
    
    % Configure clustering
    
    minExtraInitialClusterPoints = 0;
    minInitialClusterPoints = numTerms + minExtraInitialClusterPoints + addBias;
    minInitialClusterPointsForDynamicRegion = minInitialClusterPoints - 0;
    clusterMinTol = 0.02;
    clusterMaxIter = 6;
    %distance = clusterPropScale(1)
    %DirDotDir= clusterPropScale(2)
    %DirDotDisp= clusterPropScale(3)
    %clusterPropScale = [4.6 1.55 0.7]; %[4.6 1.55 0.7];
    propScale = [3.8 1.5 1];
    
    % 260 edges
    %[16 16 10 15];
    % 852 edges
    % times sizes with cuberoot(852/260) = 1.48
    %numDirDotDirClusters = clusterSizes(1); 12
    %numDirDotDispClusters = clusterSizes(2); 12
    %numPreThreshDistClusters = clusterSizes(3);6 
    %numPostThreshDistClusters = clusterSizes(4);10
    clusterSizeScale = Const.MLMoMClusterSizeScale;
    clusterSizes = round(clusterSizeScale* [16 2 10 15] * (numEdges/260).^(1/3)); %best [16 16 10 15] 
    clusterSizes(2) = 2;
    %clusterSizes = [20,15,9,14];
    for k = 1:numel(clusterSizes)
        if (clusterSizes(k) < 2)
            clusterSizes(k) = 2;
        end
    end
    %clusterSizes = [22 22 13 20]; %[16 16 10 15];
    
    % Cluster geometric properties
    
    message_fc(Const, sprintf('  Calculating cluster means '));
    tic;
    [clusterMeans, numInitClusterPoints, numInitUnclusteredPoints] = initClusterMeans(nonSingZmnProp,clusterSizes, threshDist, minInitialClusterPoints, minInitialClusterPointsForDynamicRegion);
    
    [clusterMeans, clusterInd,clusterCounts,maxClusterError, totalClusterError, clusterNumIter, clusterActualTol] =...
        calcClusterMeans(nonSingZmnProp,clusterMeans, propScale, clusterMaxIter, clusterMinTol);
    [numClusters , ~] = size(clusterMeans);
    clusterCalcTime = toc;
    
    clusterMaxEdgeLength = zeros(numClusters, 2); % m n
    
    for k = 1:numClusters
        ind = find(clusterInd(:,1) == k);
        labels = nonSingEdgeLabels(ind,:);
        lengths = zeros(numel(ind), 2);
        lengths(: , 1) = edgeLengths(labels(:, 1));
        lengths(: , 2) = edgeLengths(labels(:, 2));
        
        [~, maxLengthSumInd] = max(lengths(: , 1) + lengths(: , 2));
        maxLengthM = lengths(maxLengthSumInd , 1);
        maxLengthN = lengths(maxLengthSumInd , 2);
        clusterMaxEdgeLength(k, :) = [maxLengthM maxLengthN];
    end

    
    % ------------ REAL AND IMAG MODEL FOR EACH FREQ  ------------ 
    message_fc(Const, sprintf('  Calculating regression models '));
    tic;
    
    for i = 1:numFreq
       % message_fc(Const, sprintf('    Processing frequency %d of %d  ',i,numFreq))
        
        weightModels{i, 1} = calcWeightModel(refSelfZmn(:, i), refTriZmn(:, i), refNonSingZmn(:, i),...
    selfZmnTerms(:,:, i), triZmnTerms(:,:, i), nonSingZmnTerms(:,:, i), nonSingZmnProp,clusterInd,nonSingEdgeLabels,edgeLengths, clusterMaxEdgeLength,singDataThresh, addBias,useProjectedEdges,minPercentImprov, 1);

        weightModels{i, 2} = calcWeightModel(refSelfZmn(:, i), refTriZmn(:, i), refNonSingZmn(:, i),...
    selfZmnTerms(:,:, i), triZmnTerms(:,:, i), nonSingZmnTerms(:,:, i), nonSingZmnProp,clusterInd,nonSingEdgeLabels,edgeLengths, clusterMaxEdgeLength,singDataThresh, addBias,useProjectedEdges,minPercentImprov, 0);

    end
    regressionCalcTime = toc;
    message_fc(Const, sprintf('  Finished ML-MOM training '));
    
    % ------------ UPDATE STRUCT (EXCL. TIMING DATA)  ------------ 
    
    %cells
    mlmom.weightModels =weightModels;
    
    %matrices
    mlmom.singInd = singInd;
    mlmom.edgeLengths = edgeLengths;
    mlmom.nonSingZmnProp = nonSingZmnProp;
    mlmom.nonSingEdgeLabels = nonSingEdgeLabels;
    mlmom.maxClusterError = maxClusterError;
    mlmom.clusterMaxEdgeLength = clusterMaxEdgeLength;
    mlmom.clusterMeans = clusterMeans;
    mlmom.clusterInd = clusterInd;
    mlmom.clusterCounts = clusterCounts;
    mlmom.propScale = propScale;
    mlmom.clusterSizes = clusterSizes;

    %scalars  
    mlmom.useProjectedEdges =useProjectedEdges;
    mlmom.numInitClusterPoints = numInitClusterPoints;
    mlmom.numInitUnclusteredPoints = numInitUnclusteredPoints;
    mlmom.numFreq = numFreq;
    mlmom.totalClusterError = totalClusterError;
    mlmom.clusterNumIter = clusterNumIter;
    mlmom.clusterActualTol = clusterActualTol;
    mlmom.numClusters = numClusters;
    mlmom.maxDist = maxDist;
    mlmom.threshDist = threshDist;
    mlmom.quadPts = Const.QUAD_PTS;
    
    % ------------ CALCULATE Z MATRICES  ------------
    message_fc(Const, sprintf('  Predicting z-Matrices'));
    tic;
    includeRealCalc = 1;
    returnUnity = 1;
    [predZmn, unityZmn] = predictExtractedTerms( mlmom, selfZmnTerms, triZmnTerms,...
        nonSingZmnTerms, nonSingZmnProp, singInd,edgeLengths, includeRealCalc,returnUnity);
    mlmom.predZmn = predZmn;
    mlmom.unityZmn = unityZmn;
    predictCalcTime = toc;
    
    %real imag complex
    predFrobNorms = zeros(numFreq, 3);
    unityFrobNorms = zeros(numFreq, 3);
    refFrobNorms = zeros(numFreq, 3);
    for i = 1:numFreq
        predFrobNorms(numFreq,1) = calcFrobNorm(real(predZmn(:,:,numFreq)));
        predFrobNorms(numFreq,2) = calcFrobNorm(imag(predZmn(:,:,numFreq)));
        predFrobNorms(numFreq,3) = calcFrobNorm(predZmn(:,:,numFreq));
        
        unityFrobNorms(numFreq,1) = calcFrobNorm(real(unityZmn(:,:,numFreq)));
        unityFrobNorms(numFreq,2) = calcFrobNorm(imag(unityZmn(:,:,numFreq)));
        unityFrobNorms(numFreq,3) = calcFrobNorm(unityZmn(:,:,numFreq));
        
        refFrobNorms(numFreq,1) = calcFrobNorm(real(zMatrices.values(:,:,numFreq)));
        refFrobNorms(numFreq,2) = calcFrobNorm(imag(zMatrices.values(:,:,numFreq)));
        refFrobNorms(numFreq,3) = calcFrobNorm(zMatrices.values(:,:,numFreq));
    end
    
    mlmom.predFrobNorms = predFrobNorms;
    mlmom.unityFrobNorms = unityFrobNorms;
    mlmom.refFrobNorms = refFrobNorms;
    
    % ------------ SAVE TIMING DATA ------------
    
    mlmom.unityZmnCalcTime = unityZmnCalcTime;
    mlmom.clusterCalcTime = clusterCalcTime;
    mlmom.regressionCalcTime = regressionCalcTime;
    mlmom.predictCalcTime = predictCalcTime;
    
    % ------------ DISPLAY TIMING DATA  ------------
    message_fc(Const, sprintf('  Input z-Matrices calculation time: %.3f s',unityZmnCalcTime));
    message_fc(Const, sprintf('  Clustering calculation time: %.3f s',clusterCalcTime));
    message_fc(Const, sprintf('  Regression calculation time: %.3f s',regressionCalcTime));
    message_fc(Const, sprintf('  Prediction calculation time: %.3f s',predictCalcTime));
    message_fc(Const, sprintf('  Input z-Matrices + Prediction calculation time: %.3f s',unityZmnCalcTime + predictCalcTime));
    
    
    message_fc(Const,...
        '------------------------------------------------------------------------------------');
    
end

function frobNorm = calcFrobNorm(zMatrix)
    zMatrix = zMatrix(:);
    sum = 0;
    for k = 1:numel(zMatrix)
        sum = sum + abs(zMatrix(k))^2;
    end
    frobNorm = sqrt(sum);
end

% Threshold regression
%                 in = [0.1 0.005; 0.1 0.01; 0.1 0.015 ; 0.2 0.01 ; 0.2 0.02; 0.2 0.03;...
%                     0.4 0.02; 0.4 0.04; 0.4 0.06; 0.8 0.04; 0.8 0.08; 0.8 0.12; 1.6 0.08; 1.6 0.16; 1.6 0.24];
%                 out = [0.015; 0.025; 0.03; 0.04; 0.058; 0.068; 0.06; 0.1;0.131;0.12;0.19;0.23;0.2;0.38;0.5];
%                 w = ( (in' *in) \ in')* out ;
%                 w = [0.043777696318019  1.769794721407624 0.006263888888889];
%                 OR w = [0.049472140762463 1.769794721407624]
%                 Divide w(1) by sqrt(2)
    

    