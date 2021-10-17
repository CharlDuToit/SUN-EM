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
   
    %=====================INITIALIZE CONSTANTS ===========================
    
    clusterSizeScale = Const.MLMoMClusterSizeScale;
    minPercentImprov = Const.MLMoMMinPercentImprov;
    includeRealCalc = Const.MLMoMIncludeRealCalc;
    numFreq = zMatrices.numFreq;
    edgeLengths = Solver_setup.rwg_basis_functions_length_m;
    numEdges = Solver_setup.num_mom_basis_functions;
    

    addBias = 0;
    singDataThresh = 1e-44;
    useProjectedEdges = 0;
    sizeConst = 60; %50
    errorCode = 1;
    minExtraInitialClusterPoints = 0;
    clusterMinTol = 0.02;
    clusterMaxIter = 6;

    %=====================CALCULATE UNITY TERMS ===========================
    tic;
    [allTerms, allProperties, selfZmnTerms, triZmnTerms, nonSingZmnTerms, nonSingZmnProp,nonSingEdgeLabels, singInd] = extractZmnInfo(Const, Solver_setup);
    unityZmnCalcTime = toc;
    [refSelfZmn, refTriZmn, refNonSingZmn] = extractRefZmn(zMatrices.values, singInd );
    [~ , numTerms] = size(nonSingZmnTerms);
    
    %=====================INITIALIZE CLUSTER MEANS ===========================
    message_fc(Const, sprintf('  Calculating cluster means '));
    tic;

    minInitialClusterPoints = numTerms + minExtraInitialClusterPoints + addBias;
    
    avgEdgeLength = sum(edgeLengths)/numel(edgeLengths);
    maxDist = max(nonSingZmnProp(:,1));
    threshDist = 0.03584296*maxDist + 1.769794721407624*avgEdgeLength;
    threshPropInd = 1;
    
    avgSize = clusterSizeScale* numEdges.^2 ./(sizeConst * log(numEdges));
    numClusters = numel(nonSingZmnProp(:,1))./avgSize;
    s = ceil(sqrt(numClusters));
    clusterSizes = ceil([1.2*s s 1]);
    clusterIntervals = createIntervals(nonSingZmnProp, clusterSizes, threshDist, threshPropInd);

    [clusterMeans, numInitClusterPoints, numInitUnclusteredPoints] = initClusterMeansDynamic(nonSingZmnProp,clusterIntervals,  minInitialClusterPoints);
    
    %=====================CALCULATE CLUSTER MEANS ===========================
    [clusterMeans, clusterInd,clusterCounts,maxClusterError, totalClusterError, clusterNumIter, clusterActualTol] =...
        calcClusterMeans(nonSingZmnProp,clusterMeans, clusterMaxIter, clusterMinTol, errorCode);
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
        if (numel(maxLengthSumInd) == 2)
            maxLengthM = lengths(maxLengthSumInd , 1);
            maxLengthN = lengths(maxLengthSumInd , 2);
            clusterMaxEdgeLength(k, :) = [maxLengthM maxLengthN];
        end
    end

    %=====================CALCULATE WEIGHTS ===========================
    
    message_fc(Const, sprintf('  Calculating regression models '));
    tic;
    weightModels = cell(numFreq, 2); % real and imaginary models
    for f = 1:numFreq
        if (includeRealCalc)
            weightModels{f, 1} = calcWeightModel(real(refSelfZmn(:, f)), real(refTriZmn(:, f)), real(refNonSingZmn(:, f)),...
                real(selfZmnTerms(:,:, f)), real(triZmnTerms(:,:, f)), real(nonSingZmnTerms(:,:, f)), nonSingZmnProp,clusterInd,nonSingEdgeLabels,edgeLengths, clusterMaxEdgeLength,singDataThresh, addBias,useProjectedEdges,minPercentImprov);
        end
        weightModels{f, 2} = calcWeightModel(imag(refSelfZmn(:, f)), imag(refTriZmn(:, f)), imag(refNonSingZmn(:, f)),...
            imag(selfZmnTerms(:,:, f)), imag(triZmnTerms(:,:, f)), imag(nonSingZmnTerms(:,:, f)), nonSingZmnProp,clusterInd,nonSingEdgeLabels,edgeLengths, clusterMaxEdgeLength,singDataThresh, addBias,useProjectedEdges,minPercentImprov);

    end
    regressionCalcTime = toc;
    message_fc(Const, sprintf('  Finished ML-MOM training '));
    
    % ===================== UPDATE STRUCT (EXCL. TIMING DATA) ===================== 
    mlmom = [];
    %cells
    mlmom.weightModels =weightModels;
    
    %matrices
    mlmom.allTerms = allTerms;
    mlmom.allProperties = allProperties;
    mlmom.singInd = singInd;
    mlmom.edgeLengths = edgeLengths;
    mlmom.nonSingZmnProp = nonSingZmnProp;
    mlmom.nonSingEdgeLabels = nonSingEdgeLabels;
    mlmom.maxClusterError = maxClusterError;
    mlmom.clusterMaxEdgeLength = clusterMaxEdgeLength;
    mlmom.clusterMeans = clusterMeans;
    mlmom.clusterInd = clusterInd;
    mlmom.clusterCounts = clusterCounts;
    mlmom.clusterSizes = clusterSizes;

    %scalars
    mlmom.includeRealCalc = includeRealCalc;
    mlmom.clusterSizeScale = clusterSizeScale;
    mlmom.minPercentImprov = minPercentImprov;
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
    
    % ===================== PREDICT Z MATRICES  =====================
    message_fc(Const, sprintf('  Predicting z-Matrices'));
    tic;
    returnUnity = 1;
    [predZmn, unityZmn] = predictExtractedTerms( mlmom, selfZmnTerms, triZmnTerms,...
        nonSingZmnTerms, nonSingZmnProp, singInd,edgeLengths,returnUnity, errorCode);
    mlmom.predZmn = predZmn;
    mlmom.unityZmn = unityZmn;
    refZmn = zMatrices.values;
    mlmom.refZmn =refZmn;
    predictCalcTime = toc;
    
    %real imag complex
    predFrobNorms = zeros(numFreq, 3);
    unityFrobNorms = zeros(numFreq, 3);
    refFrobNorms = zeros(numFreq, 3);
    predRelNormPercentError = zeros(numFreq, 3);
    
    for f = 1:numFreq
        [~, predRelNormPercentError(f,1)] = calcError(real(refZmn(:, :, f)), real(predZmn(:, :, f)));
        [~, predRelNormPercentError(f,2)] = calcError(imag(refZmn(:, :, f)), imag(predZmn(:, :, f)));
        [~, predRelNormPercentError(f,3)] = calcError(refZmn(:, :, f), predZmn(:, :, f));
        
        predFrobNorms(f,1) = calcFrobNorm(real(predZmn(:,:,f)));
        predFrobNorms(f,2) = calcFrobNorm(imag(predZmn(:,:,f)));
        predFrobNorms(f,3) = calcFrobNorm(predZmn(:,:,f));
        
        unityFrobNorms(f,1) = calcFrobNorm(real(unityZmn(:,:,f)));
        unityFrobNorms(f,2) = calcFrobNorm(imag(unityZmn(:,:,f)));
        unityFrobNorms(f,3) = calcFrobNorm(unityZmn(:,:,f));
        
        refFrobNorms(f,1) = calcFrobNorm(real(zMatrices.values(:,:,f)));
        refFrobNorms(f,2) = calcFrobNorm(imag(zMatrices.values(:,:,f)));
        refFrobNorms(f,3) = calcFrobNorm(zMatrices.values(:,:,f));
    end
    
    mlmom.predFrobNorms = predFrobNorms;
    mlmom.unityFrobNorms = unityFrobNorms;
    mlmom.refFrobNorms = refFrobNorms;
    mlmom.predRelNormPercentError = predRelNormPercentError;
    
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

function clusterIntervals = createIntervals(prop, sizes, threshDist, threshPropInd)
    clusterIntervals = cell(numel(sizes),1);
    for k = 1:numel(sizes)
        ind = find(threshPropInd == k);
        
        if (numel(ind) > 0)
            td = threshDist(ind);
            minDist = min(prop(:,k));
            maxDist = max(prop(:,k));
            s = max(2, round(sizes(k)/2));
            distInterval = linspace(minDist, td, s );
            distInterval(s:(2*s -1)) = linspace(td, maxDist, s);
            clusterIntervals{k} = distInterval;
        else
           clusterIntervals{k} = linspace( min(prop(:,k)) , max(prop(:,k)) , sizes(k) + 1); 
        end
        
    end
end

% Threshold regression
%                 in = [0.1 0.005; 0.1 0.01; 0.1 0.015 ; 0.2 0.01 ; 0.2 0.02; 0.2 0.03;...
%                     0.4 0.02; 0.4 0.04; 0.4 0.06; 0.8 0.04; 0.8 0.08; 0.8 0.12; 1.6 0.08; 1.6 0.16; 1.6 0.24];
%                 out = [0.015; 0.025; 0.03; 0.04; 0.058; 0.068; 0.06; 0.1;0.131;0.12;0.19;0.23;0.2;0.38;0.5];
%                 w = ( (in' *in) \ in')* out ;
%                 w = [0.043777696318019  1.769794721407624 0.006263888888889];
%                 OR w = [0.049472140762463 1.769794721407624]
%                 Divide w(1) by sqrt(2)
    

    