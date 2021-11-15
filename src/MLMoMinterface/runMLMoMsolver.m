function [mlmom] = runMLMoMsolver(Const, Solver_setup, zMatrices, yVectors, xVectors)
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
    
    % %===================== VARIABLE CONVENTION =====================
    
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
    constMeshSize = Const.MLMoMConstMeshSize;
    numFreq = zMatrices.numFreq;
    edgeLengths = Solver_setup.rwg_basis_functions_length_m;
    numEdges = Solver_setup.num_mom_basis_functions;
    
    addBias = 0;
    singDataThresh = 1e-44;
    %useProjectedEdges = 0;
    sizeConst = 60; %50
    minExtraInitialClusterPoints = 12;
    clusterMinTol = 0.02;
    clusterMaxIter = 6;
    
    %=====================DETERMINE MESH SIZE ===========================
    minEdgeLength = min(edgeLengths);
    maxEdgeLength = max(edgeLengths);
    numLengthClusters = round(1.08*log2(maxEdgeLength/minEdgeLength));
    %numLengthClusters = 2;
    if (numLengthClusters <= 1 || constMeshSize)
        constMeshSize = 1;
    else
        constMeshSize = 0;
    end
    %constMeshSize = 0;
    %numLengthClusters = 2;
    %=====================CALCULATE UNITY TERMS =========================== 
    [termStruct, propStruct,indicesStruct, singInd ] = extractZmnInfo(Const, Solver_setup, constMeshSize);
    termAndPropCalcTime = termStruct.calcTermTime + termStruct.assignAndSwapTime + propStruct.centreDistanceTime + propStruct.edgeCentreTime;
    tic;
    [termStruct.refTwoUniqueZmn, termStruct.refThreeUniqueZmn, termStruct.refFourUniqueZmn] = extractRefZmn(zMatrices.values, singInd );
    [~ , numTerms] = size(termStruct.fourUniqueTerms);
    refZmnExtractTime = toc;
    
    %=======
    %indicesStruct = Const.indicesStruct;
    %clusterStruct = Const.clusterStruct;
    
    %=====================CALCULATE CLUSTER MEANS ===========================
    message_fc(Const, sprintf('  Calculating cluster means '));
    tic;
    
    %threshold
    avgEdgeLength = sum(edgeLengths)/numel(edgeLengths);     
    maxDist = max(propStruct.fourUniqueProp(:,1));
    threshDist = 0.03584296*maxDist + 1.769794721407624*avgEdgeLength;
    %cluster size
    avgClusterSize = clusterSizeScale* numEdges.^2 ./(sizeConst * log(numEdges));
    minClusterSize = numTerms + minExtraInitialClusterPoints + addBias;
    %initialise clusters
    clusterStruct = initClusterStruct(propStruct,threshDist, avgClusterSize, minClusterSize,numLengthClusters, constMeshSize);
    %calculate clusters
    [clusterStruct, indicesStruct] = calcClusterStruct(clusterStruct, propStruct,indicesStruct, clusterMaxIter, clusterMinTol,constMeshSize);
    
    clusterCalcTime = toc;

    %=====================CALCULATE WEIGHTS ===========================
    
    message_fc(Const, sprintf('  Calculating regression models '));
    tic;
    weightModels = loopWeightModels( termStruct,indicesStruct,propStruct, singDataThresh, includeRealCalc,addBias,minPercentImprov, numFreq);
    regressionCalcTime = toc;
    message_fc(Const, sprintf('  Finished ML-MOM training '));

    
    % ===================== PREDICT Z MATRICES  =====================
    message_fc(Const, sprintf('  Predicting z-Matrices'));
    tic;
    unityZmn = zeros(numEdges,numEdges,numFreq );
    for f = 1:numFreq
        for mm = 1:numEdges
            for nn= 1:numEdges
                unityZmn(mm,nn,f) = sum(termStruct.allTerms(mm,nn,:,f) );
            end
        end
    end
    unityZmnTime = toc;
    %[predZmn, indicesStruct] = loopPrediction_(weightModels, clusterStruct, propStruct, termStruct,indicesStruct,unityZmn, includeRealCalc, numFreq);
    [predZmn, indicesStruct, assignTime, multiplyTime] = loopPrediction(weightModels, clusterStruct, propStruct, termStruct,indicesStruct,unityZmn, includeRealCalc, numFreq);
    %predictCalcTime = toc;
    
    % ===================== CALCULATE ERROR, xVECTORS  =====================
    refZmn = zMatrices.values;
    refX = xVectors.Isol;
    predX = zeros(numEdges, numFreq);
    unityX = zeros(numEdges, numFreq);

    %real imag complex
    predRowRelNormPercentError = zeros(numEdges, numFreq);
    unityRowRelNormPercentError = zeros(numEdges, numFreq);
    
    predFrobNorms = zeros(numFreq, 3);
    unityFrobNorms = zeros(numFreq, 3);
    refFrobNorms = zeros(numFreq, 3);
    
    predRelNormPercentError = zeros(numFreq, 3);
    unityRelNormPercentError = zeros(numFreq, 3);
    
    predAngleRelNormPercentError = zeros(numFreq, 1);
    unityAngleRelNormPercentError = zeros(numFreq, 1);
    
    unityXError = zeros(numFreq, 3);
    predXError = zeros(numFreq, 3);
    
    for f = 1:numFreq
        
        unityX(:,f) = unityZmn(:,:,f)\yVectors.values(:,f);
        predX(:,f) = predZmn(:,:,f)\yVectors.values(:,f);
              
        [~,unityXError(f,1)] = calcError(real(refX(:,f)),real(unityX(:,f)));
        [~,predXError(f,1)] = calcError(real(refX(:,f)),real(predX(:,f)));
        [~,unityXError(f,2)] = calcError(imag(refX(:,f)),imag(unityX(:,f)));
        [~,predXError(f,2)] = calcError(imag(refX(:,f)),imag(predX(:,f)));
        [~,unityXError(f,3)] = calcError(refX(:,f),unityX(:,f));
        [~,predXError(f,3)] = calcError(refX(:,f),predX(:,f));
        
        [~, predRelNormPercentError(f,1)] = calcError(real(refZmn(:, :, f)), real(predZmn(:, :, f)));
        [~, predRelNormPercentError(f,2)] = calcError(imag(refZmn(:, :, f)), imag(predZmn(:, :, f)));
        [~, predRelNormPercentError(f,3)] = calcError(refZmn(:, :, f), predZmn(:, :, f));
        
        [~, unityRelNormPercentError(f,1)] = calcError(real(refZmn(:, :, f)), real(unityZmn(:, :, f)));
        [~, unityRelNormPercentError(f,2)] = calcError(imag(refZmn(:, :, f)), imag(unityZmn(:, :, f)));
        [~, unityRelNormPercentError(f,3)] = calcError(refZmn(:, :, f), unityZmn(:, :, f));
        
        %============
        [~, predAngleRelNormPercentError(f,1)] = calcError(angle(refZmn(:, :, f)), angle(predZmn(:, :, f)));       
        [~, unityAngleRelNormPercentError(f,1)] = calcError(angle(refZmn(:, :, f)), angle(unityZmn(:, :, f)));
        
        predFrobNorms(f,1) = calcFrobNorm(real(predZmn(:,:,f)));
        predFrobNorms(f,2) = calcFrobNorm(imag(predZmn(:,:,f)));
        predFrobNorms(f,3) = calcFrobNorm(predZmn(:,:,f));
        
        unityFrobNorms(f,1) = calcFrobNorm(real(unityZmn(:,:,f)));
        unityFrobNorms(f,2) = calcFrobNorm(imag(unityZmn(:,:,f)));
        unityFrobNorms(f,3) = calcFrobNorm(unityZmn(:,:,f));
        
        refFrobNorms(f,1) = calcFrobNorm(real(zMatrices.values(:,:,f)));
        refFrobNorms(f,2) = calcFrobNorm(imag(zMatrices.values(:,:,f)));
        refFrobNorms(f,3) = calcFrobNorm(zMatrices.values(:,:,f));
        
        %===========
        for mm = 1:numEdges
            [~, predRowRelNormPercentError(mm)] = calcError(refZmn(mm, :, f), predZmn(mm, :, f));
            [~, unityRowRelNormPercentError(mm)] = calcError(refZmn(mm, :, f), unityZmn(mm, :, f));
        end
    end
    

    % ===================== UPDATE STRUCT ===================== 
    mlmom = [];
    mlmom.name = "mlmom";
    mlmom.numSols = numFreq;
    mlmom.freqSamples = Solver_setup.frequencies.samples;
    mlmom.Isol = predX;
    
    %cells
    mlmom.weightModels =weightModels;
    
    %structures
    mlmom.clusterStruct = clusterStruct;
    mlmom.indicesStruct = indicesStruct;
    mlmom.propStruct = propStruct;
    mlmom.termStruct = termStruct;
    
    %matrices
    mlmom.singInd = singInd;
    mlmom.refZmn =refZmn;
    mlmom.predZmn = predZmn;
    mlmom.unityZmn = unityZmn;
    mlmom.refX = refX;
    mlmom.unityX = unityX;
    

    %scalars
    mlmom.constMeshSize = constMeshSize;
    mlmom.includeRealCalc = includeRealCalc;
    mlmom.minPercentImprov = minPercentImprov;
    mlmom.clusterSizeScale = clusterSizeScale;
    mlmom.sizeConst = sizeConst;
    mlmom.avgClusterSize = avgClusterSize;
    mlmom.numFreq = numFreq;
    mlmom.maxDist = maxDist;
    mlmom.minEdgeLength = minEdgeLength;
    mlmom.maxEdgeLength = maxEdgeLength;
    mlmom.avgEdgeLength = avgEdgeLength;
    mlmom.threshDist = threshDist;
    mlmom.singDataThresh = singDataThresh;
    mlmom.quadPts = Const.QUAD_PTS;

    mlmom.clusterMinTol = clusterMinTol;
    mlmom.clusterMaxIter = clusterMaxIter;
    
    % timing
    mlmom.calcTermTime = termStruct.calcTermTime;
    mlmom.assignAndSwapTime = termStruct.assignAndSwapTime;
    mlmom.centreDistanceTime = propStruct.centreDistanceTime;
    mlmom.edgeCentreTime = propStruct.edgeCentreTime;
    
    mlmom.termAndPropCalcTime = termAndPropCalcTime;
    mlmom.unityZmnTime = unityZmnTime;
    mlmom.refZmnExtractTime = refZmnExtractTime;
    mlmom.clusterCalcTime = clusterCalcTime;
    mlmom.regressionCalcTime = regressionCalcTime;
    trainingCalcTime = termAndPropCalcTime + refZmnExtractTime + clusterCalcTime + regressionCalcTime;
    mlmom.trainingCalcTime = trainingCalcTime;
    
    mlmom.assignTime = assignTime;
    mlmom.multiplyTime = multiplyTime;
    predictCalcTime = termAndPropCalcTime + unityZmnTime + assignTime + multiplyTime;
    mlmom.predictCalcTime = predictCalcTime;
    
    
    
    %norm and error
    mlmom.unityRowRelNormPercentError = unityRowRelNormPercentError;
    mlmom.predRowRelNormPercentError = predRowRelNormPercentError;
    mlmom.predFrobNorms = predFrobNorms;
    mlmom.unityFrobNorms = unityFrobNorms;
    mlmom.refFrobNorms = refFrobNorms;
    mlmom.unityRelNormPercentError = unityRelNormPercentError;
    mlmom.predRelNormPercentError = predRelNormPercentError;
    mlmom.unityXError = unityXError;
    mlmom.predXError = predXError;
    mlmom.predAngleRelNormPercentError = predAngleRelNormPercentError;
    mlmom.unityAngleRelNormPercentError =unityAngleRelNormPercentError;
    
    % ------------ DISPLAY TIMING DATA  ------------
    message_fc(Const, sprintf('  Term and properties calculation time: %.3f s',termAndPropCalcTime));
    message_fc(Const, sprintf('  Clustering calculation time: %.3f s',clusterCalcTime));
    %message_fc(Const, sprintf('  Regression calculation time: %.3f s',regressionCalcTime));
    message_fc(Const, sprintf('  Training calculation time: %.3f s',trainingCalcTime));
    message_fc(Const, sprintf('  Prediction calculation time: %.3f s',predictCalcTime));
    %message_fc(Const, sprintf('  Input z-Matrices + Prediction calculation time: %.3f s',termAndPropCalcTime + predictCalcTime));
    
    
    % Write the MLMoM solution to a ASCII str file, so that it can be read
    % again by FEKO (for plotting in POSTFEKO) - only if requested (i.e. if the filename is defined)
    if (~isempty(Const.SUNEMmlmomstrfilename))
        writeSolToFile(Const, mlmom);
    end%if
    
    message_fc(Const,...
        '------------------------------------------------------------------------------------');
    
end

function clusterStruct = initClusterStruct(propStruct,threshDist, avgSize, minSize,numLengthClusters, constMeshSize)

    clusterStruct = [];
    %clusterStruct.constMeshSize = propStruct.constMeshSize;
    
    clusterStruct.twoUniqueMeans = [];
    clusterStruct.twoUniqueCounts = [];
    clusterStruct.twoUniqueMaxError = [];
    clusterStruct.twoUniqueSizes = [];
    clusterStruct.twoUniqueErrorCode = propStruct.twoUniqueErrorCode;
    
    clusterStruct.threeUniqueMeans = [];
    clusterStruct.threeUniqueCounts = [];
    clusterStruct.threeUniqueMaxError = [];
    clusterStruct.threeUniqueSizes = [];
    clusterStruct.threeUniqueErrorCode = propStruct.threeUniqueErrorCode;
    
    clusterStruct.fourUniqueMeans = [];
    clusterStruct.fourUniqueCounts = [];
    clusterStruct.fourUniqueMaxError = [];
    clusterStruct.fourUniqueSizes = [];
    clusterStruct.fourUniqueErrorCode = propStruct.fourUniqueErrorCode;
    

    
    if ( constMeshSize)
        numClusters = numel(propStruct.fourUniqueProp(:,1))/avgSize;
        s = ceil(sqrt(numClusters));
        %clusterSizes = ceil([1.2*s s 1]);
        clusterSizes = ceil([1.2*s s 1]);
        clusterIntervals = createIntervals(propStruct.fourUniqueProp, clusterSizes, threshDist, 1);
        [clusterStruct.fourUniqueMeans, numInitClusterPoints, numInitUnclusteredPoints] = initClusterMeans(propStruct.fourUniqueProp,clusterIntervals,  minSize);
        clusterStruct.fourUniqueSizes = clusterSizes;
    else
        %2 unique
        numClusters = numel(propStruct.twoUniqueProp(:,1))/avgSize;
        s = ceil(numClusters);
        clusterSizes = ceil(numLengthClusters);
        clusterIntervals = createIntervals(propStruct.twoUniqueProp, clusterSizes, 0, 0);
        [clusterStruct.twoUniqueMeans, numInitClusterPoints, numInitUnclusteredPoints] = initClusterMeans(propStruct.twoUniqueProp,clusterIntervals,  minSize);
        clusterStruct.twoUniqueSizes = clusterSizes;
        clusterStruct.twoUniqueErrorCode = propStruct.twoUniqueErrorCode;
        
        %3 unique
        numClusters = numel(propStruct.threeUniqueProp(:,1))/avgSize;
        s = ceil(numClusters);
        clusterSizes = ceil([numLengthClusters numLengthClusters]);
        clusterIntervals = createIntervals(propStruct.threeUniqueProp, clusterSizes, 0, 0);
        [clusterStruct.threeUniqueMeans, numInitClusterPoints, numInitUnclusteredPoints] = initClusterMeans(propStruct.threeUniqueProp,clusterIntervals,  minSize);
        clusterStruct.threeUniqueSizes = clusterSizes;
        clusterStruct.threeUniqueErrorCode = propStruct.threeUniqueErrorCode;
        
        %4 unique
        numClusters = numel(propStruct.fourUniqueProp(:,1))/avgSize;
        s = ceil(sqrt(numClusters)/numLengthClusters);
        clusterSizes = ceil([1.2*s s 1 numLengthClusters numLengthClusters]); 
        %clusterSizes = ceil([0.8485*s 0.707*s 1 numLengthClusters numLengthClusters]);  
        clusterIntervals = createIntervals(propStruct.fourUniqueProp, clusterSizes, threshDist, 1);  
        [clusterStruct.fourUniqueMeans, numInitClusterPoints, numInitUnclusteredPoints] = initClusterMeans(propStruct.fourUniqueProp,clusterIntervals,  minSize);
        clusterStruct.fourUniqueSizes = clusterSizes;
        clusterStruct.fourUniqueErrorCode = propStruct.fourUniqueErrorCode;
    end

    %clusterIntervals = createIntervals(nonSingZmnProp, clusterSizes, threshDist, threshPropInd);

    %[clusterMeans, numInitClusterPoints, numInitUnclusteredPoints] = initClusterMeansDynamic(nonSingZmnProp,clusterIntervals,  minClusterSize);
end

function [clusterStruct, indicesStruct] = calcClusterStruct(clusterStruct, propStruct,indicesStruct, clusterMaxIter, clusterMinTol, constMeshSize)
    
    if ( ~constMeshSize)
        %2 unique
        [clusterStruct.twoUniqueMeans, indicesStruct.twoUniqueIndices(:,3),clusterStruct.twoUniqueCounts,clusterStruct.twoUniqueMaxError, totalClusterError, clusterNumIter, clusterActualTol] =...
            calcClusterMeans(propStruct.twoUniqueProp,clusterStruct.twoUniqueMeans, clusterMaxIter, clusterMinTol, clusterStruct.twoUniqueErrorCode);
        %3 unique
        [clusterStruct.threeUniqueMeans, indicesStruct.threeUniqueIndices(:,3),clusterStruct.threeUniqueCounts,clusterStruct.threeUniqueMaxError, totalClusterError, clusterNumIter, clusterActualTol] =...
            calcClusterMeans(propStruct.threeUniqueProp,clusterStruct.threeUniqueMeans, clusterMaxIter, clusterMinTol, clusterStruct.threeUniqueErrorCode);
    end
    %4 unique
    [clusterStruct.fourUniqueMeans, indicesStruct.fourUniqueIndices(:,3) ,clusterStruct.fourUniqueCounts,clusterStruct.fourUniqueMaxError, totalClusterError, clusterNumIter, clusterActualTol] =...
        calcClusterMeans(propStruct.fourUniqueProp,clusterStruct.fourUniqueMeans, clusterMaxIter, clusterMinTol, clusterStruct.fourUniqueErrorCode);


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

function weightModels = loopWeightModels( termStruct,indicesStruct,propStruct, singDataThresh, includeRealCalc,addBias,minPercentImprov, numFreq)
    weightModels = cell(numFreq, 2); % real and imaginary models
    for f = 1:numFreq
        if (includeRealCalc)
            %weightModels{f, 1} = calcWeightModel(extractTermPart(termStruct, 1),indicesStruct, singDataThresh, addBias,minPercentImprov);
            weightModels{f, 1} = calcWeightModel_temp(extractTermPart(termStruct, 1,f),indicesStruct,propStruct, singDataThresh, addBias,minPercentImprov);
        end
        %weightModels{f, 2} = calcWeightModel(extractTermPart(termStruct, 0),indicesStruct, singDataThresh, addBias,minPercentImprov);
        weightModels{f, 2} = calcWeightModel_temp(extractTermPart(termStruct, 0,f),indicesStruct, propStruct, singDataThresh, addBias,minPercentImprov);
    end
end

function [predZmn, indicesStruct, assignTimeTot, multiplyTimeTot] = loopPrediction(weightModels, clusterStruct,propStruct, termStruct,indicesStruct,unityZmn ,includeRealCalc, numFreq)
    predZmn = complex(zeros(termStruct.numEdges, termStruct.numEdges, numFreq));
    hasCalcInd = 0;
    assignTimeTot = 0;
    multiplyTimeTot = 0;
    assignTimeReal = 0;
    multiplyTimeReal = 0;
    for f = 1:numFreq
        if (includeRealCalc)
            %[predZmnReal, indicesStruct] = predictTerms( weightModels{f, 1}, clusterStruct,propStruct, extractTermPart(termStruct,1), indicesStruct, hasCalcInd);
            [predZmnReal, indicesStruct, assignTimeReal, multiplyTimeReal] = predictTerms_temp( weightModels{f, 1}, clusterStruct,propStruct, extractTermPart(termStruct,1, f), indicesStruct, hasCalcInd);
            hasCalcInd = 1;
        else
            predZmnReal = real(unityZmn(:,:,f));
        end
        %[predZmnImag, indicesStruct] = predictTerms( weightModels{f, 2}, clusterStruct,propStruct, extractTermPart(termStruct,0), indicesStruct, hasCalcInd);
        [predZmnImag, indicesStruct, assignTimeImag, multiplyTimeImag] = predictTerms_temp( weightModels{f, 2}, clusterStruct,propStruct, extractTermPart(termStruct,0, f), indicesStruct, hasCalcInd);
        hasCalcInd = 1;
        predZmn(:,:,f) = predZmnReal + 1i* predZmnImag;
        
        assignTimeTot = assignTimeTot + assignTimeImag + assignTimeReal;
        multiplyTimeTot = multiplyTimeTot +  multiplyTimeReal + multiplyTimeImag;
        
    end
end

function termStruct = extractTermPart(termStruct, useReal, f)

    if (useReal)
        termStruct.twoUniqueTerms = real(termStruct.twoUniqueTerms(:,:,f) ) ;
        termStruct.threeUniqueTerms = real(termStruct.threeUniqueTerms(:,:,f) );
        termStruct.fourUniqueTerms = real(termStruct.fourUniqueTerms(:,:,f) );
        termStruct.refTwoUniqueZmn = real(termStruct.refTwoUniqueZmn(:,f) );
        termStruct.refThreeUniqueZmn = real(termStruct.refThreeUniqueZmn(:,f) );
        termStruct.refFourUniqueZmn = real(termStruct.refFourUniqueZmn(:,f) );
    else
        termStruct.twoUniqueTerms = imag(termStruct.twoUniqueTerms(:,:,f) ) ;
        termStruct.threeUniqueTerms = imag(termStruct.threeUniqueTerms(:,:,f) );
        termStruct.fourUniqueTerms = imag(termStruct.fourUniqueTerms(:,:,f) );
        termStruct.refTwoUniqueZmn = imag(termStruct.refTwoUniqueZmn(:,f) );
        termStruct.refThreeUniqueZmn = imag(termStruct.refThreeUniqueZmn(:,f) );
        termStruct.refFourUniqueZmn = imag(termStruct.refFourUniqueZmn(:,f) );
    end
    
end

% function weightModels = loopWeightModels_( termStruct,indicesStruct,propStruct, singDataThresh, includeRealCalc,addBias,minPercentImprov, numFreq)
% weightModels = cell(numFreq, 1); % real and imaginary models
% for f = 1:numFreq
%     weightModels{f, 1} = calcWeightModel_temp(termStruct,indicesStruct, propStruct, singDataThresh, addBias,minPercentImprov);
% end
% end
% 
% function [predZmn, indicesStruct] = loopPrediction_(weightModels, clusterStruct,propStruct, termStruct,indicesStruct,unityZmn ,includeRealCalc, numFreq)
% predZmn = complex(zeros(termStruct.numEdges, termStruct.numEdges, numFreq));
% hasCalcInd = 0;
% for f = 1:numFreq
%     %[predZmnImag, indicesStruct] = predictTerms( weightModels{f, 2}, clusterStruct,propStruct, extractTermPart(termStruct,0), indicesStruct, hasCalcInd);
%     [predZmn(:,:,f), indicesStruct] = predictTerms_temp( weightModels{f, 1}, clusterStruct,propStruct, termStruct, indicesStruct, hasCalcInd);
%     hasCalcInd = 1;
% end
% end
% Threshold regression
%                 in = [0.1 0.005; 0.1 0.01; 0.1 0.015 ; 0.2 0.01 ; 0.2 0.02; 0.2 0.03;...
%                     0.4 0.02; 0.4 0.04; 0.4 0.06; 0.8 0.04; 0.8 0.08; 0.8 0.12; 1.6 0.08; 1.6 0.16; 1.6 0.24];
%                 out = [0.015; 0.025; 0.03; 0.04; 0.058; 0.068; 0.06; 0.1;0.131;0.12;0.19;0.23;0.2;0.38;0.5];
%                 w = ( (in' *in) \ in')* out ;
%                 w = [0.043777696318019  1.769794721407624 0.006263888888889];
%                 OR w = [0.049472140762463 1.769794721407624]
%                 Divide w(1) by sqrt(2)
    

    