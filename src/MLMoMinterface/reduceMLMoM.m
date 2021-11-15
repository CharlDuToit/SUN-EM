function reducedMLMoM = reduceMLMoM(mlmom)
    %only keeps mlmom member needed for prediction
    
    reducedMLMoM = [];

    reducedMLMoM.weightModels = cell(mlmom.numFreq , 2);
    for f = 1:mlmom.numFreq
        if (mlmom.includeRealCalc)
            reducedMLMoM.weightModels{f, 1}.twoUniqueWeights = mlmom.weightModels{f, 1}.twoUniqueWeights;
            reducedMLMoM.weightModels{f, 1}.predTwoUniqueRelNormPercentError = mlmom.weightModels{f, 1}.predTwoUniqueRelNormPercentError;
            reducedMLMoM.weightModels{f, 1}.unityTwoUniqueRelNormPercentError = mlmom.weightModels{f, 1}.unityTwoUniqueRelNormPercentError;
            
            reducedMLMoM.weightModels{f, 1}.threeUniqueWeightsPos = mlmom.weightModels{f, 1}.threeUniqueWeightsPos;
            reducedMLMoM.weightModels{f, 1}.threeUniqueWeightsPos = mlmom.weightModels{f, 1}.threeUniqueWeightsPos;
            reducedMLMoM.weightModels{f, 1}.predThreeUniqueRelNormPercentErrorPos = mlmom.weightModels{f, 1}.predThreeUniqueRelNormPercentErrorPos;
            reducedMLMoM.weightModels{f, 1}.unityThreeUniqueRelNormPercentErrorPos = mlmom.weightModels{f, 1}.unityThreeUniqueRelNormPercentErrorPos;
            
            reducedMLMoM.weightModels{f, 1}.fourUniqueWeightsPos = mlmom.weightModels{f, 1}.fourUniqueWeightsPos;
            reducedMLMoM.weightModels{f, 1}.fourUniqueWeightsNeg = mlmom.weightModels{f, 1}.fourUniqueWeightsNeg;
            reducedMLMoM.weightModels{f, 1}.predFourUniqueRelNormPercentErrorPos = mlmom.weightModels{f, 1}.predFourUniqueRelNormPercentErrorPos;
            reducedMLMoM.weightModels{f, 1}.unityFourUniqueRelNormPercentErrorPos = mlmom.weightModels{f, 1}.unityFourUniqueRelNormPercentErrorPos; 
            reducedMLMoM.weightModels{f, 1}.predFourUniqueRelNormPercentErrorNeg = mlmom.weightModels{f, 1}.predFourUniqueRelNormPercentErrorNeg;
            reducedMLMoM.weightModels{f, 1}.unityFourUniqueRelNormPercentErrorNeg = mlmom.weightModels{f, 1}.unityFourUniqueRelNormPercentErrorNeg; 
        end
            reducedMLMoM.weightModels{f, 2}.twoUniqueWeights = mlmom.weightModels{f, 2}.twoUniqueWeights;
            reducedMLMoM.weightModels{f, 2}.predTwoUniqueRelNormPercentError = mlmom.weightModels{f, 2}.predTwoUniqueRelNormPercentError;
            reducedMLMoM.weightModels{f, 2}.unityTwoUniqueRelNormPercentError = mlmom.weightModels{f, 2}.unityTwoUniqueRelNormPercentError;
            
            reducedMLMoM.weightModels{f, 2}.threeUniqueWeightsPos = mlmom.weightModels{f, 2}.threeUniqueWeightsPos;
            reducedMLMoM.weightModels{f, 2}.predThreeUniqueRelNormPercentErrorPos = mlmom.weightModels{f, 2}.predThreeUniqueRelNormPercentErrorPos;
            reducedMLMoM.weightModels{f, 2}.unityThreeUniqueRelNormPercentErrorPos = mlmom.weightModels{f, 2}.unityThreeUniqueRelNormPercentErrorPos;
            
            reducedMLMoM.weightModels{f, 2}.fourUniqueWeightsPos = mlmom.weightModels{f, 2}.fourUniqueWeightsPos;
            reducedMLMoM.weightModels{f, 2}.fourUniqueWeightsNeg = mlmom.weightModels{f, 2}.fourUniqueWeightsNeg;
            reducedMLMoM.weightModels{f, 2}.predFourUniqueRelNormPercentErrorPos = mlmom.weightModels{f, 2}.predFourUniqueRelNormPercentErrorPos;
            reducedMLMoM.weightModels{f, 2}.unityFourUniqueRelNormPercentErrorPos = mlmom.weightModels{f, 2}.unityFourUniqueRelNormPercentErrorPos; 
            reducedMLMoM.weightModels{f, 2}.predFourUniqueRelNormPercentErrorNeg = mlmom.weightModels{f, 2}.predFourUniqueRelNormPercentErrorNeg;
            reducedMLMoM.weightModels{f, 2}.unityFourUniqueRelNormPercentErrorNeg = mlmom.weightModels{f, 2}.unityFourUniqueRelNormPercentErrorNeg; 
    end
    reducedMLMoM.clusterStruct = [];
    reducedMLMoM.clusterStruct.twoUniqueMeans = mlmom.clusterStruct.twoUniqueMeans;
    reducedMLMoM.clusterStruct.twoUniqueMaxError = mlmom.clusterStruct.twoUniqueMaxError;
    reducedMLMoM.clusterStruct.twoUniqueErrorCode = mlmom.clusterStruct.twoUniqueErrorCode;
    reducedMLMoM.clusterStruct.twoUniqueSizes = mlmom.clusterStruct.twoUniqueSizes;
    reducedMLMoM.clusterStruct.threeUniqueMeans = mlmom.clusterStruct.threeUniqueMeans;
    reducedMLMoM.clusterStruct.threeUniqueMaxError = mlmom.clusterStruct.threeUniqueMaxError;
    reducedMLMoM.clusterStruct.threeUniqueErrorCode = mlmom.clusterStruct.threeUniqueErrorCode;
    reducedMLMoM.clusterStruct.threeUniqueSizes = mlmom.clusterStruct.threeUniqueSizes;
    reducedMLMoM.clusterStruct.fourUniqueMeans = mlmom.clusterStruct.fourUniqueMeans;
    reducedMLMoM.clusterStruct.fourUniqueMaxError = mlmom.clusterStruct.fourUniqueMaxError;
    reducedMLMoM.clusterStruct.fourUniqueErrorCode = mlmom.clusterStruct.fourUniqueErrorCode;
    reducedMLMoM.clusterStruct.fourUniqueSizes = mlmom.clusterStruct.fourUniqueSizes;
    
    %-----prediction paramenters
    reducedMLMoM.constMeshSize =mlmom.constMeshSize;
    reducedMLMoM.includeRealCalc =mlmom.includeRealCalc;
    reducedMLMoM.numFreq = mlmom.numFreq;
    reducedMLMoM.freqSamples = mlmom.freqSamples;
    reducedMLMoM.quadPts = mlmom.quadPts;
    
    % -----training parameters
    %inputs
    [numEdges, ~, ~] = size(mlmom.predZmn);
    reducedMLMoM.trainNumEdges = numEdges;
    %reducedMLMoM.trainMinPercentImprov = mlmom.minPercentImprov;
    reducedMLMoM.trainClusterSizeScale =mlmom.clusterSizeScale;
    reducedMLMoM.trainSizeConst = mlmom.sizeConst;
    reducedMLMoM.trainAvgClusterSize =mlmom.avgClusterSize;
    reducedMLMoM.trainMaxDist =mlmom.maxDist;
    reducedMLMoM.trainMinEdgeLength =mlmom.minEdgeLength;
    reducedMLMoM.trainMaxEdgeLength =mlmom.maxEdgeLength;
    reducedMLMoM.trainAvgEdgeLength =mlmom.avgEdgeLength;
    reducedMLMoM.trainThreshDist =mlmom.threshDist;
    reducedMLMoM.trainSingDataThresh =mlmom.singDataThresh;
    %results
    reducedMLMoM.trainCalcTime =mlmom.trainingCalcTime;
    reducedMLMoM.trainPredictCalcTime =mlmom.predictCalcTime;
    
    reducedMLMoM.trainCalcTermTime = mlmom.calcTermTime;
    reducedMLMoM.trainAssignAndSwapTime = mlmom.assignAndSwapTime;
    reducedMLMoM.trainCentreDistanceTime = mlmom.centreDistanceTime;
    reducedMLMoM.trainEdgeCentreTime = mlmom.edgeCentreTime;
    reducedMLMoM.trainUnityZmnTime = mlmom.unityZmnTime;
    reducedMLMoM.trainRrefZmnExtractTime = mlmom.refZmnExtractTime;
    reducedMLMoM.trainClusterCalcTime = mlmom.clusterCalcTime;
    reducedMLMoM.trainRegressionCalcTime = mlmom.regressionCalcTime;
    
    reducedMLMoM.trainAssignTime =mlmom.assignTime;
    reducedMLMoM.trainMultiplyTime = mlmom.multiplyTime;
    reducedMLMoM.trainPredRelNormPercentError =mlmom.predRelNormPercentError;
    reducedMLMoM.trainUnityRelNormPercentError= mlmom.unityRelNormPercentError;
    reducedMLMoM.trainPredXError= mlmom.predXError;
    reducedMLMoM.trainUnityXError= mlmom.unityXError;
    
end