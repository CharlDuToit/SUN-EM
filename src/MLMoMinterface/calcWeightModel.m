function [weightModel] = calcWeightModel(refSelfZmn, refTriZmn, refNonSingZmn,...
    selfZmnTerms, triZmnTerms, nonSingZmnTerms, nonSingZmnProp,clusterInd,nonSingEdgeLabels,...
    edgeLengths,clusterMaxEdgeLength,singDataThresh,  addBias,useProjectedEdges,minPercentImprov)
    
    emptyWeightsForUnityData = 0;
    
    % ------------ SELF TERMS MLR  ------------

    
    [selfWeights, predSelfZmn, unitySelfZmn, unitySelfMSE, predSelfMSE, unitySelfRelNormPercentError, predSelfRelNormPercentError] =...
        MLR(selfZmnTerms, refSelfZmn, singDataThresh, addBias, emptyWeightsForUnityData);
        
    % ------------ 3 UNIQUE TRIANGLES TERMS MLR ------------

    [triWeights, predTriZmn, unityTriZmn, unityTriMSE, predTriMSE, unityTriRelNormPercentError, predTriRelNormPercentError] =...
        MLR(triZmnTerms, refTriZmn, singDataThresh, addBias, emptyWeightsForUnityData);

    % ------------ NON SINGULAR TERMS CLUSTERING MLR  ------------

    [nonSingWeights, predNonSingZmn,unityNonSingZmn, nonSingWeightsInd, unityNonSingMSE, predNonSingMSE,...
    unityNonSingRelNormPercentError, predNonSingRelNormPercentError] =...
        calcClusterWeights(nonSingZmnTerms, refNonSingZmn,nonSingZmnProp, clusterInd,nonSingEdgeLabels,...
        edgeLengths,clusterMaxEdgeLength, singDataThresh, addBias,minPercentImprov, useProjectedEdges);    
    
    
    % ------------ UPDATE WEIGHTMODEL ------------
    weightModel = [];
    
    % ------------
    % -----Matrices
    
    % Non singular
    weightModel.nonSingZmnTerms = nonSingZmnTerms;  
    weightModel.refNonSingZmn = refNonSingZmn;
    weightModel.unityNonSingZmn = unityNonSingZmn;
    weightModel.predNonSingZmn = predNonSingZmn; 
    weightModel.nonSingWeights = nonSingWeights;
    weightModel.nonSingWeightsInd = nonSingWeightsInd;
    
    % Self
    weightModel.selfZmnTerms = selfZmnTerms;  
    weightModel.refSelfZmn = refSelfZmn;
    weightModel.unitySelfZmn = unitySelfZmn;
    weightModel.predSelfZmn = predSelfZmn; 
    weightModel.selfWeights = selfWeights;
    
    % 3 unique triangles
    weightModel.triZmnTerms = triZmnTerms;  
    weightModel.refTriZmn = refTriZmn;
    weightModel.unityTriZmn = unityTriZmn;
    weightModel.predTriZmn = predTriZmn;
    weightModel.triWeights = triWeights;
    
    % ------------
    % -----Scalars
    
   % weightModel.numNonSing = numNonSing;
    
    % Non singular
    weightModel.predNonSingRelNormPercentError = predNonSingRelNormPercentError;
    weightModel.unityNonSingRelNormPercentError = unityNonSingRelNormPercentError;
    weightModel.predNonSingMSE = predNonSingMSE;  
    weightModel.unityNonSingMSE = unityNonSingMSE;

    % Self
    weightModel.predSelfRelNormPercentError = predSelfRelNormPercentError;
    weightModel.unitySelfRelNormPercentError = unitySelfRelNormPercentError;
    weightModel.predSelfMSE = predSelfMSE;    
    weightModel.unitySelfMSE = unitySelfMSE;
    
    % 3 unique triangles
    weightModel.predTriRelNormPercentError = predTriRelNormPercentError;
    weightModel.unityTriRelNormPercentError = unityTriRelNormPercentError;
    weightModel.predTriMSE= predTriMSE;     
    weightModel.unityTriMSE = unityTriMSE;
    
end
