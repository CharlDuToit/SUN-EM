function [weightModel] = calcWeightModel(refSelfZmn, refTriZmn, refNonSingZmn,...
    selfZmnTerms, triZmnTerms, nonSingZmnTerms, nonSingZmnProp,clusterInd,nonSingEdgeLabels,edgeLengths,clusterMaxEdgeLength,singDataThresh,  addBias,useProjectedEdges,minPercentImprov, useReal)

    % Variable convention:
    % dataset + singularism + ?
    %
    % dataset : ref, pred, unity
    % singularism : nonSing, sing ( 2 or 3 triangles), self, tri
    % ? : Context should be clear, first 2 parts could be missing as well
    
    weightModel = [];
    emptyWeightsForUnityData = 0;
    
 % ------------ SELF TERMS MLR  ------------
    
    [numSelf , ~] = size(selfZmnTerms);
    if (useReal)
        selfZmnTerms = real(selfZmnTerms);
        refSelfZmn = real(refSelfZmn);
    else
        selfZmnTerms = imag(selfZmnTerms);
        refSelfZmn = imag(refSelfZmn);
    end

    unitySelfZmn = sum(selfZmnTerms, 2);
    [selfWeights, predSelfZmn] =MLR(selfZmnTerms, refSelfZmn, singDataThresh, addBias, emptyWeightsForUnityData);
    
    % ------------ 3 UNIQUE TRIANGLES TERMS MLR ------------
    
    [numTri , ~] = size(triZmnTerms);
    if (useReal)
        triZmnTerms = real(triZmnTerms);
        refTriZmn = real(refTriZmn);
    else
        triZmnTerms = imag(triZmnTerms);
        refTriZmn = imag(refTriZmn);
    end

    unityTriZmn = sum(triZmnTerms, 2);
    [triWeights, predTriZmn] =MLR(triZmnTerms, refTriZmn, singDataThresh, addBias, emptyWeightsForUnityData);

    % ------------ NON SINGULAR TERMS CLUSTERING MLR  ------------
    
    [numNonSing, ~] = size(nonSingZmnProp);
    if (useReal)
        nonSingZmnTerms = real(nonSingZmnTerms);
        refNonSingZmn = real(refNonSingZmn);
    else
        nonSingZmnTerms = imag(nonSingZmnTerms);
        refNonSingZmn = imag(refNonSingZmn);
    end

    unityNonSingZmn  = sum(nonSingZmnTerms, 2); % Same estimation as internal solver
    %minPercentImprov = 4;
    [nonSingWeights, predNonSingZmn,nonSingWeightsInd] =...
        calcClusterWeights(nonSingZmnTerms, refNonSingZmn,nonSingZmnProp, clusterInd,nonSingEdgeLabels,edgeLengths,clusterMaxEdgeLength, singDataThresh, addBias,minPercentImprov, useProjectedEdges);    
    
    % ------------ ERROR ------------
    
    % self terms
   
    predSelfDiff = predSelfZmn - refSelfZmn;
    predSelfRelNormPercentError = 100 * sqrt(sum(predSelfDiff.^2) / sum(refSelfZmn.^2));
    predSelfMSE = predSelfDiff'*predSelfDiff /numSelf;
    
    unitySelfDiff =  unitySelfZmn - refSelfZmn;
    unitySelfRelNormPercentError = 100 * sqrt(sum(unitySelfDiff.^2) / sum(refSelfZmn.^2));
    unitySelfMSE = unitySelfDiff'*unitySelfDiff ./numSelf;
    
    % 3 unique triangles
   
    predTriDiff = predTriZmn - refTriZmn;
    predTriRelNormPercentError = 100 * sqrt(sum(predTriDiff.^2) / sum(refTriZmn.^2));
    predTriMSE = predTriDiff'*predTriDiff /numTri;
    
    unityTriDiff =  unityTriZmn - refTriZmn;
    unityTriRelNormPercentError = 100 * sqrt(sum(unityTriDiff.^2) / sum(refTriZmn.^2));
    unityTriMSE = unityTriDiff'*unityTriDiff ./numTri;
    
    % Non singular
    
    predNonSingDiff = predNonSingZmn -refNonSingZmn;
    predNonSingRelNormPercentError = 100 * sqrt(sum(predNonSingDiff.^2) / sum(refNonSingZmn.^2));
    predNonSingMSE = predNonSingDiff'*predNonSingDiff /numNonSing;
    
    unityNonSingDiff =  unityNonSingZmn - refNonSingZmn;
    unityNonSingRelNormPercentError = 100 * sqrt(sum(unityNonSingDiff.^2) / sum(refNonSingZmn.^2));
    unityNonSingMSE = unityNonSingDiff'*unityNonSingDiff ./numNonSing;
    
    % ------------ UPDATE WEIGHTMODEL ------------
    % ------------
    % Matrices
    
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
    % Scalars
    
    %mlmom.totsetupTime = 0.0;
    %mlmom.totsolTime = 0.0;
    weightModel.numNonSing = numNonSing;
    
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
