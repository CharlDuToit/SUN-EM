function [weightModel] = calcWeightModel(termStruct,indicesStruct, singDataThresh, addBias,minPercentImprov)
%function [weightModel] = calcWeightModel(refSelfZmn, refTriZmn, refNonSingZmn,...
    %selfZmnTerms, triZmnTerms, nonSingZmnTerms,clusterInd,...
    %singDataThresh,  addBias,minPercentImprov)
    
    emptyWeightsForUnityData = 0;
    
    %twoUniqueWeights = [];
    %threeUniqueWeights = [];
    %fourUniqueWeights = [];
    %fourUniqueWeightsIndices = [];
    twoUniqueWeightsIndices = [];
    threeUniqueWeightsIndices = [];
    unityTwoUniqueClusterRelNormPercentError = [];
    predTwoUniqueClusterRelNormPercentError = [];
    unityThreeUniqueClusterRelNormPercentError = [];
    predThreeUniqueClusterRelNormPercentError = [];
    
    % ------------ TWO UNIQUE MLR  ------------
    [~, numCol] = size(indicesStruct.twoUniqueIndices);
    if ( numCol == 3)
        minImprov = 0;
        [twoUniqueWeights, predTwoUniqueZmn,unityTwoUniqueZmn, twoUniqueWeightsIndices, unityTwoUniqueMSE, predTwoUniqueMSE,...
            unityTwoUniqueRelNormPercentError, predTwoUniqueRelNormPercentError, unityTwoUniqueClusterRelNormPercentError,predTwoUniqueClusterRelNormPercentError] =...
            calcClusterWeights(termStruct.twoUniqueTerms, termStruct.refTwoUniqueZmn, indicesStruct.twoUniqueIndices(:,3), singDataThresh, addBias,minImprov);
    else
        [twoUniqueWeights, predTwoUniqueZmn, unityTwoUniqueZmn, unityTwoUniqueMSE, predTwoUniqueMSE, unityTwoUniqueRelNormPercentError, predTwoUniqueRelNormPercentError] =...
        MLR(termStruct.twoUniqueTerms, termStruct.refTwoUniqueZmn, singDataThresh, addBias, emptyWeightsForUnityData);           
%        [twoUniqueWeights, predTwoUniqueZmn, unityTwoUniqueZmn, unityTwoUniqueMSE, predTwoUniqueMSE, unityTwoUniqueRelNormPercentError, predTwoUniqueRelNormPercentError] =...
%            MLR(termStruct.twoUniqueTerms, termStruct.refTwoUniqueZmn, 1e50, addBias, emptyWeightsForUnityData);
        %twoUniqueWeights = ones(1,8);
    end
        
    % ------------ THREE UNIQUE MLR  ------------
    [~, numCol] = size(indicesStruct.threeUniqueIndices);
    if ( numCol == 3)
        minImprov = 0;
        [threeUniqueWeights, predThreeUniqueZmn,unityThreeUniqueZmn, threeUniqueWeightsIndices, unityThreeUniqueMSE, predThreeUniqueMSE,...
            unityThreeUniqueRelNormPercentError, predThreeUniqueRelNormPercentError, unityThreeUniqueClusterRelNormPercentError,predThreeUniqueClusterRelNormPercentError] =...
            calcClusterWeights(termStruct.threeUniqueTerms, termStruct.refThreeUniqueZmn, indicesStruct.threeUniqueIndices(:,3), singDataThresh, addBias,minImprov);
    else
         [threeUniqueWeights, predThreeUniqueZmn, unityThreeUniqueZmn, unityThreeUniqueMSE, predThreeUniqueMSE, unityThreeUniqueRelNormPercentError, predThreeUniqueRelNormPercentError] =...
         MLR(termStruct.threeUniqueTerms, termStruct.refThreeUniqueZmn, singDataThresh, addBias, emptyWeightsForUnityData);
%             [threeUniqueWeights, predThreeUniqueZmn, unityThreeUniqueZmn, unityThreeUniqueMSE, predThreeUniqueMSE, unityThreeUniqueRelNormPercentError, predThreeUniqueRelNormPercentError] =...
%         MLR(termStruct.threeUniqueTerms, termStruct.refThreeUniqueZmn, 1e50, addBias, emptyWeightsForUnityData);
        %threeUniqueWeights = ones(1,8);
    
    end
    % ------------ FOUR UNIQUE MLR  -----------

    [fourUniqueWeights, predFourUniqueZmn,unityFourUniqueZmn, fourUniqueWeightsIndices, unityFourUniqueMSE, predFourUniqueMSE,...
    unityFourUniqueRelNormPercentError, predFourUniqueRelNormPercentError, unityFourUniqueClusterRelNormPercentError, predFourUniqueClusterRelNormPercentError] =...
        calcClusterWeights(termStruct.fourUniqueTerms, termStruct.refFourUniqueZmn, indicesStruct.fourUniqueIndices(:,3), singDataThresh, addBias,minPercentImprov);    
    
    
    % ------------ UPDATE WEIGHTMODEL ------------
    weightModel = [];
    
    % ------------
    % -----Matrices
    
    % 2 unique
    %weightModel.selfZmnTerms = selfZmnTerms;  
    %weightModel.refSelfZmn = refSelfZmn;
    %weightModel.unitySelfZmn = unityTwoUniqueZmn;
    %weightModel.predSelfZmn = predTwoUniqueZmn; 
    weightModel.twoUniqueWeights = twoUniqueWeights;
    weightModel.twoUniqueWeightsIndices = twoUniqueWeightsIndices;
    weightModel.unityTwoUniqueClusterRelNormPercentError = unityTwoUniqueClusterRelNormPercentError;
    weightModel.predTwoUniqueClusterRelNormPercentError = predTwoUniqueClusterRelNormPercentError;
    
    % 3 unique triangles
    %weightModel.triZmnTerms = triZmnTerms;  
    %weightModel.refTriZmn = refTriZmn;
    %weightModel.unityTriZmn = unityThreeUniqueZmn;
    %weightModel.predTriZmn = predThreeUniqueZmn;
    weightModel.threeUniqueWeights = threeUniqueWeights;
    weightModel.threeUniqueWeightsIndices = threeUniqueWeightsIndices;
    weightModel.unityThreeUniqueClusterRelNormPercentError = unityThreeUniqueClusterRelNormPercentError;
    weightModel.predThreeUniqueClusterRelNormPercentError = predThreeUniqueClusterRelNormPercentError;
    
    % 4 unique
    %weightModel.nonSingZmnTerms = nonSingZmnTerms;  
    %weightModel.refNonSingZmn = refNonSingZmn;
    %weightModel.unityFourUniqueZmn = unityFourUniqueZmn;
    %weightModel.predFourUniqueZmn = predFourUniqueZmn; 
    weightModel.fourUniqueWeights = fourUniqueWeights;
    weightModel.fourUniqueWeightsIndices = fourUniqueWeightsIndices;
    weightModel.unityFourUniqueClusterRelNormPercentError = unityFourUniqueClusterRelNormPercentError;
    weightModel.predFourUniqueClusterRelNormPercentError = predFourUniqueClusterRelNormPercentError;
    
    % ------------
    % -----Scalars
    
   % weightModel.numNonSing = numNonSing;

    % 2 Unique
    weightModel.predTwoUniqueRelNormPercentError = predTwoUniqueRelNormPercentError;
    weightModel.unityTwoUniqueRelNormPercentError = unityTwoUniqueRelNormPercentError;
    weightModel.predTwoUniqueMSE = predTwoUniqueMSE;    
    weightModel.unityTwoUniqueMSE = unityTwoUniqueMSE;
    
    % 3 unique triangles
    weightModel.predThreeUniqueRelNormPercentError = predThreeUniqueRelNormPercentError;
    weightModel.unityThreeUniqueRelNormPercentError = unityThreeUniqueRelNormPercentError;
    weightModel.predThreeUniqueMSE= predThreeUniqueMSE;     
    weightModel.unityThreeUniqueMSE = unityThreeUniqueMSE;
    
    % 4 Unique
    weightModel.predFourUniqueRelNormPercentError = predFourUniqueRelNormPercentError;
    weightModel.unityFourUniqueRelNormPercentError = unityFourUniqueRelNormPercentError;
    weightModel.predFourUniqueMSE = predFourUniqueMSE;  
    weightModel.unityFourUniqueMSE = unityFourUniqueMSE;
    
end

function  [clusterWeights, pred,unity, clusterWeightsInd, unityMSE, predMSE,...
    unityRelNormPercentError, predRelNormPercentError, unityClusterRelNormPercentError, predClusterRelNormPercentError] =...
    calcClusterWeights(terms, refVals, clusterInd, singThresh, addBias, minPercentImprov)
    %edgeLabels only used if useProjectedEdges =1
    
    %[predNonSingMSE, predNonSingRelNormPercentError] = calcError(refVals, pred);
    %[unityNonSingMSE, unityNonSingRelNormPercentError] = calcError(refVals, unityNonSingZmn);
    
    [numObs, numTerms] = size(terms);
    %[numObs, 1] = size(ClusterInd);
    %[numObs, 1] = size(refVals)
    [numClusters ,~] = max(clusterInd);
    clusterWeights = ones(numClusters,numTerms );
    clusterWeightsInd = zeros(numClusters,1 );
    pred = zeros(numObs,1);
    returnEmptyWeightsIfUnity = 0;%1
    
    predClusterRelNormPercentError = zeros(numClusters,1);
    unityClusterRelNormPercentError = zeros(numClusters,1);
   % addBias = 0;
   % 
   count= 0;
    if (addBias)
        clusterWeights(:, numTerms +1 ) = zeros(numClusters, 1);
    end
    
    nonUnityClusterWeightsCount = 0;
    for k = 1:numClusters
        
        ind = find(clusterInd(:,1) == k);
        t = terms(ind,:);       
        unity  = sum(t, 2);
        
        ref = refVals(ind);
        [weights, pred(ind,1),~ ]=MLR(t,ref, singThresh, addBias, returnEmptyWeightsIfUnity);

        if (eq(size(weights), [0 0]))
            pred(ind,1) = unity;
        else

          urnpe = 100* sqrt(((unity - ref)' * (unity - ref))/ sum(ref.^2) ) ;
          prnpe = 100* sqrt(((pred(ind,1) - ref)' * (pred(ind,1) - ref))/ sum(ref.^2) ) ;
          
          %did not correct error enough
          if (urnpe - prnpe < minPercentImprov)
              pred(ind,1) = unity;
          else
              nonUnityClusterWeightsCount = nonUnityClusterWeightsCount + 1;
              count = count + numel(ref);
              clusterWeights(nonUnityClusterWeightsCount,:) = weights;
              clusterWeightsInd(nonUnityClusterWeightsCount) = k;
              predClusterRelNormPercentError(nonUnityClusterWeightsCount) = prnpe;
              unityClusterRelNormPercentError(nonUnityClusterWeightsCount) = urnpe;
          end
            
        end %if (eq(size(weights), [0 0]))
        
    end% for k
    clusterWeights = clusterWeights(1:nonUnityClusterWeightsCount, :);
    clusterWeightsInd = clusterWeightsInd(1:nonUnityClusterWeightsCount,:);
    predClusterRelNormPercentError = predClusterRelNormPercentError(1:nonUnityClusterWeightsCount, :);
    unityClusterRelNormPercentError = unityClusterRelNormPercentError(1:nonUnityClusterWeightsCount, :);
    
    %error
    unity  = sum(terms, 2);
    [predMSE, predRelNormPercentError] = calcError(refVals, pred);
    [unityMSE, unityRelNormPercentError] = calcError(refVals, unity);

end
