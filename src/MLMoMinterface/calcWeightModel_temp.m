function [weightModel] = calcWeightModel_temp(termStruct,indicesStruct,propStruct, singDataThresh, addBias,minPercentImprov)
%function [weightModel] = calcWeightModel(refSelfZmn, refTriZmn, refNonSingZmn,...
    %selfZmnTerms, triZmnTerms, nonSingZmnTerms,clusterInd,...
    %singDataThresh,  addBias,minPercentImprov)
    
    emptyWeightsForUnityData = 0;
    
    %twoUniqueWeights = [];
    %threeUniqueWeights = [];
    %fourUniqueWeights = [];
    %fourUniqueWeightsIndices = [];
    %twoUniqueWeightsIndices = [];
    %threeUniqueWeightsIndices = [];
    
    unityTwoUniqueClusterRelNormPercentError = [];
    predTwoUniqueClusterRelNormPercentError = [];
    
    unityThreeUniqueClusterRelNormPercentErrorPos = [];
    predThreeUniqueClusterRelNormPercentErrorPos = [];
    unityThreeUniqueClusterRelNormPercentErrorNeg = [];
    predThreeUniqueClusterRelNormPercentErrorNeg = [];
    
    % ------------ TWO UNIQUE MLR  ------------
    [~, numCol] = size(indicesStruct.twoUniqueIndices);
    if ( numCol == 3)
        minImprov = 0;
        [twoUniqueWeights, predTwoUniqueZmn,unityTwoUniqueZmn, unityTwoUniqueMSE, predTwoUniqueMSE,...
            unityTwoUniqueRelNormPercentError, predTwoUniqueRelNormPercentError, unityTwoUniqueClusterRelNormPercentError,predTwoUniqueClusterRelNormPercentError] =...
            calcClusterWeights(termStruct.twoUniqueTerms, termStruct.refTwoUniqueZmn, indicesStruct.twoUniqueIndices(:,3), singDataThresh, addBias);
    else
        [twoUniqueWeights, predTwoUniqueZmn, unityTwoUniqueZmn, unityTwoUniqueMSE, predTwoUniqueMSE, unityTwoUniqueRelNormPercentError, predTwoUniqueRelNormPercentError] =...
        MLR(termStruct.twoUniqueTerms, termStruct.refTwoUniqueZmn, singDataThresh, addBias, emptyWeightsForUnityData);           
%        [twoUniqueWeights, predTwoUniqueZmn, unityTwoUniqueZmn, unityTwoUniqueMSE, predTwoUniqueMSE, unityTwoUniqueRelNormPercentError, predTwoUniqueRelNormPercentError] =...
%            MLR(termStruct.twoUniqueTerms, termStruct.refTwoUniqueZmn, 1e50, addBias, emptyWeightsForUnityData);
        %twoUniqueWeights = ones(1,8);
    end
        
    % ------------ THREE UNIQUE MLR  ------------
    %positive group
%     ind = find(propStruct.threeUniqueSwapGroups == 0);  
%     threeUniqueTermsPos = termStruct.threeUniqueTerms(ind,:);
%     refThreeUniqueZmnPos = termStruct.refThreeUniqueZmn(ind,:);
%     
%     %negative group
%     ind = find(propStruct.threeUniqueSwapGroups == 1);
%     threeUniqueTermsNeg = termStruct.threeUniqueTerms(ind,:);
%     refThreeUniqueZmnNeg = termStruct.refThreeUniqueZmn(ind,:);
%     
%     [~, numCol] = size(indicesStruct.threeUniqueIndices);
%     if ( numCol == 3)
%         threeUniqueIndicesPos = indicesStruct.threeUniqueIndices(:,3);
%         ind = find(propStruct.threeUniqueSwapGroups == 0); 
%         threeUniqueIndicesPos = threeUniqueIndicesPos(ind,:);
%         [threeUniqueWeightsPos, ~,~, unityThreeUniqueMSEPos, predThreeUniqueMSEPos,...
%             unityThreeUniqueRelNormPercentErrorPos, predThreeUniqueRelNormPercentErrorPos, unityThreeUniqueClusterRelNormPercentErrorPos,predThreeUniqueClusterRelNormPercentErrorPos] =...
%             calcClusterWeights(threeUniqueTermsPos, refThreeUniqueZmnPos, threeUniqueIndicesPos, singDataThresh, addBias);
%         
%         threeUniqueIndicesNeg = indicesStruct.threeUniqueIndices(:,3);
%         ind = find(propStruct.threeUniqueSwapGroups == 1);
%         threeUniqueIndicesNeg = threeUniqueIndicesNeg(ind,:);
%         [threeUniqueWeightsNeg, ~,~, unityThreeUniqueMSENeg, predThreeUniqueMSENeg,...
%             unityThreeUniqueRelNormPercentErrorNeg, predThreeUniqueRelNormPercentErrorNeg, unityThreeUniqueClusterRelNormPercentErrorNeg,predThreeUniqueClusterRelNormPercentErrorNeg] =...
%             calcClusterWeights(threeUniqueTermsNeg, refThreeUniqueZmnNeg, threeUniqueIndicesNeg, singDataThresh, addBias);
%     else
%          [threeUniqueWeightsPos, ~, ~, unityThreeUniqueMSEPos, predThreeUniqueMSEPos, unityThreeUniqueRelNormPercentErrorPos, predThreeUniqueRelNormPercentErrorPos] =...
%          MLR(threeUniqueTermsPos, refThreeUniqueZmnPos, singDataThresh, addBias, emptyWeightsForUnityData);
%      
%          [threeUniqueWeightsNeg, ~, ~, unityThreeUniqueMSENeg, predThreeUniqueMSENeg, unityThreeUniqueRelNormPercentErrorNeg, predThreeUniqueRelNormPercentErrorNeg] =...
%          MLR(threeUniqueTermsNeg, refThreeUniqueZmnNeg, singDataThresh, addBias, emptyWeightsForUnityData);
%     
%     end
    %=========================
    threeUniqueTermsPos = termStruct.threeUniqueTerms;
    refThreeUniqueZmnPos = termStruct.refThreeUniqueZmn;
    
    [~, numCol] = size(indicesStruct.threeUniqueIndices);
    if ( numCol == 3)
        threeUniqueIndicesPos = indicesStruct.threeUniqueIndices(:,3);
        [threeUniqueWeightsPos, ~,~, unityThreeUniqueMSEPos, predThreeUniqueMSEPos,...
            unityThreeUniqueRelNormPercentErrorPos, predThreeUniqueRelNormPercentErrorPos, unityThreeUniqueClusterRelNormPercentErrorPos,predThreeUniqueClusterRelNormPercentErrorPos] =...
            calcClusterWeights(threeUniqueTermsPos, refThreeUniqueZmnPos, threeUniqueIndicesPos, singDataThresh, addBias);   
    else
         [threeUniqueWeightsPos, ~, ~, unityThreeUniqueMSEPos, predThreeUniqueMSEPos, unityThreeUniqueRelNormPercentErrorPos, predThreeUniqueRelNormPercentErrorPos] =...
         MLR(threeUniqueTermsPos, refThreeUniqueZmnPos, singDataThresh, addBias, emptyWeightsForUnityData);
     
    end
    % ------------ FOUR UNIQUE MLR  -----------
    %positive group
    ind = find(propStruct.fourUniqueSwapGroups == 0);
    fourUniqueIndicesPos = indicesStruct.fourUniqueIndices(:,3);
    fourUniqueIndicesPos = fourUniqueIndicesPos(ind,:);
    fourUniqueTermsPos = termStruct.fourUniqueTerms(ind,:);
    refFourUniqueZmnPos = termStruct.refFourUniqueZmn(ind,:);
    [fourUniqueWeightsPos, ~,~, unityFourUniqueMSEPos, predFourUniqueMSEPos,...
    unityFourUniqueRelNormPercentErrorPos, predFourUniqueRelNormPercentErrorPos, unityFourUniqueClusterRelNormPercentErrorPos, predFourUniqueClusterRelNormPercentErrorPos] =...
        calcClusterWeights(fourUniqueTermsPos, refFourUniqueZmnPos, fourUniqueIndicesPos, singDataThresh, addBias);    
    
    %z = getZElements(mlmom.indicesStruct.fourUniqueIndices,zMatrices.values,13);
    
    %negative group
    ind = find(propStruct.fourUniqueSwapGroups == 1);
    fourUniqueIndicesNeg = indicesStruct.fourUniqueIndices(:,3);
    fourUniqueIndicesNeg = fourUniqueIndicesNeg(ind,:);
    fourUniqueTermsNeg = termStruct.fourUniqueTerms(ind,:);
    refFourUniqueZmnNeg = termStruct.refFourUniqueZmn(ind,:);
    [fourUniqueWeightsNeg, ~,~, unityFourUniqueMSENeg, predFourUniqueMSENeg,...
    unityFourUniqueRelNormPercentErrorNeg, predFourUniqueRelNormPercentErrorNeg, unityFourUniqueClusterRelNormPercentErrorNeg, predFourUniqueClusterRelNormPercentErrorNeg] =...
        calcClusterWeights(fourUniqueTermsNeg, refFourUniqueZmnNeg, fourUniqueIndicesNeg, singDataThresh, addBias); 
    
    %ind = find(fourUniqueIndicesNeg == 13);
    %z = refFourUniqueZmnNeg(ind);
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
    %weightModel.twoUniqueWeightsIndices = twoUniqueWeightsIndices;
    weightModel.unityTwoUniqueClusterRelNormPercentError = unityTwoUniqueClusterRelNormPercentError;
    weightModel.predTwoUniqueClusterRelNormPercentError = predTwoUniqueClusterRelNormPercentError;
    
    % 3 unique triangles
    %weightModel.triZmnTerms = triZmnTerms;  
    %weightModel.refTriZmn = refTriZmn;
    %weightModel.unityTriZmn = unityThreeUniqueZmn;
    %weightModel.predTriZmn = predThreeUniqueZmn;
    weightModel.threeUniqueWeightsPos = threeUniqueWeightsPos;
    %weightModel.threeUniqueWeightsIndices = threeUniqueWeightsIndices;
    weightModel.unityThreeUniqueClusterRelNormPercentErrorPos = unityThreeUniqueClusterRelNormPercentErrorPos;
    weightModel.predThreeUniqueClusterRelNormPercentErrorPos = predThreeUniqueClusterRelNormPercentErrorPos;
    
    %==========
    %weightModel.threeUniqueWeightsNeg = threeUniqueWeightsNeg;
   % weightModel.unityThreeUniqueClusterRelNormPercentErrorNeg = unityThreeUniqueClusterRelNormPercentErrorNeg;
    %weightModel.predThreeUniqueClusterRelNormPercentErrorNeg = predThreeUniqueClusterRelNormPercentErrorNeg;
    %==========
    
    % 4 unique
    %weightModel.nonSingZmnTerms = nonSingZmnTerms;  
    %weightModel.refNonSingZmn = refNonSingZmn;
    %weightModel.unityFourUniqueZmn = unityFourUniqueZmn;
    %weightModel.predFourUniqueZmn = predFourUniqueZmn; 
    weightModel.fourUniqueWeightsPos = fourUniqueWeightsPos;
    %weightModel.fourUniqueWeightsIndicesPos = fourUniqueWeightsIndicesPos;
    weightModel.unityFourUniqueClusterRelNormPercentErrorPos = unityFourUniqueClusterRelNormPercentErrorPos;
    weightModel.predFourUniqueClusterRelNormPercentErrorPos = predFourUniqueClusterRelNormPercentErrorPos;
    
    weightModel.fourUniqueWeightsNeg = fourUniqueWeightsNeg;
    %weightModel.fourUniqueWeightsIndicesNeg = fourUniqueWeightsIndicesNeg;
    weightModel.unityFourUniqueClusterRelNormPercentErrorNeg = unityFourUniqueClusterRelNormPercentErrorNeg;
    weightModel.predFourUniqueClusterRelNormPercentErrorNeg = predFourUniqueClusterRelNormPercentErrorNeg;
    
    % ------------
    % -----Scalars
    
   % weightModel.numNonSing = numNonSing;

    % 2 Unique
    weightModel.predTwoUniqueRelNormPercentError = predTwoUniqueRelNormPercentError;
    weightModel.unityTwoUniqueRelNormPercentError = unityTwoUniqueRelNormPercentError;
    weightModel.predTwoUniqueMSE = predTwoUniqueMSE;    
    weightModel.unityTwoUniqueMSE = unityTwoUniqueMSE;
    
    % 3 unique triangles
    weightModel.predThreeUniqueRelNormPercentErrorPos = predThreeUniqueRelNormPercentErrorPos;
    weightModel.unityThreeUniqueRelNormPercentErrorPos = unityThreeUniqueRelNormPercentErrorPos;
    weightModel.predThreeUniqueMSEPos= predThreeUniqueMSEPos;     
    weightModel.unityThreeUniqueMSEPos = unityThreeUniqueMSEPos;
    
    %==========
    %weightModel.predThreeUniqueRelNormPercentErrorNeg = predThreeUniqueRelNormPercentErrorNeg;
    %weightModel.unityThreeUniqueRelNormPercentErrorNeg = unityThreeUniqueRelNormPercentErrorNeg;
    %weightModel.predThreeUniqueMSENeg = predThreeUniqueMSENeg;     
    %weightModel.unityThreeUniqueMSENeg = unityThreeUniqueMSENeg;
    %==========
    
    % 4 Unique
    weightModel.predFourUniqueRelNormPercentErrorPos = predFourUniqueRelNormPercentErrorPos;
    weightModel.unityFourUniqueRelNormPercentErrorPos = unityFourUniqueRelNormPercentErrorPos;
    weightModel.predFourUniqueMSEPos = predFourUniqueMSEPos;  
    weightModel.unityFourUniqueMSEPos = unityFourUniqueMSEPos;
    
    weightModel.predFourUniqueRelNormPercentErrorNeg = predFourUniqueRelNormPercentErrorNeg;
    weightModel.unityFourUniqueRelNormPercentErrorNeg = unityFourUniqueRelNormPercentErrorNeg;
    weightModel.predFourUniqueMSENeg = predFourUniqueMSENeg;  
    weightModel.unityFourUniqueMSENeg = unityFourUniqueMSENeg;
    
end

function  [clusterWeights, pred,unity, unityMSE, predMSE,...
    unityRelNormPercentError, predRelNormPercentError, unityClusterRelNormPercentError, predClusterRelNormPercentError] =...
    calcClusterWeights(terms, refVals, clusterInd, singThresh, addBias)
    %edgeLabels only used if useProjectedEdges =1
    
    %[predNonSingMSE, predNonSingRelNormPercentError] = calcError(refVals, pred);
    %[unityNonSingMSE, unityNonSingRelNormPercentError] = calcError(refVals, unityNonSingZmn);
    
    [numObs, numTerms] = size(terms);
    %[numObs, 1] = size(ClusterInd);
    %[numObs, 1] = size(refVals)
    [numClusters ,~] = max(clusterInd);
    clusterWeights = ones(numClusters,numTerms );
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
          if (false)
          %if (urnpe - prnpe < minPercentImprov)
              pred(ind,1) = unity;
          else
              nonUnityClusterWeightsCount = nonUnityClusterWeightsCount + 1;
              count = count + numel(ref);
              clusterWeights(nonUnityClusterWeightsCount,:) = weights;
              predClusterRelNormPercentError(nonUnityClusterWeightsCount) = prnpe;
              unityClusterRelNormPercentError(nonUnityClusterWeightsCount) = urnpe;
          end
            
        end %if (eq(size(weights), [0 0]))
        
    end% for k
    clusterWeights = clusterWeights(1:nonUnityClusterWeightsCount, :);
    predClusterRelNormPercentError = predClusterRelNormPercentError(1:nonUnityClusterWeightsCount, :);
    unityClusterRelNormPercentError = unityClusterRelNormPercentError(1:nonUnityClusterWeightsCount, :);
    
    %error
    unity  = sum(terms, 2);
    [predMSE, predRelNormPercentError] = calcError(refVals, pred);
    [unityMSE, unityRelNormPercentError] = calcError(refVals, unity);

end
