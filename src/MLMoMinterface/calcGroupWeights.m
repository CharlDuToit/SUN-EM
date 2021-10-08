function [groupWeights, projZmn, predZmn] = calcGroupWeights(groupIndices, projTerms, refZmn, singDataThresh)
    %groupIndices: struct members are arrays of size [N,3]    
    %refZmn: [numNewEdges, numNewEdges] 
    %projTerms : [numNewEdges, numNewEdges, numTerms]
    
    %------------initialise------------
    
    groupWeights = [];
    returnEmptyWeightsIfUnity = 0;
    addOnes = 0;
    [numNewEdges, ~, ~] = size(projTerms);
    predZmn = zeros(numNewEdges,numNewEdges);
    projZmn = zeros(numNewEdges,numNewEdges);
    
    %------------2 unique------------
    %----- twoUnique_int
    
    %-----twoUnique_pos
    
    %-----twoUnique_neg
    
    %-----twoUnique_ext
    
    %------------3 unique------------
    %-----threeUnique_intInt
    
    %-----threeUnique_intPos
    
    %-----threeUnique_intNeg
    
    %-----threeUnique_posInt
    
    %-----threeUnique_negInt
    
    %-----threeUnique_extNotExt
    
    %-----threeUnique_notExtExt
    
    %-----threeUnique_parIntParInt
    
    %-----threeUnique_extExt
    
    %------------4 unique------------
    %-----fourUnique_intInt
    indices = groupIndices.fourUnique_intInt_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.fourUnique_intInt_weights = clustersWeights;
    groupWeights.fourUnique_intInt_projMSE = clustersProjMSE;
    groupWeights.fourUnique_intInt_predMSE = clustersPredMSE;
    groupWeights.fourUnique_intInt_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.fourUnique_intInt_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.fourUnique_intInt_totalProjMSE = totalProjMSE;
    groupWeights.fourUnique_intInt_totalPredMSE = totalPredMSE;
    groupWeights.fourUnique_intInt_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.fourUnique_intInt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
        %-----fourUnique_intPos
    
    %-----fourUnique_intNeg
    
    %-----fourUnique_posInt
    
    %-----fourUnique_negInt
    
    %-----fourUnique_posPos
    
    %-----fourUnique_posNeg
    
    %-----fourUnique_negPos
    
    %-----fourUnique_negNeg
    
    %-----fourUnique_notExtNotExt
    
    %-----fourUnique_notExtExt
    
    %-----fourUnique_extNotExt
    
    %-----fourUnique_extExt
    
    
             
end

function [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity )
    
    mmInd = indices(:,1);
    nnInd = indices(:,2);
    numClusters = max(indices(:,3));
    [~, ~, numTerms] = size(projTerms);
    clustersWeights = zeros(numClusters, numTerms);
    clustersProjMSE = zeros(numClusters, 1);
    clustersPredMSE = zeros(numClusters, 1);
    clustersProjRelNormPercentError = zeros(numClusters, 1);
    clustersPredRelNormPercentError = zeros(numClusters, 1);
    clustersTerms = zeros(numel(mmInd), numTerms);
    clustersRef =zeros(numel(mmInd), 1);
    clustersPred = zeros(numel(mmInd), 1);
    clustersProj = zeros(numel(mmInd), 1);
    for k = 1:numel(mmInd)
        clustersTerms(k,:) = projTerms(mmInd(k), nnInd(k), :);
        clustersRef(k) = refZmn(mmInd(k), nnInd(k));
    end
    
    totalProjSquaredError = 0;
    totalPredSquaredError= 0;
    for k = 1:numClusters
        ind = find( indices(:,3) == k );
        terms = clustersTerms(ind, :);
        ref = clustersRef(ind);

        [weights, pred, proj, projMSE, predMSE, projRelNormPercentError, predRelNormPercentError] =...
            MLR(terms,ref, singDataThresh, addOnes, returnEmptyWeightsIfUnity);
        clustersWeights(k, :) = weights;
        clustersPred(ind) = pred;
        clustersProj(ind) = proj;
        clustersProjMSE(k) = projMSE;
        clustersPredMSE(k) = predMSE;
        clustersProjRelNormPercentError(k) = projRelNormPercentError;
        clustersPredRelNormPercentError(k) = predRelNormPercentError;
        
        totalProjSquaredError = totalProjSquaredError + projMSE * numel(ind);
        totalPredSquaredError = totalPredSquaredError + predMSE * numel(ind);
        %pred ./ ref;
        %proj ./ ref;
    end
    totalProjRelNormPercentError = 100 * sqrt(totalProjSquaredError /sum(ref.^2) );
    totalPredRelNormPercentError = 100 * sqrt(totalPredSquaredError /sum(ref.^2) );
    
    totalProjMSE = totalProjSquaredError ./ numel(mmInd);
    totalPredMSE = totalPredSquaredError ./ numel(mmInd);
    
    %100* sqrt(( (clustersPred - clustersRef)' * (clustersPred - clustersRef) )/ sum(ref.^2));

end