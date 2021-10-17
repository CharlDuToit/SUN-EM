function [groupErrors, projZmn, predZmn] = applyGroupWeights(weightModel, groupIndices, projTerms, refZmn, includeError)
    %groupIndices: struct members are arrays of size [N,3]    
    %refZmn: [numNewEdges, numNewEdges] 
    %projTerms : [numNewEdges, numNewEdges, numTerms]
    
    %------------initialise------------
    
    groupErrors = [];
    [numNewEdges, ~, ~] = size(projTerms);
    predZmn = zeros(numNewEdges,numNewEdges);
    projZmn = zeros(numNewEdges,numNewEdges);
    
    %------------2 unique------------
    %----- twoUnique_int
    indices = groupIndices.twoUnique_int_indices; 
    weights = weightModel.twoUnique_int_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    if (includeError)
        groupErrors.twoUnique_int_projMSE = clustersProjMSE;
        groupErrors.twoUnique_int_predMSE = clustersPredMSE;
        groupErrors.twoUnique_int_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.twoUnique_int_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.twoUnique_int_totalProjMSE = totalProjMSE;
        groupErrors.twoUnique_int_totalPredMSE = totalPredMSE;
        groupErrors.twoUnique_int_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.twoUnique_int_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----twoUnique_pos
    indices = groupIndices.twoUnique_pos_indices;  
    weights = weightModel.twoUnique_pos_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    if (includeError)
        groupErrors.twoUnique_pos_projMSE = clustersProjMSE;
        groupErrors.twoUnique_pos_predMSE = clustersPredMSE;
        groupErrors.twoUnique_pos_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.twoUnique_pos_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.twoUnique_pos_totalProjMSE = totalProjMSE;
        groupErrors.twoUnique_pos_totalPredMSE = totalPredMSE;
        groupErrors.twoUnique_pos_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.twoUnique_pos_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----twoUnique_neg
    indices = groupIndices.twoUnique_neg_indices;
    weights = weightModel.twoUnique_neg_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    if (includeError)
        groupErrors.twoUnique_neg_projMSE = clustersProjMSE;
        groupErrors.twoUnique_neg_predMSE = clustersPredMSE;
        groupErrors.twoUnique_neg_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.twoUnique_neg_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.twoUnique_neg_totalProjMSE = totalProjMSE;
        groupErrors.twoUnique_neg_totalPredMSE = totalPredMSE;
        groupErrors.twoUnique_neg_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.twoUnique_neg_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----twoUnique_ext
    indices = groupIndices.twoUnique_ext_indices; 
    weights = weightModel.twoUnique_ext_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    if (includeError)
        groupErrors.twoUnique_ext_projMSE = clustersProjMSE;
        groupErrors.twoUnique_ext_predMSE = clustersPredMSE;
        groupErrors.twoUnique_ext_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.twoUnique_ext_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.twoUnique_ext_totalProjMSE = totalProjMSE;
        groupErrors.twoUnique_ext_totalPredMSE = totalPredMSE;
        groupErrors.twoUnique_ext_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.twoUnique_ext_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %------------3 unique------------
    %-----threeUnique_intInt
    indices = groupIndices.threeUnique_intInt_indices; 
    weights = weightModel.threeUnique_intInt_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    if (includeError)
        groupErrors.threeUnique_intInt_projMSE = clustersProjMSE;
        groupErrors.threeUnique_intInt_predMSE = clustersPredMSE;
        groupErrors.threeUnique_intInt_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.threeUnique_intInt_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.threeUnique_intInt_totalProjMSE = totalProjMSE;
        groupErrors.threeUnique_intInt_totalPredMSE = totalPredMSE;
        groupErrors.threeUnique_intInt_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.threeUnique_intInt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----threeUnique_intPos
    indices = groupIndices.threeUnique_intPos_indices;
    weights = weightModel.threeUnique_intPos_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    if (includeError)
        groupErrors.threeUnique_intPos_projMSE = clustersProjMSE;
        groupErrors.threeUnique_intPos_predMSE = clustersPredMSE;
        groupErrors.threeUnique_intPos_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.threeUnique_intPos_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.threeUnique_intPos_totalProjMSE = totalProjMSE;
        groupErrors.threeUnique_intPos_totalPredMSE = totalPredMSE;
        groupErrors.threeUnique_intPos_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.threeUnique_intPos_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----threeUnique_intNeg
    indices = groupIndices.threeUnique_intNeg_indices;  
    weights = weightModel.threeUnique_intNeg_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    if (includeError)
        groupErrors.threeUnique_intNeg_projMSE = clustersProjMSE;
        groupErrors.threeUnique_intNeg_predMSE = clustersPredMSE;
        groupErrors.threeUnique_intNeg_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.threeUnique_intNeg_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.threeUnique_intNeg_totalProjMSE = totalProjMSE;
        groupErrors.threeUnique_intNeg_totalPredMSE = totalPredMSE;
        groupErrors.threeUnique_intNeg_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.threeUnique_intNeg_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    %-----threeUnique_posInt
    indices = groupIndices.threeUnique_posInt_indices; 
    weights = weightModel.threeUnique_posInt_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    if (includeError)
        groupErrors.threeUnique_posInt_projMSE = clustersProjMSE;
        groupErrors.threeUnique_posInt_predMSE = clustersPredMSE;
        groupErrors.threeUnique_posInt_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.threeUnique_posInt_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.threeUnique_posInt_totalProjMSE = totalProjMSE;
        groupErrors.threeUnique_posInt_totalPredMSE = totalPredMSE;
        groupErrors.threeUnique_posInt_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.threeUnique_posInt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----threeUnique_negInt
    indices = groupIndices.threeUnique_negInt_indices;  
    weights = weightModel.threeUnique_negInt_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    if (includeError)
        groupErrors.threeUnique_negInt_projMSE = clustersProjMSE;
        groupErrors.threeUnique_negInt_predMSE = clustersPredMSE;
        groupErrors.threeUnique_negInt_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.threeUnique_negInt_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.threeUnique_negInt_totalProjMSE = totalProjMSE;
        groupErrors.threeUnique_negInt_totalPredMSE = totalPredMSE;
        groupErrors.threeUnique_negInt_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.threeUnique_negInt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----threeUnique_extNotExt
    indices = groupIndices.threeUnique_extNotExt_indices;  
    weights = weightModel.threeUnique_extNotExt_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    if (includeError)
        groupErrors.threeUnique_extNotExt_projMSE = clustersProjMSE;
        groupErrors.threeUnique_extNotExt_predMSE = clustersPredMSE;
        groupErrors.threeUnique_extNotExt_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.threeUnique_extNotExt_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.threeUnique_extNotExt_totalProjMSE = totalProjMSE;
        groupErrors.threeUnique_extNotExt_totalPredMSE = totalPredMSE;
        groupErrors.threeUnique_extNotExt_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.threeUnique_extNotExt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----threeUnique_notExtExt
    indices = groupIndices.threeUnique_notExtExt_indices;
    weights = weightModel.threeUnique_notExtExt_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    if (includeError)
        groupErrors.threeUnique_notExtExt_projMSE = clustersProjMSE;
        groupErrors.threeUnique_notExtExt_predMSE = clustersPredMSE;
        groupErrors.threeUnique_notExtExt_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.threeUnique_notExtExt_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.threeUnique_notExtExt_totalProjMSE = totalProjMSE;
        groupErrors.threeUnique_notExtExt_totalPredMSE = totalPredMSE;
        groupErrors.threeUnique_notExtExt_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.threeUnique_notExtExt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    %-----threeUnique_parIntParInt
    indices = groupIndices.threeUnique_parIntParInt_indices; 
    weights = weightModel.threeUnique_parIntParInt_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    if (includeError)
        groupErrors.threeUnique_parIntParInt_projMSE = clustersProjMSE;
        groupErrors.threeUnique_parIntParInt_predMSE = clustersPredMSE;
        groupErrors.threeUnique_parIntParInt_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.threeUnique_parIntParInt_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.threeUnique_parIntParInt_totalProjMSE = totalProjMSE;
        groupErrors.threeUnique_parIntParInt_totalPredMSE = totalPredMSE;
        groupErrors.threeUnique_parIntParInt_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.threeUnique_parIntParInt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----threeUnique_extExt
    if (numel(groupIndices.threeUnique_extExt_indices(:,1)) > 0)
        indices = groupIndices.threeUnique_extExt_indices;
        weights = weightModel.threeUnique_extExt_weights;
        [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
            clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
            applyWeights(weights,indices,projTerms, refZmn, includeError );
        if (includeError)
            groupErrors.threeUnique_extExt_projMSE = clustersProjMSE;
            groupErrors.threeUnique_extExt_predMSE = clustersPredMSE;
            groupErrors.threeUnique_extExt_projRelNormPercentError = clustersProjRelNormPercentError;
            groupErrors.threeUnique_extExt_predRelNormPercentError = clustersPredRelNormPercentError;
            groupErrors.threeUnique_extExt_totalProjMSE = totalProjMSE;
            groupErrors.threeUnique_extExt_totalPredMSE = totalPredMSE;
            groupErrors.threeUnique_extExt_totalProjRelNormPercentError= totalProjRelNormPercentError;
            groupErrors.threeUnique_extExt_totalPredRelNormPercentError= totalPredRelNormPercentError;
        end
        for k = 1:numel(indices(:,3))
            projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
            predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
        end
    end
    
    
    %------------4 unique------------
    %-----fourUnique_intInt
    indices = groupIndices.fourUnique_intInt_indices;  
    weights = weightModel.fourUnique_intInt_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    if (includeError)
        groupErrors.fourUnique_intInt_projMSE = clustersProjMSE;
        groupErrors.fourUnique_intInt_predMSE = clustersPredMSE;
        groupErrors.fourUnique_intInt_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.fourUnique_intInt_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.fourUnique_intInt_totalProjMSE = totalProjMSE;
        groupErrors.fourUnique_intInt_totalPredMSE = totalPredMSE;
        groupErrors.fourUnique_intInt_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.fourUnique_intInt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----fourUnique_intPos
    indices = groupIndices.fourUnique_intPos_indices;
    weights = weightModel.fourUnique_intPos_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    if (includeError)
        groupErrors.fourUnique_intPos_projMSE = clustersProjMSE;
        groupErrors.fourUnique_intPos_predMSE = clustersPredMSE;
        groupErrors.fourUnique_intPos_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.fourUnique_intPos_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.fourUnique_intPos_totalProjMSE = totalProjMSE;
        groupErrors.fourUnique_intPos_totalPredMSE = totalPredMSE;
        groupErrors.fourUnique_intPos_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.fourUnique_intPos_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----fourUnique_intNeg
    indices = groupIndices.fourUnique_intNeg_indices;  
    weights = weightModel.fourUnique_intNeg_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    if (includeError)
        groupErrors.fourUnique_intNeg_projMSE = clustersProjMSE;
        groupErrors.fourUnique_intNeg_predMSE = clustersPredMSE;
        groupErrors.fourUnique_intNeg_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.fourUnique_intNeg_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.fourUnique_intNeg_totalProjMSE = totalProjMSE;
        groupErrors.fourUnique_intNeg_totalPredMSE = totalPredMSE;
        groupErrors.fourUnique_intNeg_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.fourUnique_intNeg_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----fourUnique_posInt
    indices = groupIndices.fourUnique_posInt_indices;
    weights = weightModel.fourUnique_posInt_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    if (includeError)
        groupErrors.fourUnique_posInt_projMSE = clustersProjMSE;
        groupErrors.fourUnique_posInt_predMSE = clustersPredMSE;
        groupErrors.fourUnique_posInt_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.fourUnique_posInt_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.fourUnique_posInt_totalProjMSE = totalProjMSE;
        groupErrors.fourUnique_posInt_totalPredMSE = totalPredMSE;
        groupErrors.fourUnique_posInt_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.fourUnique_posInt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end


    %-----fourUnique_negInt
    indices = groupIndices.fourUnique_negInt_indices;
    weights = weightModel.fourUnique_negInt_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    if (includeError)
        groupErrors.fourUnique_negInt_projMSE = clustersProjMSE;
        groupErrors.fourUnique_negInt_predMSE = clustersPredMSE;
        groupErrors.fourUnique_negInt_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.fourUnique_negInt_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.fourUnique_negInt_totalProjMSE = totalProjMSE;
        groupErrors.fourUnique_negInt_totalPredMSE = totalPredMSE;
        groupErrors.fourUnique_negInt_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.fourUnique_negInt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----fourUnique_posPos
    indices = groupIndices.fourUnique_posPos_indices;  
    weights = weightModel.fourUnique_posPos_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    if (includeError)
        groupErrors.fourUnique_posPos_projMSE = clustersProjMSE;
        groupErrors.fourUnique_posPos_predMSE = clustersPredMSE;
        groupErrors.fourUnique_posPos_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.fourUnique_posPos_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.fourUnique_posPos_totalProjMSE = totalProjMSE;
        groupErrors.fourUnique_posPos_totalPredMSE = totalPredMSE;
        groupErrors.fourUnique_posPos_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.fourUnique_posPos_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----fourUnique_posNeg
    indices = groupIndices.fourUnique_posNeg_indices;
    weights = weightModel.fourUnique_posNeg_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    if (includeError)
        groupErrors.fourUnique_posNeg_projMSE = clustersProjMSE;
        groupErrors.fourUnique_posNeg_predMSE = clustersPredMSE;
        groupErrors.fourUnique_posNeg_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.fourUnique_posNeg_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.fourUnique_posNeg_totalProjMSE = totalProjMSE;
        groupErrors.fourUnique_posNeg_totalPredMSE = totalPredMSE;
        groupErrors.fourUnique_posNeg_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.fourUnique_posNeg_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----fourUnique_negPos
    indices = groupIndices.fourUnique_negPos_indices;
    weights = weightModel.fourUnique_negPos_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    if (includeError)
        groupErrors.fourUnique_negPos_projMSE = clustersProjMSE;
        groupErrors.fourUnique_negPos_predMSE = clustersPredMSE;
        groupErrors.fourUnique_negPos_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.fourUnique_negPos_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.fourUnique_negPos_totalProjMSE = totalProjMSE;
        groupErrors.fourUnique_negPos_totalPredMSE = totalPredMSE;
        groupErrors.fourUnique_negPos_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.fourUnique_negPos_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----fourUnique_negNeg
    indices = groupIndices.fourUnique_negNeg_indices;
    weights = weightModel.fourUnique_negNeg_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    if (includeError)
        groupErrors.fourUnique_negNeg_projMSE = clustersProjMSE;
        groupErrors.fourUnique_negNeg_predMSE = clustersPredMSE;
        groupErrors.fourUnique_negNeg_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.fourUnique_negNeg_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.fourUnique_negNeg_totalProjMSE = totalProjMSE;
        groupErrors.fourUnique_negNeg_totalPredMSE = totalPredMSE;
        groupErrors.fourUnique_negNeg_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.fourUnique_negNeg_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----fourUnique_notExtNotExt
    indices = groupIndices.fourUnique_notExtNotExt_indices;
    weights = weightModel.fourUnique_notExtNotExt_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    if (includeError)
        groupErrors.fourUnique_notExtNotExt_projMSE = clustersProjMSE;
        groupErrors.fourUnique_notExtNotExt_predMSE = clustersPredMSE;
        groupErrors.fourUnique_notExtNotExt_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.fourUnique_notExtNotExt_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.fourUnique_notExtNotExt_totalProjMSE = totalProjMSE;
        groupErrors.fourUnique_notExtNotExt_totalPredMSE = totalPredMSE;
        groupErrors.fourUnique_notExtNotExt_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.fourUnique_notExtNotExt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----fourUnique_notExtExt
    indices = groupIndices.fourUnique_notExtExt_indices;
    weights = weightModel.fourUnique_notExtExt_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    if (includeError)
        groupErrors.fourUnique_notExtExt_projMSE = clustersProjMSE;
        groupErrors.fourUnique_notExtExt_predMSE = clustersPredMSE;
        groupErrors.fourUnique_notExtExt_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.fourUnique_notExtExt_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.fourUnique_notExtExt_totalProjMSE = totalProjMSE;
        groupErrors.fourUnique_notExtExt_totalPredMSE = totalPredMSE;
        groupErrors.fourUnique_notExtExt_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.fourUnique_notExtExt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----fourUnique_extNotExt
    indices = groupIndices.fourUnique_extNotExt_indices;
    weights = weightModel.fourUnique_extNotExt_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
        clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
        applyWeights(weights,indices,projTerms, refZmn, includeError );
    
    if (includeError)
        groupErrors.fourUnique_extNotExt_projMSE = clustersProjMSE;
        groupErrors.fourUnique_extNotExt_predMSE = clustersPredMSE;
        groupErrors.fourUnique_extNotExt_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.fourUnique_extNotExt_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.fourUnique_extNotExt_totalProjMSE = totalProjMSE;
        groupErrors.fourUnique_extNotExt_totalPredMSE = totalPredMSE;
        groupErrors.fourUnique_extNotExt_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.fourUnique_extNotExt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----fourUnique_extExt
    indices = groupIndices.fourUnique_extExt_indices;
    weights = weightModel.fourUnique_extExt_weights;
    [ clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    applyWeights(weights,indices,projTerms, refZmn, includeError );
      
    if (includeError)
        groupErrors.fourUnique_extExt_projMSE = clustersProjMSE;
        groupErrors.fourUnique_extExt_predMSE = clustersPredMSE;
        groupErrors.fourUnique_extExt_projRelNormPercentError = clustersProjRelNormPercentError;
        groupErrors.fourUnique_extExt_predRelNormPercentError = clustersPredRelNormPercentError;
        groupErrors.fourUnique_extExt_totalProjMSE = totalProjMSE;
        groupErrors.fourUnique_extExt_totalPredMSE = totalPredMSE;
        groupErrors.fourUnique_extExt_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupErrors.fourUnique_extExt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    end
    
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    
             
end

function [clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    applyWeights(clustersWeights, indices,projTerms, refZmn, includeError)
    
    mmInd = indices(:,1);
    nnInd = indices(:,2);
    numClusters = max(indices(:,3));
    [~, ~, numTerms] = size(projTerms);

    clustersTerms = zeros(numel(mmInd), numTerms);
    clustersPred = zeros(numel(mmInd), 1);
    clustersProj = zeros(numel(mmInd), 1);

    if (includeError)
        clustersRef =zeros(numel(mmInd), 1);
        clustersProjMSE = zeros(numClusters, 1);
        clustersPredMSE = zeros(numClusters, 1);
        clustersProjRelNormPercentError = zeros(numClusters, 1);
        clustersPredRelNormPercentError = zeros(numClusters, 1);
    else
        clustersRef = [];
        clustersProjMSE = [];
        clustersPredMSE = [];
        clustersProjRelNormPercentError = [];
        clustersPredRelNormPercentError = [];
    end

    for k = 1:numel(mmInd)
        clustersTerms(k,:) = projTerms(mmInd(k), nnInd(k), :);
        if (includeError)
            clustersRef(k) = refZmn(mmInd(k), nnInd(k));
        end
    end
    
    totalProjSquaredError = 0;
    totalPredSquaredError= 0;
    for k = 1:numClusters
        ind = find( indices(:,3) == k );
        clustersPred(ind) = clustersTerms(ind, :) * clustersWeights(k, :)';
        clustersProj(ind) = sum(clustersTerms(ind, :), 2);
        
        if (includeError)
            ref = clustersRef(ind);
            
            projMSE = (clustersProj(ind) - ref)' * (clustersProj(ind) - ref);
            projRelNormPercentError = 100* sqrt(projMSE/ sum(ref.^2) ) ;
            totalProjSquaredError = totalProjSquaredError + projMSE;
            projMSE = projMSE /numel(ind);
            
            predMSE = (clustersPred(ind) - ref)' * (clustersPred(ind) - ref);
            predRelNormPercentError = 100* sqrt(predMSE/ sum(ref.^2) ) ;
            totalPredSquaredError = totalPredSquaredError + predMSE ;
            predMSE = predMSE /numel(ind);
            
            clustersProjMSE(k) = projMSE;
            clustersPredMSE(k) = predMSE;
            clustersProjRelNormPercentError(k) = projRelNormPercentError;
            clustersPredRelNormPercentError(k) = predRelNormPercentError;
            
        end
        
    end %for k = 1:numClusters
    
    totalProjMSE =0;
    totalPredMSE=0;
    totalProjRelNormPercentError = 0; 
    totalPredRelNormPercentError = 0;
    if (includeError)
        totalProjRelNormPercentError = 100 * sqrt(totalProjSquaredError /sum(clustersRef.^2) );
        totalPredRelNormPercentError = 100 * sqrt(totalPredSquaredError /sum(clustersRef.^2) );
        
        totalProjMSE = totalProjSquaredError ./ numel(mmInd);
        totalPredMSE = totalPredSquaredError ./ numel(mmInd);
    end

end