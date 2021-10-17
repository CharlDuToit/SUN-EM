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
    indices = groupIndices.twoUnique_int_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.twoUnique_int_weights = clustersWeights;
    groupWeights.twoUnique_int_projMSE = clustersProjMSE;
    groupWeights.twoUnique_int_predMSE = clustersPredMSE;
    groupWeights.twoUnique_int_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.twoUnique_int_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.twoUnique_int_totalProjMSE = totalProjMSE;
    groupWeights.twoUnique_int_totalPredMSE = totalPredMSE;
    groupWeights.twoUnique_int_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.twoUnique_int_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----twoUnique_pos
    indices = groupIndices.twoUnique_pos_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.twoUnique_pos_weights = clustersWeights;
    groupWeights.twoUnique_pos_projMSE = clustersProjMSE;
    groupWeights.twoUnique_pos_predMSE = clustersPredMSE;
    groupWeights.twoUnique_pos_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.twoUnique_pos_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.twoUnique_pos_totalProjMSE = totalProjMSE;
    groupWeights.twoUnique_pos_totalPredMSE = totalPredMSE;
    groupWeights.twoUnique_pos_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.twoUnique_pos_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----twoUnique_neg
    indices = groupIndices.twoUnique_neg_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.twoUnique_neg_weights = clustersWeights;
    groupWeights.twoUnique_neg_projMSE = clustersProjMSE;
    groupWeights.twoUnique_neg_predMSE = clustersPredMSE;
    groupWeights.twoUnique_neg_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.twoUnique_neg_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.twoUnique_neg_totalProjMSE = totalProjMSE;
    groupWeights.twoUnique_neg_totalPredMSE = totalPredMSE;
    groupWeights.twoUnique_neg_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.twoUnique_neg_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----twoUnique_ext
    indices = groupIndices.twoUnique_ext_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.twoUnique_ext_weights = clustersWeights;
    groupWeights.twoUnique_ext_projMSE = clustersProjMSE;
    groupWeights.twoUnique_ext_predMSE = clustersPredMSE;
    groupWeights.twoUnique_ext_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.twoUnique_ext_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.twoUnique_ext_totalProjMSE = totalProjMSE;
    groupWeights.twoUnique_ext_totalPredMSE = totalPredMSE;
    groupWeights.twoUnique_ext_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.twoUnique_ext_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %------------3 unique------------
    %-----threeUnique_intInt
    indices = groupIndices.threeUnique_intInt_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.threeUnique_intInt_weights = clustersWeights;
    groupWeights.threeUnique_intInt_projMSE = clustersProjMSE;
    groupWeights.threeUnique_intInt_predMSE = clustersPredMSE;
    groupWeights.threeUnique_intInt_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.threeUnique_intInt_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.threeUnique_intInt_totalProjMSE = totalProjMSE;
    groupWeights.threeUnique_intInt_totalPredMSE = totalPredMSE;
    groupWeights.threeUnique_intInt_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.threeUnique_intInt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----threeUnique_intPos
    indices = groupIndices.threeUnique_intPos_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.threeUnique_intPos_weights = clustersWeights;
    groupWeights.threeUnique_intPos_projMSE = clustersProjMSE;
    groupWeights.threeUnique_intPos_predMSE = clustersPredMSE;
    groupWeights.threeUnique_intPos_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.threeUnique_intPos_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.threeUnique_intPos_totalProjMSE = totalProjMSE;
    groupWeights.threeUnique_intPos_totalPredMSE = totalPredMSE;
    groupWeights.threeUnique_intPos_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.threeUnique_intPos_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----threeUnique_intNeg
    indices = groupIndices.threeUnique_intNeg_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.threeUnique_intNeg_weights = clustersWeights;
    groupWeights.threeUnique_intNeg_projMSE = clustersProjMSE;
    groupWeights.threeUnique_intNeg_predMSE = clustersPredMSE;
    groupWeights.threeUnique_intNeg_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.threeUnique_intNeg_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.threeUnique_intNeg_totalProjMSE = totalProjMSE;
    groupWeights.threeUnique_intNeg_totalPredMSE = totalPredMSE;
    groupWeights.threeUnique_intNeg_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.threeUnique_intNeg_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    %-----threeUnique_posInt
    indices = groupIndices.threeUnique_posInt_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.threeUnique_posInt_weights = clustersWeights;
    groupWeights.threeUnique_posInt_projMSE = clustersProjMSE;
    groupWeights.threeUnique_posInt_predMSE = clustersPredMSE;
    groupWeights.threeUnique_posInt_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.threeUnique_posInt_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.threeUnique_posInt_totalProjMSE = totalProjMSE;
    groupWeights.threeUnique_posInt_totalPredMSE = totalPredMSE;
    groupWeights.threeUnique_posInt_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.threeUnique_posInt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----threeUnique_negInt
    indices = groupIndices.threeUnique_negInt_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.threeUnique_negInt_weights = clustersWeights;
    groupWeights.threeUnique_negInt_projMSE = clustersProjMSE;
    groupWeights.threeUnique_negInt_predMSE = clustersPredMSE;
    groupWeights.threeUnique_negInt_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.threeUnique_negInt_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.threeUnique_negInt_totalProjMSE = totalProjMSE;
    groupWeights.threeUnique_negInt_totalPredMSE = totalPredMSE;
    groupWeights.threeUnique_negInt_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.threeUnique_negInt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----threeUnique_extNotExt
    indices = groupIndices.threeUnique_extNotExt_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.threeUnique_extNotExt_weights = clustersWeights;
    groupWeights.threeUnique_extNotExt_projMSE = clustersProjMSE;
    groupWeights.threeUnique_extNotExt_predMSE = clustersPredMSE;
    groupWeights.threeUnique_extNotExt_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.threeUnique_extNotExt_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.threeUnique_extNotExt_totalProjMSE = totalProjMSE;
    groupWeights.threeUnique_extNotExt_totalPredMSE = totalPredMSE;
    groupWeights.threeUnique_extNotExt_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.threeUnique_extNotExt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----threeUnique_notExtExt
    indices = groupIndices.threeUnique_notExtExt_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.threeUnique_notExtExt_weights = clustersWeights;
    groupWeights.threeUnique_notExtExt_projMSE = clustersProjMSE;
    groupWeights.threeUnique_notExtExt_predMSE = clustersPredMSE;
    groupWeights.threeUnique_notExtExt_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.threeUnique_notExtExt_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.threeUnique_notExtExt_totalProjMSE = totalProjMSE;
    groupWeights.threeUnique_notExtExt_totalPredMSE = totalPredMSE;
    groupWeights.threeUnique_notExtExt_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.threeUnique_notExtExt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    %-----threeUnique_parIntParInt
    indices = groupIndices.threeUnique_parIntParInt_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.threeUnique_parIntParInt_weights = clustersWeights;
    groupWeights.threeUnique_parIntParInt_projMSE = clustersProjMSE;
    groupWeights.threeUnique_parIntParInt_predMSE = clustersPredMSE;
    groupWeights.threeUnique_parIntParInt_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.threeUnique_parIntParInt_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.threeUnique_parIntParInt_totalProjMSE = totalProjMSE;
    groupWeights.threeUnique_parIntParInt_totalPredMSE = totalPredMSE;
    groupWeights.threeUnique_parIntParInt_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.threeUnique_parIntParInt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----threeUnique_extExt
    if (numel(groupIndices.threeUnique_extExt_indices(:,1)) > 0)
        indices = groupIndices.threeUnique_extExt_indices;
        [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
            clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
            clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
        
        groupWeights.threeUnique_extExt_weights = clustersWeights;
        groupWeights.threeUnique_extExt_projMSE = clustersProjMSE;
        groupWeights.threeUnique_extExt_predMSE = clustersPredMSE;
        groupWeights.threeUnique_extExt_projRelNormPercentError = clustersProjRelNormPercentError;
        groupWeights.threeUnique_extExt_predRelNormPercentError = clustersPredRelNormPercentError;
        groupWeights.threeUnique_extExt_totalProjMSE = totalProjMSE;
        groupWeights.threeUnique_extExt_totalPredMSE = totalPredMSE;
        groupWeights.threeUnique_extExt_totalProjRelNormPercentError= totalProjRelNormPercentError;
        groupWeights.threeUnique_extExt_totalPredRelNormPercentError= totalPredRelNormPercentError;
        for k = 1:numel(indices(:,3))
            projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
            predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
        end
        
    else
        groupWeights.threeUnique_extExt_weights = [];
    end
    
    
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
    indices = groupIndices.fourUnique_intPos_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.fourUnique_intPos_weights = clustersWeights;
    groupWeights.fourUnique_intPos_projMSE = clustersProjMSE;
    groupWeights.fourUnique_intPos_predMSE = clustersPredMSE;
    groupWeights.fourUnique_intPos_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.fourUnique_intPos_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.fourUnique_intPos_totalProjMSE = totalProjMSE;
    groupWeights.fourUnique_intPos_totalPredMSE = totalPredMSE;
    groupWeights.fourUnique_intPos_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.fourUnique_intPos_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----fourUnique_intNeg
    indices = groupIndices.fourUnique_intNeg_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.fourUnique_intNeg_weights = clustersWeights;
    groupWeights.fourUnique_intNeg_projMSE = clustersProjMSE;
    groupWeights.fourUnique_intNeg_predMSE = clustersPredMSE;
    groupWeights.fourUnique_intNeg_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.fourUnique_intNeg_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.fourUnique_intNeg_totalProjMSE = totalProjMSE;
    groupWeights.fourUnique_intNeg_totalPredMSE = totalPredMSE;
    groupWeights.fourUnique_intNeg_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.fourUnique_intNeg_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----fourUnique_posInt
    indices = groupIndices.fourUnique_posInt_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.fourUnique_posInt_weights = clustersWeights;
    groupWeights.fourUnique_posInt_projMSE = clustersProjMSE;
    groupWeights.fourUnique_posInt_predMSE = clustersPredMSE;
    groupWeights.fourUnique_posInt_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.fourUnique_posInt_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.fourUnique_posInt_totalProjMSE = totalProjMSE;
    groupWeights.fourUnique_posInt_totalPredMSE = totalPredMSE;
    groupWeights.fourUnique_posInt_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.fourUnique_posInt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----fourUnique_negInt
    indices = groupIndices.fourUnique_negInt_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.fourUnique_negInt_weights = clustersWeights;
    groupWeights.fourUnique_negInt_projMSE = clustersProjMSE;
    groupWeights.fourUnique_negInt_predMSE = clustersPredMSE;
    groupWeights.fourUnique_negInt_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.fourUnique_negInt_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.fourUnique_negInt_totalProjMSE = totalProjMSE;
    groupWeights.fourUnique_negInt_totalPredMSE = totalPredMSE;
    groupWeights.fourUnique_negInt_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.fourUnique_negInt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----fourUnique_posPos
    indices = groupIndices.fourUnique_posPos_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.fourUnique_posPos_weights = clustersWeights;
    groupWeights.fourUnique_posPos_projMSE = clustersProjMSE;
    groupWeights.fourUnique_posPos_predMSE = clustersPredMSE;
    groupWeights.fourUnique_posPos_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.fourUnique_posPos_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.fourUnique_posPos_totalProjMSE = totalProjMSE;
    groupWeights.fourUnique_posPos_totalPredMSE = totalPredMSE;
    groupWeights.fourUnique_posPos_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.fourUnique_posPos_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----fourUnique_posNeg
    indices = groupIndices.fourUnique_posNeg_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.fourUnique_posNeg_weights = clustersWeights;
    groupWeights.fourUnique_posNeg_projMSE = clustersProjMSE;
    groupWeights.fourUnique_posNeg_predMSE = clustersPredMSE;
    groupWeights.fourUnique_posNeg_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.fourUnique_posNeg_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.fourUnique_posNeg_totalProjMSE = totalProjMSE;
    groupWeights.fourUnique_posNeg_totalPredMSE = totalPredMSE;
    groupWeights.fourUnique_posNeg_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.fourUnique_posNeg_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----fourUnique_negPos
    indices = groupIndices.fourUnique_negPos_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.fourUnique_negPos_weights = clustersWeights;
    groupWeights.fourUnique_negPos_projMSE = clustersProjMSE;
    groupWeights.fourUnique_negPos_predMSE = clustersPredMSE;
    groupWeights.fourUnique_negPos_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.fourUnique_negPos_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.fourUnique_negPos_totalProjMSE = totalProjMSE;
    groupWeights.fourUnique_negPos_totalPredMSE = totalPredMSE;
    groupWeights.fourUnique_negPos_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.fourUnique_negPos_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----fourUnique_negNeg
    indices = groupIndices.fourUnique_negNeg_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.fourUnique_negNeg_weights = clustersWeights;
    groupWeights.fourUnique_negNeg_projMSE = clustersProjMSE;
    groupWeights.fourUnique_negNeg_predMSE = clustersPredMSE;
    groupWeights.fourUnique_negNeg_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.fourUnique_negNeg_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.fourUnique_negNeg_totalProjMSE = totalProjMSE;
    groupWeights.fourUnique_negNeg_totalPredMSE = totalPredMSE;
    groupWeights.fourUnique_negNeg_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.fourUnique_negNeg_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----fourUnique_notExtNotExt
    indices = groupIndices.fourUnique_notExtNotExt_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.fourUnique_notExtNotExt_weights = clustersWeights;
    groupWeights.fourUnique_notExtNotExt_projMSE = clustersProjMSE;
    groupWeights.fourUnique_notExtNotExt_predMSE = clustersPredMSE;
    groupWeights.fourUnique_notExtNotExt_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.fourUnique_notExtNotExt_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.fourUnique_notExtNotExt_totalProjMSE = totalProjMSE;
    groupWeights.fourUnique_notExtNotExt_totalPredMSE = totalPredMSE;
    groupWeights.fourUnique_notExtNotExt_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.fourUnique_notExtNotExt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----fourUnique_notExtExt
    indices = groupIndices.fourUnique_notExtExt_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.fourUnique_notExtExt_weights = clustersWeights;
    groupWeights.fourUnique_notExtExt_projMSE = clustersProjMSE;
    groupWeights.fourUnique_notExtExt_predMSE = clustersPredMSE;
    groupWeights.fourUnique_notExtExt_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.fourUnique_notExtExt_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.fourUnique_notExtExt_totalProjMSE = totalProjMSE;
    groupWeights.fourUnique_notExtExt_totalPredMSE = totalPredMSE;
    groupWeights.fourUnique_notExtExt_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.fourUnique_notExtExt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----fourUnique_extNotExt
    indices = groupIndices.fourUnique_extNotExt_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.fourUnique_extNotExt_weights = clustersWeights;
    groupWeights.fourUnique_extNotExt_projMSE = clustersProjMSE;
    groupWeights.fourUnique_extNotExt_predMSE = clustersPredMSE;
    groupWeights.fourUnique_extNotExt_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.fourUnique_extNotExt_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.fourUnique_extNotExt_totalProjMSE = totalProjMSE;
    groupWeights.fourUnique_extNotExt_totalPredMSE = totalPredMSE;
    groupWeights.fourUnique_extNotExt_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.fourUnique_extNotExt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
    %-----fourUnique_extExt
    indices = groupIndices.fourUnique_extExt_indices;  
    [clustersWeights, clustersProj,clustersPred , clustersProjMSE, clustersPredMSE,clustersProjRelNormPercentError,...
    clustersPredRelNormPercentError, totalProjMSE, totalPredMSE, totalProjRelNormPercentError, totalPredRelNormPercentError ] =...
    clusterMLR(indices,projTerms, refZmn, singDataThresh, addOnes, returnEmptyWeightsIfUnity );
      
    groupWeights.fourUnique_extExt_weights = clustersWeights;
    groupWeights.fourUnique_extExt_projMSE = clustersProjMSE;
    groupWeights.fourUnique_extExt_predMSE = clustersPredMSE;
    groupWeights.fourUnique_extExt_projRelNormPercentError = clustersProjRelNormPercentError;
    groupWeights.fourUnique_extExt_predRelNormPercentError = clustersPredRelNormPercentError;
    groupWeights.fourUnique_extExt_totalProjMSE = totalProjMSE;
    groupWeights.fourUnique_extExt_totalPredMSE = totalPredMSE;
    groupWeights.fourUnique_extExt_totalProjRelNormPercentError= totalProjRelNormPercentError;
    groupWeights.fourUnique_extExt_totalPredRelNormPercentError= totalPredRelNormPercentError;
    for k = 1:numel(indices(:,3))
        projZmn(indices(k,1),  indices(k,2) ) = clustersProj(k);
        predZmn(indices(k,1),  indices(k,2) ) = clustersPred(k);
    end
    
             
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

    end % for k = 1:numClusters
    totalProjRelNormPercentError = 100 * sqrt(totalProjSquaredError /sum(clustersRef.^2) );
    totalPredRelNormPercentError = 100 * sqrt(totalPredSquaredError /sum(clustersRef.^2) );
    
    totalProjMSE = totalProjSquaredError ./ numel(mmInd);
    totalPredMSE = totalPredSquaredError ./ numel(mmInd);
    
    %100* sqrt(( (clustersPred - clustersRef)' * (clustersPred - clustersRef) )/ sum(ref.^2));

end