function [groupMeans, groupIndices] = calcGroupMeans(groupMeans, groupIndices, maxIter, minTol)
    
    %------------2 unique------------
    %----- twoUnique_int
    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.twoUnique_int_properties, groupMeans.twoUnique_int_means,...
        maxIter, minTol, groupMeans.twoUnique_int_errorCode);
    groupMeans.twoUnique_int_means = clusterMeans;
    groupMeans.twoUnique_int_counts = clusterCounts;
    groupIndices.twoUnique_int_indices(:,3) = clusterInd;
    
    %-----twoUnique_pos
    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.twoUnique_pos_properties, groupMeans.twoUnique_pos_means,...
        maxIter, minTol, groupMeans.twoUnique_pos_errorCode);
    groupMeans.twoUnique_pos_means = clusterMeans;
    groupMeans.twoUnique_pos_counts = clusterCounts;
    groupIndices.twoUnique_pos_indices(:,3) = clusterInd;
    
    %-----twoUnique_neg
    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.twoUnique_neg_properties, groupMeans.twoUnique_neg_means,...
        maxIter, minTol, groupMeans.twoUnique_neg_errorCode);
    groupMeans.twoUnique_neg_means = clusterMeans;
    groupMeans.twoUnique_neg_counts = clusterCounts;
    groupIndices.twoUnique_neg_indices(:,3) = clusterInd;
    
    %-----twoUnique_ext
    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.twoUnique_ext_properties, groupMeans.twoUnique_ext_means,...
        maxIter, minTol, groupMeans.twoUnique_ext_errorCode);
    groupMeans.twoUnique_ext_means = clusterMeans;
    groupMeans.twoUnique_ext_counts = clusterCounts;
    groupIndices.twoUnique_ext_indices(:,3) = clusterInd;
    
    %------------3 unique------------
    %-----threeUnique_intInt
    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.threeUnique_intInt_properties, groupMeans.threeUnique_intInt_means,...
        maxIter, minTol, groupMeans.threeUnique_intInt_errorCode);  
    groupMeans.threeUnique_intInt_means = clusterMeans;
    groupMeans.threeUnique_intInt_counts = clusterCounts;
    groupIndices.threeUnique_intInt_indices(:,3) = clusterInd;
    
    %-----threeUnique_intPos
    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.threeUnique_intPos_properties, groupMeans.threeUnique_intPos_means,...
        maxIter, minTol, groupMeans.threeUnique_intPos_errorCode);
    groupMeans.threeUnique_intPos_means = clusterMeans;
    groupMeans.threeUnique_intPos_counts = clusterCounts;
    groupIndices.threeUnique_intPos_indices(:,3) = clusterInd;
    
    %-----threeUnique_intNeg
    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.threeUnique_intNeg_properties, groupMeans.threeUnique_intNeg_means,...
        maxIter, minTol, groupMeans.threeUnique_intNeg_errorCode);
    groupMeans.threeUnique_intNeg_means = clusterMeans;
    groupMeans.threeUnique_intNeg_counts = clusterCounts;
    groupIndices.threeUnique_intNeg_indices(:,3) = clusterInd;
    
    %-----threeUnique_posInt
    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.threeUnique_posInt_properties, groupMeans.threeUnique_posInt_means,...
        maxIter, minTol, groupMeans.threeUnique_posInt_errorCode);
    groupMeans.threeUnique_posInt_means = clusterMeans;
    groupMeans.threeUnique_posInt_counts = clusterCounts;
    groupIndices.threeUnique_posInt_indices(:,3) = clusterInd;
    
    %-----threeUnique_negInt
    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.threeUnique_negInt_properties, groupMeans.threeUnique_negInt_means,...
        maxIter, minTol, groupMeans.threeUnique_negInt_errorCode);
    groupMeans.threeUnique_negInt_means = clusterMeans;
    groupMeans.threeUnique_negInt_counts = clusterCounts;
    groupIndices.threeUnique_negInt_indices(:,3) = clusterInd;
    
    %-----threeUnique_extNotExt
    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.threeUnique_extNotExt_properties, groupMeans.threeUnique_extNotExt_means,...
        maxIter, minTol, groupMeans.threeUnique_extNotExt_errorCode);
    groupMeans.threeUnique_extNotExt_means = clusterMeans;
    groupMeans.threeUnique_extNotExt_counts = clusterCounts;
    groupIndices.threeUnique_extNotExt_indices(:,3) = clusterInd;
    
    %-----threeUnique_notExtExt
    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.threeUnique_notExtExt_properties, groupMeans.threeUnique_notExtExt_means,...
        maxIter, minTol, groupMeans.threeUnique_notExtExt_errorCode);
    groupMeans.threeUnique_notExtExt_means = clusterMeans;
    groupMeans.threeUnique_notExtExt_counts = clusterCounts;
    groupIndices.threeUnique_notExtExt_indices(:,3) = clusterInd;
    
    %-----threeUnique_parIntParInt
    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.threeUnique_parIntParInt_properties, groupMeans.threeUnique_parIntParInt_means,...
        maxIter, minTol, groupMeans.threeUnique_parIntParInt_errorCode);
    groupMeans.threeUnique_parIntParInt_means = clusterMeans;
    groupMeans.threeUnique_parIntParInt_counts = clusterCounts;
    groupIndices.threeUnique_parIntParInt_indices(:,3) = clusterInd;
    
    %-----threeUnique_extExt
    if (numel(groupIndices.threeUnique_extExt_indices(:,1)) > 0)
        [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
            calcClusterMeans(groupMeans.threeUnique_extExt_properties, groupMeans.threeUnique_extExt_means,...
            maxIter, minTol, groupMeans.threeUnique_extExt_errorCode);
        groupMeans.threeUnique_extExt_means = clusterMeans;
        groupMeans.threeUnique_extExt_counts = clusterCounts;
        groupIndices.threeUnique_extExt_indices(:,3) = clusterInd;
    end
    
    %------------4 unique------------
    %-----fourUnique_intInt

     [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
         calcClusterMeans(groupMeans.fourUnique_intInt_properties, groupMeans.fourUnique_intInt_means,...
         maxIter, minTol, groupMeans.fourUnique_intInt_errorCode);
    groupMeans.fourUnique_intInt_means = clusterMeans;
    groupMeans.fourUnique_intInt_counts = clusterCounts;
    groupIndices.fourUnique_intInt_indices(:,3) = clusterInd;
    
    %-----fourUnique_intPos
    
    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.fourUnique_intPos_properties, groupMeans.fourUnique_intPos_means,...
        maxIter, minTol, groupMeans.fourUnique_intPos_errorCode);
    groupMeans.fourUnique_intPos_means = clusterMeans;
    groupMeans.fourUnique_intPos_counts = clusterCounts;
    groupIndices.fourUnique_intPos_indices(:,3) = clusterInd;
    
    %-----fourUnique_intNeg
    
    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.fourUnique_intNeg_properties, groupMeans.fourUnique_intNeg_means,...
        maxIter, minTol, groupMeans.fourUnique_intNeg_errorCode);
    groupMeans.fourUnique_intNeg_means = clusterMeans;
    groupMeans.fourUnique_intNeg_counts = clusterCounts;
    groupIndices.fourUnique_intNeg_indices(:,3) = clusterInd;
    
    %-----fourUnique_posInt
    
    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.fourUnique_posInt_properties, groupMeans.fourUnique_posInt_means,...
        maxIter, minTol, groupMeans.fourUnique_posInt_errorCode);
    groupMeans.fourUnique_posInt_means = clusterMeans;
    groupMeans.fourUnique_posInt_counts = clusterCounts;
    groupIndices.fourUnique_posInt_indices(:,3) = clusterInd;
    
    %-----fourUnique_negInt

    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.fourUnique_negInt_properties, groupMeans.fourUnique_negInt_means,...
        maxIter, minTol, groupMeans.fourUnique_negInt_errorCode);
    groupMeans.fourUnique_negInt_means = clusterMeans;
    groupMeans.fourUnique_negInt_counts = clusterCounts;
    groupIndices.fourUnique_negInt_indices(:,3) = clusterInd;
    
    %-----fourUnique_posPos
    
    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.fourUnique_posPos_properties, groupMeans.fourUnique_posPos_means,...
        maxIter, minTol, groupMeans.fourUnique_posPos_errorCode);
    groupMeans.fourUnique_posPos_means = clusterMeans;
    groupMeans.fourUnique_posPos_counts = clusterCounts;
    groupIndices.fourUnique_posPos_indices(:,3) = clusterInd;
    
    %-----fourUnique_posNeg
        
    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.fourUnique_posNeg_properties, groupMeans.fourUnique_posNeg_means,...
        maxIter, minTol, groupMeans.fourUnique_posNeg_errorCode);
    groupMeans.fourUnique_posNeg_means = clusterMeans;
    groupMeans.fourUnique_posNeg_counts = clusterCounts;
    groupIndices.fourUnique_posNeg_indices(:,3) = clusterInd;
    
    %-----fourUnique_negPos
    
    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.fourUnique_negPos_properties, groupMeans.fourUnique_negPos_means,...
        maxIter, minTol, groupMeans.fourUnique_negPos_errorCode);
    groupMeans.fourUnique_negPos_means = clusterMeans;
    groupMeans.fourUnique_negPos_counts = clusterCounts;
    groupIndices.fourUnique_negPos_indices(:,3) = clusterInd;
    
    %-----fourUnique_negNeg
        
    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.fourUnique_negNeg_properties, groupMeans.fourUnique_negNeg_means,...
        maxIter, minTol, groupMeans.fourUnique_negNeg_errorCode);
    groupMeans.fourUnique_negNeg_means = clusterMeans;
    groupMeans.fourUnique_negNeg_counts = clusterCounts;
    groupIndices.fourUnique_negNeg_indices(:,3) = clusterInd;
    
    %-----fourUnique_notExtNotExt
    
    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.fourUnique_notExtNotExt_properties, groupMeans.fourUnique_notExtNotExt_means,...
        maxIter, minTol, groupMeans.fourUnique_notExtNotExt_errorCode);
    groupMeans.fourUnique_notExtNotExt_means = clusterMeans;
    groupMeans.fourUnique_notExtNotExt_counts = clusterCounts;
    groupIndices.fourUnique_notExtNotExt_indices(:,3) = clusterInd;
    
    %-----fourUnique_notExtExt
        
    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.fourUnique_notExtExt_properties, groupMeans.fourUnique_notExtExt_means,...
        maxIter, minTol, groupMeans.fourUnique_notExtExt_errorCode);
    groupMeans.fourUnique_notExtExt_means = clusterMeans;
    groupMeans.fourUnique_notExtExt_counts = clusterCounts;
    groupIndices.fourUnique_notExtExt_indices(:,3) = clusterInd;
    
    %-----fourUnique_extNotExt
            
    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.fourUnique_extNotExt_properties, groupMeans.fourUnique_extNotExt_means,...
        maxIter, minTol, groupMeans.fourUnique_extNotExt_errorCode);
    groupMeans.fourUnique_extNotExt_means = clusterMeans;
    groupMeans.fourUnique_extNotExt_counts = clusterCounts;
    groupIndices.fourUnique_extNotExt_indices(:,3) = clusterInd;
    
    %-----fourUnique_extExt
                
    [clusterMeans, clusterInd, clusterCounts,~, ~, ~, ~] =...
        calcClusterMeans(groupMeans.fourUnique_extExt_properties, groupMeans.fourUnique_extExt_means,...
        maxIter, minTol, groupMeans.fourUnique_extExt_errorCode);
    groupMeans.fourUnique_extExt_means = clusterMeans;
    groupMeans.fourUnique_extExt_counts = clusterCounts;
    groupIndices.fourUnique_extExt_indices(:,3) = clusterInd;
    
end