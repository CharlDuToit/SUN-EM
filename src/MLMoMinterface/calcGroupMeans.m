function [groupMeans, groupIndices] = calcGroupMeans(groupMeans, groupIndices, maxIter, minTol)
    
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

     [clusterMeans, clusterInd, clusterCounts,maxClusterError, totError, numIter, tol] =...
         calcClusterMeans(groupMeans.fourUnique_intInt_properties, groupMeans.fourUnique_intInt_means,...
         groupMeans.fourUnique_intInt_propScale, maxIter, minTol, groupMeans.fourUnique_intInt_errorCode);
    
    groupMeans.fourUnique_intInt_means = clusterMeans;
    groupMeans.fourUnique_intInt_counts = clusterCounts;
    groupIndices.fourUnique_intInt_indices(:,3) = clusterInd;
    
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