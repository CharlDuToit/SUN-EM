function groupIndices = assignGroupClusters(groupMeans, groupIndices,newLinkOld, newProperties, oldProperties,newEdgeLengths, oldEdgeLengths, newRhoProperties, oldRhoProperties)
    %groupMeans from trained model
    %groupIndices from different solver setup
    
    %------------2 unique------------
    %----- twoUnique_int
    prop = createProp(groupIndices.twoUnique_int_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.twoUnique_int_errorCode);
    [~, clusterInd, ~,~, ~] = assignClusters(prop, groupMeans.twoUnique_int_means,...
        groupMeans.twoUnique_int_errorCode);
    groupIndices.twoUnique_int_indices(:,3) = clusterInd;
    
    %-----twoUnique_pos
    prop = createProp(groupIndices.twoUnique_pos_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.twoUnique_pos_errorCode);
    [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.twoUnique_pos_means,...
        groupMeans.twoUnique_pos_errorCode);
    groupIndices.twoUnique_pos_indices(:,3) = clusterInd;
    
    %-----twoUnique_neg
     prop = createProp(groupIndices.twoUnique_neg_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.twoUnique_neg_errorCode);
    [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.twoUnique_neg_means,...
        groupMeans.twoUnique_neg_errorCode);
    groupIndices.twoUnique_neg_indices(:,3) = clusterInd;
    
    %-----twoUnique_ext
     prop = createProp(groupIndices.twoUnique_ext_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.twoUnique_ext_errorCode);
    [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.twoUnique_ext_means,...
        groupMeans.twoUnique_ext_errorCode);
    groupIndices.twoUnique_ext_indices(:,3) = clusterInd;
    
    %------------3 unique------------
    %-----threeUnique_intInt
     prop = createProp(groupIndices.threeUnique_intInt_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.threeUnique_intInt_errorCode);
    [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.threeUnique_intInt_means,...
        groupMeans.threeUnique_intInt_errorCode);  
    groupIndices.threeUnique_intInt_indices(:,3) = clusterInd;
    
    %-----threeUnique_intPos
    prop = createProp(groupIndices.threeUnique_intPos_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.threeUnique_intPos_errorCode);
    [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.threeUnique_intPos_means,...
        groupMeans.threeUnique_intPos_errorCode);
    groupIndices.threeUnique_intPos_indices(:,3) = clusterInd;
    
    %-----threeUnique_intNeg
    prop = createProp(groupIndices.threeUnique_intNeg_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.threeUnique_intNeg_errorCode);
    [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.threeUnique_intNeg_means,...
        groupMeans.threeUnique_intNeg_errorCode);
    groupIndices.threeUnique_intNeg_indices(:,3) = clusterInd;
    
    %-----threeUnique_posInt
    prop = createProp(groupIndices.threeUnique_posInt_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.threeUnique_posInt_errorCode);
    [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.threeUnique_posInt_means,...
        groupMeans.threeUnique_posInt_errorCode);
    groupIndices.threeUnique_posInt_indices(:,3) = clusterInd;
    
    %-----threeUnique_negInt
    prop = createProp(groupIndices.threeUnique_negInt_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.threeUnique_negInt_errorCode);
    [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.threeUnique_negInt_means,...
        groupMeans.threeUnique_negInt_errorCode);
    groupIndices.threeUnique_negInt_indices(:,3) = clusterInd;
    
    %-----threeUnique_extNotExt
    prop = createProp(groupIndices.threeUnique_extNotExt_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.threeUnique_extNotExt_errorCode);
    [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.threeUnique_extNotExt_means,...
        groupMeans.threeUnique_extNotExt_errorCode);
    groupIndices.threeUnique_extNotExt_indices(:,3) = clusterInd;
    
    %-----threeUnique_notExtExt
    prop = createProp(groupIndices.threeUnique_notExtExt_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.threeUnique_notExtExt_errorCode);
    [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.threeUnique_notExtExt_means,...
        groupMeans.threeUnique_notExtExt_errorCode);
    groupIndices.threeUnique_notExtExt_indices(:,3) = clusterInd;
    
    %-----threeUnique_parIntParInt
    prop = createProp(groupIndices.threeUnique_parIntParInt_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.threeUnique_parIntParInt_errorCode);
    [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.threeUnique_parIntParInt_means,...
        groupMeans.threeUnique_parIntParInt_errorCode);
    groupIndices.threeUnique_parIntParInt_indices(:,3) = clusterInd;
    
    %-----threeUnique_extExt
    prop = createProp(groupIndices.threeUnique_extExt_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.threeUnique_extExt_errorCode);
    if (numel(groupIndices.threeUnique_extExt_indices(:,1)) > 0)
        [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.threeUnique_extExt_means,...
            groupMeans.threeUnique_extExt_errorCode);
        groupIndices.threeUnique_extExt_indices(:,3) = clusterInd;
    end
    
    %------------4 unique------------
    %-----fourUnique_intInt
    prop = createProp(groupIndices.fourUnique_intInt_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.fourUnique_intInt_errorCode);
     [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.fourUnique_intInt_means,...
         groupMeans.fourUnique_intInt_errorCode);
    groupIndices.fourUnique_intInt_indices(:,3) = clusterInd;
    
    %-----fourUnique_intPos
    
    prop = createProp(groupIndices.fourUnique_intPos_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.fourUnique_intPos_errorCode);
    [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.fourUnique_intPos_means,...
        groupMeans.fourUnique_intPos_errorCode);
    groupIndices.fourUnique_intPos_indices(:,3) = clusterInd;
    
    %-----fourUnique_intNeg
    
    prop = createProp(groupIndices.fourUnique_intNeg_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.fourUnique_intNeg_errorCode);
    [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.fourUnique_intNeg_means,...
        groupMeans.fourUnique_intNeg_errorCode);
    groupIndices.fourUnique_intNeg_indices(:,3) = clusterInd;
    
    %-----fourUnique_posInt
    prop = createProp(groupIndices.fourUnique_posInt_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.fourUnique_posInt_errorCode);
    [~, clusterInd, ~,~, ~] = assignClusters(prop, groupMeans.fourUnique_posInt_means,...
        groupMeans.fourUnique_posInt_errorCode);
    groupIndices.fourUnique_posInt_indices(:,3) = clusterInd;
    
    %-----fourUnique_negInt
    prop = createProp(groupIndices.fourUnique_negInt_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.fourUnique_negInt_errorCode);
    [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.fourUnique_negInt_means,...
        groupMeans.fourUnique_negInt_errorCode);
    groupIndices.fourUnique_negInt_indices(:,3) = clusterInd;
    
    %-----fourUnique_posPos
    prop = createProp(groupIndices.fourUnique_posPos_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.fourUnique_posPos_errorCode);
    [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.fourUnique_posPos_means,...
        groupMeans.fourUnique_posPos_errorCode);
    groupIndices.fourUnique_posPos_indices(:,3) = clusterInd;
    
    %-----fourUnique_posNeg
    prop = createProp(groupIndices.fourUnique_posNeg_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.fourUnique_posNeg_errorCode);
    [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.fourUnique_posNeg_means,...
        groupMeans.fourUnique_posNeg_errorCode);
    groupIndices.fourUnique_posNeg_indices(:,3) = clusterInd;
    
    %-----fourUnique_negPos
    prop = createProp(groupIndices.fourUnique_negPos_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.fourUnique_negPos_errorCode);
    [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.fourUnique_negPos_means,...
        groupMeans.fourUnique_negPos_errorCode);
    groupIndices.fourUnique_negPos_indices(:,3) = clusterInd;
    
    %-----fourUnique_negNeg
     prop = createProp(groupIndices.fourUnique_negNeg_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.fourUnique_negNeg_errorCode);
    [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.fourUnique_negNeg_means,...
        groupMeans.fourUnique_negNeg_errorCode);
    groupIndices.fourUnique_negNeg_indices(:,3) = clusterInd;
    
    %-----fourUnique_notExtNotExt
    prop = createProp(groupIndices.fourUnique_notExtNotExt_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.fourUnique_notExtNotExt_errorCode);
    [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.fourUnique_notExtNotExt_means,...
        groupMeans.fourUnique_notExtNotExt_errorCode);
    groupIndices.fourUnique_notExtNotExt_indices(:,3) = clusterInd;
    
    %-----fourUnique_notExtExt
    prop = createProp(groupIndices.fourUnique_notExtExt_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.fourUnique_notExtExt_errorCode);    
    [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.fourUnique_notExtExt_means,...
        groupMeans.fourUnique_notExtExt_errorCode);
    groupIndices.fourUnique_notExtExt_indices(:,3) = clusterInd;
    
    %-----fourUnique_extNotExt
    prop = createProp(groupIndices.fourUnique_extNotExt_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.fourUnique_extNotExt_errorCode); 
    [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.fourUnique_extNotExt_means,...
        groupMeans.fourUnique_extNotExt_errorCode);
    groupIndices.fourUnique_extNotExt_indices(:,3) = clusterInd;
    
    %-----fourUnique_extExt
    prop = createProp(groupIndices.fourUnique_extExt_indices , newLinkOld, newProperties,...
        oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, groupMeans.fourUnique_extExt_errorCode); 
    [~, clusterInd, ~,~, ~] =assignClusters(prop, groupMeans.fourUnique_extExt_means,...
        groupMeans.fourUnique_extExt_errorCode);
    groupIndices.fourUnique_extExt_indices(:,3) = clusterInd;
    
end

function prop = createProp(indices, newLinkOld, newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties, oldRhoProperties, errorCode)
    mmInd = indices(:,1);
    nnInd = indices(:,2);
    switch errorCode
        case 1 %3 properties, new dist, new orientation, new radial
            prop = zeros(numel(mmInd), 3);
            for k = 1:numel(mmInd)
                prop(k,:) = newProperties(mmInd(k), nnInd(k), :);
            end
        case 2 %4 properties, new dist, new radial, old dist, old radial
            %fourUnique_posInt, intInt
            prop = zeros(numel(mmInd), 4);
            for k = 1:numel(mmInd)
                prop(k,1:2) = newProperties(mmInd(k), nnInd(k), [1 3]);
                prop(k,3:4) = oldProperties(newLinkOld(mmInd(k), nnInd(k), 1), newLinkOld(mmInd(k), nnInd(k), 2), [1 3]);
            end
        case 3 %3 properties new dist, new radial, old radial
            %threeUnique posInt, intInt
            prop = zeros(numel(mmInd), 3);
            for k = 1:numel(mmInd)
                prop(k,1:2) = newProperties(mmInd(k), nnInd(k), [1 3]);
                prop(k,3) = oldProperties(newLinkOld(mmInd(k), nnInd(k), 1), newLinkOld(mmInd(k), nnInd(k), 2), 3);
            end
        case 4 %6 properties, 3 new, 3 old
            %same as 7
            %fourUnique_extNotExt, notExtnotExt
            prop = zeros(numel(mmInd), 6);
            for k = 1:numel(mmInd)
                prop(k,1:3) = newProperties(mmInd(k), nnInd(k), :);
                prop(k,4:6) = oldProperties(newLinkOld(mmInd(k), nnInd(k), 1), newLinkOld(mmInd(k), nnInd(k), 2), :);
            end
        case 5 %1 property, new edge length
            %twoUnique int
            %numClusters = numel(mmInd)./avgSize;
            prop = newEdgeLengths(mmInd);
        case 6 %4 properties, new dist, new radial, old dist, old radial
            %threeUnique_parIntParInt, extNotExt, extExt
            prop = zeros(numel(mmInd), 4);
            for k = 1:numel(mmInd)
                prop(k,1:2) = newProperties(mmInd(k), nnInd(k), [1 3]);
                prop(k,3:4) = oldProperties(newLinkOld(mmInd(k), nnInd(k), 1), newLinkOld(mmInd(k), nnInd(k), 2), [1 3]);
            end
         case 7 %6 properties, 3 new, 3 old
            %fourUnique_extNotExt, notExtnotExt
            prop = zeros(numel(mmInd), 6);
            for k = 1:numel(mmInd)
                prop(k,1:3) = newProperties(mmInd(k), nnInd(k), :);
                prop(k,4:6) = oldProperties(newLinkOld(mmInd(k), nnInd(k), 1), newLinkOld(mmInd(k), nnInd(k), 2), :);
            end
         case 8 %1 property, new edge length,new posneg , old posneg
            %twoUnique pos neg int
            prop = zeros(numel(mmInd), 3);
            for k = 1:numel(mmInd)
                prop(k,1) = newEdgeLengths(k);
                prop(k,2) = newRhoProperties(mmInd(k), nnInd(k), 2);
                prop(k,3) = oldRhoProperties(newLinkOld(mmInd(k), nnInd(k), 1), newLinkOld(mmInd(k), nnInd(k), 2), 2);
            end
         case 9%3 property, new edge length,new posneg , old posneg, old pospos
            %twoUnique pos neg ext
            prop = zeros(numel(mmInd), 5);
            for k = 1:numel(mmInd)
                prop(k,1) = newEdgeLengths(k);
                prop(k,2) = newRhoProperties(mmInd(k), nnInd(k), 1);
                prop(k,3) = newRhoProperties(mmInd(k), nnInd(k), 2);
                prop(k,4) = oldRhoProperties(newLinkOld(mmInd(k), nnInd(k), 1), newLinkOld(mmInd(k), nnInd(k), 2), 1);
                prop(k,5) = oldRhoProperties(newLinkOld(mmInd(k), nnInd(k), 1), newLinkOld(mmInd(k), nnInd(k), 2), 2);
            end
    end
    
end