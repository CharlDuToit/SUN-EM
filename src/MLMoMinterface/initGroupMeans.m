function groupMeans = initGroupMeans(newLinkOld, newProperties, oldProperties,newEdgeLengths, oldEdgeLengths, groupIndices, minSize)
    %newProperties: [numNewEdges, numNewEdges, 3]

    groupMeans = [];
    [~,~,numNewProp] = size(newProperties);
    [~,~,numOldProp] = size(oldProperties);
    
    %------------ 2 unique--------------
    numEdgeIntervals = 1;
    clusterIntervals = cell(1,1);
    %----- twoUnique_int
    prop = newEdgeLengths(groupIndices.twoUnique_int_indices(:,1));
    clusterIntervals{1} = linspace( min(prop(:,1)) , max(prop(:,1)) , numEdgeIntervals + 1);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.twoUnique_int_means = clusterMeans;
    groupMeans.twoUnique_int_properties = prop;
    groupMeans.twoUnique_int_numClass = numClass;
    groupMeans.twoUnique_int_numNoClass = numNoClass;
    groupMeans.twoUnique_int_counts = 0;
    
    %-----twoUnique_pos
    prop = newEdgeLengths(groupIndices.twoUnique_pos_indices(:,1));
    clusterIntervals{1} = linspace( min(prop(:,1)) , max(prop(:,1)) , numEdgeIntervals + 1);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.twoUnique_pos_means = clusterMeans;
    groupMeans.twoUnique_pos_properties = prop;
    groupMeans.twoUnique_pos_numClass = numClass;
    groupMeans.twoUnique_pos_numNoClass = numNoClass;
    groupMeans.twoUnique_pos_counts = 0;
    
    %-----twoUnique_neg
    prop = newEdgeLengths(groupIndices.twoUnique_neg_indices(:,1));
    clusterIntervals{1} = linspace( min(prop(:,1)) , max(prop(:,1)) , numEdgeIntervals + 1);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.twoUnique_neg_means = clusterMeans;
    groupMeans.twoUnique_neg_properties = prop;
    groupMeans.twoUnique_neg_numClass = numClass;
    groupMeans.twoUnique_neg_numNoClass = numNoClass;
    groupMeans.twoUnique_neg_counts = 0;
    
    %-----twoUnique_ext
    prop = newEdgeLengths(groupIndices.twoUnique_ext_indices(:,1));
    clusterIntervals{1} = linspace( min(prop(:,1)) , max(prop(:,1)) , numEdgeIntervals + 1);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.twoUnique_ext_means = clusterMeans;
    groupMeans.twoUnique_ext_properties = prop;
    groupMeans.twoUnique_ext_numClass = numClass;
    groupMeans.twoUnique_ext_numNoClass = numNoClass;
    groupMeans.twoUnique_ext_counts = 0;
    
    %------------3 unique--------------
    %clusterIntervals = cell(3,1);
    
    %numDistIntervals = 1;
    %numOrientationIntervals = 1;
    %numRadialIntervals = 1;
    
    %-----threeUnique_intInt
    clusterSizes = [1 1 1];
    propScale = [1 1 1]; % used by calcGroupMeans
    errorCode = 1; % used by calcGroupMeans
    
    mmInd = groupIndices.threeUnique_intInt_indices(:,1);
    nnInd = groupIndices.threeUnique_intInt_indices(:,2);
    prop = zeros(numel(mmInd), numNewProp);
    for k = 1:numel(mmInd)
        prop(k,:) = newProperties(mmInd(k), nnInd(k), :);
    end
    clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.threeUnique_intInt_means = clusterMeans;
    groupMeans.threeUnique_intInt_properties = prop;
    groupMeans.threeUnique_intInt_propScale = propScale;
    groupMeans.threeUnique_intInt_errorCode = errorCode;
    groupMeans.threeUnique_intInt_numClass = numClass;
    groupMeans.threeUnique_intInt_numNoClass = numNoClass;
    groupMeans.threeUnique_intInt_counts = 0;
    
    %-----threeUnique_intPos
    
    %-----threeUnique_intNeg
    
    %-----threeUnique_posInt
    
    %-----threeUnique_negInt
    
    %-----threeUnique_extNotExt
    
    %-----threeUnique_notExtExt
    
    %-----threeUnique_parIntParInt
    
    %-----threeUnique_extExt
    
    %------------4 unique--------------
    %clusterIntervals = cell(3,1);
    %numDistIntervals = 2;
    %numOrientationIntervals = 1;
    %numRadialIntervals = 1;
    
    %-----fourUnique_intInt
    % no need to cluster by orientation
    clusterSizes = [2 1 2 1];
    propScale = [3 1 3 1];
    errorCode = 2;
    
    mmInd = groupIndices.fourUnique_intInt_indices(:,1);
    nnInd = groupIndices.fourUnique_intInt_indices(:,2);
    prop = zeros(numel(mmInd), numNewProp-1 + numOldProp-1);
    for k = 1:numel(mmInd)
        prop(k,1:2) = newProperties(mmInd(k), nnInd(k), [1 3]);
        prop(k,3:4) = oldProperties(newLinkOld(mmInd(k), nnInd(k), 1), newLinkOld(mmInd(k), nnInd(k), 2), [1 3]);
    end
    clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.fourUnique_intInt_means = clusterMeans;
    groupMeans.fourUnique_intInt_properties = prop;
    groupMeans.fourUnique_intInt_propScale = propScale;
    groupMeans.fourUnique_intInt_errorCode = errorCode;
    groupMeans.fourUnique_intInt_numClass = numClass;
    groupMeans.fourUnique_intInt_numNoClass = numNoClass;
    groupMeans.fourUnique_intInt_counts = 0;
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

function clusterIntervals = createIntervals(prop, sizes)
    clusterIntervals = cell(numel(sizes),1);
    for k = 1:numel(sizes)
        clusterIntervals{k} = linspace( min(prop(:,k)) , max(prop(:,k)) , sizes(k) + 1);
    end
end