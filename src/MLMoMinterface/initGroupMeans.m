function groupMeans = initGroupMeans(groupIndices, newLinkOld, newProperties, oldProperties,newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties, minSize)
    %newProperties: [numNewEdges, numNewEdges, 3]

    groupMeans = [];
    numNewEdges = numel(newEdgeLengths(:, 1));
    sizeConst = 60; %50, 
    avgSize = numNewEdges.^2 ./(sizeConst * log(numNewEdges).^2);
    
    %avgEdgeLength = sum(newEdgeLengths)/numel(newEdgeLengths);
    %maxDist = max(newProperties(:,:,1));
    %minDist = min(newProperties(:,:,1));
    %threshDist = 0.03584296*maxDist + 1.769794721407624*avgEdgeLength;
    %------------ 2 unique--------------
    %----- twoUnique_int
    errorCode = 5;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.twoUnique_int_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
    %clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.twoUnique_int_means = clusterMeans;
    groupMeans.twoUnique_int_properties = prop;
    groupMeans.twoUnique_int_errorCode = errorCode;
%     groupMeans.twoUnique_int_numClass = numClass;
%     groupMeans.twoUnique_int_numNoClass = numNoClass;
    groupMeans.twoUnique_int_counts = 0;
    
    %-----twoUnique_pos
    errorCode = 9;
    minSizeTwoUniqueNegOrPosOrExt = 4;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.twoUnique_pos_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
    %clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSizeTwoUniqueNegOrPosOrExt);
    groupMeans.twoUnique_pos_means = clusterMeans;
    groupMeans.twoUnique_pos_properties = prop;
    groupMeans.twoUnique_pos_errorCode = errorCode;
%     groupMeans.twoUnique_pos_numClass = numClass;
%     groupMeans.twoUnique_pos_numNoClass = numNoClass;
    groupMeans.twoUnique_pos_counts = 0;
    
    %-----twoUnique_neg
    errorCode = 9;
    
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.twoUnique_neg_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
    %clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSizeTwoUniqueNegOrPosOrExt);
    groupMeans.twoUnique_neg_means = clusterMeans;
    groupMeans.twoUnique_neg_properties = prop;
    groupMeans.twoUnique_neg_errorCode = errorCode;
%     groupMeans.twoUnique_neg_numClass = numClass;
%     groupMeans.twoUnique_neg_numNoClass = numNoClass;
    groupMeans.twoUnique_neg_counts = 0;
    
    %-----twoUnique_ext
    errorCode = 9;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.twoUnique_ext_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
    %clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSizeTwoUniqueNegOrPosOrExt);
    groupMeans.twoUnique_ext_means = clusterMeans;
    groupMeans.twoUnique_ext_properties = prop;
    groupMeans.twoUnique_ext_errorCode = errorCode;
%     groupMeans.twoUnique_ext_numClass = numClass;
%     groupMeans.twoUnique_ext_numNoClass = numNoClass;
    groupMeans.twoUnique_ext_counts = 0;
    
    %------------3 unique--------------
    %clusterIntervals = cell(3,1);
    
    %numDistIntervals = 1;
    %numOrientationIntervals = 1;
    %numRadialIntervals = 1;
    
    %-----threeUnique_intInt
    errorCode = 3;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.threeUnique_intInt_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
    %clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSizeTwoUniqueNegOrPosOrExt);
    groupMeans.threeUnique_intInt_means = clusterMeans;
    groupMeans.threeUnique_intInt_properties = prop;
    groupMeans.threeUnique_intInt_errorCode = errorCode;
%     groupMeans.threeUnique_intInt_numClass = numClass;
%     groupMeans.threeUnique_intInt_numNoClass = numNoClass;
    groupMeans.threeUnique_intInt_counts = 0;
    
    %-----threeUnique_intPos
    errorCode = 3;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.threeUnique_intPos_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
    %clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.threeUnique_intPos_means = clusterMeans;
    groupMeans.threeUnique_intPos_properties = prop;
    groupMeans.threeUnique_intPos_errorCode = errorCode;
%     groupMeans.threeUnique_intPos_numClass = numClass;
%     groupMeans.threeUnique_intPos_numNoClass = numNoClass;
    groupMeans.threeUnique_intPos_counts = 0;
    
    %-----threeUnique_intNeg
    errorCode = 3;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.threeUnique_intNeg_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
   % clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.threeUnique_intNeg_means = clusterMeans;
    groupMeans.threeUnique_intNeg_properties = prop;
    groupMeans.threeUnique_intNeg_errorCode = errorCode;
%     groupMeans.threeUnique_intNeg_numClass = numClass;
%     groupMeans.threeUnique_intNeg_numNoClass = numNoClass;
    groupMeans.threeUnique_intNeg_counts = 0;
    
    %-----threeUnique_posInt
    errorCode = 3;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.threeUnique_posInt_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
    %clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.threeUnique_posInt_means = clusterMeans;
    groupMeans.threeUnique_posInt_properties = prop;
    groupMeans.threeUnique_posInt_errorCode = errorCode;
%     groupMeans.threeUnique_posInt_numClass = numClass;
%     groupMeans.threeUnique_posInt_numNoClass = numNoClass;
    groupMeans.threeUnique_posInt_counts = 0;
    
    %-----threeUnique_negInt
    errorCode = 3;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.threeUnique_negInt_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
    %clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.threeUnique_negInt_means = clusterMeans;
    groupMeans.threeUnique_negInt_properties = prop;
    groupMeans.threeUnique_negInt_errorCode = errorCode;
%     groupMeans.threeUnique_negInt_numClass = numClass;
%     groupMeans.threeUnique_negInt_numNoClass = numNoClass;
    groupMeans.threeUnique_negInt_counts = 0;
    
    %-----threeUnique_extNotExt
    errorCode = 6;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.threeUnique_extNotExt_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
    %clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.threeUnique_extNotExt_means = clusterMeans;
    groupMeans.threeUnique_extNotExt_properties = prop;
    groupMeans.threeUnique_extNotExt_errorCode = errorCode;
%     groupMeans.threeUnique_extNotExt_numClass = numClass;
%     groupMeans.threeUnique_extNotExt_numNoClass = numNoClass;
    groupMeans.threeUnique_extNotExt_counts = 0;
    
    %-----threeUnique_notExtExt
    errorCode = 6;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.threeUnique_notExtExt_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
    %clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.threeUnique_notExtExt_means = clusterMeans;
    groupMeans.threeUnique_notExtExt_properties = prop;
    groupMeans.threeUnique_notExtExt_errorCode = errorCode;
%     groupMeans.threeUnique_notExtExt_numClass = numClass;
%     groupMeans.threeUnique_notExtExt_numNoClass = numNoClass;
    groupMeans.threeUnique_notExtExt_counts = 0;
    
    %-----threeUnique_parIntParInt
    errorCode = 6;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.threeUnique_parIntParInt_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
    %clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.threeUnique_parIntParInt_means = clusterMeans;
    groupMeans.threeUnique_parIntParInt_properties = prop;
    groupMeans.threeUnique_parIntParInt_errorCode = errorCode;
%     groupMeans.threeUnique_parIntParInt_numClass = numClass;
%     groupMeans.threeUnique_parIntParInt_numNoClass = numNoClass;
    groupMeans.threeUnique_parIntParInt_counts = 0;
    
    %-----threeUnique_extExt
    errorCode = 6;
    if (numel(groupIndices.threeUnique_extExt_indices(:,1)) > 0)
        [prop, clusterIntervals] =createPropAndIntervals(groupIndices.threeUnique_extExt_indices, newLinkOld,...
            newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
        %clusterIntervals = createIntervals(prop, clusterSizes);
        [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    else
        clusterMeans = [];
        prop = [];
        numClass = 0;
        numNoClass = 0;
    end

    groupMeans.threeUnique_extExt_means = clusterMeans;
    groupMeans.threeUnique_extExt_properties = prop;
    groupMeans.threeUnique_extExt_errorCode = errorCode;
%     groupMeans.threeUnique_extExt_numClass = numClass;
%     groupMeans.threeUnique_extExt_numNoClass = numNoClass;
    groupMeans.threeUnique_extExt_counts = 0;
    
    %------------4 unique--------------
    %clusterIntervals = cell(3,1);
    %numDistIntervals = 2;
    %numOrientationIntervals = 1;
    %numRadialIntervals = 1;
    
    %-----fourUnique_intInt
    errorCode = 2;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.fourUnique_intInt_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
    %clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.fourUnique_intInt_means = clusterMeans;
    groupMeans.fourUnique_intInt_properties = prop;
    groupMeans.fourUnique_intInt_errorCode = errorCode;
%     groupMeans.fourUnique_intInt_numClass = numClass;
%     groupMeans.fourUnique_intInt_numNoClass = numNoClass;
    groupMeans.fourUnique_intInt_counts = 0;
    
    %-----fourUnique_intPos
    errorCode = 2;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.fourUnique_intPos_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
   % clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.fourUnique_intPos_means = clusterMeans;
    groupMeans.fourUnique_intPos_properties = prop;
    groupMeans.fourUnique_intPos_errorCode = errorCode;
%     groupMeans.fourUnique_intPos_numClass = numClass;
%     groupMeans.fourUnique_intPos_numNoClass = numNoClass;
    groupMeans.fourUnique_intPos_counts = 0;
    
    %-----fourUnique_intNeg
    errorCode = 2;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.fourUnique_intNeg_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
    %clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.fourUnique_intNeg_means = clusterMeans;
    groupMeans.fourUnique_intNeg_properties = prop;
    groupMeans.fourUnique_intNeg_errorCode = errorCode;
%     groupMeans.fourUnique_intNeg_numClass = numClass;
%     groupMeans.fourUnique_intNeg_numNoClass = numNoClass;
    groupMeans.fourUnique_intNeg_counts = 0;
    
    %-----fourUnique_posInt
    errorCode = 2;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.fourUnique_posInt_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
    %clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.fourUnique_posInt_means = clusterMeans;
    groupMeans.fourUnique_posInt_properties = prop;
    groupMeans.fourUnique_posInt_errorCode = errorCode;
%     groupMeans.fourUnique_posInt_numClass = numClass;
%     groupMeans.fourUnique_posInt_numNoClass = numNoClass;
    groupMeans.fourUnique_posInt_counts = 0;
    
    %-----fourUnique_negInt
    errorCode = 2;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.fourUnique_negInt_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
    %clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.fourUnique_negInt_means = clusterMeans;
    groupMeans.fourUnique_negInt_properties = prop;
    groupMeans.fourUnique_negInt_errorCode = errorCode;
%     groupMeans.fourUnique_negInt_numClass = numClass;
%     groupMeans.fourUnique_negInt_numNoClass = numNoClass;
    groupMeans.fourUnique_negInt_counts = 0;
    
    %-----fourUnique_posPos
    errorCode = 2;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.fourUnique_posPos_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
    %clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.fourUnique_posPos_means = clusterMeans;
    groupMeans.fourUnique_posPos_properties = prop;
    groupMeans.fourUnique_posPos_errorCode = errorCode;
%     groupMeans.fourUnique_posPos_numClass = numClass;
%     groupMeans.fourUnique_posPos_numNoClass = numNoClass;
    groupMeans.fourUnique_posPos_counts = 0;
    
    %-----fourUnique_posNeg
    errorCode = 2;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.fourUnique_posNeg_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
%    clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.fourUnique_posNeg_means = clusterMeans;
    groupMeans.fourUnique_posNeg_properties = prop;
    groupMeans.fourUnique_posNeg_errorCode = errorCode;
%     groupMeans.fourUnique_posNeg_numClass = numClass;
%     groupMeans.fourUnique_posNeg_numNoClass = numNoClass;
    groupMeans.fourUnique_posNeg_counts = 0;
    
    %-----fourUnique_negPos
    errorCode = 2;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.fourUnique_negPos_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
    %clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.fourUnique_negPos_means = clusterMeans;
    groupMeans.fourUnique_negPos_properties = prop;
    groupMeans.fourUnique_negPos_errorCode = errorCode;
%     groupMeans.fourUnique_negPos_numClass = numClass;
%     groupMeans.fourUnique_negPos_numNoClass = numNoClass;
    groupMeans.fourUnique_negPos_counts = 0;
    
    %-----fourUnique_negNeg
    errorCode = 2;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.fourUnique_negNeg_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
    %clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.fourUnique_negNeg_means = clusterMeans;
    groupMeans.fourUnique_negNeg_properties = prop;
    groupMeans.fourUnique_negNeg_errorCode = errorCode;
%     groupMeans.fourUnique_negNeg_numClass = numClass;
%     groupMeans.fourUnique_negNeg_numNoClass = numNoClass;
    groupMeans.fourUnique_negNeg_counts = 0;
    
    %-----fourUnique_notExtNotExt
    errorCode = 7;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.fourUnique_notExtNotExt_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
    %clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.fourUnique_notExtNotExt_means = clusterMeans;
    groupMeans.fourUnique_notExtNotExt_properties = prop;
    groupMeans.fourUnique_notExtNotExt_errorCode = errorCode;
%     groupMeans.fourUnique_notExtNotExt_numClass = numClass;
%     groupMeans.fourUnique_notExtNotExt_numNoClass = numNoClass;
    groupMeans.fourUnique_notExtNotExt_counts = 0;
    
    %-----fourUnique_notExtExt
    errorCode = 4;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.fourUnique_notExtExt_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
    %clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.fourUnique_notExtExt_means = clusterMeans;
    groupMeans.fourUnique_notExtExt_properties = prop;
    groupMeans.fourUnique_notExtExt_errorCode = errorCode;
%     groupMeans.fourUnique_notExtExt_numClass = numClass;
%     groupMeans.fourUnique_notExtExt_numNoClass = numNoClass;
    groupMeans.fourUnique_notExtExt_counts = 0;
    
    %-----fourUnique_extNotExt
    errorCode = 4;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.fourUnique_extNotExt_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
    %clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.fourUnique_extNotExt_means = clusterMeans;
    groupMeans.fourUnique_extNotExt_properties = prop;
    groupMeans.fourUnique_extNotExt_errorCode = errorCode;
%     groupMeans.fourUnique_extNotExt_numClass = numClass;
%     groupMeans.fourUnique_extNotExt_numNoClass = numNoClass;
    groupMeans.fourUnique_extNotExt_counts = 0;
    
    %-----fourUnique_extExt
    errorCode = 4;
    [prop, clusterIntervals] =createPropAndIntervals(groupIndices.fourUnique_extExt_indices, newLinkOld,...
        newProperties, oldProperties, newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode);
    %clusterIntervals = createIntervals(prop, clusterSizes);
    [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(prop,clusterIntervals,  minSize);
    groupMeans.fourUnique_extExt_means = clusterMeans;
    groupMeans.fourUnique_extExt_properties = prop;
    groupMeans.fourUnique_extExt_errorCode = errorCode;
%     groupMeans.fourUnique_extExt_numClass = numClass;
%     groupMeans.fourUnique_extExt_numNoClass = numNoClass;
    groupMeans.fourUnique_extExt_counts = 0;

end

function [prop, clusterIntervals] = createPropAndIntervals(indices, newLinkOld, newProperties, oldProperties,...
    newEdgeLengths, oldEdgeLengths,newRhoProperties,oldRhoProperties,avgSize, errorCode)

    mmInd = indices(:,1);
    nnInd = indices(:,2);
    
    newAvgEdgeLength = sum(newEdgeLengths)/numel(newEdgeLengths);
    newMaxDist = max(newProperties(:,:,1));
    newMaxDist = max(newMaxDist);
    %newMinDist = min(newProperties(:,:,1));
    newThreshDist = 0.03584296*newMaxDist + 1.769794721407624*newAvgEdgeLength;
    
    oldAvgEdgeLength = sum(oldEdgeLengths)/numel(oldEdgeLengths);
    oldMaxDist = max(oldProperties(:,:,1));
    oldMaxDist = max(oldMaxDist);
    %oldMinDist = min(oldProperties(:,:,1));
    oldThreshDist = 0.03584296*oldMaxDist + 1.769794721407624*oldAvgEdgeLength;
    
    threshPropInd = [];
    threshDist = [];
    switch errorCode
        case 1 %3 properties, new dist, new orientation, new radial
           
            prop = zeros(numel(mmInd), 3);
            numClusters = numel(mmInd)./avgSize;
            s = ceil(sqrt(numClusters));
            sizes = [s s 1];
            %sizes = [18 2 11 11]; orientation, radial, preThresh, post Thresh
            %sizes = [2 1 1 2 1 1] * numel(mmInd)/ 100;
            for k = 1:numel(mmInd)
                prop(k,:) = newProperties(mmInd(k), nnInd(k), :);
            end
        case 2 %4 properties, new dist, new radial, old dist, old radial
            %fourUnique_posInt, intInt
            prop = zeros(numel(mmInd), 4);
            numClusters = 2*numel(mmInd)./avgSize; % 2
            s = ceil(sqrt(numClusters));
            sizes = [s 1 s 1];
            for k = 1:numel(mmInd)
                prop(k,1:2) = newProperties(mmInd(k), nnInd(k), [1 3]);
                prop(k,3:4) = oldProperties(newLinkOld(mmInd(k), nnInd(k), 1), newLinkOld(mmInd(k), nnInd(k), 2), [1 3]);
            end
            threshDist = [newThreshDist, oldThreshDist];
            threshPropInd = [1 3];
        case 3 %3 properties new dist, new radial, old radial
            %threeUnique posInt, intInt
            prop = zeros(numel(mmInd), 3);
            %numClusters = numel(mmInd)./avgSize
            numClusters = numel(mmInd)./sqrt(avgSize);
            s = ceil(numClusters);
            sizes = [s 1 1];
            for k = 1:numel(mmInd)
                prop(k,1:2) = newProperties(mmInd(k), nnInd(k), [1 3]);
                prop(k,3) = oldProperties(newLinkOld(mmInd(k), nnInd(k), 1), newLinkOld(mmInd(k), nnInd(k), 2), 3);
            end
        case 4 %6 properties, 3 new, 3 old
            %fourUnique_extNotExt, notExtExt
            prop = zeros(numel(mmInd), 6);
            numClusters = 3*numel(mmInd)./avgSize; % 0.5
            s = ceil(numClusters.^0.25);
            sizes = [s s 1 s s 1];
            for k = 1:numel(mmInd)
                prop(k,1:3) = newProperties(mmInd(k), nnInd(k), :);
                prop(k,4:6) = oldProperties(newLinkOld(mmInd(k), nnInd(k), 1), newLinkOld(mmInd(k), nnInd(k), 2), :);
            end
            threshDist = [newThreshDist, oldThreshDist];
            threshPropInd = [1 4];

        case 5 %1 property, new edge length
            %twoUnique
            %numClusters = numel(mmInd)./avgSize;
            numClusters = numel(mmInd)./sqrt(avgSize);
            sizes = ceil(numClusters);
            prop = newEdgeLengths(mmInd);
        case 6 %4 properties, new dist, new radial, old dist, old radial
            %threeUnique_parIntParInt, extNotExt, extExt
            prop = zeros(numel(mmInd), 4);
            numClusters = 2*numel(mmInd)./sqrt(avgSize);
            s = ceil(sqrt(numClusters));
            sizes = [s 1 s 1];
            for k = 1:numel(mmInd)
                prop(k,1:2) = newProperties(mmInd(k), nnInd(k), [1 3]);
                prop(k,3:4) = oldProperties(newLinkOld(mmInd(k), nnInd(k), 1), newLinkOld(mmInd(k), nnInd(k), 2), [1 3]);
            end
         case 7 %6 properties, 3 new, 3 old
            %fourUnique_notExtnotExt
            prop = zeros(numel(mmInd), 6);
            numClusters = 3*numel(mmInd)./avgSize; % 0.5
            s = ceil(numClusters.^0.25);
            %sizes = [s s 1 s s 1];
            sizes = [s s 1 s s 1];
            for k = 1:numel(mmInd)
                prop(k,1:3) = newProperties(mmInd(k), nnInd(k), :);
                prop(k,4:6) = oldProperties(newLinkOld(mmInd(k), nnInd(k), 1), newLinkOld(mmInd(k), nnInd(k), 2), :);
            end
         case 8 %3 property, new edge length,new posneg , old posneg
            %twoUnique pos neg ext
            %numClusters = numel(mmInd)./avgSize;
            numClusters = 2*numel(mmInd)./avgSize;
            s = ceil(numClusters);
            sizes = [1 1 3*s];
            prop = zeros(numel(mmInd), 3);
            for k = 1:numel(mmInd)
                prop(k,1) = newEdgeLengths(k);
                prop(k,2) = newRhoProperties(mmInd(k), nnInd(k), 2);
                prop(k,3) = oldRhoProperties(newLinkOld(mmInd(k), nnInd(k), 1), newLinkOld(mmInd(k), nnInd(k), 2), 2);
            end
         case 9%5 property, new edge length
            %twoUnique pos neg ext
            %numClusters = numel(mmInd)./avgSize;
            numClusters = 2*numel(mmInd)./avgSize;
            s = ceil(numClusters);
            sizes = [1 s s s s];
            prop = zeros(numel(mmInd), 5);
            for k = 1:numel(mmInd)
                prop(k,1) = newEdgeLengths(k);
                prop(k,2) = newRhoProperties(mmInd(k), nnInd(k), 1);
                prop(k,3) = newRhoProperties(mmInd(k), nnInd(k), 2);
                prop(k,4) = oldRhoProperties(newLinkOld(mmInd(k), nnInd(k), 1), newLinkOld(mmInd(k), nnInd(k), 2), 1);
                prop(k,5) = oldRhoProperties(newLinkOld(mmInd(k), nnInd(k), 1), newLinkOld(mmInd(k), nnInd(k), 2), 2);
            end
         %case 10 
             %twoUnique ext

    end % switch
    clusterIntervals = createIntervals(prop, sizes, threshDist, threshPropInd);
    
end

function clusterIntervals = createIntervals(prop, sizes, threshDist, threshPropInd)
    clusterIntervals = cell(numel(sizes),1);
    for k = 1:numel(sizes)
        ind = find(threshPropInd == k);
        
        if (numel(ind) > 0)
            td = threshDist(ind);
            minDist = min(prop(:,k));
            maxDist = max(prop(:,k));
            s = max(2, round(sizes(k)/2));
            distInterval = linspace(minDist, td, s );
            distInterval(s:(2*s -1)) = linspace(td, maxDist, s);
            clusterIntervals{k} = distInterval;
        else
           clusterIntervals{k} = linspace( min(prop(:,k)) , max(prop(:,k)) , sizes(k) + 1); 
        end
        
    end
end