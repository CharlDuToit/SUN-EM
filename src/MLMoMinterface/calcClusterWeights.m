function  [clusterWeights, pred] = calcClusterWeights(terms, refVals, clusterInd,singThresh, addBias)
    [numObs, numTerms] = size(terms);
    %[numObs, 1] = size(ClusterInd);
    %[numObs, 1] = size(refVals)
    [numClusters ,~] = max(clusterInd);
    clusterWeights = ones(numClusters,numTerms );
    pred = zeros(numObs,1);
    returnEmptyWeightsIfSingular = 0;
    
    if (addBias)
        clusterWeights(:, numTerms +1 ) = zeros(numClusters, 1);
    end
    
    for i = 1:numClusters
        ind = find(clusterInd(:,1) == i);
        t = terms(ind,:);
        ref = refVals(ind);
        [clusterWeights(i,:), pred(ind,1) ]= MLR(t,ref, singThresh, addBias, returnEmptyWeightsIfSingular);
    end

end