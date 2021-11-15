function [clusterTotals, clusterInd, clusterCounts,maxClusterError, totError, pointErrors] = assignClusters(prop,clusterMeans, errorCode )
    [numPoints, numProp] = size(prop);
    [numClusters, ~] = size(clusterMeans);
    
    maxClusterError = zeros(numClusters,1);
    clusterInd = zeros(numPoints, 1);
    pointErrors = zeros(numPoints, 1);
    clusterTotals = zeros(numClusters,numProp);
    clusterCounts = zeros(numClusters,1);
    
    totError = 0;
    
    for k = 1:numPoints
        %[err, ind] =min((clusterMeans(:,1) - prop(k,1)).^2 + (clusterMeans(:,2) - prop(k,2)).^2 + (clusterMeans(:,3) - prop(k,3)).^2);
        %if (numel(clusterMeans) == 0)
            %a= 5;
        %end
        [err, ind] = calcMinClusterError(clusterMeans, prop(k,:), errorCode);
        clusterTotals(ind, :) = clusterTotals(ind, :) + prop(k,:);
        clusterCounts(ind) = clusterCounts(ind) + 1;
        clusterInd(k) = ind;
        if (maxClusterError(ind) < err)
            maxClusterError(ind) = err;
        end
        totError = totError + err;
        pointErrors(k) = err;
    end
end