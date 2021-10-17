function [clusterMeans, clusterInd, clusterCounts,maxClusterError, totError, numIter, tol] = calcClusterMeans(prop,clusterMeans, maxIter, tol, errorCode)
   
    oldTotError = 0;
    for i = 1:maxIter
        [clusterMeans,totError] = calcNewMeans(prop,clusterMeans, errorCode );
        t = abs(oldTotError - totError)/totError;
        oldTotError = totError;
        if (t < tol )
            break;
        end
    end
    numIter = i;
    %tol = t;
    
    %cluster indicies are still assigned to the previous means
    [~, clusterInd,clusterCounts,maxClusterError, totError] = assignClusters(prop,clusterMeans, errorCode );
    tol = abs(oldTotError - totError)/totError;
    
end

function [clusterMeans, totError] = calcNewMeans(prop,clusterMeans, errorCode )

    [clusterTotals, ~ ,clusterCounts,~, totError] = assignClusters(prop,clusterMeans, errorCode );
    
    nonEmptyClusterInd = find(clusterCounts > 0);  
    numNonEmptyClusters = numel(nonEmptyClusterInd);
    
    for i = 1:numNonEmptyClusters
        ind = nonEmptyClusterInd(i);
        clusterMeans(ind,:) = clusterTotals(ind,:) / clusterCounts(ind);   
    end
    clusterMeans = clusterMeans(nonEmptyClusterInd, :);
    
end

% function [clusterTotals, clusterInd, clusterCounts,maxClusterError, totError] = assignClusters(prop,clusterMeans, errorCode )
%     [numPoints, numProp] = size(prop);
%     [numClusters, ~] = size(clusterMeans);
%     
%     maxClusterError = zeros(numClusters,1);
%     clusterInd = zeros(numPoints, 1);
%     clusterTotals = zeros(numClusters,numProp);
%     clusterCounts = zeros(numClusters,1);
%     
%     totError = 0;
%     
%     for k = 1:numPoints
%         %[err, ind] =min((clusterMeans(:,1) - prop(k,1)).^2 + (clusterMeans(:,2) - prop(k,2)).^2 + (clusterMeans(:,3) - prop(k,3)).^2);
%         [err, ind] = calcMinClusterError(clusterMeans, prop(k,:), errorCode);
%         clusterTotals(ind, :) = clusterTotals(ind, :) + prop(k,:);
%         clusterCounts(ind) = clusterCounts(ind) + 1;
%         clusterInd(k) = ind;
%         if (maxClusterError(ind) < err)
%             maxClusterError(ind) = err;
%         end
%         totError = totError + err;
%     end
% end
