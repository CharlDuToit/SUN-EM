function [clusterMeans, clusterInd, totError, numIter, tol] = calcClusterMeans(prop,clusterMeans, propScale, maxIter, tol)
    [~, numProp] = size(prop);
    
     %Scale properties
    for i =1:numProp
        clusterMeans(:, i) = clusterMeans(:, i) *propScale(i);
        prop(:,i) = prop(:,i) *propScale(i);
    end
    
    
    oldTotError = 0;
    for i = 1:maxIter
        [clusterMeans,clusterInd, totError] = calcNewMeans(prop,clusterMeans );
        t = abs(oldTotError - totError)/totError;
        oldTotError = totError;
        if (t < tol )
            break;
        end
    end
    numIter = i;
    tol = t;
    
    %Scale properties back
    for i =1:numProp
        clusterMeans(:, i) = clusterMeans(:, i) /propScale(i);
        prop(:,i) = prop(:,i) / propScale(i);
    end
    
    %[clusterMeans,propClusterInd, totError] = calcNewMeans(prop,clusterMeans );
    
end


function [clusterMeans, clusterInd, totError] = calcNewMeans(prop,clusterMeans )
    [numPoints, numProp] = size(prop);
    [numClusters, ~] = size(clusterMeans);
    
    clusterInd = zeros(numPoints, 1);
    clusterTotals = zeros(numClusters,numProp);
    clusterCounts = zeros(numClusters,1);
    totError = 0;
    
    for i = 1:numPoints
        [err, ind] =min((clusterMeans(:,1) - prop(i,1)).^2 + (clusterMeans(:,2) - prop(i,2)).^2 + (clusterMeans(:,3) - prop(i,3)).^2);
        clusterTotals(ind, :) = clusterTotals(ind, :) + prop(i,:);
        clusterCounts(ind, :) = clusterCounts(ind, :) + 1;
        clusterInd(i) = ind;
        totError = totError + err;
    end
    
    nonEmptyClusterInd = find(clusterCounts > 0);  
    numNonEmptyClusters = numel(nonEmptyClusterInd);
    
    for i = 1:numNonEmptyClusters
        ind = nonEmptyClusterInd(i);
        clusterMeans(ind,:) = clusterTotals(ind,:) / clusterCounts(ind);   
    end
    clusterMeans = clusterMeans(nonEmptyClusterInd, :);
    
end

%clusterMeans = [0 0 0; 1 1 1; 2 2 2; 3 3 3; 4 4 4; 5 5 5];
%prop = [4.5 4.6 4.6; 3.9 3.9 3.9; 0 0 1];


%(clusterMeans(:,1) - prop(:,1)).^2 + (clusterMeans(:,2) - prop(:,2)).^2 + (clusterMeans(:,3) - prop(:,3)).^2;
%(clusterMeans(:,1) - prop(1,1)).^2 + (clusterMeans(:,2) - prop(1,2)).^2 + (clusterMeans(:,3) - prop(1,3)).^2
%[val, ind] =min((clusterMeans(:,1) - prop(1,1)).^2 + (clusterMeans(:,2) - prop(1,2)).^2 + (clusterMeans(:,3) - prop(1,3)).^2)