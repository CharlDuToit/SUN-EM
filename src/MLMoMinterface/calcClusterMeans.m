function [clusterMeans, clusterInd, clusterCounts,maxClusterError, totError, numIter, tol] = calcClusterMeans(prop,clusterMeans, propScale, maxIter, tol)
    
    
     %Scale properties
%      [~, numProp] = size(prop);
%     for i =1:numProp
%         clusterMeans(:, i) = clusterMeans(:, i) *propScale(i);
%         prop(:,i) = prop(:,i) *propScale(i);
%     end
    
    
    oldTotError = 0;
    for i = 1:maxIter
        [clusterMeans,totError] = calcNewMeans(prop,clusterMeans,propScale );
        t = abs(oldTotError - totError)/totError;
        oldTotError = totError;
        if (t < tol )
            break;
        end
    end
    numIter = i;
    %tol = t;
    
    %cluster indicies are still assigned to the previous means
    [~, clusterInd,clusterCounts,maxClusterError, totError] = assignClusters(prop,clusterMeans,propScale );
    tol = abs(oldTotError - totError)/totError;
    
    %Scale properties back
%     for i =1:numProp
%         clusterMeans(:, i) = clusterMeans(:, i) /propScale(i);
%         prop(:,i) = prop(:,i) / propScale(i);
%     end
    
end


function [clusterMeans, totError] = calcNewMeans(prop,clusterMeans,propScale )

    [clusterTotals, ~,clusterCounts,~, totError] = assignClusters(prop,clusterMeans,propScale );
    
    nonEmptyClusterInd = find(clusterCounts > 0);  
    numNonEmptyClusters = numel(nonEmptyClusterInd);
    
    for i = 1:numNonEmptyClusters
        ind = nonEmptyClusterInd(i);
        clusterMeans(ind,:) = clusterTotals(ind,:) / clusterCounts(ind);   
    end
    clusterMeans = clusterMeans(nonEmptyClusterInd, :);
    
end

function [clusterTotals, clusterInd,clusterCounts,maxClusterError, totError] = assignClusters(prop,clusterMeans,propScale )
    [numPoints, numProp] = size(prop);
    [numClusters, ~] = size(clusterMeans);
    
    maxClusterError = zeros(numClusters,1);
    clusterInd = zeros(numPoints, 1);
    clusterTotals = zeros(numClusters,numProp);
    clusterCounts = zeros(numClusters,1);
    totError = 0;
    
    for k = 1:numPoints
        %[err, ind] =min((clusterMeans(:,1) - prop(k,1)).^2 + (clusterMeans(:,2) - prop(k,2)).^2 + (clusterMeans(:,3) - prop(k,3)).^2);
        [err, ind] = calcMinError(clusterMeans, prop(k,:),propScale);
        clusterTotals(ind, :) = clusterTotals(ind, :) + prop(k,:);
        clusterCounts(ind, :) = clusterCounts(ind, :) + 1;
        clusterInd(k) = ind;
        if (maxClusterError(ind) < err)
            maxClusterError(ind) = err;
        end
        totError = totError + err;
    end
end

function [err, ind] = calcMinError(clusterMeans, prop, propScale)
    %multiplying distance of cluster and prop does not affect log
   % propScale = [4 1 0.6];
    %[err, ind] =min(abs( log(clusterMeans(:,1)/prop(1) ) ) + (clusterMeans(:,2) - prop(2)).^2 + (clusterMeans(:,3) - prop(3)).^2);
    [err, ind] =min(propScale(1) * abs( log(clusterMeans(:,1)/prop(1) ) ) +...
    propScale(2)*abs(clusterMeans(:,2) - prop(2)) + ...
    propScale(3)*abs(clusterMeans(:,3) - prop(3)));
end
%clusterMeans = [0 0 0; 1 1 1; 2 2 2; 3 3 3; 4 4 4; 5 5 5];
%prop = [4.5 4.6 4.6; 3.9 3.9 3.9; 0 0 1];


%(clusterMeans(:,1) - prop(:,1)).^2 + (clusterMeans(:,2) - prop(:,2)).^2 + (clusterMeans(:,3) - prop(:,3)).^2;
%(clusterMeans(:,1) - prop(1,1)).^2 + (clusterMeans(:,2) - prop(1,2)).^2 + (clusterMeans(:,3) - prop(1,3)).^2
%[val, ind] =min((clusterMeans(:,1) - prop(1,1)).^2 + (clusterMeans(:,2) - prop(1,2)).^2 + (clusterMeans(:,3) - prop(1,3)).^2)