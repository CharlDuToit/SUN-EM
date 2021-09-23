%function  [clusterWeights, pred] = calcClusterWeights(terms, refVals, clusterInd, singThresh, addBias)
function  [clusterWeights, pred,clusterWeightsInd] = calcClusterWeights(terms, refVals,nonSingZmnProp, clusterInd,edgeLabels, edgeLengths,clusterMaxEdgeLength, singThresh, addBias, minPercentImprov, useProjectedEdges)
    [numObs, numTerms] = size(terms);
    %[numObs, 1] = size(ClusterInd);
    %[numObs, 1] = size(refVals)
    [numClusters ,~] = max(clusterInd);
    clusterWeights = ones(numClusters,numTerms );
    clusterWeightsInd = zeros(numClusters,1 );
    pred = zeros(numObs,1);
    returnEmptyWeightsIfUnity = 1;
    predClusterRelNormPercentError = zeros(numClusters,1);
    unityClusterRelNormPercentError = zeros(numClusters,1);
   % addBias = 0;
    if (addBias)
        clusterWeights(:, numTerms +1 ) = zeros(numClusters, 1);
    end
    
    nonUnityClusterWeightsCount = 0;
    AInd = [1 3 5 7];
    PhiInd =[2 4 6 8];
    for k = 1:numClusters
        
        ind = find(clusterInd(:,1) == k);
        t = terms(ind,:);       
        unity  = sum(t, 2);
        
        if (~useProjectedEdges)
            ref = refVals(ind);
            [weights, pred(ind,1),~ ]=MLR(t,ref, singThresh, addBias, returnEmptyWeightsIfUnity);
        else
            [zeroUnityInd, ~] = find(unity == 0);
            if (numel(zeroUnityInd) > 0)
                [nonZeroUnityInd, ~] = find(unity ~= 0);
                ind = ind(nonZeroUnityInd);
                t = terms(ind,:);
                unity = unity(nonZeroUnityInd);
            end
            
            ref = refVals(ind);
            numPoints = numel(ind);
            prop = nonSingZmnProp(ind,:);
            
            labels = edgeLabels(ind,:);
            lengths = zeros(numPoints, 2);
            lengths(: , 1) = edgeLengths(labels(:, 1));
            lengths(: , 2) = edgeLengths(labels(:, 2));
                        
            a = clusterMaxEdgeLength(k,1)./lengths(:,1);
            b = clusterMaxEdgeLength(k,2)./lengths(:,2);
            
%             minDist = min(prop(:,1));
%             c = prop(:,1)./minDist;
            
            projectedRef = zeros(numPoints,1);
%             projectedRef(:,1) = (sum( t(:,AInd), 2) ./ unity(:,1) ) .* ref(:,1).* a .* b.^2 .* c;
%             projectedRef(:,1) = projectedRef(:,1) + (sum( t(:,PhiInd), 2) ./ unity(:,1) ) .* ref(:,1) .* b .*c;
            projectedRef(:,1) = (sum( t(:,AInd), 2) ./ unity(:,1) ) .* ref(:,1).* a .* b.^2;
            projectedRef(:,1) = projectedRef(:,1) + (sum( t(:,PhiInd), 2) ./ unity(:,1) ) .* ref(:,1) .* b;
            projectedTerms = zeros(numPoints, numTerms);
%             projectedTerms(:,AInd) = t(:,AInd) .* a .* b.^2 .*c;
%             projectedTerms(:,PhiInd) = t(:,PhiInd) .* b .*c;
            projectedTerms(:,AInd) = t(:,AInd) .* a .* b.^2;
            projectedTerms(:,PhiInd) = t(:,PhiInd) .* b;
            
            [weights, ~,regressedProjectedTerms ]= MLR(projectedTerms,projectedRef, singThresh, addBias, returnEmptyWeightsIfUnity);
            
            regressedTerms = zeros(numPoints, numTerms);
%             regressedTerms(:,AInd) = regressedProjectedTerms(:,AInd) ./( a .* b.^2 .*c);
%             regressedTerms(:,PhiInd) = regressedProjectedTerms(:,PhiInd) ./( b .*c);
            regressedTerms(:,AInd) = regressedProjectedTerms(:,AInd) ./( a .* b.^2);
            regressedTerms(:,PhiInd) = regressedProjectedTerms(:,PhiInd) ./ b;
            pred(ind,1) = sum(regressedTerms,2);
        
        end
       

        if (eq(size(weights), [0 0]))
            pred(ind,1) = unity;
        else

          unityClusterRelNormPercentError(k) = 100* sqrt(((unity - ref)' * (unity - ref))/ sum(ref.^2) ) ;
          predClusterRelNormPercentError(k) = 100* sqrt(((pred(ind,1) - ref)' * (pred(ind,1) - ref))/ sum(ref.^2) ) ;
          
          %did not correct error enough
          if (unityClusterRelNormPercentError(k) - predClusterRelNormPercentError(k) < minPercentImprov)
              pred(ind,1) = unity;
          else
              nonUnityClusterWeightsCount = nonUnityClusterWeightsCount + 1;
              clusterWeights(nonUnityClusterWeightsCount,:) = weights;
              clusterWeightsInd(nonUnityClusterWeightsCount) = k;
          end
            

        end
        
    end% for k
    clusterWeights = clusterWeights(1:nonUnityClusterWeightsCount, :);
    clusterWeightsInd = clusterWeightsInd(1:nonUnityClusterWeightsCount,:);

end