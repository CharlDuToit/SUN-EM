function [predZmn, unityZmn] = predictExtractedTerms( mlmom, selfZmnTerms, triZmnTerms, nonSingZmnTerms, nonSingZmnProp, singInd,edgeLengths, includeRealCalc, returnUnity)
    
    %numClusters = numel(mlmom.clusterMeans(:,1));
    %maxClusterError = 10*ones(numClusters,1);
    useProjectedEdges = mlmom.useProjectedEdges;
    clusterMaxEdgeLength = mlmom.clusterMaxEdgeLength;
    maxClusterError = mlmom.maxClusterError;
    
    clusterMeans = mlmom.clusterMeans ;
    propScale = mlmom.propScale ;
    %clusterMeans = clusterMeans .* propScale;
    numFreq = mlmom.numFreq;

    [numWeights,~]= size(mlmom.weightModels{1, 2}.selfWeights);
    [ ~, numTerms ] = size(selfZmnTerms);
    [numEdges, ~] = size(singInd) ;
    [ numTri, ~ ] = size(triZmnTerms);
    [ numNonSing, ~ ] = size(nonSingZmnTerms);
    
    unityZmn = [];
    if (returnUnity)
        unityZmn = complex(zeros(numEdges,numEdges,numFreq));
    %else
       % unityZmn = zeros(numEdges,numEdges,numFreq);
    end
    
    predZmn = complex(zeros(numEdges,numEdges,numFreq));
    addBias = numWeights - numTerms;
    
    if (addBias)
        selfZmnTerms(:, numTerms + 1, :) = ones(numEdges, 1,1) + 1i* ones(numEdges, 1,1);
        triZmnTerms(:, numTerms + 1, :) = ones(numTri, 1,1)+ 1i* ones(numTri, 1,1);
        nonSingZmnTerms(:, numTerms + 1, :) = ones(numNonSing, 1,1)+ 1i* ones(numNonSing, 1,1);
    end
    
    AInd = [1 3 5 7];
    PhiInd =[2 4 6 8];
    weightModels = mlmom.weightModels;
    for freq = 1:mlmom.numFreq
        % NB NB NB NB
        % extractZmnInfo loops through data in same way
        % fillZmnTermsByEdge loops in same way
        
        nonSingCount = 0;
        triCount = 0;
        
        realWeightsInd = weightModels{freq,1}.nonSingWeightsInd;
        imagWeightsInd = weightModels{freq,2}.nonSingWeightsInd;
        
        for mm = 1:numEdges
            for nn = 1:numEdges
                if (mm == nn)
                    if (returnUnity)
                        unityZmn(mm,nn,freq) = sum(selfZmnTerms(mm, 1:numTerms, freq)) ;
                    end
                    if (~includeRealCalc )
                        if (returnUnity)
                            predZmn(mm,nn,freq) = real(unityZmn(mm,nn,freq));
                        else
                            predZmn(mm,nn,freq) = sum(real(selfZmnTerms(mm, 1:numTerms, freq))) ;
                        end
                    else
                        predZmn(mm,nn,freq) = real(selfZmnTerms(mm, :, freq)) *weightModels{freq,1}.selfWeights;
                    end
                    predZmn(mm,nn,freq) = predZmn(mm,nn,freq) + 1i * imag(selfZmnTerms(mm, :, freq)) *weightModels{freq,2}.selfWeights;
                    
                elseif (singInd(mm,nn))
                    triCount = triCount + 1;
                    if (returnUnity)
                        unityZmn(mm,nn,freq) = sum(triZmnTerms(triCount, 1:numTerms, freq)) ;
                    end
                    if ( ~includeRealCalc )
                        if (returnUnity)
                            predZmn(mm,nn,freq) = real(unityZmn(mm,nn,freq));
                        else
                            predZmn(mm,nn,freq) = sum(real(triZmnTerms(triCount, 1:numTerms, freq))) ;
                        end

                    else
                        predZmn(mm,nn,freq) = real(triZmnTerms(triCount, :, freq)) *weightModels{freq,1}.triWeights;
                    end
                    predZmn(mm,nn,freq) = predZmn(mm,nn,freq) + 1i * imag(triZmnTerms(triCount, :, freq)) *weightModels{freq,2}.triWeights;
                    
                else
                    nonSingCount = nonSingCount + 1;
                    if (returnUnity)
                        unityZmn(mm,nn,freq) = sum(nonSingZmnTerms(nonSingCount, 1:numTerms, freq)) ;
                    end
                    prop = nonSingZmnProp(nonSingCount,:);                  
                    [err, ind] = calcMinError(clusterMeans(imagWeightsInd,:), prop, propScale);
                    
                    %if (maxClusterError(ind) < err)
                    clusterIndex = imagWeightsInd(ind);
                    if (maxClusterError(clusterIndex) < err)
                        if (returnUnity)
                            predZmn(mm,nn,freq) = unityZmn(mm,nn,freq);
                        else
                            predZmn(mm,nn,freq) = sum(nonSingZmnTerms(nonSingCount, 1:numTerms, freq)) ;
                        end
                    else
                        if ( ~includeRealCalc )
                            if (returnUnity)
                                predZmn(mm,nn,freq) = real(unityZmn(mm,nn,freq));
                            else
                                predZmn(mm,nn,freq) = sum(real(nonSingZmnTerms(nonSingCount, 1:numTerms, freq))) ;
                            end
                        else
                            realIndex = find(realWeightsInd == clusterIndex);
                            if (numel(realIndex) > 0)
                                %use weights
                                if (~useProjectedEdges)
                                    predZmn(mm,nn,freq) = real(nonSingZmnTerms(nonSingCount, :, freq)) * weightModels{freq,1}.nonSingWeights(realIndex(1),:)';
                                else
                                    t = real(nonSingZmnTerms(nonSingCount, :, freq));
                                    a = clusterMaxEdgeLength(clusterIndex,1)./edgeLengths(mm);
                                    b = clusterMaxEdgeLength(clusterIndex,2)./edgeLengths(nn);
                                    
                                    projectedTerms = zeros(1, numTerms);
                                    projectedTerms(1,AInd) = t(1,AInd) .* a .* b.^2;
                                    projectedTerms(1,PhiInd) = t(1,PhiInd) .* b;
                                    
                                    weights =  weightModels{freq,1}.nonSingWeights(realIndex(1),:);
                                    
                                    regressedProjectedTerms = projectedTerms(1,:) .* weights;
                                    regressedTerms = zeros(1, numTerms);
                                    
                                    regressedTerms(1,AInd) = regressedProjectedTerms(1,AInd) ./( a .* b.^2);
                                    regressedTerms(1,PhiInd) = regressedProjectedTerms(1,PhiInd) ./ b;
                                    predZmn(mm,nn,freq) = sum(regressedTerms,2);
                                end
                                
                            else
                                if (returnUnity)
                                    predZmn(mm,nn,freq) = real(unityZmn(mm,nn,freq));
                                else
                                    predZmn(mm,nn,freq) = sum(real(nonSingZmnTerms(nonSingCount, 1:numTerms, freq))) ;
                                end
                            end
                            %predZmn(mm,nn,freq) = real(nonSingZmnTerms(nonSingCount, :, freq)) * weightModels{freq,1}.nonSingWeights(ind,:)';
                            
                        end
                        % use weights
                        if (~useProjectedEdges)
                             predZmn(mm,nn,freq) = predZmn(mm,nn,freq) + 1i * imag(nonSingZmnTerms(nonSingCount, :, freq)) *weightModels{freq,2}.nonSingWeights(ind,:)';
                        else
                            t = imag(nonSingZmnTerms(nonSingCount, :, freq));
                            a = clusterMaxEdgeLength(clusterIndex,1)./edgeLengths(mm);
                            b = clusterMaxEdgeLength(clusterIndex,2)./edgeLengths(nn);
                            
                            projectedTerms = zeros(1, numTerms);
                            projectedTerms(1,AInd) = t(1,AInd) .* a .* b.^2;
                            projectedTerms(1,PhiInd) = t(1,PhiInd) .* b;
                            
                            weights =  weightModels{freq,2}.nonSingWeights(ind,:);
                            
                            regressedProjectedTerms = projectedTerms(1,:) .* weights;
                            regressedTerms = zeros(1, numTerms);
                            
                            regressedTerms(1,AInd) = regressedProjectedTerms(1,AInd) ./( a .* b.^2);
                            regressedTerms(1,PhiInd) = regressedProjectedTerms(1,PhiInd) ./ b;
                            predZmn(mm,nn,freq) = predZmn(mm,nn,freq) + 1i*sum(regressedTerms,2);
                        end                                           
                    end
                end
            end % for nn
        end % for mm
    end % for freq

end

function [err, ind] = calcMinError(clusterMeans, prop, propScale)
    %multiplying distance of cluster and prop does not affect log
   % propScale = [4 1 0.6];
    %[err, ind] =min(abs( log(clusterMeans(:,1)/prop(1) ) ) + (clusterMeans(:,2) - prop(2)).^2 + (clusterMeans(:,3) - prop(3)).^2);
    [err, ind] =min(propScale(1) * abs( log(clusterMeans(:,1)/prop(1) ) ) +...
    propScale(2)*abs(clusterMeans(:,2) - prop(2)) + ...
    propScale(3)*abs(clusterMeans(:,3) - prop(3)));
end