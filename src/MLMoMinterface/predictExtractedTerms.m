function [zMatrices] = predictExtractedTerms( mlmom, selfZmnTerms, triZmnTerms, nonSingZmnTerms, nonSingZmnProp, singInd, includeRealCalc)

    clusterMeans = mlmom.clusterMeans ;
    clusterPropScale = mlmom.clusterPropScale ;
    clusterMeans = clusterMeans .* clusterPropScale;
    numFreq = mlmom.numFreq;

    [numWeights,~]= size(mlmom.weightModels{1, 2}.selfWeights);
    [ ~, numTerms ] = size(selfZmnTerms);
    [numEdges, ~] = size(singInd) ;
    [ numTri, ~ ] = size(triZmnTerms);
    [ numNonSing, ~ ] = size(nonSingZmnTerms);
    
    zMatrices = complex(zeros(numEdges,numEdges,numFreq));
    addBias = numWeights - numTerms;
    
    if (addBias)
        selfZmnTerms(:, numTerms + 1, :) = ones(numEdges, 1,1);
        triZmnTerms(:, numTerms + 1, :) = ones(numTri, 1,1);
        nonSingZmnTerms(:, numTerms + 1, :) = ones(numNonSing, 1,1);
    end
       
    for freq = 1:mlmom.numFreq
        % NB NB NB NB
        % extractZmnInfo loops through data in same way
        % fillZmnTermsByEdge loops in same way  

        % real and imag
        for k = 1:2  
            nonSingCount = 0;
            triCount = 0;
            weightModel = mlmom.weightModels{freq, k};
            for mm = 1:numEdges
                for nn = 1:numEdges
                    if (mm == nn)
                        
                        if (k == 1 && ~includeRealCalc )
                            zMatrices(mm,nn,freq) = sum(real(selfZmnTerms(mm, :, freq))) - addBias;
                        elseif (k ==1)
                            zMatrices(mm,nn,freq) = real(selfZmnTerms(mm, :, freq)) *weightModel.selfWeights;
                        else
                            zMatrices(mm,nn,freq) = zMatrices(mm,nn,freq) + 1i * imag(selfZmnTerms(mm, :, freq)) *weightModel.selfWeights;
                        end
                        
                    elseif (singInd(mm,nn))
                        triCount = triCount + 1;
                        
                        if (k == 1 && ~includeRealCalc )
                            zMatrices(mm,nn,freq) = sum(real(triZmnTerms(triCount, :, freq))) - addBias;
                        elseif (k ==1)
                            zMatrices(mm,nn,freq) = real(triZmnTerms(triCount, :, freq)) *weightModel.triWeights;
                        else
                            zMatrices(mm,nn,freq) = zMatrices(mm,nn,freq) + 1i * imag(triZmnTerms(triCount, :, freq)) *weightModel.triWeights;
                        end
                        
                    else
                        nonSingCount = nonSingCount + 1;
                        
                        if (k == 1 && ~includeRealCalc )
                            zMatrices(mm,nn,freq) = sum(real(nonSingZmnTerms(nonSingCount, :, freq))) - addBias;
                            continue
                        end
                         prop = nonSingZmnProp(nonSingCount,:) .* clusterPropScale;
                        % means already scaled
                        [err, ind] =min((clusterMeans(:,1) - prop(1)).^2 +...
                            (clusterMeans(:,2) - prop(2)).^2 +...
                            (clusterMeans(:,3) - prop(3)).^2);

                        if (k ==1)
                            zMatrices(mm,nn,freq) = real(nonSingZmnTerms(nonSingCount, :, freq)) * weightModel.nonSingWeights(ind,:)';
                        else
                            zMatrices(mm,nn,freq) = zMatrices(mm,nn,freq) +...
                                1i * imag(nonSingZmnTerms(nonSingCount, :, freq)) * weightModel.nonSingWeights(ind,:)';
                        end
                        
                    end
                    
                end % for nn
            end % for mm
            
            
        end % for ii
    end

end