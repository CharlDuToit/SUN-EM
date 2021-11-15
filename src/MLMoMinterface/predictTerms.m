%function [predZmn, unityZmn] = predictExtractedTerms( mlmom, selfZmnTerms, triZmnTerms, nonSingZmnTerms, nonSingZmnProp, singInd,edgeLengths, returnUnity)
function [predZmn, indicesStruct] = predictTerms( weightModel, clusterStruct,propStruct, termStruct, indicesStruct, hasCalcInd)

    %numClusters = numel(mlmom.clusterMeans(:,1));
    %maxClusterError = 10*ones(numClusters,1);
    %includeRealCalc = mlmom.includeRealCalc;
%     useProjectedEdges = 0;%mlmom.useProjectedEdges;
%     clusterMaxEdgeLength = 0;%mlmom.clusterMaxEdgeLength;
%     maxClusterError = 0;%mlmom.maxClusterError;
%     
%     clusterMeans = mlmom.clusterMeans ;
%     %clusterMeans = clusterMeans .* propScale;
%     numFreq = mlmom.numFreq;
%     errorCode = mlmom.errorCode;
% 
%     [numWeights,~]= size(mlmom.weightModels{1, 2}.selfWeights);
%     [ ~, numTerms ] = size(selfZmnTerms);
%     [numEdges, ~] = size(singInd) ;
%     [ numTri, ~ ] = size(triZmnTerms);
%     [ numNonSing, ~ ] = size(nonSingZmnTerms);
%     
%     unityZmn = [];
%     if (returnUnity)
%         unityZmn = complex(zeros(numEdges,numEdges,numFreq));
%     end
%     
%     predZmn = complex(zeros(numEdges,numEdges,numFreq));
%     addBias = numWeights - numTerms;
%     
%     if (addBias)
%         selfZmnTerms(:, numTerms + 1, :) = ones(numEdges, 1,1) + 1i* ones(numEdges, 1,1);
%         triZmnTerms(:, numTerms + 1, :) = ones(numTri, 1,1)+ 1i* ones(numTri, 1,1);
%         nonSingZmnTerms(:, numTerms + 1, :) = ones(numNonSing, 1,1)+ 1i* ones(numNonSing, 1,1);
%     end
%     
%     AInd = [1 3 5 7];
%     PhiInd =[2 4 6 8];
%     weightModels = mlmom.weightModels;
%     for freq = 1:mlmom.numFreq
%         % NB NB NB NB
%         % extractZmnInfo loops through data in same way
%         % fillZmnTermsByEdge loops in same way
%         
%         nonSingCount = 0;
%         triCount = 0;
%         
%         if (includeRealCalc)
%              realWeightsInd = weightModels{freq,1}.nonSingWeightsInd;
%         end
%         imagWeightsInd = weightModels{freq,2}.nonSingWeightsInd;
%         
%         for mm = 1:numEdges
%             for nn = 1:numEdges
%                 %=========2 Unique ============
%                 if (mm == nn)
%                     if (returnUnity)
%                         unityZmn(mm,nn,freq) = sum(selfZmnTerms(mm, 1:numTerms, freq)) ;
%                     end
%                     if (~includeRealCalc )
%                         if (returnUnity)
%                             predZmn(mm,nn,freq) = real(unityZmn(mm,nn,freq));
%                         else
%                             predZmn(mm,nn,freq) = sum(real(selfZmnTerms(mm, 1:numTerms, freq))) ;
%                         end
%                     else
%                         predZmn(mm,nn,freq) = real(selfZmnTerms(mm, :, freq)) *weightModels{freq,1}.selfWeights;
%                     end
%                     predZmn(mm,nn,freq) = predZmn(mm,nn,freq) + 1i * imag(selfZmnTerms(mm, :, freq)) *weightModels{freq,2}.selfWeights;
%                 %=========3 Unique ============    
%                 elseif (singInd(mm,nn))
%                     triCount = triCount + 1;
%                     if (returnUnity)
%                         unityZmn(mm,nn,freq) = sum(triZmnTerms(triCount, 1:numTerms, freq)) ;
%                     end
%                     if ( ~includeRealCalc )
%                         if (returnUnity)
%                             predZmn(mm,nn,freq) = real(unityZmn(mm,nn,freq));
%                         else
%                             predZmn(mm,nn,freq) = sum(real(triZmnTerms(triCount, 1:numTerms, freq))) ;
%                         end
% 
%                     else
%                         predZmn(mm,nn,freq) = real(triZmnTerms(triCount, :, freq)) *weightModels{freq,1}.triWeights;
%                     end
%                     predZmn(mm,nn,freq) = predZmn(mm,nn,freq) + 1i * imag(triZmnTerms(triCount, :, freq)) *weightModels{freq,2}.triWeights;
%                 %=========4 Unique ============   
%                 else
%                     nonSingCount = nonSingCount + 1;
%                     if (returnUnity)
%                         unityZmn(mm,nn,freq) = sum(nonSingZmnTerms(nonSingCount, 1:numTerms, freq)) ;
%                     end
%                     prop = nonSingZmnProp(nonSingCount,:);                  
%                     [err, ind] = calcMinClusterError(clusterMeans(imagWeightsInd,:), prop, errorCode);
%                     
%                     %if (maxClusterError(ind) < err)
%                     clusterIndex = imagWeightsInd(ind);
%                     if (maxClusterError(clusterIndex) < err)
%                         if (returnUnity)
%                             predZmn(mm,nn,freq) = unityZmn(mm,nn,freq);
%                         else
%                             predZmn(mm,nn,freq) = sum(nonSingZmnTerms(nonSingCount, 1:numTerms, freq)) ;
%                         end
%                     else
%                         if ( ~includeRealCalc )
%                             if (returnUnity)
%                                 predZmn(mm,nn,freq) = real(unityZmn(mm,nn,freq));
%                             else
%                                 predZmn(mm,nn,freq) = sum(real(nonSingZmnTerms(nonSingCount, 1:numTerms, freq))) ;
%                             end
%                         else
%                             realIndex = find(realWeightsInd == clusterIndex);
%                             if (numel(realIndex) > 0)
%                                 %use weights
%                                 if (~useProjectedEdges)
%                                     predZmn(mm,nn,freq) = real(nonSingZmnTerms(nonSingCount, :, freq)) * weightModels{freq,1}.nonSingWeights(realIndex(1),:)';
%                                 else
%                                     t = real(nonSingZmnTerms(nonSingCount, :, freq));
%                                     a = clusterMaxEdgeLength(clusterIndex,1)./edgeLengths(mm);
%                                     b = clusterMaxEdgeLength(clusterIndex,2)./edgeLengths(nn);
%                                     
%                                     projectedTerms = zeros(1, numTerms);
%                                     projectedTerms(1,AInd) = t(1,AInd) .* a .* b.^2;
%                                     projectedTerms(1,PhiInd) = t(1,PhiInd) .* b;
%                                     
%                                     weights =  weightModels{freq,1}.nonSingWeights(realIndex(1),:);
%                                     
%                                     regressedProjectedTerms = projectedTerms(1,:) .* weights;
%                                     regressedTerms = zeros(1, numTerms);
%                                     
%                                     regressedTerms(1,AInd) = regressedProjectedTerms(1,AInd) ./( a .* b.^2);
%                                     regressedTerms(1,PhiInd) = regressedProjectedTerms(1,PhiInd) ./ b;
%                                     predZmn(mm,nn,freq) = sum(regressedTerms,2);
%                                 end
%                                 
%                             else %if (numel(realIndex) > 0)
%                                 if (returnUnity)
%                                     predZmn(mm,nn,freq) = real(unityZmn(mm,nn,freq));
%                                 else
%                                     predZmn(mm,nn,freq) = sum(real(nonSingZmnTerms(nonSingCount, 1:numTerms, freq))) ;
%                                 end
%                             end
%                             %predZmn(mm,nn,freq) = real(nonSingZmnTerms(nonSingCount, :, freq)) * weightModels{freq,1}.nonSingWeights(ind,:)';
%                             
%                         end
%                         % use weights
%                         if (~useProjectedEdges)
%                              predZmn(mm,nn,freq) = predZmn(mm,nn,freq) + 1i * imag(nonSingZmnTerms(nonSingCount, :, freq)) *weightModels{freq,2}.nonSingWeights(ind,:)';
%                         else
%                             t = imag(nonSingZmnTerms(nonSingCount, :, freq));
%                             a = clusterMaxEdgeLength(clusterIndex,1)./edgeLengths(mm);
%                             b = clusterMaxEdgeLength(clusterIndex,2)./edgeLengths(nn);
%                             
%                             projectedTerms = zeros(1, numTerms);
%                             projectedTerms(1,AInd) = t(1,AInd) .* a .* b.^2;
%                             projectedTerms(1,PhiInd) = t(1,PhiInd) .* b;
%                             
%                             weights =  weightModels{freq,2}.nonSingWeights(ind,:);
%                             
%                             regressedProjectedTerms = projectedTerms(1,:) .* weights;
%                             regressedTerms = zeros(1, numTerms);
%                             
%                             regressedTerms(1,AInd) = regressedProjectedTerms(1,AInd) ./( a .* b.^2);
%                             regressedTerms(1,PhiInd) = regressedProjectedTerms(1,PhiInd) ./ b;
%                             predZmn(mm,nn,freq) = predZmn(mm,nn,freq) + 1i*sum(regressedTerms,2);
%                         end                                           
%                     end
%                 end
%             end % for nn
%         end % for mm
%     end % for freq
    
    %==========================
    numEdges = termStruct.numEdges;
    predZmn = zeros(numEdges,numEdges);
    

    %=====2 unique
    if ( clusterStruct.twoUniqueErrorCode ~= 0)
        means = clusterStruct.twoUniqueMeans(weightModel.twoUniqueWeightsIndices, :);
        twoUniqueMaxError = clusterStruct.twoUniqueMaxError(weightModel.twoUniqueWeightsIndices, :);
        if (~hasCalcInd)
            [~, indicesStruct.twoUniqueIndices(:,3), ~,~, ~, indicesStruct.twoUniqueIndices(:,4)] = assignClusters(propStruct.twoUniqueProp, means,...
                clusterStruct.twoUniqueErrorCode);
        end
        twoUniqueIndices = indicesStruct.twoUniqueIndices;
        %terms = termStruct.twoUniqueTerms;
        %weightModel.twoUniqueWeights
        for k = 1:numel(twoUniqueIndices(:,1))
            if (indicesStruct.twoUniqueIndices(k,4) <= twoUniqueMaxError(twoUniqueIndices(k,3) ) )
                predZmn(twoUniqueIndices(k,1), twoUniqueIndices(k,2)) = termStruct.twoUniqueTerms(k,:) * weightModel.twoUniqueWeights(twoUniqueIndices(k,3), :)';
            else
                predZmn(twoUniqueIndices(k,1), twoUniqueIndices(k,2)) = sum(termStruct.twoUniqueTerms(k,:));
            end
        end
    else
        twoUniqueIndices = indicesStruct.twoUniqueIndices;
        for k = 1:numel(twoUniqueIndices(:,1))
            predZmn(twoUniqueIndices(k,1), twoUniqueIndices(k,2)) = termStruct.twoUniqueTerms(k,:) * weightModel.twoUniqueWeights;
            %predZmn(twoUniqueIndices(k,1), twoUniqueIndices(k,2)) = sum(termStruct.twoUniqueTerms(k,:));
        end
    end
    %=====3 unique
    if ( clusterStruct.threeUniqueErrorCode ~= 0)
        means = clusterStruct.threeUniqueMeans(weightModel.threeUniqueWeightsIndices, :);
        threeUniqueMaxError = clusterStruct.threeUniqueMaxError(weightModel.threeUniqueWeightsIndices, :);
        if (~hasCalcInd)
            [~, indicesStruct.threeUniqueIndices(:,3), ~,~, ~, indicesStruct.threeUniqueIndices(:,4)] = assignClusters(propStruct.threeUniqueProp, means,...
                clusterStruct.threeUniqueErrorCode);
        end
        
        threeUniqueIndices = indicesStruct.threeUniqueIndices;
        for k = 1:numel(threeUniqueIndices(:,1))
            if (indicesStruct.threeUniqueIndices(k,4) <= threeUniqueMaxError(threeUniqueIndices(k,3) ) )
                predZmn(threeUniqueIndices(k,1), threeUniqueIndices(k,2)) = termStruct.threeUniqueTerms(k,:) * weightModel.threeUniqueWeights(threeUniqueIndices(k,3), :)';
            else
                predZmn(threeUniqueIndices(k,1), threeUniqueIndices(k,2)) = sum(termStruct.threeUniqueTerms(k,:));
            end
        end
    else
        threeUniqueIndices = indicesStruct.threeUniqueIndices;
        for k = 1:numel(threeUniqueIndices(:,1))
            predZmn(threeUniqueIndices(k,1), threeUniqueIndices(k,2)) = termStruct.threeUniqueTerms(k,:) * weightModel.threeUniqueWeights;
            %predZmn(threeUniqueIndices(k,1), threeUniqueIndices(k,2)) = sum(termStruct.threeUniqueTerms(k,:));
        end
    end
    %=====4 unique
    %tic;
    means = clusterStruct.fourUniqueMeans(weightModel.fourUniqueWeightsIndices, :);
    fourUniqueMaxErrors = clusterStruct.fourUniqueMaxError(weightModel.fourUniqueWeightsIndices, :);
    if (~hasCalcInd)
        [~, indicesStruct.fourUniqueIndices(:,3), ~,~, ~, indicesStruct.fourUniqueIndices(:,4)] = assignClusters(propStruct.fourUniqueProp, means,...
            clusterStruct.fourUniqueErrorCode);
    end
    %time = toc;
    %a = 5;
    %tic;
    fourUniqueIndices = indicesStruct.fourUniqueIndices;
    for k = 1:numel(fourUniqueIndices(:,1))
        if ( fourUniqueIndices(k,4) <= fourUniqueMaxErrors(fourUniqueIndices(k,3)) )
            predZmn(fourUniqueIndices(k,1), fourUniqueIndices(k,2)) = termStruct.fourUniqueTerms(k,:) * weightModel.fourUniqueWeights(fourUniqueIndices(k,3), :)';
            %predZmn(fourUniqueIndices(k,1), fourUniqueIndices(k,2)) = sum(termStruct.fourUniqueTerms(k,:));
        else
            predZmn(fourUniqueIndices(k,1), fourUniqueIndices(k,2)) = sum(termStruct.fourUniqueTerms(k,:));
        end
    end
    %time2 = toc;
    

end

% function [clustersUnity,clustersPred , clustersUnityMSE, clustersPredMSE,clustersUnityRelNormPercentError,...
%     clustersPredRelNormPercentError, totalUnityMSE, totalPredMSE, totalUnityRelNormPercentError, totalPredRelNormPercentError ] =...
%     applyWeights(clustersWeights, indices,terms, refZmn, includeError)
%     
%     mmInd = indices(:,1);
%     nnInd = indices(:,2);
%     numClusters = max(indices(:,3));
%     [~, ~, numTerms] = size(terms);
% 
%     clustersTerms = zeros(numel(mmInd), numTerms);
%     clustersPred = zeros(numel(mmInd), 1);
%     clustersUnity = zeros(numel(mmInd), 1);
% 
%     if (includeError)
%         clustersRef =zeros(numel(mmInd), 1);
%         clustersUnityMSE = zeros(numClusters, 1);
%         clustersPredMSE = zeros(numClusters, 1);
%         clustersUnityRelNormPercentError = zeros(numClusters, 1);
%         clustersPredRelNormPercentError = zeros(numClusters, 1);
%     else
%         clustersRef = [];
%         clustersUnityMSE = [];
%         clustersPredMSE = [];
%         clustersUnityRelNormPercentError = [];
%         clustersPredRelNormPercentError = [];
%     end
% 
%     for k = 1:numel(mmInd)
%         clustersTerms(k,:) = terms(mmInd(k), nnInd(k), :);
%         if (includeError)
%             clustersRef(k) = refZmn(mmInd(k), nnInd(k));
%         end
%     end
%     
%     totalProjSquaredError = 0;
%     totalPredSquaredError= 0;
%     for k = 1:numClusters
%         ind = find( indices(:,3) == k );
%         clustersPred(ind) = clustersTerms(ind, :) * clustersWeights(k, :)';
%         clustersUnity(ind) = sum(clustersTerms(ind, :), 2);
%         
%         if (includeError)
%             ref = clustersRef(ind);
%             
%             projMSE = (clustersUnity(ind) - ref)' * (clustersUnity(ind) - ref);
%             projRelNormPercentError = 100* sqrt(projMSE/ sum(ref.^2) ) ;
%             totalProjSquaredError = totalProjSquaredError + projMSE;
%             projMSE = projMSE /numel(ind);
%             
%             predMSE = (clustersPred(ind) - ref)' * (clustersPred(ind) - ref);
%             predRelNormPercentError = 100* sqrt(predMSE/ sum(ref.^2) ) ;
%             totalPredSquaredError = totalPredSquaredError + predMSE ;
%             predMSE = predMSE /numel(ind);
%             
%             clustersUnityMSE(k) = projMSE;
%             clustersPredMSE(k) = predMSE;
%             clustersUnityRelNormPercentError(k) = projRelNormPercentError;
%             clustersPredRelNormPercentError(k) = predRelNormPercentError;
%             
%         end
%         
%     end %for k = 1:numClusters
%     
%     totalUnityMSE =0;
%     totalPredMSE=0;
%     totalUnityRelNormPercentError = 0; 
%     totalPredRelNormPercentError = 0;
%     if (includeError)
%         totalUnityRelNormPercentError = 100 * sqrt(totalProjSquaredError /sum(clustersRef.^2) );
%         totalPredRelNormPercentError = 100 * sqrt(totalPredSquaredError /sum(clustersRef.^2) );
%         
%         totalUnityMSE = totalProjSquaredError ./ numel(mmInd);
%         totalPredMSE = totalPredSquaredError ./ numel(mmInd);
%     end

%end

