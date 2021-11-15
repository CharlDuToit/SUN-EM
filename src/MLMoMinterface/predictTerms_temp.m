%function [predZmn, unityZmn] = predictExtractedTerms( mlmom, selfZmnTerms, triZmnTerms, nonSingZmnTerms, nonSingZmnProp, singInd,edgeLengths, returnUnity)
function [predZmn, indicesStruct, assignTime, multiplyTime] = predictTerms_temp( weightModel, clusterStruct,propStruct, termStruct, indicesStruct, hasCalcInd)

    %numClusters = numel(mlmom.clusterMeans(:,1));
    %maxClusterError = 10*ones(numClusters,1);
    %includeRealCalc = mlmom.includeRealCalc;

    %==========================
    numEdges = termStruct.numEdges;
    predZmn = zeros(numEdges,numEdges);
    assignTime = 0;
    multiplyTime = 0;
    

    %=====2 unique
    if ( clusterStruct.twoUniqueErrorCode ~= 0)
        tic;
        means = clusterStruct.twoUniqueMeans;
        twoUniqueMaxError = clusterStruct.twoUniqueMaxError;
        if (~hasCalcInd)
            [~, indicesStruct.twoUniqueIndices(:,3), ~,~, ~, indicesStruct.twoUniqueIndices(:,4)] = assignClusters(propStruct.twoUniqueProp, means,...
                clusterStruct.twoUniqueErrorCode);
        end
        assignTime = toc;
        tic;
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
        multiplyTime = toc;
    else
        tic;
        twoUniqueIndices = indicesStruct.twoUniqueIndices;
        for k = 1:numel(twoUniqueIndices(:,1))
            predZmn(twoUniqueIndices(k,1), twoUniqueIndices(k,2)) = termStruct.twoUniqueTerms(k,:) * weightModel.twoUniqueWeights;
            %predZmn(twoUniqueIndices(k,1), twoUniqueIndices(k,2)) = sum(termStruct.twoUniqueTerms(k,:));
        end
        multiplyTime = toc;
    end
    %=====3 unique
    if ( clusterStruct.threeUniqueErrorCode ~= 0)
        tic;
        means = clusterStruct.threeUniqueMeans;
        threeUniqueMaxError = clusterStruct.threeUniqueMaxError;
        if (~hasCalcInd)
            [~, indicesStruct.threeUniqueIndices(:,3), ~,~, ~, indicesStruct.threeUniqueIndices(:,4)] = assignClusters(propStruct.threeUniqueProp, means,...
                clusterStruct.threeUniqueErrorCode);
        end
        assignTime = assignTime + toc;
        tic;
        threeUniqueIndices = indicesStruct.threeUniqueIndices;
        for k = 1:numel(threeUniqueIndices(:,1))
            if (indicesStruct.threeUniqueIndices(k,4) <= threeUniqueMaxError(threeUniqueIndices(k,3) ) )
                if (true) %pos
                %if (propStruct.threeUniqueSwapGroups(k) == 0) %pos
                    predZmn(threeUniqueIndices(k,1), threeUniqueIndices(k,2)) = termStruct.threeUniqueTerms(k,:) * weightModel.threeUniqueWeightsPos(threeUniqueIndices(k,3), :)';
                else
                    predZmn(threeUniqueIndices(k,1), threeUniqueIndices(k,2)) = termStruct.threeUniqueTerms(k,:) * weightModel.threeUniqueWeightsNeg(threeUniqueIndices(k,3), :)';
                end
            else
                predZmn(threeUniqueIndices(k,1), threeUniqueIndices(k,2)) = sum(termStruct.threeUniqueTerms(k,:));
            end
        end
        multiplyTime = multiplyTime + toc;
    else
        tic;
        threeUniqueIndices = indicesStruct.threeUniqueIndices;
        for k = 1:numel(threeUniqueIndices(:,1))
            if (true)
            %if (propStruct.threeUniqueSwapGroups(k) == 0) %pos
                predZmn(threeUniqueIndices(k,1), threeUniqueIndices(k,2)) = termStruct.threeUniqueTerms(k,:) * weightModel.threeUniqueWeightsPos;
                %predZmn(threeUniqueIndices(k,1), threeUniqueIndices(k,2)) = sum(termStruct.threeUniqueTerms(k,:));
            else
                predZmn(threeUniqueIndices(k,1), threeUniqueIndices(k,2)) = termStruct.threeUniqueTerms(k,:) * weightModel.threeUniqueWeightsNeg;
            end
            
            %predZmn(threeUniqueIndices(k,1), threeUniqueIndices(k,2)) = sum(termStruct.threeUniqueTerms(k,:));
        end
        multiplyTime = multiplyTime + toc;
    end
    %=====4 unique
    tic;
    means = clusterStruct.fourUniqueMeans;
    fourUniqueMaxErrors = clusterStruct.fourUniqueMaxError;
    if (~hasCalcInd)
        [~, indicesStruct.fourUniqueIndices(:,3), ~,~, ~, indicesStruct.fourUniqueIndices(:,4)] = assignClusters(propStruct.fourUniqueProp, means,...
            clusterStruct.fourUniqueErrorCode);
    end
    assignTime = assignTime + toc;
    tic;
    fourUniqueIndices = indicesStruct.fourUniqueIndices;
    for k = 1:numel(fourUniqueIndices(:,1))
        if ( fourUniqueIndices(k,4) <= fourUniqueMaxErrors(fourUniqueIndices(k,3)) )
            if (propStruct.fourUniqueSwapGroups(k) == 0) %pos
                predZmn(fourUniqueIndices(k,1), fourUniqueIndices(k,2)) = termStruct.fourUniqueTerms(k,:) * weightModel.fourUniqueWeightsPos(fourUniqueIndices(k,3), :)';
            else
                predZmn(fourUniqueIndices(k,1), fourUniqueIndices(k,2)) = termStruct.fourUniqueTerms(k,:) * weightModel.fourUniqueWeightsNeg(fourUniqueIndices(k,3), :)';
            end
            
            %predZmn(fourUniqueIndices(k,1), fourUniqueIndices(k,2)) = termStruct.fourUniqueTerms(k,:) * weightModel.fourUniqueWeights(fourUniqueIndices(k,3), :)';
            %predZmn(fourUniqueIndices(k,1), fourUniqueIndices(k,2)) = sum(termStruct.fourUniqueTerms(k,:));
        else
            predZmn(fourUniqueIndices(k,1), fourUniqueIndices(k,2)) = sum(termStruct.fourUniqueTerms(k,:));
        end
    end
    multiplyTime = multiplyTime + toc;

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

