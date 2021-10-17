%function [weights, pred ]= MLR(terms,ref, singThresh, addOnes, returnEmptyWeightsIfSingular)
%function [weights, pred,newTerms, predRelNormPercentError, unityRelNormPercentError ]= MLR(terms,ref, singThresh, addOnes, returnEmptyWeightsIfSingular)
%function [weights, pred,newTerms  ]= MLR(terms,ref, singThresh, addOnes, returnEmptyWeightsIfUnity)
function [weights, pred, unity, unityMSE, predMSE, unityRelNormPercentError, predRelNormPercentError] =...
    MLR(terms,ref, singThresh, addOnes, returnEmptyWeightsIfUnity)
  
    %clusterTerms =  nonSingZmnTerms(ind, :);
    unity = sum(terms, 2);
    
    [N,d] = size(terms);
    newTerms = zeros(N,d);
    
    [unityMSE, unityRelNormPercentError] = calcError(ref, unity);
    predRelNormPercentError = unityRelNormPercentError;
    predMSE = unityMSE;
    
    if (N < d + addOnes)
        if (returnEmptyWeightsIfUnity)
            pred =unity;
            weights = [];
        else
            weights = ones(d, 1);
            if (addOnes)
                weights(d +1,1) = 0;
            end
            w = sum(ref ./ unity )/numel(unity);
            weights = weights * w;
            pred = terms * weights;
            for k = 1:N
                newTerms(k,:) = terms(k,1:d) .* weights(1:d)';
            end
            [predMSE, predRelNormPercentError] = calcError(ref, pred);
        end
        
        return;
    end
    
    if (addOnes)
        terms(:,d+1) = ones(N,1);
    end
    a = terms' *terms;
    if (abs(cond(a) * det(a)) < singThresh)
        %pred = unity;
        if (returnEmptyWeightsIfUnity)
            pred = unity;
            weights = [];
        else
            weights = ones(d, 1);
            if (addOnes)
                weights(d +1,1) = 0;
            end
            w = sum(ref ./ unity )/numel(unity);
            weights = weights * w;
            pred = terms * weights;
            for k = 1:N
                newTerms(k,:) = terms(k,1:d) .* weights(1:d)';
            end
            [predMSE, predRelNormPercentError] = calcError(ref, pred);
        end
        
    else
        weights = ( a \ terms')* ref ;
        pred = terms * weights;
        for k = 1:N
            newTerms(k,:) = terms(k,1:d) .* weights(1:d)';
        end
        [predMSE, predRelNormPercentError] = calcError(ref, pred);

    end

end