%function [weights, pred ]= MLR(terms,ref, singThresh, addOnes, returnEmptyWeightsIfSingular)
%function [weights, pred,newTerms, predRelNormPercentError, unityRelNormPercentError ]= MLR(terms,ref, singThresh, addOnes, returnEmptyWeightsIfSingular)
%function [weights, pred,newTerms  ]= MLR(terms,ref, singThresh, addOnes, returnEmptyWeightsIfUnity)
function [weights, pred, unity, unityMSE, predMSE, unityRelNormPercentError, predRelNormPercentError] =...
    MLR(terms,ref, singThresh, addOnes, returnEmptyWeightsIfUnity)
  
    %clusterTerms =  nonSingZmnTerms(ind, :);
    unity = sum(terms, 2);
    
    [N,d] = size(terms);
    newTerms = zeros(N,d);
    
    unityMSE = (unity - ref)' * (unity - ref);

    unityRelNormPercentError = 100* sqrt(unityMSE/ sum(ref.^2) ) ;
    predRelNormPercentError = unityRelNormPercentError;
    unityMSE = unityMSE /N;
    predMSE = unityMSE;
    
    if (N < d + addOnes)
        %pred = sum(terms, 2);
        pred =unity;
        if (returnEmptyWeightsIfUnity)
            weights = [];
        else
            weights = ones(d, 1);
            if (addOnes)
                weights(d +1,1) = 0;
            end
        end
        
        return;
    end
    
    if (addOnes)
        terms(:,d+1) = ones(N,1);
    end
    a = terms' *terms;
    if (abs(cond(a) * det(a)) < singThresh) % -50 1.8752 982
        pred = unity;
        if (returnEmptyWeightsIfUnity)
            weights = [];
        else
            weights = ones(d, 1);
            if (addOnes)
                weights(d +1,1) = 0;
            end
        end
        

    else
        weights = ( a \ terms')* ref ;
        pred = terms * weights;
        for k = 1:N
            newTerms(k,:) = terms(k,1:d) .* weights(1:d)';
        end
        predMSE = (pred - ref)' * (pred - ref);
        predRelNormPercentError = 100* sqrt(predMSE/ sum(ref.^2) ) ;
        predMSE = predMSE /N;

    end

end