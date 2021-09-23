%function [weights, pred ]= MLR(terms,ref, singThresh, addOnes, returnEmptyWeightsIfSingular)
%function [weights, pred,newTerms, predRelNormPercentError, unityRelNormPercentError ]= MLR(terms,ref, singThresh, addOnes, returnEmptyWeightsIfSingular)
function [weights, pred,newTerms ]= MLR(terms,ref, singThresh, addOnes, returnEmptyWeightsIfUnity)
  
    %clusterTerms =  nonSingZmnTerms(ind, :);
    unity = sum(terms, 2);
    
    [N,d] = size(terms);
    newTerms = zeros(N,d);
    
    %unityMSE = (unity - ref)' * (unity - ref) /N;
    %unityRelNormPercentError = 100* sqrt(((unity - ref)' * (unity - ref))/ sum(ref.^2) ) ;
    %predRelNormPercentError = unityRelNormPercentError;
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
    if (cond(a) * det(a) < singThresh) % -50 1.8752 982
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
        %predMSE = (pred - ref)' * (pred - ref) /N;
        %predRelNormPercentError = 100* sqrt(((pred - ref)' * (pred - ref))/ sum(ref.^2) ) ;
%         if (unityRelNormPercentError -predRelNormPercentError < 0 )
%             pred = unity;
%             weights = ones(d, 1);
%             if (addOnes)
%                 weights(d +1,1) = 0;
%             end
%         end
    end

end