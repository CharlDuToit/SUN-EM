function [weights, pred ]= MLR(terms,ref, singThresh, addOnes, returnEmptyWeightsIfSingular)
  
    %clusterTerms =  nonSingZmnTerms(ind, :);
    [N,d] = size(terms);
    if (N < d + addOnes)
        pred = sum(terms, 2);
        weights = ones(d, 1);
        if (addOnes)
            weights(d +1,1) = 0;
        end
        return;
    end
    
    if (addOnes)
        terms(:,d+1) = ones(N,1);
    end
    a = terms' *terms;
    if (cond(a) * det(a) < singThresh) % -50 1.8752 982
        if (returnEmptyWeightsIfSingular)
            weights = [];
        else
            weights = ones(d, 1);
            if (addOnes)
                weights(d +1,1) = 0;
            end
        end
        
        if (addOnes)
            pred = sum(terms(:,1:d), 2);
        else
            pred = sum(terms, 2);
        end     

    else
        weights = ( a \ terms')* ref ;
        pred = terms * weights;
    end

end