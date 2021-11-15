function [twoUniqueZmn, threeUniqueZmn, fourUniqueZmn] = extractRefZmn(zMatrices, singInd)

    [numEdges , ~, numFreq] = size(zMatrices);
    
    fourUniqueZmn = zeros(numEdges^2 - numEdges, numFreq);
    twoUniqueZmn = zeros(numEdges,  numFreq);
    threeUniqueZmn = zeros(numEdges * 8, numFreq  );
    
    nonSingCount = 0;
    triCount = 0;
    for mm = 1:numEdges
        for nn = 1:numEdges           
            if (mm == nn)
                twoUniqueZmn(mm, :) = zMatrices(mm,mm, : );
            elseif (singInd(mm,nn))
                triCount = triCount + 1;
                threeUniqueZmn(triCount, :) = zMatrices(mm,nn, : );
            else
                nonSingCount = nonSingCount + 1;               
                fourUniqueZmn(nonSingCount, :) = zMatrices(mm,nn, : );
            end
        end      
    end
    threeUniqueZmn = threeUniqueZmn(1:triCount, :);
    fourUniqueZmn = fourUniqueZmn(1:nonSingCount, :);
end