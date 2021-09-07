function [selfZmn, triZmn, nonSingZmn] = extractRefZmn(zMatrices, singInd)

    [numEdges , ~, numFreq] = size(zMatrices);
    
    nonSingZmn = zeros(numEdges^2 - numEdges, numFreq);
    selfZmn = zeros(numEdges,  numFreq);
    triZmn = zeros(numEdges * 8, numFreq  );
    
    nonSingCount = 0;
    triCount = 0;
    for mm = 1:numEdges
        for nn = 1:numEdges           
            if (mm == nn)
                selfZmn(mm, :) = zMatrices(mm,mm, : );
            elseif (singInd(mm,nn))
                triCount = triCount + 1;
                triZmn(triCount, :) = zMatrices(mm,nn, : );
            else
                nonSingCount = nonSingCount + 1;               
                nonSingZmn(nonSingCount, :) = zMatrices(mm,nn, : );
            end
        end      
    end
    triZmn = triZmn(1:triCount, :);
    nonSingZmn = nonSingZmn(1:nonSingCount, :);
end