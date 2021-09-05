function [nonSingZmn] = extractNonSingZmn(refZmn, singIndices)
    [numEdges , ~] = size(refZmn);
    numObs = numEdges^2 - numEdges;
    %nonSingZmn = zeros(numObs, 2);
    nonSingZmn = zeros(numObs, 1);
    i = 0;
    for mm = 1:numEdges
        for nn = 1:numEdges
            if (mm == nn)
                continue
            end
            if (singIndices(mm,nn))
                continue
            end
            i = i + 1;
            %nonSingZmn(i, 1) = real(refZmn(mm,nn));
            %nonSingZmn(i, 2) = imag(refZmn(mm,nn));
            nonSingZmn(i, 1) = imag(refZmn(mm,nn)); 
        end      
    end
    nonSingZmn = nonSingZmn(1:i, :);
end