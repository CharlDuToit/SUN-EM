function [ZmnRows] = extractZmnRows(Zmn, ind)
    % ind is the same size as Zmn
    [numEdges , ~] = size(Zmn);
    numObs = numEdges^2 - numEdges;
    %nonSingZmn = zeros(numObs, 2);
    ZmnRows = zeros(numObs, 1);
    i = 0;
    for mm = 1:numEdges
        for nn = 1:numEdges
            if (mm == nn)
                continue
            end
            if (ind(mm,nn))
                continue
            end
            i = i + 1;
            %nonSingZmn(i, 1) = real(refZmn(mm,nn));
            %nonSingZmn(i, 2) = imag(refZmn(mm,nn));
            ZmnRows(i, 1) = imag(Zmn(mm,nn)); 
        end      
    end
    ZmnRows = ZmnRows(1:i, :);
end