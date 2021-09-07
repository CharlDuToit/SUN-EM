function [selfZmnTerms, triZmnTerms, nonSingZmnTerms, nonSingZmnProp, singInd ] = extractZmnInfo(Const, Solver_setup)

    [terms,singInd, dist,edge_mm_dir_dot_edge_nn_dir, edge_mm_dir_dot_edge_nn_disp] = fillZmnTermsByEdge(Const,Solver_setup);
    [numEdges, ~, numTerms, numFreq] = size(terms);
    
    nonSingZmnTerms = zeros(numEdges^2 - numEdges, numTerms, numFreq);
    selfZmnTerms = zeros(numEdges , numTerms , numFreq);
    triZmnTerms = zeros(numEdges * 8,numTerms,numFreq  );
    
    % geometric properties not dependant on frequency
    nonSingZmnProp = zeros(numEdges^2 - numEdges, 3);
    
    nonSingCount = 0;
    triCount = 0;
    for mm = 1:numEdges
        for nn = 1:numEdges
            if (mm == nn)
                selfZmnTerms(mm, :, :) = terms(mm,mm, :, : );
            elseif (singInd(mm,nn))
                triCount = triCount + 1;
                triZmnTerms(triCount, :, :) = terms(mm,nn, :, : );
            else
                nonSingCount = nonSingCount + 1;
                nonSingZmnProp(nonSingCount, 1) = dist(mm,nn);
                nonSingZmnProp(nonSingCount, 2) = edge_mm_dir_dot_edge_nn_dir(mm,nn);
                nonSingZmnProp(nonSingCount, 3) = edge_mm_dir_dot_edge_nn_disp(mm,nn);
                
                nonSingZmnTerms(nonSingCount, :, :) = terms(mm,nn, :, : );
            end
        end % for nn         
    end % for mm
    nonSingZmnTerms = nonSingZmnTerms(1:nonSingCount, :, :);
    nonSingZmnProp = nonSingZmnProp(1:nonSingCount, :, :);
    triZmnTerms = triZmnTerms(1:triCount, :, :);
  
end