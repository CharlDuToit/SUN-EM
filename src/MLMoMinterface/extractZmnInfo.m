function [allTerms, allProperties, selfZmnTerms, triZmnTerms, nonSingZmnTerms, nonSingZmnProp,nonSingEdgeLabels, singInd ] = extractZmnInfo(Const, Solver_setup)
    
    %[allTerms,singInd, dist,edge_mm_dir_dot_edge_nn_dir, edge_mm_dir_dot_edge_nn_disp] = fillZmnTermsByEdge(Const,Solver_setup);
    [allTerms,singInd, allProperties] = fillZmnTermsByEdge(Const,Solver_setup);
    [numEdges, ~, numTerms, numFreq] = size(allTerms);
    
    nonSingZmnTerms = zeros(numEdges^2 - numEdges, numTerms, numFreq);
    selfZmnTerms = zeros(numEdges , numTerms , numFreq);
    triZmnTerms = zeros(numEdges * 8,numTerms,numFreq  );
    
    % geometric properties not dependant on frequency
    nonSingZmnProp = zeros(numEdges^2 - numEdges, 3);
    nonSingEdgeLabels = zeros(numEdges^2 - numEdges, 2);
    
    nonSingCount = 0;
    triCount = 0;
    for mm = 1:numEdges
        for nn = 1:numEdges
            if (mm == nn)
                selfZmnTerms(mm, :, :) = allTerms(mm,mm, :, : );
            elseif (singInd(mm,nn))
                triCount = triCount + 1;
                triZmnTerms(triCount, :, :) = allTerms(mm,nn, :, : );
            else
                nonSingCount = nonSingCount + 1;
%                 nonSingZmnProp(nonSingCount, 1) = dist(mm,nn);
%                 nonSingZmnProp(nonSingCount, 2) = edge_mm_dir_dot_edge_nn_dir(mm,nn);
%                 nonSingZmnProp(nonSingCount, 3) = edge_mm_dir_dot_edge_nn_disp(mm,nn);
                
                nonSingZmnProp(nonSingCount, 1) = allProperties(mm,nn,1);
                nonSingZmnProp(nonSingCount, 2) = allProperties(mm,nn,2);
                nonSingZmnProp(nonSingCount, 3) = allProperties(mm,nn,3);
                
                nonSingZmnTerms(nonSingCount, :, :) = allTerms(mm,nn, :, : );
                
                nonSingEdgeLabels(nonSingCount, :) = [mm,nn];
            end
        end % for nn         
    end % for mm
    nonSingZmnTerms = nonSingZmnTerms(1:nonSingCount, :, :);
    nonSingZmnProp = nonSingZmnProp(1:nonSingCount, :, :);
    nonSingEdgeLabels = nonSingEdgeLabels(1:nonSingCount, : );
    triZmnTerms = triZmnTerms(1:triCount, :, :);
    
  
end