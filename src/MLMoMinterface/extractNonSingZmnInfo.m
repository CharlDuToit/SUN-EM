function [nonSingZmnTerms, nonSingZmnProp, singInd] = extractNonSingZmnInfo(Const, Solver_setup)
    %x : [N , d]
    %t : [N ,2 ]
    [Phi_mn_mns_source_pls_matrix, Phi_mn_mns_source_mns_matrix,...
    Phi_mn_pls_source_pls_matrix, Phi_mn_pls_source_mns_matrix,...
    Amn_pls_source_pls_matrix, Amn_pls_source_mns_matrix,...
    Amn_mns_source_pls_matrix, Amn_mns_source_mns_matrix, ...
    dist,edge_mm_dir_dot_edge_nn_dir, edge_mm_dir_dot_edge_nn_disp, singInd] = fillNonSingZmnTermsByEdge(Const,Solver_setup);

    numEdges =Solver_setup.num_mom_basis_functions;
   % [numObs, ~] = size(singIndices);
    %inputTrainingData = zeros(numObs, 16);
    nonSingZmnTerms = zeros(numEdges^2 - numEdges, 8);
    nonSingZmnProp = zeros(numEdges^2 - numEdges, 3);
    i = 0;
    for mm = 1:numEdges
        for nn = 1:numEdges
            if (mm == nn)
                continue
            end
            if (singInd(mm,nn))
                continue
            end
            i = i + 1;
          
%             inputTrainingData(i, 1) = real(Phi_mn_mns_source_pls_matrix(mm,nn));
%             inputTrainingData(i, 2) = real(Phi_mn_mns_source_mns_matrix(mm,nn));
%             inputTrainingData(i, 3) = real(Phi_mn_pls_source_pls_matrix(mm,nn));
%             inputTrainingData(i, 4) = real(Phi_mn_pls_source_mns_matrix(mm,nn));
%             inputTrainingData(i, 5) = real(Amn_pls_source_pls_matrix(mm,nn));
%             inputTrainingData(i, 6) = real(Amn_pls_source_mns_matrix(mm,nn));
%             inputTrainingData(i, 7) = real(Amn_mns_source_pls_matrix(mm,nn));
%             inputTrainingData(i, 8) = real(Amn_mns_source_mns_matrix(mm,nn));
%             
%             inputTrainingData(i, 9) = imag(Phi_mn_mns_source_pls_matrix(mm,nn));
%             inputTrainingData(i, 10) = imag(Phi_mn_mns_source_mns_matrix(mm,nn));
%             inputTrainingData(i, 11) = imag(Phi_mn_pls_source_pls_matrix(mm,nn));
%             inputTrainingData(i, 12) = imag(Phi_mn_pls_source_mns_matrix(mm,nn));
%             inputTrainingData(i, 13) = imag(Amn_pls_source_pls_matrix(mm,nn));
%             inputTrainingData(i, 14) = imag(Amn_pls_source_mns_matrix(mm,nn));
%             inputTrainingData(i, 15) = imag(Amn_mns_source_pls_matrix(mm,nn));
%             inputTrainingData(i, 16) = imag(Amn_mns_source_mns_matrix(mm,nn));
            nonSingZmnProp(i, 1) = dist(mm,nn);
            nonSingZmnProp(i, 2) = edge_mm_dir_dot_edge_nn_dir(mm,nn);
            nonSingZmnProp(i, 3) = edge_mm_dir_dot_edge_nn_disp(mm,nn);

            nonSingZmnTerms(i, 1) = imag(Phi_mn_mns_source_pls_matrix(mm,nn));
            nonSingZmnTerms(i, 2) = imag(Phi_mn_mns_source_mns_matrix(mm,nn));
            nonSingZmnTerms(i, 3) = imag(Phi_mn_pls_source_pls_matrix(mm,nn));
            nonSingZmnTerms(i, 4) = imag(Phi_mn_pls_source_mns_matrix(mm,nn));
            nonSingZmnTerms(i, 5) = imag(Amn_pls_source_pls_matrix(mm,nn));
            nonSingZmnTerms(i, 6) = imag(Amn_pls_source_mns_matrix(mm,nn));
            nonSingZmnTerms(i, 7) = imag(Amn_mns_source_pls_matrix(mm,nn));
            nonSingZmnTerms(i, 8) = imag(Amn_mns_source_mns_matrix(mm,nn));
        end % for
        
    end % for
    nonSingZmnTerms = nonSingZmnTerms(1:i, :);
    nonSingZmnProp = nonSingZmnProp(1:i, :);
  
end