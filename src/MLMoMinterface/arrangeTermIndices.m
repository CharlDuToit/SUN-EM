function [indices] = arrangeTermIndices(mmPlusCentre, mmMinusCentre, nnPlusCentre, nnMinusCentre)
    % Arrange so that plus centres are nearest
    % Centres are [1, 3] xyz
    indices = linspace(1,8,8);
    %signs = ones(8);
    
    mmPlus_nnPlus = norm(mmPlusCentre-nnPlusCentre);
    mmPlus_nnMinus = norm(mmPlusCentre-nnMinusCentre);
    mmMinus_nnPlus= norm(mmMinusCentre-nnPlusCentre);
    mmMinus_nnMinus = norm(mmMinusCentre-nnMinusCentre);
    
    %terms(:,:,1) = A_m_pls_n_pls
    %terms(:,:,2) = Phi_m_pls_n_pls
    %terms(:,:,3) = A_m_pls_n_mns
    %terms(:,:,4) = Phi_m_pls_n_mns
    %terms(:,:,5) = A_m_mns_n_pls
    %terms(:,:,6) = Phi_m_mns_n_pls
    %terms(:,:,7) = A_m_mns_n_mns
    %terms(:,:,8) = Phi_m_mns_n_mns
    
    [~, ind] = sort([mmPlus_nnPlus,mmPlus_nnMinus,mmMinus_nnPlus,mmMinus_nnMinus]);
    %sorted from smallest to largest
    % min
    switch ind(1)
        case 1 % mmPlus_nnPlus
        case 2  %mmPlus_nnMinus
            indices(1:8) = [3,4, 1, 2, 7, 8, 5, 6];
        case 3 % mmMinus_nnPlus
            indices(1:8) = [5,6, 7,8, 1, 2, 3, 4];
        case 4 % mmMinus_nnMinus
            indices(1:8) = [1,2,5,6,3,4,7,8];
    end
end