%function [allTerms, allEdgeCentreProperties, twoUniqueTerms, threeUniqueTerms, fourUniqueTerms, fourUniqueProp,fourUniqueEdgeLengths, threeUniqueEdgeLengths, singInd ] = extractZmnInfo(Const, Solver_setup,constMeshSize)
function [termStruct, propStruct,indicesStruct, singInd ] = extractZmnInfo(Const, Solver_setup,constMeshSize)
    
    tic;
    centreDistances= createCentreDistances(Solver_setup );
    centreDistanceTime = toc;
    tic;
    [allTerms,singInd] = fillZmnTermsByEdge(Const,Solver_setup, centreDistances);%9.3, 9.2 9.5 9.6 9.37 9.98 9.27
    %[allTerms,singInd] = fillZmnTermsByEdge(Const,Solver_setup);% 9.5 10.1 10 9.9 9.4
    calcTermTime = toc;
    tic;
    allEdgeCentreProperties = createEdgeCentreProperties(Solver_setup); %0.31
    [numEdges, ~, numTerms, numFreq] = size(allTerms);
    edgeCentreTime = toc;
    tic;
    %centreDistances= createCentreDistances(Solver_setup );
    
    fourUniqueTerms = zeros(numEdges^2 - numEdges, numTerms, numFreq);
    twoUniqueTerms = zeros(numEdges , numTerms , numFreq);
    threeUniqueTerms = zeros(numEdges * 8,numTerms,numFreq  );
    
    % geometric properties not dependant on frequency
    fourUniqueProp = zeros(numEdges^2 - numEdges, 3);
    fourUniqueEdgeLengths = zeros(numEdges^2 - numEdges, 2);
    threeUniqueEdgeLengths = zeros(numEdges * 8, 2);
    
    %indices
    twoUniqueIndices = zeros(numEdges, 2);
    threeUniqueIndices = zeros(numEdges * 8, 2);
    fourUniqueIndices = zeros(numEdges^2 - numEdges, 2);
    
    %swap group
    threeUniqueSwapGroups = zeros(numEdges * 8, 1);
    fourUniqueSwapGroups = zeros(numEdges^2 - numEdges, 1);
    
    fourUniqueCount = 0;
    threeUniqueCount = 0;
    for mm = 1:numEdges
        for nn = 1:numEdges
            if (mm == nn)
                twoUniqueTerms(mm, :, :) = allTerms(mm,mm, :, : );
                twoUniqueIndices(mm,:) = [mm ,mm];
            elseif (singInd(mm,nn))
%                 [indices] = arrangeTermIndices(Solver_setup.triangle_centre_point(Solver_setup.rwg_basis_functions_trianglePlus(mm),:), ...
%                     Solver_setup.triangle_centre_point(Solver_setup.rwg_basis_functions_triangleMinus(mm),:), ...
%                     Solver_setup.triangle_centre_point(Solver_setup.rwg_basis_functions_trianglePlus(nn),:), ...
%                     Solver_setup.triangle_centre_point(Solver_setup.rwg_basis_functions_triangleMinus(nn),:));
%                  allTerms(mm,nn,:,:) = allTerms(mm,nn,indices,:);

                 [allTerms(mm,nn,:,:), swapGroup] = arrangeTerms(allTerms(mm,nn,:,:), Solver_setup.rwg_basis_functions_trianglePlus(mm),...
                     Solver_setup.rwg_basis_functions_triangleMinus(mm),...
                     Solver_setup.rwg_basis_functions_trianglePlus(nn),...
                     Solver_setup.rwg_basis_functions_triangleMinus(nn),...
                     centreDistances );
%                  
%                  [allTerms(mm,nn,:,:)] = arrangeTerms(allTerms(mm,nn,:,:),...
%                      centreDistances(Solver_setup.rwg_basis_functions_trianglePlus(mm), Solver_setup.rwg_basis_functions_trianglePlus(nn)),...
%                      centreDistances(Solver_setup.rwg_basis_functions_trianglePlus(mm), Solver_setup.rwg_basis_functions_triangleMinus(nn) ),...
%                      centreDistances(Solver_setup.rwg_basis_functions_triangleMinus(mm), Solver_setup.rwg_basis_functions_trianglePlus(nn)),...
%                      centreDistances(Solver_setup.rwg_basis_functions_triangleMinus(mm), Solver_setup.rwg_basis_functions_triangleMinus(nn)) );
                 
                threeUniqueCount = threeUniqueCount + 1;
                threeUniqueSwapGroups(threeUniqueCount, 1) = swapGroup;
                threeUniqueTerms(threeUniqueCount, :, :) = allTerms(mm,nn, :, : );
                threeUniqueEdgeLengths(threeUniqueCount, :) = [Solver_setup.rwg_basis_functions_length_m(mm),...
                    Solver_setup.rwg_basis_functions_length_m(nn)];
                threeUniqueIndices(threeUniqueCount,:) = [mm ,nn];
            else
%                 [indices] = arrangeTermIndices(Solver_setup.triangle_centre_point(Solver_setup.rwg_basis_functions_trianglePlus(mm),:), ...
%                     Solver_setup.triangle_centre_point(Solver_setup.rwg_basis_functions_triangleMinus(mm),:), ...
%                     Solver_setup.triangle_centre_point(Solver_setup.rwg_basis_functions_trianglePlus(nn),:), ...
%                     Solver_setup.triangle_centre_point(Solver_setup.rwg_basis_functions_triangleMinus(nn),:));
%                 allTerms(mm,nn,:,:) = allTerms(mm,nn,indices,:);
                
                [allTerms(mm,nn,:,:), swapGroup] = arrangeTerms(allTerms(mm,nn,:,:), Solver_setup.rwg_basis_functions_trianglePlus(mm),...
                    Solver_setup.rwg_basis_functions_triangleMinus(mm),...
                    Solver_setup.rwg_basis_functions_trianglePlus(nn),...
                    Solver_setup.rwg_basis_functions_triangleMinus(nn),...
                    centreDistances );

                 
%                  [allTerms(mm,nn,:,:)] = arrangeTerms(allTerms(mm,nn,:,:),...
%                      centreDistances(Solver_setup.rwg_basis_functions_trianglePlus(mm), Solver_setup.rwg_basis_functions_trianglePlus(nn)),...
%                      centreDistances(Solver_setup.rwg_basis_functions_trianglePlus(mm), Solver_setup.rwg_basis_functions_triangleMinus(nn) ),...
%                      centreDistances(Solver_setup.rwg_basis_functions_triangleMinus(mm), Solver_setup.rwg_basis_functions_trianglePlus(nn)),...
%                      centreDistances(Solver_setup.rwg_basis_functions_triangleMinus(mm), Solver_setup.rwg_basis_functions_triangleMinus(nn)) );
                
                fourUniqueCount = fourUniqueCount + 1;
                
                fourUniqueSwapGroups(fourUniqueCount, 1) = swapGroup;
                
                fourUniqueProp(fourUniqueCount, 1) = allEdgeCentreProperties(mm,nn,1);
                fourUniqueProp(fourUniqueCount, 2) = allEdgeCentreProperties(mm,nn,2);
                fourUniqueProp(fourUniqueCount, 3) = allEdgeCentreProperties(mm,nn,3);
                
                fourUniqueTerms(fourUniqueCount, :, :) = allTerms(mm,nn, :, : );
                
                fourUniqueEdgeLengths(fourUniqueCount, :) = [Solver_setup.rwg_basis_functions_length_m(mm),...
                    Solver_setup.rwg_basis_functions_length_m(nn)];
                fourUniqueIndices(fourUniqueCount,:) = [mm ,nn];
            end
        end % for nn         
    end % for mm
    assignAndSwapTime = toc;
    
    fourUniqueTerms = fourUniqueTerms(1:fourUniqueCount, :, :);
    fourUniqueProp = fourUniqueProp(1:fourUniqueCount, :, :);
    fourUniqueEdgeLengths = fourUniqueEdgeLengths(1:fourUniqueCount, : );
    fourUniqueIndices = fourUniqueIndices(1:fourUniqueCount,:);
    fourUniqueSwapGroups = fourUniqueSwapGroups(1:fourUniqueCount,:);
     
    threeUniqueTerms = threeUniqueTerms(1:threeUniqueCount, :, :);
    threeUniqueEdgeLengths = threeUniqueEdgeLengths(1:threeUniqueCount, :);
    threeUniqueIndices = threeUniqueIndices(1:threeUniqueCount,:);
    threeUniqueSwapGroups = threeUniqueSwapGroups(1:threeUniqueCount,:);

    twoUniqueProp = [];
    threeUniqueProp = [];
    twoUniqueErrorCode = 0;
    threeUniqueErrorCode = 0;
    fourUniqueErrorCode = 1;
    %property allocation, errorCode must match calcMinClusterError 
    if (~constMeshSize)
        twoUniqueProp = Solver_setup.rwg_basis_functions_length_m;
        threeUniqueProp = threeUniqueEdgeLengths;
        fourUniqueProp(:, 4:5) = fourUniqueEdgeLengths;
        twoUniqueErrorCode = 5;
        threeUniqueErrorCode = 12;
        fourUniqueErrorCode = 11;
    end
   
    termStruct = [];
    termStruct.calcTermTime = calcTermTime;
    termStruct.assignAndSwapTime = assignAndSwapTime;
    termStruct.allTerms = allTerms;
    termStruct.numEdges = numEdges;
    termStruct.twoUniqueTerms = twoUniqueTerms;
    termStruct.threeUniqueTerms = threeUniqueTerms;
    termStruct.fourUniqueTerms = fourUniqueTerms;
    
    propStruct = [];
    %propStruct.constMeshSize = constMeshSize;
   % propStruct.allEdgeCentreProperties = allEdgeCentreProperties; 
    propStruct.centreDistanceTime = centreDistanceTime;
    propStruct.edgeCentreTime = edgeCentreTime;
    propStruct.twoUniqueProp = twoUniqueProp;
    propStruct.twoUniqueErrorCode = twoUniqueErrorCode;
    
    propStruct.threeUniqueProp = threeUniqueProp;
    propStruct.threeUniqueErrorCode = threeUniqueErrorCode;
    propStruct.threeUniqueSwapGroups = threeUniqueSwapGroups;
    
    propStruct.fourUniqueProp = fourUniqueProp;
    propStruct.fourUniqueErrorCode = fourUniqueErrorCode;
    propStruct.fourUniqueSwapGroups = fourUniqueSwapGroups;
    
    indicesStruct = [];
    indicesStruct.twoUniqueIndices = twoUniqueIndices;
    indicesStruct.threeUniqueIndices = threeUniqueIndices;
    indicesStruct.fourUniqueIndices = fourUniqueIndices;
    

end

% function [indices] = arrangeTermIndices(mmPlusCentre, mmMinusCentre, nnPlusCentre, nnMinusCentre)
%     % Arrange so that plus centres are nearest
%     % Centres are [1, 3] xyz
%     indices = linspace(1,8,8);
%     %signs = ones(8);
%     
%     mmPlus_nnPlus = norm(mmPlusCentre-nnPlusCentre);
%     mmPlus_nnMinus = norm(mmPlusCentre-nnMinusCentre);
%     mmMinus_nnPlus= norm(mmMinusCentre-nnPlusCentre);
%     mmMinus_nnMinus = norm(mmMinusCentre-nnMinusCentre);
%     
%     %terms(:,:,1) = A_m_pls_n_pls
%     %terms(:,:,2) = Phi_m_pls_n_pls
%     %terms(:,:,3) = A_m_pls_n_mns
%     %terms(:,:,4) = Phi_m_pls_n_mns
%     %terms(:,:,5) = A_m_mns_n_pls
%     %terms(:,:,6) = Phi_m_mns_n_pls
%     %terms(:,:,7) = A_m_mns_n_mns
%     %terms(:,:,8) = Phi_m_mns_n_mns
%     
%     [~, ind] = sort([mmPlus_nnPlus,mmPlus_nnMinus,mmMinus_nnPlus,mmMinus_nnMinus]);
%     %sorted from smallest to largest
%     % min
%     switch ind(1)
%         case 1 % mmPlus_nnPlus
%         case 2  %mmPlus_nnMinus
%             indices(1:8) = [3,4, 1, 2, 7, 8, 5, 6];
%         case 3 % mmMinus_nnPlus
%             indices(1:8) = [5,6, 7,8, 1, 2, 3, 4];
%         case 4 % mmMinus_nnMinus
%             %indices(1:8) = [1,2,5,6,3,4,7,8];
%             indices(1:8) = [7,8,5,6,3,4,1,2];
%     end
% end

function [terms, swapGroup] = arrangeTerms(terms, mmPlus,mmMinus,nnPlus,nnMinus, centreDistances)
%function [terms] = arrangeTerms(terms, mmPlus_nnPlus,mmPlus_nnMinus,mmMinus_nnPlus,mmMinus_nnMinus)
    % Arrange so that plus centres are nearest
    % Centres are [1, 3] xyz
    %indices = linspace(1,8,8);
    %signs = ones(8);
    
    % swapType :0, ++ or -- closest
    %          :1, +- or -+ closest
    
    %mmPlus_nnPlus = norm(mmPlusCentre-nnPlusCentre);
   % mmPlus_nnMinus = norm(mmPlusCentre-nnMinusCentre);
    %mmMinus_nnPlus= norm(mmMinusCentre-nnPlusCentre);
    %mmMinus_nnMinus = norm(mmMinusCentre-nnMinusCentre);
    
    %terms(:,:,1) = A_m_pls_n_pls
    %terms(:,:,2) = Phi_m_pls_n_pls
    %terms(:,:,3) = A_m_pls_n_mns
    %terms(:,:,4) = Phi_m_pls_n_mns
    %terms(:,:,5) = A_m_mns_n_pls
    %terms(:,:,6) = Phi_m_mns_n_pls
    %terms(:,:,7) = A_m_mns_n_mns
    %terms(:,:,8) = Phi_m_mns_n_mns
    
%     [indices] = arrangeTermIndices(Solver_setup.triangle_centre_point(Solver_setup.rwg_basis_functions_trianglePlus(mm),:), ...
%         Solver_setup.triangle_centre_point(Solver_setup.rwg_basis_functions_triangleMinus(mm),:), ...
%         Solver_setup.triangle_centre_point(Solver_setup.rwg_basis_functions_trianglePlus(nn),:), ...
%         Solver_setup.triangle_centre_point(Solver_setup.rwg_basis_functions_triangleMinus(nn),:));
    
    %[~, ind] = sort([mmPlus_nnPlus,mmPlus_nnMinus,mmMinus_nnPlus,mmMinus_nnMinus]);
    [~, ind] = sort([centreDistances(mmPlus,nnPlus),centreDistances(mmPlus,nnMinus),centreDistances(mmMinus,nnPlus),centreDistances(mmMinus,nnMinus)]);
    %sorted from smallest to largest
    % min
    swapGroup = 0;
    switch ind(1)
        case 1 % mmPlus_nnPlus
            swapGroup = 0;
        case 2  %mmPlus_nnMinus
            terms(1,1,:,:) = terms(1,1,[3,4, 1, 2, 7, 8, 5, 6],:);
            swapGroup = 1;
        case 3 % mmMinus_nnPlus
            terms(1,1,:,:) = terms(1,1,[5,6, 7,8, 1, 2, 3, 4],:);
            swapGroup = 1;
        case 4 % mmMinus_nnMinus
            %terms(1,1,:,:) = terms(1,1,[1,2,5,6,3,4,7,8],:);
            terms(1,1,:,:) = terms(1,1,[7,8,5,6,3,4,1,2],:);
            swapGroup = 0;
    end
end

function centreDistances= createCentreDistances(Solver_setup )
    
    %numEdges = Solver_setup.num_mom_basis_functions;
    numTri = Solver_setup.num_metallic_triangles;
    centreDistances = zeros(numTri,numTri);
    for pp = 1:numTri
        for qq = 1:numTri
            centreDistances(pp,qq) = norm(Solver_setup.triangle_centre_point(pp,:) -  Solver_setup.triangle_centre_point(qq,:)); 
        end
    end
    
end