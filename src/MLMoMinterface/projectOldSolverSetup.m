function [newAllTerms, newSingInd, newLinkOld, groupIndices] =...
    projectOldSolverSetup(new_solver_setup, old_solver_setup,newEdgeLinkOldEdge, newEdgeParallelExternalEdgeLinkOldInternalEdge, oldAllTerms, oldSingInd, newCentreDistances, oldCentreDistances  )
    %1 project based on properties
    %2 multiply terms weights calculated in old mlmom
    
    [numOldEdges, ~, numTerms, numFreq] = size(oldAllTerms);
    %[~, ~, numProp] = size(oldAllProperties);
    
    numNewEdges = new_solver_setup.num_mom_basis_functions;
    
    newAllTerms = zeros(numNewEdges,numNewEdges,numTerms, numFreq);
    %newAllProperties = zeros(numNewEdges,numNewEdges,numProp);
    newSingInd = zeros(numNewEdges,numNewEdges);
    
    %for testing purposes
    %calculatedNonSingInd = zeros(numNewEdges,numNewEdges);
    %calculatedTriInd = zeros(numNewEdges,numNewEdges);
    %calculatedSelfInd = zeros(numNewEdges,numNewEdges);
    
    %link new [mm,nn] to old [mm,nn]
    newLinkOld = zeros(numNewEdges,numNewEdges,2);
    
    %Indices

    groupIndices = [];
    % 2 unique
    twoUnique_int_indices = zeros(numOldEdges *2,2);
    twoUnique_pos_indices = zeros(numOldEdges ,2);
    twoUnique_neg_indices = zeros(numOldEdges ,2);
    twoUnique_ext_indices = zeros(numOldEdges ,2);
    
    %3 unique
    threeUnique_intInt_indices= zeros(numOldEdges*4,2);
    threeUnique_intPos_indices= zeros(numOldEdges*2,2);
    threeUnique_intNeg_indices= zeros(numOldEdges*2,2);
    threeUnique_posInt_indices= zeros(numOldEdges*2,2);
    threeUnique_negInt_indices= zeros(numOldEdges*2,2);
    
    threeUnique_extNotExt_indices= zeros(numOldEdges*4,2);
    threeUnique_notExtExt_indices= zeros(numOldEdges*4,2);
    threeUnique_parIntParInt_indices= zeros(numOldEdges * 4,2);
    threeUnique_extExt_indices = zeros(numOldEdges,2);
    
    %4 unique
    fourUnique_intInt_indices = zeros(numNewEdges.^2,2);
    fourUnique_intPos_indices = zeros(numNewEdges.^2 ,2);
    fourUnique_intNeg_indices = zeros(numNewEdges.^2,2);
    fourUnique_posInt_indices = zeros(numNewEdges.^2,2);
    fourUnique_negInt_indices = zeros(numNewEdges.^2,2);
    
    fourUnique_posPos_indices = zeros(numNewEdges.^2,2);
    fourUnique_posNeg_indices = zeros(numNewEdges.^2,2);
    fourUnique_negPos_indices = zeros(numNewEdges.^2,2);
    fourUnique_negNeg_indices = zeros(numNewEdges.^2,2);
    
    fourUnique_notExtNotExt_indices = zeros(numNewEdges.^2 ,2);
    fourUnique_notExtExt_indices = zeros(numNewEdges.^2 ,2);
    fourUnique_extNotExt_indices = zeros(numNewEdges.^2 ,2);
    fourUnique_extExt_indices = zeros(numNewEdges.^2 ,2);
    
    
    %counts
    % 2 unique
    twoUnique_int_count = 0;
    twoUnique_pos_count = 0;
    twoUnique_neg_count = 0;
    twoUnique_ext_count = 0;
    
    %3 unique
    threeUnique_intInt_count= 0;
    threeUnique_intPos_count= 0;
    threeUnique_intNeg_count= 0;
    threeUnique_posInt_count= 0;
    threeUnique_negInt_count= 0;
    
    threeUnique_extNotExt_count= 0;
    threeUnique_notExtExt_count= 0;
    threeUnique_parIntParInt_count= 0;
    threeUnique_extExt_count = 0;
    
    %4 unique
    fourUnique_intInt_count  = 0;
    fourUnique_intPos_count  = 0;
    fourUnique_intNeg_count  = 0;
    fourUnique_posInt_count  = 0;
    fourUnique_negInt_count  = 0;
    
    fourUnique_posPos_count  = 0;
    fourUnique_posNeg_count  = 0;
    fourUnique_negPos_count  = 0;
    fourUnique_negNeg_count  = 0;
    
    fourUnique_notExtNotExt_count  = 0;
    fourUnique_notExtExt_count  = 0;
    fourUnique_extNotExt_count  = 0;
    fourUnique_extExt_count  = 0;
    
    
    %terms(:,:,1) = A_m_pls_n_pls
    %terms(:,:,2) = Phi_m_pls_n_pls
    %terms(:,:,3) = A_m_pls_n_mns
    %terms(:,:,4) = Phi_m_pls_n_mns
    %terms(:,:,5) = A_m_mns_n_pls
    %terms(:,:,6) = Phi_m_mns_n_pls
    %terms(:,:,7) = A_m_mns_n_mns
    %terms(:,:,8) = Phi_m_mns_n_mns
    
    % if edge length of source (mm) is multiplied by scalar 'a'
    % and edge length of observation (nn) is multiplied by scalar 'b'
    % then the new terms will transformed:
    % A_new = a*b*b*A_old
    % Phi_new = b*Phi_old
    
    % if the distance is multiplied by a scalar 'c'
    % then the new terms will be transformed:
    % A_new = A_old/c;
    % Phi_new = Phi_old/c;
    
    % if a = 0.5, b = 0.5, c = 0.5 then:
    % A_new = 0.25*A_old
    % Phi_new = Phi_old
    
    % if a = 0.5, b = 0.5, c = 1 then:
    % A_new = 0.125*A_old
    % Phi_new = 0.5*Phi_old
    
    %k = 2pi
   % newEdgeParallelExternalEdgeCount = 0;
    AInd = [1 3 5 7];
    PhiInd =[2 4 6 8];
    newTriPlus = new_solver_setup.rwg_basis_functions_trianglePlus;
    newTriMinus = new_solver_setup.rwg_basis_functions_triangleMinus;
    %newEdgeNodes = new_solver_setup.rwg_basis_functions_shared_edge_nodes;
    
    
    %========= ESTIMATED FROM OLD SOLVER SETUP =======
    for freq = 1:numFreq
        for mm =1:numNewEdges
            mmLinkOldEdge = newEdgeLinkOldEdge(mm,1);
            mmLinkType = newEdgeLinkOldEdge(mm,2);
            for nn= 1:numNewEdges
                nnLinkOldEdge = newEdgeLinkOldEdge(nn,1);
                nnLinkType = newEdgeLinkOldEdge(nn,2);
                
                if (mm == nn)
                    %============= 2 unique triangles =============
                    
                    newSingInd(mm,nn) = 1;
                    %calculatedSelfInd(mm,nn) = 1;
                    
                    if (mmLinkType == 0) % lies on old edge
                        newAllTerms(mm,nn,:,freq) = oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,:,freq);
                        twoUnique_int_count = twoUnique_int_count  + 1;
                        twoUnique_int_indices(twoUnique_int_count, :) = [mm,nn];
                        newLinkOld(mm,nn, :) = [mmLinkOldEdge, mmLinkOldEdge];
                    elseif (mmLinkType == 1 ) % parallel to old edge, lies in trianlge plus
                        % 2 new triangles will have some symmetry
                        % symmetry of: terms([1 2]) = terms([7 8]) 
                        %              terms([3 4]) = terms([5 6]) 
                        % triangle plus of old edge has same shape as 
                        % triangle plus of new edge
                        newAllTerms(mm,nn,[1 2],freq) = oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,[1 2],freq);
                        newAllTerms(mm,nn,[7 8],freq) = newAllTerms(mm,nn,[1 2],freq);
                        % reuse old terms for nonSing terms (+-, -+)
                        % will have errors if old trianlge minus has
                        % different shape than new triangle minus
                        newAllTerms(mm,nn,[3 4 5 6],freq) = oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,[3 4 5 6],freq);
                        
                        twoUnique_pos_count =twoUnique_pos_count  + 1;
                        twoUnique_pos_indices(twoUnique_pos_count, :) = [mm,nn];
                        newLinkOld(mm,nn, :) = [mmLinkOldEdge, mmLinkOldEdge];
                        
                    elseif (mmLinkType == -1 )% parallel to old edge, lies in triangle minus
                        newAllTerms(mm,nn,[7 8],freq) = oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,[7 8],freq);
                        newAllTerms(mm,nn,[1 2],freq) = newAllTerms(mm,nn,[7 8],freq);
                        newAllTerms(mm,nn,[3 4 5 6],freq) = oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,[3 4 5 6],freq);
                        
                        twoUnique_neg_count =twoUnique_neg_count  + 1;
                        twoUnique_neg_indices(twoUnique_neg_count, :) = [mm,nn];
                        newLinkOld(mm,nn, :) = [mmLinkOldEdge, mmLinkOldEdge];

                    else % parallel to an external edge
                        % reuse any old edge that the new edge contained
                        % within one of its triangles
                        % even more proned to error
                        % will be 1 or 2 old edges to select from
                        
                        ind = find( newEdgeParallelExternalEdgeLinkOldInternalEdge(:, 1) == mm);
                        oldEdge = newEdgeParallelExternalEdgeLinkOldInternalEdge(ind(1), 2);        
                        newAllTerms(mm,nn,:,freq) = oldAllTerms(oldEdge,oldEdge,:,freq);
                        
                        twoUnique_ext_count =twoUnique_ext_count  + 1;
                        twoUnique_ext_indices(twoUnique_ext_count, :) = [mm,nn];
                        newLinkOld(mm,nn, :) = [oldEdge, oldEdge];

                    end %if (mmLinkType == 0)
                    
                    newAllTerms(mm,nn,AInd,freq) = 0.125* newAllTerms(mm,nn,AInd,freq);
                    newAllTerms(mm,nn,PhiInd,freq) = 0.5* newAllTerms(mm,nn,PhiInd,freq);
                    
                elseif (newTriPlus(mm) == newTriPlus(nn) || newTriPlus(mm) == newTriMinus(nn) ||...
                        newTriMinus(mm) == newTriPlus(nn) || newTriMinus(mm) == newTriMinus(nn))  %if mm == nn

                    %================ 3 unique triangles ================
                    newSingInd(mm,nn) = 1;
                    %calculatedTriInd(mm,nn) = 1;
                    % There are 4 link types: -1, 0, 1, 2
                    % 16 combinations
                    
                    % MAIN CASE 1
                    if (mmLinkType ~= 2 && nnLinkType ~= 2)
                        %------- both on old edge ----
                        if (mmLinkType == 0 && nnLinkType == 0) %1
                            threeUnique_intInt_count = threeUnique_intInt_count + 1;
                            threeUnique_intInt_indices(threeUnique_intInt_count, :) = [mm,nn];
                        %------- 1 edge on old edge, 1 edge parallel old edge ----
                            %  4 cases
                        elseif (mmLinkType == 0 && nnLinkType == 1  ) %2
                            %mm lies on old edge, nn lies in positive triangle
                            % mPlus_nPlus and mMinus_nPlus terms will have same
                            % shape
                            %same shape[1 2 5 6], reuse shape [3 4 7 8]
                            threeUnique_intPos_count = threeUnique_intPos_count + 1;
                            threeUnique_intPos_indices(threeUnique_intPos_count, :) = [mm,nn];
                        elseif (mmLinkType == 0 && nnLinkType == -1  ) %3
                            % mPlus_nMinus and mMinus_nMinus
                            %same shape[3 4 7 8], reuse shape [1 2 5 6]
                            threeUnique_intNeg_count = threeUnique_intNeg_count + 1;
                            threeUnique_intNeg_indices(threeUnique_intNeg_count, :) = [mm,nn];
                        elseif (mmLinkType == 1 && nnLinkType == 0  ) %4
                            % mPlus_nPlus and mPlus_nMinus
                            %same shape[1 2 3 4], reuse shape [5 6 7 8]
                            threeUnique_posInt_count = threeUnique_posInt_count + 1;
                            threeUnique_posInt_indices(threeUnique_posInt_count, :) = [mm,nn];
                        elseif (mmLinkType == -1 && nnLinkType == 0  ) %5
                            % mMinus_nPlus and mMinus_nMinus
                            %same shape[5 6 7 8], reuse shape [1 2 3 4]                           
                            threeUnique_negInt_count = threeUnique_negInt_count + 1;
                            threeUnique_negInt_indices(threeUnique_negInt_count, :) = [mm,nn];
                         %------- both parallel to old edge ----
                         % 4 cases
                        elseif ((mmLinkType == -1 && nnLinkType == -1 )|| ...
                                (mmLinkType == -1 && nnLinkType == 1) || ...
                                (mmLinkType == 1 && nnLinkType == -1) || ...
                                (mmLinkType == 1 && nnLinkType == 1 )) %6 7 8 9 
                            
                            %newAllTerms(mm,nn,:,freq) = swapTerms(oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,:,freq), 3);
                            threeUnique_parIntParInt_count = threeUnique_parIntParInt_count + 1;
                            threeUnique_parIntParInt_indices(threeUnique_parIntParInt_count, :) = [mm,nn];
                        end
                        newLinkOld(mm,nn, :) = [mmLinkOldEdge, nnLinkOldEdge];
                        newAllTerms(mm,nn,:,freq) = oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,:,freq);
                        
                        if (threeUnique_parIntParInt_count > 0)
                            if ( eq(threeUnique_parIntParInt_indices(threeUnique_parIntParInt_count, :), [mm,nn])  )
                                newAllTerms(mm,nn,:,freq) =  swapTerms(newAllTerms(mm,nn,:,freq), 3);
                            end
                        end
                        
                    %------- 1 edge NOT parralel external, 1 edge parallel external ----
                    % 6 cases
                    % MAIN CASE 2
                    elseif (mmLinkType ~= 2 && nnLinkType == 2  ) %10 11 12
                        
                        ind = find( newEdgeParallelExternalEdgeLinkOldInternalEdge(:, 1) == nn & newEdgeParallelExternalEdgeLinkOldInternalEdge(:, 2) ~= mmLinkOldEdge);
                        if (numel(ind) > 0)
                            oldEdgeNN = newEdgeParallelExternalEdgeLinkOldInternalEdge(ind(1), 2);
                        else
                            ind = find (oldSingInd(mmLinkOldEdge, :) == 1);
                            ind2 = find ((ind(:) ~= mmLinkOldEdge));
                            oldEdgeNN = ind(ind2(1));
                        end
                        newAllTerms(mm,nn,:,freq) = compareAndSwapTerms_3unique(oldAllTerms(mmLinkOldEdge, oldEdgeNN, : ,freq),...
                            newCentreDistances, oldCentreDistances, mm, nn, mmLinkOldEdge, oldEdgeNN  );
                        newLinkOld(mm,nn, :) = [mmLinkOldEdge, oldEdgeNN];
                        threeUnique_notExtExt_count = threeUnique_notExtExt_count + 1;
                        threeUnique_notExtExt_indices(threeUnique_notExtExt_count, :) = [mm,nn];
                        
                    % MAIN CASE 3
                    elseif (mmLinkType == 2 && nnLinkType ~= 2  ) %13 14 15
                        
                        ind = find( newEdgeParallelExternalEdgeLinkOldInternalEdge(:, 1) == mm & newEdgeParallelExternalEdgeLinkOldInternalEdge(:, 2) ~= nnLinkOldEdge);
                        if (numel(ind) > 0)
                            oldEdgeMM = newEdgeParallelExternalEdgeLinkOldInternalEdge(ind(1), 2);
                        else
                            ind = find (oldSingInd(:, nnLinkOldEdge) == 1);
                            ind2 = find ((ind(:) ~= nnLinkOldEdge));
                            oldEdgeMM = ind(ind2(1));
                        end
                        newAllTerms(mm,nn,:,freq) = compareAndSwapTerms_3unique(oldAllTerms(oldEdgeMM, nnLinkOldEdge, : ,freq),...
                            newCentreDistances, oldCentreDistances, mm, nn,oldEdgeMM, nnLinkOldEdge );
                        newLinkOld(mm,nn, :) = [oldEdgeMM, nnLinkOldEdge];
                        threeUnique_extNotExt_count = threeUnique_extNotExt_count + 1;
                        threeUnique_extNotExt_indices(threeUnique_extNotExt_count, :) = [mm,nn];
                         
                    %------- both parallel external ----
                    % MAIN CASE 5
                    elseif (mmLinkType == 2 && nnLinkType == 2  ) %16
                        %only possible if both edges in same old triangle
                        %parallelExternal_parallelExternalTriInd(mm,nn) = 1;
                        
                        ind = find( newEdgeParallelExternalEdgeLinkOldInternalEdge(:, 1) == nn);
                        oldEdgeNN = newEdgeParallelExternalEdgeLinkOldInternalEdge(ind, 2);
                        ind = find (oldSingInd(:, oldEdgeNN) == 1);
                        ind2 = find ((ind(:) ~= oldEdgeNN));
                        oldEdgeMM = ind(ind2(1));
                        newAllTerms(mm,nn,:,freq) = compareAndSwapTerms_3unique(oldAllTerms(oldEdgeMM, oldEdgeNN, : ,freq),...
                            newCentreDistances, oldCentreDistances, mm, nn,oldEdgeMM, oldEdgeNN );
                        
                        newLinkOld(mm,nn, :) = [oldEdgeMM, oldEdgeNN];
                        threeUnique_extExt_count  = threeUnique_extExt_count  + 1;
                        threeUnique_extExt_indices (threeUnique_extExt_count , :) = [mm,nn];
                    end
                    newAllTerms(mm,nn,AInd,freq) = 0.125* newAllTerms(mm,nn,AInd,freq);
                    newAllTerms(mm,nn,PhiInd,freq) = 0.5* newAllTerms(mm,nn,PhiInd,freq);
                    
                else % if mm == nn
                    
                    %============= 4 unique triangles ==============
                    %calculatedNonSingInd(mm,nn) = 1;
                    %%=== Main case 1 ====
                    % old relationship is non singular
                    % new edges are not parallel external edges
                    if (mmLinkType ~= 2 && nnLinkType ~= 2 && ~oldSingInd(mmLinkOldEdge,nnLinkOldEdge))
                        % There are 4 link types: -1, 0, 1
                        % 9 cases
                        
                        oldMM = mmLinkOldEdge;
                        oldNN = nnLinkOldEdge;
                        %------- both on old edge ----
                        if (mmLinkType == 0 && nnLinkType == 0) %1 
                            %calculatedNonSingInd(mm,nn) = 1;
                            fourUnique_intInt_count = fourUnique_intInt_count + 1;
                            fourUnique_intInt_indices(fourUnique_intInt_count , :) = [mm,nn];
                            
                            %------- 1 edge on old edge, 1 edge parallel old edge ----
                            %  4 cases
                        elseif (mmLinkType == 0 && nnLinkType == 1  ) %2
                            %mm lies on old edge, nn lies in positive triangle
                            % mPlus_nPlus and mMinus_nPlus terms will have
                            % same shape
                             %same shape terms [1 2 5 6], reuse [3 4 7 8]
                             fourUnique_intPos_count = fourUnique_intPos_count + 1;
                             fourUnique_intPos_indices(fourUnique_intPos_count, :) = [mm,nn];
                        elseif (mmLinkType == 0 && nnLinkType == -1  ) %3
                            % mPlus_nMinus and mMinus_nMinus
                             %same shape terms [3 4 7 8],reuse [1 2 5 6]
                             fourUnique_intNeg_count = fourUnique_intNeg_count + 1;
                             fourUnique_intNeg_indices(fourUnique_intNeg_count, :) = [mm,nn];
                             
                        elseif (mmLinkType == 1 && nnLinkType == 0  ) %4
                            % mPlus_nPlus and mPlus_nMinus
                            %same shape terms [1 2 3 4],reuse [5 6 7 8]
                            fourUnique_posInt_count  = fourUnique_posInt_count + 1;
                            fourUnique_posInt_indices(fourUnique_posInt_count, :) = [mm,nn];
  
                        elseif (mmLinkType == -1 && nnLinkType == 0  ) %5
                            % mMinus_nPlus and mMinus_nMinus
                            %same shape terms [5 6 7 8], reuse [1 2 3 4]
                            fourUnique_negInt_count  = fourUnique_negInt_count + 1;
                            fourUnique_negInt_indices(fourUnique_negInt_count, :) = [mm,nn];
                            
                            %------- both parallel to old edge ----
                            % 4 cases
                        elseif (mmLinkType == -1 && nnLinkType == -1 ) %6
                             %same shape [ 7 8],reuse [1 2 3 4 5 6]
                            fourUnique_negNeg_count = fourUnique_negNeg_count + 1;
                            fourUnique_negNeg_indices(fourUnique_negNeg_count, :) = [mm,nn];
                            
                        elseif (mmLinkType == -1 && nnLinkType == 1) %7
                             %same shape [ 5 6], reuse [1 2 3 4 7 8]
                             fourUnique_negPos_count = fourUnique_negPos_count + 1;
                             fourUnique_negPos_indices(fourUnique_negPos_count, :) = [mm,nn];
                            
                        elseif (mmLinkType == 1 && nnLinkType == -1) %8
                             %same terms [ 3 4], reuse [1 2 5 6 7 8]
                             fourUnique_posNeg_count = fourUnique_posNeg_count + 1;
                             fourUnique_posNeg_indices(fourUnique_posNeg_count, :) = [mm,nn];
                             
                        elseif (mmLinkType == 1 && nnLinkType == 1 ) %9
                             %same shape [ 1 2 ], reuse [3 4 5 6 7 8]
                             fourUnique_posPos_count = fourUnique_posPos_count + 1;
                             fourUnique_posPos_indices(fourUnique_posPos_count, :) = [mm,nn];

                        end % MAIN CASE 1: if (mmLinkType == 0 && nnLinkType == 0) %1 
                        
                    
                    %%=== Main case 2 ====
                    elseif (mmLinkType ~= 2 && nnLinkType ~= 2 && oldSingInd(mmLinkOldEdge,nnLinkOldEdge))
                        % old terms involve 2 or 3 unique triangles
                        
                        oldEdgeNN = nnLinkOldEdge;
                        
                        ind = find (oldSingInd(:, oldEdgeNN) == 1);
                        ind2 = find (  (ind ~= mmLinkOldEdge) & (ind ~= oldEdgeNN) );
                        oldEdgeNNList = ind(ind2);
                        
                        for k = 1:numel(oldEdgeNNList)
                            if (oldSingInd(mmLinkOldEdge,oldEdgeNNList(k)) == 0)
                                oldEdgeNN = oldEdgeNNList(k);
                                break;
                            end
                        end
                        
                        if (oldSingInd(mmLinkOldEdge,oldEdgeNN) == 1 )
                            for k = 1:numel(oldEdgeNNList)
                                if (oldSingInd(mmLinkOldEdge,oldEdgeNN) == 0 )
                                    break
                                end
                                ind = find (oldSingInd(:, oldEdgeNNList(k)) == 1);
                                if ( numel(oldEdgeNNList) == 1)
                                    ind2 = find (  (ind ~= mmLinkOldEdge) & (ind ~= oldEdgeNN) & (ind ~= oldEdgeNNList(1))  );
                                elseif ( numel(oldEdgeNNList) == 2)
                                    ind2 = find (  (ind ~= mmLinkOldEdge) & (ind ~= oldEdgeNN) & (ind ~= oldEdgeNNList(1)) & (ind ~= oldEdgeNNList(2)) );
                                elseif ( numel(oldEdgeNNList) == 3)
                                    ind2 = find (  (ind ~= mmLinkOldEdge) & (ind ~= oldEdgeNN) & (ind ~= oldEdgeNNList(1)) & (ind ~= oldEdgeNNList(2)) & (ind ~= oldEdgeNNList(3)) );
                                elseif ( numel(oldEdgeNNList) == 4)
                                    ind2 = find (  (ind ~= mmLinkOldEdge) & (ind ~= oldEdgeNN) & (ind ~= oldEdgeNNList(1)) & (ind ~= oldEdgeNNList(2)) & (ind ~= oldEdgeNNList(3)) & (ind ~= oldEdgeNNList(4)));
                                elseif ( numel(oldEdgeNNList) == 5)
                                    ind2 = find (  (ind ~= mmLinkOldEdge) & (ind ~= oldEdgeNN) & (ind ~= oldEdgeNNList(1)) & (ind ~= oldEdgeNNList(2)) & (ind ~= oldEdgeNNList(3)) & (ind ~= oldEdgeNNList(4)) & (ind ~= oldEdgeNNList(5)));
                                end
                                oldEdgeNNList2 = ind(ind2);
                                for j = 1:numel(oldEdgeNNList2)
                                    if (oldSingInd(mmLinkOldEdge,oldEdgeNNList2(j)) == 0)
                                        oldEdgeNN = oldEdgeNNList2(j);
                                        break;
                                    end
                                end
                            end
                        end
                        fourUnique_notExtNotExt_count = fourUnique_notExtNotExt_count + 1;
                        fourUnique_notExtNotExt_indices(fourUnique_notExtNotExt_count, :) = [mm,nn];
                        oldMM = mmLinkOldEdge;
                        oldNN = oldEdgeNN;  

                    %%=== Main case 3 ====    
                    elseif ( mmLinkType ~= 2 && nnLinkType == 2)
                        % nn: parallel to external edge
                        
                        fourUnique_notExtExt_count = fourUnique_notExtExt_count+ 1;
                        
                        %oldEdgeNN is the nearest old edge to external nn
                        ind = find( newEdgeParallelExternalEdgeLinkOldInternalEdge(:, 1) == nn);
                        
                        oldEdgeNN = newEdgeParallelExternalEdgeLinkOldInternalEdge(ind, 2);
                        for k = 1:numel(oldEdgeNN)
                            if (~oldSingInd(mmLinkOldEdge,oldEdgeNN(k)))
                                %calculatedNonSingInd(mm,nn) = 1;
                                fourUnique_notExtExt_indices(fourUnique_notExtExt_count, :) = [mm,nn];
                                oldMM = mmLinkOldEdge;
                                oldNN = oldEdgeNN(k);
                                break;
                            end
                        end
                        
                        %if (calculatedNonSingInd(mm,nn) == 0)
                        if (fourUnique_notExtExt_indices(fourUnique_notExtExt_count, 1) == 0)
                            for k = 1:numel(oldEdgeNN)
                                if(fourUnique_notExtExt_indices(fourUnique_notExtExt_count, 1) ~= 0)
                                    break
                                end
                                if (oldEdgeNN(k) ~= mmLinkOldEdge)
                                    % old group: 3 unique
                                    ind = find (oldSingInd(mmLinkOldEdge, :) == 1);
                                    ind2 = find ((ind(:) ~= mmLinkOldEdge) & (ind(:) ~= oldEdgeNN(k)));
                                    oldEdgeMM = ind(ind2);
                                    for j = 1:numel(oldEdgeMM)
                                        if (oldSingInd(oldEdgeMM(j), oldEdgeNN(k)) == 0)
                                            fourUnique_notExtExt_indices(fourUnique_notExtExt_count, :) = [mm,nn];
                                            oldMM = oldEdgeMM(j);
                                            oldNN = oldEdgeNN(k);
                                            break
                                        end
                                    end
                                end
                            end
                        end %if (fourUnique_notExtExt_indices(fourUnique_notExtExt_count, 1) == 0)
                        
                        %if (calculatedNonSingInd(mm,nn) == 0)
                        if (fourUnique_notExtExt_indices(fourUnique_notExtExt_count, 1) == 0)
                            % old group: 2 unique triangles (eg 21, 4)
                            % oldEdgeNN can only have 1 element
                            %calculatedNonSingInd(mm,nn) = 1;
                            ind = find (oldSingInd(mmLinkOldEdge, :) == 1);
                            ind2 = find ((ind(:) ~= mmLinkOldEdge));
                            otherOldMM = ind(ind2);
                            
                            ind = find (oldSingInd(otherOldMM(1), :) == 1);
                            ind2 = find (  (ind(:) ~= otherOldMM(1)) & (ind(:) ~= otherOldMM(2)) & (ind(:) ~= mmLinkOldEdge)  );
                            otherOldMM = ind(ind2(1));
                            
                            fourUnique_notExtExt_indices(fourUnique_notExtExt_count, :) = [mm,nn];
                            oldMM = otherOldMM;
                            oldNN = oldEdgeNN(1);
                            
                        end
                        
                        
                    %%=== Main case 4 ====    
                    elseif (mmLinkType == 2 && nnLinkType ~= 2)
                        % mm: parallel to external edge
                        
                        fourUnique_extNotExt_count = fourUnique_extNotExt_count+ 1;
                        
                        ind = find( newEdgeParallelExternalEdgeLinkOldInternalEdge(:, 1) == mm);
                        
                        oldEdgeMM = newEdgeParallelExternalEdgeLinkOldInternalEdge(ind, 2);
                        for k = 1:numel(oldEdgeMM)
                            if (~oldSingInd(oldEdgeMM(k),nnLinkOldEdge ))
                                
                                fourUnique_extNotExt_indices(fourUnique_extNotExt_count, :) = [mm,nn];
                                oldMM = oldEdgeMM(k);
                                oldNN = nnLinkOldEdge;
                                break;
                            end
                        end
                        
                        %if (calculatedNonSingInd(mm,nn) == 0)
                        if (fourUnique_extNotExt_indices(fourUnique_extNotExt_count, 1) == 0)
                            for k = 1:numel(oldEdgeMM)
                                if(fourUnique_extNotExt_indices(fourUnique_extNotExt_count, 1)  ~= 0)
                                    break
                                end
                                if (oldEdgeMM(k) ~= nnLinkOldEdge)
                                    % old group: 3 unique
                                    ind = find (oldSingInd(:, nnLinkOldEdge) == 1);
                                    ind2 = find ((ind(:) ~= nnLinkOldEdge) & (ind(:) ~= oldEdgeMM(k)));
                                    oldEdgeNN = ind(ind2);
                                    for j = 1:numel(oldEdgeNN)
                                        if (oldSingInd(oldEdgeMM(k), oldEdgeNN(j)) == 0)
                                            fourUnique_extNotExt_indices(fourUnique_extNotExt_count, :) = [mm,nn];
                                            oldMM = oldEdgeMM(k);
                                            oldNN = oldEdgeNN(j);
                                            break
                                        end
                                    end
                                end
                            end
                        end %if (fourUnique_extNotExt_indices(fourUnique_extNotExt_count, 1) == 0)
                        
                        %if (calculatedNonSingInd(mm,nn) == 0)
                        if (fourUnique_extNotExt_indices(fourUnique_extNotExt_count, 1) == 0)
                            % old group: 2 unique triangles (eg 21, 4)
                            % oldEdgeNN can only have 1 element
                            ind = find (oldSingInd(:, nnLinkOldEdge) == 1);
                            ind2 = find ((ind(:) ~= nnLinkOldEdge));
                            otherOldNN = ind(ind2);
                            
                            ind = find (oldSingInd(:,otherOldNN(1)) == 1);
                            ind2 = find (  (ind(:) ~= otherOldNN(1)) & (ind(:) ~= otherOldNN(2)) & (ind(:) ~= nnLinkOldEdge)  );
                            otherOldNN = ind(ind2(1));
                            
                            fourUnique_extNotExt_indices(fourUnique_extNotExt_count, :) = [mm,nn];
                            oldMM = oldEdgeMM(1);
                            oldNN = otherOldNN;
                        end

                    %%=== Main case 5 ====
                    elseif (mmLinkType == 2 && nnLinkType == 2)
                         
                        fourUnique_extExt_count = fourUnique_extExt_count + 1;
                        
                        ind = find( newEdgeParallelExternalEdgeLinkOldInternalEdge(:, 1) == mm);
                        oldEdgeMM = newEdgeParallelExternalEdgeLinkOldInternalEdge(ind, 2);
                        
                        ind = find( newEdgeParallelExternalEdgeLinkOldInternalEdge(:, 1) == nn);
                        oldEdgeNN = newEdgeParallelExternalEdgeLinkOldInternalEdge(ind, 2);
                        
                        for k = 1: numel(oldEdgeMM)
                            if (fourUnique_extExt_indices(fourUnique_extExt_count, 1) ~= 0)
                                break
                            end
                            for j = 1:numel(oldEdgeNN)
                                if (~oldSingInd(oldEdgeMM(k),oldEdgeNN(j)))
                                    fourUnique_extExt_indices(fourUnique_extExt_count, :) = [mm,nn];
                                    oldMM = oldEdgeMM(k);
                                    oldNN = oldEdgeNN(j);
                                    break
                                end
                            end
                        end     

                    end % if (mmLinkType ~= 2 && nnLinkType ~= 2 && ~oldSingInd(mmLinkOldEdge,nnLinkOldEdge))
                
                    newLinkOld(mm,nn, :) = [oldMM, oldNN];
                    
                    [newAllTerms(mm,nn,:,freq)] = projectCentreDistances(oldAllTerms(oldMM,oldNN,:,freq),...
                        newCentreDistances(mm,nn,:), oldCentreDistances(oldMM,oldNN,:));
                    
                    newAllTerms(mm,nn,AInd,freq) = 0.125* newAllTerms(mm,nn,AInd,freq); %0.125 works best
                    newAllTerms(mm,nn,PhiInd,freq) = 0.5* newAllTerms(mm,nn,PhiInd,freq); % 0.25 works best
                    
                end % if (mm == nn)
            end %for nn= 1:numNewEdges
            
        end %for mm =1:numNewEdges
    end %for freq = 1:numFreq
    
    %2 unique
    twoUnique_int_indices = twoUnique_int_indices(1:twoUnique_int_count, :);
    twoUnique_pos_indices = twoUnique_pos_indices(1:twoUnique_pos_count, :);
    twoUnique_neg_indices = twoUnique_neg_indices(1:twoUnique_neg_count, :);
    twoUnique_ext_indices = twoUnique_ext_indices(1:twoUnique_ext_count, :);

    %3 unique
    threeUnique_intInt_indices = threeUnique_intInt_indices(1:threeUnique_intInt_count, :);
    threeUnique_intPos_indices = threeUnique_intPos_indices(1:threeUnique_intPos_count, :);
    threeUnique_intNeg_indices = threeUnique_intNeg_indices(1:threeUnique_intNeg_count, :);
    threeUnique_posInt_indices = threeUnique_posInt_indices(1:threeUnique_posInt_count, :);
    threeUnique_negInt_indices = threeUnique_negInt_indices(1:threeUnique_negInt_count, :);
    
    threeUnique_extNotExt_indices = threeUnique_extNotExt_indices(1:threeUnique_extNotExt_count, :);
    threeUnique_notExtExt_indices = threeUnique_notExtExt_indices(1:threeUnique_notExtExt_count, :);
    threeUnique_parIntParInt_indices = threeUnique_parIntParInt_indices(1:threeUnique_parIntParInt_count, :);
    threeUnique_extExt_indices = threeUnique_extExt_indices(1:threeUnique_extExt_count, :);
    
    %4 unique
    
    fourUnique_intInt_indices = fourUnique_intInt_indices(1:fourUnique_intInt_count, :);
    fourUnique_intPos_indices = fourUnique_intPos_indices(1:fourUnique_intPos_count, :);
    fourUnique_intNeg_indices = fourUnique_intNeg_indices(1:fourUnique_intNeg_count, :);
    fourUnique_posInt_indices = fourUnique_posInt_indices(1:fourUnique_posInt_count, :);
    fourUnique_negInt_indices = fourUnique_negInt_indices(1:fourUnique_negInt_count, :);
    
    fourUnique_posPos_indices = fourUnique_posPos_indices(1:fourUnique_posPos_count, :);
    fourUnique_posNeg_indices = fourUnique_posNeg_indices(1:fourUnique_posNeg_count, :);
    fourUnique_negPos_indices = fourUnique_negPos_indices(1:fourUnique_negPos_count, :);
    fourUnique_negNeg_indices = fourUnique_negNeg_indices(1:fourUnique_negNeg_count, :);
    
    fourUnique_notExtNotExt_indices = fourUnique_notExtNotExt_indices(1:fourUnique_notExtNotExt_count, :);
    fourUnique_notExtExt_indices = fourUnique_notExtExt_indices(1:fourUnique_notExtExt_count, :);
    fourUnique_extNotExt_indices = fourUnique_extNotExt_indices(1:fourUnique_extNotExt_count, :);
    fourUnique_extExt_indices = fourUnique_extExt_indices(1:fourUnique_extExt_count, :);
    
    %copy to groupIndices struct
    
    %2 unique
    groupIndices.twoUnique_int_indices =twoUnique_int_indices ;
    groupIndices.twoUnique_pos_indices = twoUnique_pos_indices;
    groupIndices.twoUnique_neg_indices = twoUnique_neg_indices;
    groupIndices.twoUnique_ext_indices = twoUnique_ext_indices;
    
    %3 unique
    groupIndices.threeUnique_intInt_indices = threeUnique_intInt_indices;
    groupIndices.threeUnique_intPos_indices = threeUnique_intPos_indices;
    groupIndices.threeUnique_intNeg_indices = threeUnique_intNeg_indices;
    groupIndices.threeUnique_posInt_indices = threeUnique_posInt_indices;
    groupIndices.threeUnique_negInt_indices = threeUnique_negInt_indices;
    
    groupIndices.threeUnique_extNotExt_indices = threeUnique_extNotExt_indices;
    groupIndices.threeUnique_notExtExt_indices =threeUnique_notExtExt_indices ;
    groupIndices.threeUnique_parIntParInt_indices = threeUnique_parIntParInt_indices;
    groupIndices.threeUnique_extExt_indices = threeUnique_extExt_indices;
    
    %4 unique
    groupIndices.fourUnique_intInt_indices = fourUnique_intInt_indices ;
    groupIndices.fourUnique_intPos_indices = fourUnique_intPos_indices ;
    groupIndices.fourUnique_intNeg_indices = fourUnique_intNeg_indices;
    groupIndices.fourUnique_posInt_indices = fourUnique_posInt_indices ;
    groupIndices.fourUnique_negInt_indices = fourUnique_negInt_indices ;
    
    groupIndices.fourUnique_posPos_indices = fourUnique_posPos_indices ;
    groupIndices.fourUnique_posNeg_indices = fourUnique_posNeg_indices;
    groupIndices.fourUnique_negPos_indices = fourUnique_negPos_indices;
    groupIndices.fourUnique_negNeg_indices = fourUnique_negNeg_indices;
    
    groupIndices.fourUnique_notExtNotExt_indices = fourUnique_notExtNotExt_indices;
    groupIndices.fourUnique_notExtExt_indices = fourUnique_notExtExt_indices ;
    groupIndices.fourUnique_extNotExt_indices = fourUnique_extNotExt_indices;
    groupIndices.fourUnique_extExt_indices = fourUnique_extExt_indices;
    

end

function [swappedTerms] = swapTerms(terms, swapCode )
    % There are 2 terms associated with every triangle pair between 4 triangles
    %terms(1) = A_m_pls_n_pls
    %terms(2) = Phi_m_pls_n_pls
    %terms(3) = A_m_pls_n_mns
    %terms(4) = Phi_m_pls_n_mns
    %terms(5) = A_m_mns_n_pls
    %terms(6) = Phi_m_mns_n_pls
    %terms(7) = A_m_mns_n_mns
    %terms(8) = Phi_m_mns_n_mns
    
    %If 4 terms between 2 triangle pairs are swapped then the other 4 terms
    % are also swapped in the same manner
    %So there are only 3 swap codes
    %1: mPlus_nPlus > mPlus_nMinus
    %2: mPlus_nPlus > mMinus_nPlus
    %3: mPlus_nPlus > mMinus_nMinus
    
    switch swapCode
        case 1
            swappedTerms = terms([3 4 1 2 7 8 5 6]);
            %swappedTerms(1,1,[1 3 5 7]) = -1 *swappedTerms(1,1,[1 3 5 7]) ;
        case 2
            swappedTerms = terms([5 6 7 8 1 2 3 4]);
            %swappedTerms(1,1,[1 3 5 7]) = -1 *swappedTerms(1,1,[1 3 5 7]) ;
        case 3
            swappedTerms = terms([7 8 5 6 3 4 1 2]);
        otherwise
            swappedTerms = terms;    
    end


end

function [swappedTerms] = swapTermCodes(terms, termCodeOne, termCodeTwo )
    %termCode:
    % 1: mPlus_nPlus
    % 2: mPlus_nMinus
    % 3: mMinus_nPlus
    % 4: mMinus_nMinus
    
    % termCodeOne = termCodeTwo will have no effect
    switch termCodeOne
        case 1
            switch termCodeTwo
                case 1
                    swappedTerms = terms;
                case 2
                    swappedTerms = swapTerms(terms, 1 );
                case 3
                    swappedTerms = swapTerms(terms, 2 );
                case 4
                    swappedTerms = swapTerms(terms, 3 );
                otherwise
                    swappedTerms = terms;
            end
        case 2
            switch termCodeTwo
                case 1
                    swappedTerms = swapTerms(terms, 1 );
                case 2
                    swappedTerms = terms;
                case 3
                    swappedTerms = swapTerms(terms, 3 );
                case 4
                    swappedTerms = swapTerms(terms, 2 );
                otherwise
                    swappedTerms = terms;
            end
            
        case 3
            switch termCodeTwo
                case 1
                    swappedTerms = swapTerms(terms, 2 );
                case 2
                    swappedTerms = swapTerms(terms, 3 );
                case 3
                    swappedTerms = terms;
                case 4
                    swappedTerms = swapTerms(terms, 1 );
                otherwise
                    swappedTerms = terms;
            end
            
        case 4
            switch termCodeTwo
                case 1
                    swappedTerms = swapTerms(terms, 3 );
                case 2
                    swappedTerms = swapTerms(terms, 2 );
                case 3
                    swappedTerms = swapTerms(terms, 1 );
                case 4
                    swappedTerms = terms;
                otherwise
                    swappedTerms = terms;
            end
        otherwise 
            swappedTerms = terms;
    end
    
end

%newCentreDistances, oldCentreDistances 
%function [swappedTerms] = compareAndSwapTerms_3unique(oldTerms, new_solver_setup, old_solver_setup, oldMM, oldNN, mm, nn )
function [swappedTerms] = compareAndSwapTerms_3unique(oldTerms,newCentreDistances, oldCentreDistances,mm, nn, oldMM, oldNN  )
%     oldMPlus = old_solver_setup.rwg_basis_functions_trianglePlus(oldMM);
%     oldMMinus = old_solver_setup.rwg_basis_functions_triangleMinus(oldMM);
%     oldNPlus = old_solver_setup.rwg_basis_functions_trianglePlus(oldNN);
%     oldNMinus = old_solver_setup.rwg_basis_functions_triangleMinus(oldNN);
    
    [~, ind] = sort(oldCentreDistances(oldMM,oldNN, :));
    %sorted from smallest to largest
    % min
    ind = ind(:);
    switch ind(1)
        case 1 % mmPlus_nnPlus
            oldCode = 1;
        case 2  %mmPlus_nnMinus
            oldCode = 2;
        case 3 % mmMinus_nnPlus
            oldCode = 3;
        case 4 % mmMinus_nnMinus
            oldCode = 4;
        otherwise
            swappedTerms = oldTerms;
            return;
    end
    
    %identify terms associated with middle triangle
%     if (oldMPlus == oldNPlus)
%         oldCode = 1;
%     elseif (oldMPlus == oldNMinus)
%         oldCode = 2;
%     elseif (oldMMinus == oldNPlus)
%         oldCode = 3;
%     elseif (oldMMinus == oldNMinus)
%         oldCode = 4;
%     else
%         swappedTerms = oldTerms;
%         return;
%     end
    
    %-----
%     mPlus = new_solver_setup.rwg_basis_functions_trianglePlus(mm);
%     mMinus = new_solver_setup.rwg_basis_functions_triangleMinus(mm);
%     nPlus = new_solver_setup.rwg_basis_functions_trianglePlus(nn);
%     nMinus = new_solver_setup.rwg_basis_functions_triangleMinus(nn);
    
    [~, ind] = sort(newCentreDistances(mm,nn, :));
    %identify terms associated with middle triangle
    ind = ind(:);
    switch ind(1)
        case 1 % mmPlus_nnPlus
            newCode = 1;
        case 2  %mmPlus_nnMinus
            newCode = 2;
        case 3 % mmMinus_nnPlus
            newCode = 3;
        case 4 % mmMinus_nnMinus
            newCode = 4;
        otherwise
            swappedTerms = oldTerms;
            return;
    end
    
%     if (mPlus == nPlus)
%         newCode = 1;
%     elseif (mPlus == nMinus)
%         newCode = 2;
%     elseif (mMinus == nPlus)
%         newCode = 3;
%     elseif (mMinus == nMinus)
%         newCode = 4;
%     else
%         swappedTerms = oldTerms;
%         return;
%     end
    
    swappedTerms = swapTermCodes(oldTerms, oldCode, newCode );
end

function [newTerms] = projectCentreDistances(oldTerms, newCentreDistances, oldCentreDistances )
    
    newTerms = zeros(1,1,8);

%     mmPlusCentre = new_solver_setup.triangle_centre_point(new_solver_setup.rwg_basis_functions_trianglePlus(mm),:);
%     mmMinusCentre = new_solver_setup.triangle_centre_point(new_solver_setup.rwg_basis_functions_triangleMinus(mm),:);
%     nnPlusCentre =  new_solver_setup.triangle_centre_point(new_solver_setup.rwg_basis_functions_trianglePlus(nn),:);
%     nnMinusCentre = new_solver_setup.triangle_centre_point(new_solver_setup.rwg_basis_functions_triangleMinus(nn),:);
%     
%     oldmmPlusCentre = old_solver_setup.triangle_centre_point(old_solver_setup.rwg_basis_functions_trianglePlus(oldMM),:);
%     oldmmMinusCentre = old_solver_setup.triangle_centre_point(old_solver_setup.rwg_basis_functions_triangleMinus(oldMM),:);
%     oldnnPlusCentre =  old_solver_setup.triangle_centre_point(old_solver_setup.rwg_basis_functions_trianglePlus(oldNN),:);
%     oldnnMinusCentre = old_solver_setup.triangle_centre_point(old_solver_setup.rwg_basis_functions_triangleMinus(oldNN),:);
%     
%     mmPlus_nnPlus = norm(mmPlusCentre-nnPlusCentre);
%     mmPlus_nnMinus = norm(mmPlusCentre-nnMinusCentre);
%     mmMinus_nnPlus= norm(mmMinusCentre-nnPlusCentre);
%     mmMinus_nnMinus = norm(mmMinusCentre-nnMinusCentre);
%     
%     oldmmPlus_nnPlus = norm(oldmmPlusCentre-oldnnPlusCentre);
%     oldmmPlus_nnMinus = norm(oldmmPlusCentre-oldnnMinusCentre);
%     oldmmMinus_nnPlus= norm(oldmmMinusCentre-oldnnPlusCentre);
%     oldmmMinus_nnMinus = norm(oldmmMinusCentre-oldnnMinusCentre);
    
    % NB NB 1 INDICES DEPENDANT ON ORDER OF TERMS
    % NB NB 2 Triangle plus/minus labels are already alligned
%     newTerms(1) = oldmmPlus_nnPlus/mmPlus_nnPlus;
%     newTerms(3) = oldmmPlus_nnMinus/mmPlus_nnMinus;
%     newTerms(5) = oldmmMinus_nnPlus/mmMinus_nnPlus;
%     newTerms(7) = oldmmMinus_nnMinus/mmMinus_nnMinus; 

%     newTerms(1) = oldCentreDistances(oldMM,oldNN,1)/newCentreDistances(mm,nn,1);
%     newTerms(3) = oldCentreDistances(oldMM,oldNN,2)/newCentreDistances(mm,nn,2);
%     newTerms(5) = oldCentreDistances(oldMM,oldNN,3)/newCentreDistances(mm,nn,3);
%     newTerms(7) = oldCentreDistances(oldMM,oldNN,4)/newCentreDistances(mm,nn,4);
    
    newTerms(1) = oldCentreDistances(1,1,1)/newCentreDistances(1,1,1);
    newTerms(3) = oldCentreDistances(1,1,2)/newCentreDistances(1,1,2);
    newTerms(5) = oldCentreDistances(1,1,3)/newCentreDistances(1,1,3);
    newTerms(7) = oldCentreDistances(1,1,4)/newCentreDistances(1,1,4);
    newTerms([2 4 6 8]) = newTerms([1 3 5 7]);
    
    newTerms = oldTerms .* newTerms;
    
%     [~, ind] = sort([mmPlus_nnPlus,mmPlus_nnMinus,mmMinus_nnPlus,mmMinus_nnMinus]);
%     %sorted from smallest to largest
%     % min
%     switch ind(1)
%         case 1 % mmPlus_nnPlus
%         case 2  %mmPlus_nnMinus
%              newTerms = swapTerms(newTerms, 1 );
%         case 3 % mmMinus_nnPlus
%             newTerms = swapTerms(newTerms, 2 );
%         case 4 % mmMinus_nnMinus
%             newTerms = swapTerms(newTerms, 3 );
%     end
    
end