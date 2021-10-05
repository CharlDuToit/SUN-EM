function [newAllTerms, newAllProperties, newSingInd, calculatedNonSingInd, calculatedTriInd, calculatedSelfInd] = projectOldSolverSetup(new_solver_setup, old_solver_setup,newEdgeLinkOldEdge, newEdgeParallelExternalEdgeLinkOldInternalEdge, oldAllTerms, oldAllProperties, oldSingInd )
    %1 project based on properties
    %2 multiply terms weights calculated in old mlmom
    
    [numOldEdges, ~, numTerms, numFreq] = size(oldAllTerms);
    [~, ~, numProp] = size(oldAllProperties);
    
    numNewEdges = new_solver_setup.num_mom_basis_functions;
    
    newAllTerms = zeros(numNewEdges,numNewEdges,numTerms, numFreq);
    newAllProperties = zeros(numNewEdges,numNewEdges,numProp);
    newSingInd = zeros(numNewEdges,numNewEdges);
    
    %for testing purposes
    calculatedNonSingInd = zeros(numNewEdges,numNewEdges);
    calculatedTriInd = zeros(numNewEdges,numNewEdges);
    calculatedSelfInd = zeros(numNewEdges,numNewEdges);
    
    %calculate after first main for loop
    onOldEdge_ParallelExternalTriInd = zeros(numNewEdges,numNewEdges);
    onOldEdge_OnOldEdgeTriInd = zeros(numNewEdges,numNewEdges);
    parallelInternal_parallelExternalTriInd = zeros(numNewEdges,numNewEdges);
    parallelExternal_parallelExternalTriInd = zeros(numNewEdges,numNewEdges);
    
    %onOldEdge_parallelExternalNonSingInd = zeros(numNewEdges,numNewEdges);
    
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
    
    tic;
    %=========FIRST LOOP - ESTIMATED FROM OLD SOLVER SETUP =======
    for freq = 1:numFreq
        for mm =1:numNewEdges
            mmLinkOldEdge = newEdgeLinkOldEdge(mm,1);
            mmLinkType = newEdgeLinkOldEdge(mm,2);
            for nn= 1:numNewEdges
                nnLinkOldEdge = newEdgeLinkOldEdge(nn,1);
                nnLinkType = newEdgeLinkOldEdge(nn,2);
                
                if (mm == nn)
                    %============= 2 unique triangles =============
                    
                    newAllProperties(mm,nn,1:3) = [0, 1, 0];
                    newSingInd(mm,nn) = 1;
                    calculatedSelfInd(mm,nn) = 1; % move this to one of the cases to double check only those cases
                    if (mmLinkType == 0) % lies on old edge
                        newAllTerms(mm,nn,:,freq) = oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,:,freq);

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
                        
                    elseif (mmLinkType == -1 )% parallel to old edge, lies in triangle minus
                        newAllTerms(mm,nn,[7 8],freq) = oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,[7 8],freq);
                        newAllTerms(mm,nn,[1 2],freq) = newAllTerms(mm,nn,[7 8],freq);
                        newAllTerms(mm,nn,[3 4 5 6],freq) = oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,[3 4 5 6],freq);

                    else % parallel to an external edge
                        % reuse any old edge that the new edge contained
                        % within one of its triangles
                        % even more proned to error
                        % will be 1 or 2 old edges to select from
                        
                        ind = find( newEdgeParallelExternalEdgeLinkOldInternalEdge(:, 1) == mm);
                        %newEdgeParallelExternalEdgeCount = newEdgeParallelExternalEdgeCount + 1;
                        %oldEdge = newEdgeParallelExternalEdgeLinkOldInternalEdge(newEdgeParallelExternalEdgeCount, 2);
                        oldEdge = newEdgeParallelExternalEdgeLinkOldInternalEdge(ind(1), 2);
                        
%                         if (newEdgeParallelExternalEdgeLinkOldInternalEdge(newEdgeParallelExternalEdgeCount, 1) ~= mm) % to test if indexing always ligns up 
%                             indexingiswrongSetBreakpointHere = 1;
%                         end
                        newAllTerms(mm,nn,:,freq) = oldAllTerms(oldEdge,oldEdge,:,freq);

                    end %if (mmLinkType == 0)
                    newAllTerms(mm,nn,AInd,freq) = 0.125* newAllTerms(mm,nn,AInd,freq);
                    newAllTerms(mm,nn,PhiInd,freq) = 0.5* newAllTerms(mm,nn,PhiInd,freq);
                elseif (newTriPlus(mm) == newTriPlus(nn) || newTriPlus(mm) == newTriMinus(nn) ||...
                        newTriMinus(mm) == newTriPlus(nn) || newTriMinus(mm) == newTriMinus(nn))  %if mm == nn
                    
                    %================ 3 unique triangles ================
                    newSingInd(mm,nn) = 1;
                    % There are 4 link types: -1, 0, 1, 2
                    % 16 combinations
                    %------- both on old edge ----
                    if (mmLinkType == 0 && nnLinkType == 0) %1
                        calculatedTriInd(mm,nn) = 1;
                        newAllTerms(mm,nn,:,freq) = oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,:,freq);
                        onOldEdge_OnOldEdgeTriInd(mm,nn) =1 ;
                        
                    %------- 1 edge on old edge, 1 edge parallel old edge ----
                    %  4 cases 
                    elseif (mmLinkType == 0 && nnLinkType == 1  ) %2
                        %mm lies on old edge, nn lies in positive triangle
                        % mPlus_nPlus and mMinus_nPlus terms will have same
                        % shape
                        calculatedTriInd(mm,nn) = 1;
                        newAllTerms(mm,nn,[1 2 5 6],freq) = oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,[1 2 5 6],freq);
                        % simply reuse the rest - proned to error if
                        % different shape
                        newAllTerms(mm,nn,[3 4 7 8],freq) = oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,[3 4 7 8],freq);
                    elseif (mmLinkType == 0 && nnLinkType == -1  ) %3
                       % mPlus_nMinus and mMinus_nMinus
                        calculatedTriInd(mm,nn) = 1;
                        newAllTerms(mm,nn,[3 4 7 8],freq) = oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,[3 4 7 8],freq);
                        % simply reuse rest
                        newAllTerms(mm,nn,[1 2 5 6],freq) = oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,[1 2 5 6],freq);
                    elseif (mmLinkType == 1 && nnLinkType == 0  ) %4
                        % mPlus_nPlus and mPlus_nMinus
                        calculatedTriInd(mm,nn) = 1;
                        newAllTerms(mm,nn,[1 2 3 4],freq) = oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,[1 2 3 4],freq);
                        % simply reuse rest
                        newAllTerms(mm,nn,[5 6 7 8],freq) = oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,[5 6 7 8],freq);
                    elseif (mmLinkType == -1 && nnLinkType == 0  ) %5
                        % mMinus_nPlus and mMinus_nMinus
                        calculatedTriInd(mm,nn) = 1;
                        newAllTerms(mm,nn,[5 6 7 8],freq) = oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,[5 6 7 8],freq);
                        % simply reuse rest
                        newAllTerms(mm,nn,[1 2 3 4],freq) = oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,[1 2 3 4],freq);
                        
                    %------- 1 edge on old edge, 1 edge parallel external ----
                    % 2 cases
                    elseif (mmLinkType == 0 && nnLinkType == 2  ) %6
                        %m lies on old edge, n parallel to external edge
                        %addTriangles ensures that new triangle plus of n
                        % is never central (4th) triangle
                        onOldEdge_ParallelExternalTriInd(mm,nn) = 1;
                        calculatedTriInd(mm,nn) = 1;
                    elseif (mmLinkType == 2 && nnLinkType == 0  ) %7
                        onOldEdge_ParallelExternalTriInd(mm,nn) = 1;
                        calculatedTriInd(mm,nn) = 1;
                    %------- both parallel to old edge ----
                    % 4 cases
                    elseif ((mmLinkType == -1 && nnLinkType == -1 )|| ...
                            (mmLinkType == -1 && nnLinkType == 1) || ...
                            (mmLinkType == 1 && nnLinkType == -1) || ...
                            (mmLinkType == 1 && nnLinkType == 1 )) %8
                        calculatedTriInd(mm,nn) = 1;
                        newAllTerms(mm,nn,:,freq) = swapTerms(oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,:,freq), 3);
                    %------- 1 parallel to old edge, 1 parallel external ----
                    % 4 cases
                    elseif (mmLinkType == -1 && nnLinkType == 2  ) %12
                        parallelInternal_parallelExternalTriInd(mm,nn) = 1;
                        calculatedTriInd(mm,nn) = 1;
                    elseif (mmLinkType == 1 && nnLinkType == 2  ) %13
                        parallelInternal_parallelExternalTriInd(mm,nn) = 1;
                        calculatedTriInd(mm,nn) = 1;
                    elseif (mmLinkType == 2 && nnLinkType == -1  ) %14
                        parallelInternal_parallelExternalTriInd(mm,nn) = 1;
                        calculatedTriInd(mm,nn) = 1;
                    elseif (mmLinkType == 2 && nnLinkType == 1  ) %15
                        parallelInternal_parallelExternalTriInd(mm,nn) = 1;
                        calculatedTriInd(mm,nn) = 1;
                    %------- both parallel external ----
                    elseif (mmLinkType == 2 && nnLinkType == 2  ) %16
                        parallelExternal_parallelExternalTriInd(mm,nn) = 1;
                        calculatedTriInd(mm,nn) = 1;
                    end
                    newAllTerms(mm,nn,AInd,freq) = 0.125* newAllTerms(mm,nn,AInd,freq);
                    newAllTerms(mm,nn,PhiInd,freq) = 0.5* newAllTerms(mm,nn,PhiInd,freq);
                    
                else % if mm == nn
                    
                    %============= 4 unique triangles ==============
                    %%=== Main case 1 ====
                    % old relationship is non singular
                    % new edges are not parallel external edges
                    if (mmLinkType ~= 2 && nnLinkType ~= 2 && ~oldSingInd(mmLinkOldEdge,nnLinkOldEdge))
                        % There are 4 link types: -1, 0, 1
                        % 9 cases
                        calculatedNonSingInd(mm,nn) = 1;
                        %------- both on old edge ----
                        if (mmLinkType == 0 && nnLinkType == 0) %1 
                            %calculatedNonSingInd(mm,nn) = 1;
                            [newAllTerms(mm,nn,:,freq)] = projectCentreDistances(oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,:,freq),...
                                new_solver_setup, old_solver_setup, mm, nn ,mmLinkOldEdge,nnLinkOldEdge );
                            
                            %------- 1 edge on old edge, 1 edge parallel old edge ----
                            %  4 cases
                        elseif (mmLinkType == 0 && nnLinkType == 1  ) %2
                            %mm lies on old edge, nn lies in positive triangle
                            % mPlus_nPlus and mMinus_nPlus terms will have same
                            % shape
                             %same shape terms [1 2 5 6]
                             %reuse rest [3 4 7 8]
                            %calculatedNonSingInd(mm,nn) = 1;
                            [newAllTerms(mm,nn,:,freq)] = projectCentreDistances(oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,:,freq),...
                                new_solver_setup, old_solver_setup, mm, nn ,mmLinkOldEdge,nnLinkOldEdge );
                        elseif (mmLinkType == 0 && nnLinkType == -1  ) %3
                            % mPlus_nMinus and mMinus_nMinus
                             %same shape terms [3 4 7 8]
                             %reuse rest [1 2 5 6]
                            %calculatedNonSingInd(mm,nn) = 1;
                            [newAllTerms(mm,nn,:,freq)] = projectCentreDistances(oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,:,freq),...
                                new_solver_setup, old_solver_setup, mm, nn ,mmLinkOldEdge,nnLinkOldEdge );
                        elseif (mmLinkType == 1 && nnLinkType == 0  ) %4
                            % mPlus_nPlus and mPlus_nMinus
                             %same shape terms [1 2 3 4]
                             %reuse rest [5 6 7 8]
                            %calculatedNonSingInd(mm,nn) = 1;
                            [newAllTerms(mm,nn,:,freq)] = projectCentreDistances(oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,:,freq),...
                                new_solver_setup, old_solver_setup, mm, nn ,mmLinkOldEdge,nnLinkOldEdge );
                        elseif (mmLinkType == -1 && nnLinkType == 0  ) %5
                            % mMinus_nPlus and mMinus_nMinus
                             %same shape terms [5 6 7 8]
                             %reuse rest [1 2 3 4]
                            %calculatedNonSingInd(mm,nn) = 1;
                            [newAllTerms(mm,nn,:,freq)] = projectCentreDistances(oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,:,freq),...
                                new_solver_setup, old_solver_setup, mm, nn ,mmLinkOldEdge,nnLinkOldEdge );
                            
                            %------- both parallel to old edge ----
                            % 4 cases
                        elseif (mmLinkType == -1 && nnLinkType == -1 ) %6
                             %same shape terms [ 7 8]
                             %reuse rest [1 2 3 4 5 6]
                            %calculatedNonSingInd(mm,nn) = 1;
                            [newAllTerms(mm,nn,:,freq)] = projectCentreDistances(oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,:,freq),...
                                new_solver_setup, old_solver_setup, mm, nn ,mmLinkOldEdge,nnLinkOldEdge );
                        elseif (mmLinkType == -1 && nnLinkType == 1) %7
                             %same shape terms [ 5 6]
                             %reuse rest [1 2 3 4 7 8]
                            %calculatedNonSingInd(mm,nn) = 1;
                            [newAllTerms(mm,nn,:,freq)] = projectCentreDistances(oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,:,freq),...
                                new_solver_setup, old_solver_setup, mm, nn ,mmLinkOldEdge,nnLinkOldEdge );
                        elseif (mmLinkType == 1 && nnLinkType == -1) %8
                             %same shape terms [ 3 4]
                             %reuse rest [1 2 5 6 7 8]
                            %calculatedNonSingInd(mm,nn) = 1;
                            [newAllTerms(mm,nn,:,freq)] = projectCentreDistances(oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,:,freq),...
                                new_solver_setup, old_solver_setup, mm, nn ,mmLinkOldEdge,nnLinkOldEdge );
                        elseif (mmLinkType == 1 && nnLinkType == 1 ) %9
                             %same shape terms [ 1 2 ]
                             %reuse rest [3 4 5 6 7 8]
                            %calculatedNonSingInd(mm,nn) = 1;
                            [newAllTerms(mm,nn,:,freq)] = projectCentreDistances(oldAllTerms(mmLinkOldEdge,nnLinkOldEdge,:,freq),...
                                new_solver_setup, old_solver_setup, mm, nn ,mmLinkOldEdge,nnLinkOldEdge );
                        end % MAIN CASE 1: if (mmLinkType == 0 && nnLinkType == 0) %1 
                        
                    elseif (mmLinkType ~= 2 && nnLinkType ~= 2 && oldSingInd(mmLinkOldEdge,nnLinkOldEdge))
                        % old terms involve 2 or 3 unique triangles
                        
                    elseif ((mmLinkType == -1 || mmLinkType == 0 || mmLinkType == 1)&& nnLinkType == 2)
                        % nn: parallel to external edge
                        
                        % 1 : external, 2 : on old edge
                        %onOldEdge_parallelExternalNonSingInd(mm,nn) = 1;
                        
                        ind = find( newEdgeParallelExternalEdgeLinkOldInternalEdge(:, 1) == nn);
                        
                        %externalNNlinkedOldEdgeNN is the nearest old edge
                        %to nn
                        externalNNlinkedOldEdgeNN = newEdgeParallelExternalEdgeLinkOldInternalEdge(ind(1), 2);
                        
                        if (~oldSingInd(mmLinkOldEdge,externalNNlinkedOldEdgeNN))
                            
                            calculatedNonSingInd(mm,nn) = 1;
                            %newAllProperties(mm,nn, :) = calcProp(new_solver_setup, mm, nn );
                            [newAllTerms(mm,nn,:,freq)] = projectCentreDistances(oldAllTerms(mmLinkOldEdge,externalNNlinkedOldEdgeNN,:,freq),...
                                new_solver_setup, old_solver_setup, mm, nn ,mmLinkOldEdge,externalNNlinkedOldEdgeNN );
                            % need to also do some rotations
%                             newAllTerms(mm,nn,AInd,freq) = newAllTerms(mm,nn,AInd,freq) * ...
%                             cos(newAllProperties(mm,nn, 2))/ ...
%                             cos(oldAllProperties(mmLinkOldEdge,externalNNlinkedOldEdgeNN,2));
                        elseif (numel(ind) > 1)
                            prev_externalNNlinkedOldEdgeNN = externalNNlinkedOldEdgeNN;
                            externalNNlinkedOldEdgeNN = newEdgeParallelExternalEdgeLinkOldInternalEdge(ind(2), 2);
                            if (~oldSingInd(mmLinkOldEdge,externalNNlinkedOldEdgeNN))
                                calculatedNonSingInd(mm,nn) = 1;
                                [newAllTerms(mm,nn,:,freq)] = projectCentreDistances(oldAllTerms(mmLinkOldEdge,externalNNlinkedOldEdgeNN,:,freq),...
                                new_solver_setup, old_solver_setup, mm, nn ,mmLinkOldEdge,externalNNlinkedOldEdgeNN );
                            else
                                if (prev_externalNNlinkedOldEdgeNN ~= mmLinkOldEdge)
                                    % old group: 3 unique triangles
                                    calculatedNonSingInd(mm,nn) = 1;
                                    ind = find (oldSingInd(:, prev_externalNNlinkedOldEdgeNN) == 1);     
                                    ind2 = find ((ind(:) ~= mmLinkOldEdge) & (ind(:) ~= prev_externalNNlinkedOldEdgeNN));
                                    otherOldMM = ind(ind2(1));
                                    [newAllTerms(mm,nn,:,freq)] = projectCentreDistances(oldAllTerms(otherOldMM,prev_externalNNlinkedOldEdgeNN,:,freq),...
                                        new_solver_setup, old_solver_setup, mm, nn ,otherOldMM,prev_externalNNlinkedOldEdgeNN );
                                  
                                elseif (externalNNlinkedOldEdgeNN ~= mmLinkOldEdge)
                                    % old group: 3 unique triangles
                                    calculatedNonSingInd(mm,nn) = 1;
                                    ind = find (oldSingInd(:, externalNNlinkedOldEdgeNN) == 1);     
                                    ind2 = find ((ind(:) ~= mmLinkOldEdge) & (ind(:) ~= externalNNlinkedOldEdgeNN));
                                    otherOldMM = ind(ind2(1));
                                    [newAllTerms(mm,nn,:,freq)] = projectCentreDistances(oldAllTerms(otherOldMM,externalNNlinkedOldEdgeNN,:,freq),...
                                        new_solver_setup, old_solver_setup, mm, nn ,otherOldMM,externalNNlinkedOldEdgeNN );
                                else
                                   % b = 5; THIS LINE WILL NEVER EXECUTE
                                end
                            end
                            
                        else
                            if (externalNNlinkedOldEdgeNN ~= mmLinkOldEdge)
                                % old group: 3 unique triangles
                                calculatedNonSingInd(mm,nn) = 1;
                                ind = find (oldSingInd(:, externalNNlinkedOldEdgeNN) == 1);
                                ind2 = find ((ind(:) ~= mmLinkOldEdge) & (ind(:) ~= externalNNlinkedOldEdgeNN));
                                otherOldMM = ind(ind2(1));
                                [newAllTerms(mm,nn,:,freq)] = projectCentreDistances(oldAllTerms(otherOldMM,externalNNlinkedOldEdgeNN,:,freq),...
                                    new_solver_setup, old_solver_setup, mm, nn ,otherOldMM,externalNNlinkedOldEdgeNN );
                            else
                                % old group: 2 unique triangles
                                calculatedNonSingInd(mm,nn) = 1;
                                ind = find (oldSingInd(mmLinkOldEdge, :) == 1);
                                ind2 = find ((ind(:) ~= mmLinkOldEdge));
                                otherOldMM = ind(ind2);
                                
                                ind = find (oldSingInd(otherOldMM(1), :) == 1);
                                ind2 = find (  (ind(:) ~= otherOldMM(1)) & (ind(:) ~= otherOldMM(2)) & (ind(:) ~= mmLinkOldEdge)  );
                                otherOldMM = ind(ind2(1));
                                
                                [newAllTerms(mm,nn,:,freq)] = projectCentreDistances(oldAllTerms(otherOldMM,externalNNlinkedOldEdgeNN,:,freq),...
                                    new_solver_setup, old_solver_setup, mm, nn ,otherOldMM,externalNNlinkedOldEdgeNN );
                            end
                            
                        end
                        
                    elseif (mmLinkType == 2 && (nnLinkType == -1 || nnLinkType == 0 || nnLinkType == 1) )
                        % mm: parallel to external edge
                        
                        ind = find( newEdgeParallelExternalEdgeLinkOldInternalEdge(:, 1) == mm);
                        
                        %externalMMlinkedOldEdgeMM is the nearest old edge
                        %to mm
                        externalMMlinkedOldEdgeMM = newEdgeParallelExternalEdgeLinkOldInternalEdge(ind(1), 2);
                        
                        if (~oldSingInd(externalMMlinkedOldEdgeMM,nnLinkOldEdge))  
                            calculatedNonSingInd(mm,nn) = 1;
                            %newAllProperties(mm,nn, :) = calcProp(new_solver_setup, mm, nn );
                            [newAllTerms(mm,nn,:,freq)] = projectCentreDistances(oldAllTerms(externalMMlinkedOldEdgeMM,nnLinkOldEdge,:,freq),...
                                new_solver_setup, old_solver_setup, mm, nn ,externalMMlinkedOldEdgeMM,nnLinkOldEdge );
                        elseif (numel(ind) > 1)
                            prev_externalMMlinkedOldEdgeMM = externalMMlinkedOldEdgeMM;
                            externalMMlinkedOldEdgeMM = newEdgeParallelExternalEdgeLinkOldInternalEdge(ind(2), 2);
                            if (~oldSingInd(externalMMlinkedOldEdgeMM,nnLinkOldEdge))  
                                calculatedNonSingInd(mm,nn) = 1;
                                [newAllTerms(mm,nn,:,freq)] = projectCentreDistances(oldAllTerms(externalMMlinkedOldEdgeMM,nnLinkOldEdge,:,freq),...
                                    new_solver_setup, old_solver_setup, mm, nn ,externalMMlinkedOldEdgeMM,nnLinkOldEdge );
                            else
                                if (prev_externalMMlinkedOldEdgeMM ~= nnLinkOldEdge)
                                    % old group: 3 unique triangles
                                    calculatedNonSingInd(mm,nn) = 1;
                                    ind = find (oldSingInd(prev_externalMMlinkedOldEdgeMM, :) == 1);     
                                    ind2 = find ((ind(:) ~= nnLinkOldEdge) & (ind(:) ~= prev_externalMMlinkedOldEdgeMM));
                                    otherOldNN = ind(ind2(1));
                                    [newAllTerms(mm,nn,:,freq)] = projectCentreDistances(oldAllTerms(prev_externalMMlinkedOldEdgeMM,otherOldNN,:,freq),...
                                        new_solver_setup, old_solver_setup, mm, nn ,prev_externalMMlinkedOldEdgeMM,otherOldNN );
                                  
                                elseif (externalMMlinkedOldEdgeMM ~= nnLinkOldEdge)
                                    % old group: 3 unique triangles
                                    calculatedNonSingInd(mm,nn) = 1;
                                    ind = find (oldSingInd(externalMMlinkedOldEdgeMM, :) == 1);     
                                    ind2 = find ((ind(:) ~= nnLinkOldEdge) & (ind(:) ~= externalMMlinkedOldEdgeMM));
                                    otherOldNN = ind(ind2(1));
                                    [newAllTerms(mm,nn,:,freq)] = projectCentreDistances(oldAllTerms(externalMMlinkedOldEdgeMM,otherOldNN,:,freq),...
                                        new_solver_setup, old_solver_setup, mm, nn ,externalMMlinkedOldEdgeMM,otherOldNN );
                                else
                                   % b = 5; THIS LINE WILL NEVER EXECUTE
                                end
                            end
                        else
                            if (externalMMlinkedOldEdgeMM ~= nnLinkOldEdge)
                                % old group: 3 unique triangles
                                calculatedNonSingInd(mm,nn) = 1;
                                ind = find (oldSingInd(externalMMlinkedOldEdgeMM, :) == 1);
                                ind2 = find ((ind(:) ~= nnLinkOldEdge) & (ind(:) ~= externalMMlinkedOldEdgeMM));
                                otherOldNN = ind(ind2(1));
                                [newAllTerms(mm,nn,:,freq)] = projectCentreDistances(oldAllTerms(externalMMlinkedOldEdgeMM,otherOldNN,:,freq),...
                                    new_solver_setup, old_solver_setup, mm, nn ,externalMMlinkedOldEdgeMM,otherOldNN );
                            else
                                % old group: 2 unique triangles
                                calculatedNonSingInd(mm,nn) = 1;
                                ind = find (oldSingInd(:, nnLinkOldEdge) == 1);
                                ind2 = find ((ind(:) ~= nnLinkOldEdge));
                                otherOldNN = ind(ind2);
                                
                                ind = find (oldSingInd(:, otherOldNN(1)) == 1);
                                ind2 = find (  (ind(:) ~= otherOldNN(1)) & (ind(:) ~= otherOldNN(2)) & (ind(:) ~= nnLinkOldEdge)  );
                                otherOldNN = ind(ind2(1));
                                
                                [newAllTerms(mm,nn,:,freq)] = projectCentreDistances(oldAllTerms(externalMMlinkedOldEdgeMM,otherOldNN,:,freq),...
                                    new_solver_setup, old_solver_setup, mm, nn ,externalMMlinkedOldEdgeMM,otherOldNN );
                            end
                        end
                        
                    elseif (mmLinkType == 2 && nnLinkType == 2)
                        
                        ind = find( newEdgeParallelExternalEdgeLinkOldInternalEdge(:, 1) == mm);
                        %oldEdgeMM1 = newEdgeParallelExternalEdgeLinkOldInternalEdge(ind(1), 2);
                        oldEdgeMM = newEdgeParallelExternalEdgeLinkOldInternalEdge(ind, 2);
                        
                        ind = find( newEdgeParallelExternalEdgeLinkOldInternalEdge(:, 1) == nn);
                        oldEdgeNN = newEdgeParallelExternalEdgeLinkOldInternalEdge(ind, 2);
                        
                        for k = 1: numel(oldEdgeMM)
                            for j = 1:numel(oldEdgeNN)
                                if (calculatedNonSingInd(mm,nn) == 1)
                                    break
                                end
                                if (~oldSingInd(oldEdgeMM(k),oldEdgeNN(j)))
                                    calculatedNonSingInd(mm,nn) = 1;
                                    [newAllTerms(mm,nn,:,freq)] = projectCentreDistances(oldAllTerms(oldEdgeMM(k),oldEdgeNN(j),:,freq),...
                                        new_solver_setup, old_solver_setup, mm, nn ,oldEdgeMM(k),oldEdgeNN(j) );
                                    break
                                end
                            end
                        end
%                         oldEdgeMM2 = 0;
%                         if (numel(ind) > 1)
%                             oldEdgeMM2 = newEdgeParallelExternalEdgeLinkOldInternalEdge(ind(2), 2);
%                         end
                            

                    end % if (mmLinkType ~= 2 && nnLinkType ~= 2 && ~oldSingInd(mmLinkOldEdge,nnLinkOldEdge))
                
                    newAllTerms(mm,nn,AInd,freq) = 0.125* newAllTerms(mm,nn,AInd,freq); %0.125 works best
                    newAllTerms(mm,nn,PhiInd,freq) = 0.5* newAllTerms(mm,nn,PhiInd,freq); % 0.25 works best
                    
                end % if (mm == nn)
            end %for nn= 1:numNewEdges
            
        end %for mm =1:numNewEdges
    end %for freq = 1:numFreq

    %=========SECOND LOOP - ESTIMATED FROM FIRST LOOP =======
    %reuse new terms generated in first loop
    for freq = 1:numFreq
        for mm =1:numNewEdges
            %mmLinkOldEdge = newEdgeLinkOldEdge(mm,1);
            %mmLinkType = newEdgeLinkOldEdge(mm,2);
            for nn= 1:numNewEdges
               % nnLinkOldEdge = newEdgeLinkOldEdge(nn,1);
               % nnLinkType = newEdgeLinkOldEdge(nn,2);
                
                if (onOldEdge_ParallelExternalTriInd(mm,nn) == 1)
                    %addTriangles method ensures that new triangle plus of
                    %a new edge which is parallel to an external edge
                    % is never central (4th) triangle
                    % i.e. the free vertex of triangle plus is an old node,
                    % not a new node in the center of an external edge
                    ind = find(onOldEdge_OnOldEdgeTriInd(mm, :) == 1);
                    
                    if (numel(ind)> 0)
                        %mm is on old edge, nn is parallel external edge
                        onOldEdgeNN = ind(1);
                        newAllTerms(mm,nn,:,freq) = compareAndSwapTerms_3unique(new_solver_setup, newAllTerms(mm,onOldEdgeNN,:,freq),mm,onOldEdgeNN, mm, nn );
                    else
                        ind = find(onOldEdge_OnOldEdgeTriInd(:, nn) == 1);
                        %mm is parallel to external edge, nn is on old edge
                        onOldEdgeMM = ind(1);
                        newAllTerms(mm,nn,:,freq) = compareAndSwapTerms_3unique(new_solver_setup, newAllTerms(onOldEdgeMM,nn,:,freq), onOldEdgeMM,nn, mm, nn );
                    end
                end % if (onOldEdgeParallelExternalTriInd(mm,nn) == 1)
                
                if (parallelInternal_parallelExternalTriInd(mm,nn) == 1)
                    
                    %on old edge nn indices
                    ind = find(onOldEdge_ParallelExternalTriInd(mm, :) == 1);
                    %ind = find(onOldEdgeOnOldEdgeTriInd(mm, :) == 1);
                    
                    if (numel(ind)> 0)
                        %mm is parallel external edge,
                        %nn is parallel internal edge
                        onOldEdgeNN = ind(1);
                        ind = find(onOldEdge_OnOldEdgeTriInd(:, onOldEdgeNN) == 1);
                        onOldEdgeMM = ind(1);
                        newAllTerms(mm,nn,:,freq) = compareAndSwapTerms_3unique(new_solver_setup, newAllTerms(onOldEdgeMM,onOldEdgeNN,:,freq),onOldEdgeMM,onOldEdgeNN, mm, nn );
                    else
                        ind = find(onOldEdge_ParallelExternalTriInd(:, nn) == 1);
                        %mm is parallel internal edge,
                        %nn is parallel external edge
                        onOldEdgeMM = ind(1);
                        ind = find(onOldEdge_OnOldEdgeTriInd(onOldEdgeMM, :) == 1);
                        onOldEdgeNN = ind(1);
                        newAllTerms(mm,nn,:,freq) = compareAndSwapTerms_3unique(new_solver_setup, newAllTerms(onOldEdgeMM,onOldEdgeNN,:,freq), onOldEdgeMM,onOldEdgeNN, mm, nn );
                    end
                end % if (parallelInternalparallelExternalTriInd(mm,nn) == 1)
                
                if (parallelExternal_parallelExternalTriInd(mm,nn) == 1)
                    
                    ind = find(onOldEdge_ParallelExternalTriInd(mm, :) == 1);
                    %ind = find(onOldEdgeOnOldEdgeTriInd(mm, :) == 1);
                    
                    if (numel(ind)> 0)
                        %mm is parallel external edge,
                        %nn is parallel external edge
                        onOldEdgeNN = ind(1);
                        ind = find(onOldEdge_OnOldEdgeTriInd(:, onOldEdgeNN) == 1);
                        onOldEdgeMM = ind(1);
                        newAllTerms(mm,nn,:,freq) = compareAndSwapTerms_3unique(new_solver_setup, newAllTerms(onOldEdgeMM,onOldEdgeNN,:,freq),onOldEdgeMM,onOldEdgeNN, mm, nn );
                    else
                        %This code below will never execute? to be tested
                        %mm is parallel external edge,
                        %nn is parallel external edge
                        ind = find(onOldEdge_ParallelExternalTriInd(:, nn) == 1);
                        onOldEdgeMM = ind(1);
                        ind = find(onOldEdge_OnOldEdgeTriInd(onOldEdgeMM, :) == 1);
                        onOldEdgeNN = ind(1);
                        newAllTerms(mm,nn,:,freq) = compareAndSwapTerms_3unique(new_solver_setup, newAllTerms(onOldEdgeMM,onOldEdgeNN,:,freq), onOldEdgeMM,onOldEdgeNN, mm, nn );
                    end
                end % if (parallelExternal_parallelExternalTriInd(mm,nn) == 1)
                

            end %for nn= 1:numNewEdges
            
        end %for mm =1:numNewEdges
    end %for freq = 1:numFreq

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
        case 2
            swappedTerms = terms([5 6 7 8 1 2 3 4]);
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

function [swappedTerms] = compareAndSwapTerms_3unique(new_solver_setup, terms, termsM, termsN, m, n )
    termsMPlus = new_solver_setup.rwg_basis_functions_trianglePlus(termsM);
    termsMMinus = new_solver_setup.rwg_basis_functions_triangleMinus(termsM);
    termsNPlus = new_solver_setup.rwg_basis_functions_trianglePlus(termsN);
    termsNMinus = new_solver_setup.rwg_basis_functions_triangleMinus(termsN);
    
    %identify terms associated with middle triangle
    if (termsMPlus == termsNPlus)
        termCodeOne = 1;
    elseif (termsMPlus == termsNMinus)
        termCodeOne = 2;
    elseif (termsMMinus == termsNPlus)
        termCodeOne = 3;
    elseif (termsMMinus == termsNMinus)
        termCodeOne = 4;
    else
        swappedTerms = terms;
        return;
    end
    
    %-----
    mPlus = new_solver_setup.rwg_basis_functions_trianglePlus(m);
    mMinus = new_solver_setup.rwg_basis_functions_triangleMinus(m);
    nPlus = new_solver_setup.rwg_basis_functions_trianglePlus(n);
    nMinus = new_solver_setup.rwg_basis_functions_triangleMinus(n);
    
    %identify terms associated with middle triangle
    if (mPlus == nPlus)
        codeTwo = 1;
    elseif (mPlus == nMinus)
        codeTwo = 2;
    elseif (mMinus == nPlus)
        codeTwo = 3;
    elseif (mMinus == nMinus)
        codeTwo = 4;
    else
        swappedTerms = terms;
        return;
    end
    
    swappedTerms = swapTermCodes(terms, termCodeOne, codeTwo );
end

function [newTerms] = projectCentreDistances(oldTerms, new_solver_setup, old_solver_setup, mm, nn ,mmLinkOldEdge,nnLinkOldEdge )
    
    newTerms = zeros(1,1,8);

    mmPlusCentre = new_solver_setup.triangle_centre_point(new_solver_setup.rwg_basis_functions_trianglePlus(mm),:);
    mmMinusCentre = new_solver_setup.triangle_centre_point(new_solver_setup.rwg_basis_functions_triangleMinus(mm),:);
    nnPlusCentre =  new_solver_setup.triangle_centre_point(new_solver_setup.rwg_basis_functions_trianglePlus(nn),:);
    nnMinusCentre = new_solver_setup.triangle_centre_point(new_solver_setup.rwg_basis_functions_triangleMinus(nn),:);
    
    oldmmPlusCentre = old_solver_setup.triangle_centre_point(old_solver_setup.rwg_basis_functions_trianglePlus(mmLinkOldEdge),:);
    oldmmMinusCentre = old_solver_setup.triangle_centre_point(old_solver_setup.rwg_basis_functions_triangleMinus(mmLinkOldEdge),:);
    oldnnPlusCentre =  old_solver_setup.triangle_centre_point(old_solver_setup.rwg_basis_functions_trianglePlus(nnLinkOldEdge),:);
    oldnnMinusCentre = old_solver_setup.triangle_centre_point(old_solver_setup.rwg_basis_functions_triangleMinus(nnLinkOldEdge),:);
    
    mmPlus_nnPlus = norm(mmPlusCentre-nnPlusCentre);
    mmPlus_nnMinus = norm(mmPlusCentre-nnMinusCentre);
    mmMinus_nnPlus= norm(mmMinusCentre-nnPlusCentre);
    mmMinus_nnMinus = norm(mmMinusCentre-nnMinusCentre);
    
    oldmmPlus_nnPlus = norm(oldmmPlusCentre-oldnnPlusCentre);
    oldmmPlus_nnMinus = norm(oldmmPlusCentre-oldnnMinusCentre);
    oldmmMinus_nnPlus= norm(oldmmMinusCentre-oldnnPlusCentre);
    oldmmMinus_nnMinus = norm(oldmmMinusCentre-oldnnMinusCentre);
    
    % NB NB 1 INDICES DEPENDANT ON ORDER OF TERMS
    % NB NB 2 Triangle plus/minus labels are already alligned
    newTerms(1) = oldmmPlus_nnPlus/mmPlus_nnPlus;
    newTerms(3) = oldmmPlus_nnMinus/mmPlus_nnMinus;
    newTerms(5) = oldmmMinus_nnPlus/mmMinus_nnPlus;
    newTerms(7) = oldmmMinus_nnMinus/mmMinus_nnMinus; 
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