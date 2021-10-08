%function prop = calcProp(Solver_setup, mmDir,nnDir,mmEdgeCenter, nnEdgeCenter )
function properties = calcProperties(Solver_setup)

    numEdges = Solver_setup.num_mom_basis_functions;
    
    edgeDir = Solver_setup.nodes_xyz(Solver_setup.rwg_basis_functions_shared_edge_nodes(:,1), :)...
        - Solver_setup.nodes_xyz(Solver_setup.rwg_basis_functions_shared_edge_nodes(:,2), :);
    edgeDirMag = sqrt(edgeDir(:,1).^2 + edgeDir(:,2).^2 + edgeDir(:,3).^2);
    edgeDir = edgeDir ./edgeDirMag;
    
    properties = zeros(numEdges,numEdges,3);
    
    for mm = 1:numEdges
        %for nn = 1:numEdges
        for nn = 1:numEdges
            
            disp = Solver_setup.rwg_basis_functions_shared_edge_centre(nn, :) -...
                Solver_setup.rwg_basis_functions_shared_edge_centre(mm, :);
            
            properties(mm,nn,1) =  norm(disp);
            properties(mm,nn,2) = calcAngle(edgeDir(mm,:),edgeDir(nn,:) );
            
            %1st  symmetry
            if (properties(mm,nn,2) < 0)
                properties(mm,nn,2) = properties(mm,nn,2) + pi;
            end
            %avoid divide by zero
            if (properties(mm,nn,1) ~= 0)
                properties(mm,nn,3) = calcAngle(edgeDir(mm,:),disp );
                % 1st symmetry
                if (properties(mm,nn,3)< 0)
                    properties(mm,nn,3) = properties(mm,nn,3) + pi;
                end
                % 2nd symmetry
                if (properties(mm,nn,3) > pi/2)
                    properties(mm,nn,2) = pi - properties(mm,nn,2) ;
                    properties(mm,nn,3) = pi - properties(mm,nn,3);
                end
            end
        end
    end
    % prop: 1 dist, 2 orientation, 3 radial
    
   
end