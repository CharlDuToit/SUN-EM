%function prop = calcProp(Solver_setup, mmDir,nnDir,mmEdgeCenter, nnEdgeCenter )
function [edgeCentreProperties] = calcProperties(Solver_setup)

    numEdges = Solver_setup.num_mom_basis_functions;
    
    edgeDir = Solver_setup.nodes_xyz(Solver_setup.rwg_basis_functions_shared_edge_nodes(:,1), :)...
        - Solver_setup.nodes_xyz(Solver_setup.rwg_basis_functions_shared_edge_nodes(:,2), :);
    %edgeDirMag = sqrt(edgeDir(:,1).^2 + edgeDir(:,2).^2 + edgeDir(:,3).^2);
    %edgeDir = edgeDir ./edgeDirMag;
    
    %properties = zeros(numEdges,numEdges,3);
    edgeCentreProperties = zeros(numEdges,numEdges,3);
    
    for mm = 1:numEdges
        %for nn = 1:numEdges
        for nn = 1:numEdges
            
            disp = Solver_setup.rwg_basis_functions_shared_edge_centre(nn, :) -...
                Solver_setup.rwg_basis_functions_shared_edge_centre(mm, :);
            
            edgeCentreProperties(mm,nn,1) =  norm(disp);
            edgeCentreProperties(mm,nn,2) = calcAngle(edgeDir(mm,:),edgeDir(nn,:) );
            
            %1st  symmetry
            if (edgeCentreProperties(mm,nn,2) < 0)
                edgeCentreProperties(mm,nn,2) = edgeCentreProperties(mm,nn,2) + pi;
            end
            %avoid divide by zero
            if (edgeCentreProperties(mm,nn,1) ~= 0)
                edgeCentreProperties(mm,nn,3) = calcAngle(edgeDir(mm,:),disp );
                % 1st symmetry
                if (edgeCentreProperties(mm,nn,3)< 0)
                    edgeCentreProperties(mm,nn,3) = edgeCentreProperties(mm,nn,3) + pi;
                end
                % 2nd symmetry
                if (edgeCentreProperties(mm,nn,3) > pi/2)
                    edgeCentreProperties(mm,nn,2) = pi - edgeCentreProperties(mm,nn,2) ;
                    edgeCentreProperties(mm,nn,3) = pi - edgeCentreProperties(mm,nn,3);
                end
            end
        end
    end
    % prop: 1 dist, 2 orientation, 3 radial
    
   
end