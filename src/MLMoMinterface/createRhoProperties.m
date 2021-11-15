%function prop = calcProp(Solver_setup, mmDir,nnDir,mmEdgeCenter, nnEdgeCenter )
function [rhoProperties] = createRhoProperties(Solver_setup)

    numEdges = Solver_setup.num_mom_basis_functions;
    
    %----
    edgeDir = Solver_setup.nodes_xyz(Solver_setup.rwg_basis_functions_shared_edge_nodes(:,1), :)...
        - Solver_setup.nodes_xyz(Solver_setup.rwg_basis_functions_shared_edge_nodes(:,2), :);
    %----
    %properties = zeros(numEdges,numEdges,3);
    rhoProperties = zeros(numEdges,numEdges,4);
    rhoCPls = Solver_setup.rho_c_pls;
    rhoCMns = Solver_setup.rho_c_mns;
    for mm = 1:numEdges
        %for nn = 1:numEdges
        for nn = 1:numEdges
            %----
%             rhoProperties(mm,nn,1) = calcAngle(rhoCPls(mm,:),rhoCPls(nn,:) );
%             rhoProperties(mm,nn,2) = calcAngle(rhoCPls(mm,:),rhoCMns(nn,:) );
%             rhoProperties(mm,nn,3) = calcAngle(rhoCMns(mm,:),rhoCPls(nn,:) );
%             rhoProperties(mm,nn,4) = calcAngle(rhoCMns(mm,:),rhoCMns(nn,:) );
            %----
            rhoProperties(mm,nn,1) = calcAngle(edgeDir(mm,:),rhoCPls(mm,:) );
            rhoProperties(mm,nn,2) = calcAngle(edgeDir(mm,:),rhoCMns(mm,:) );
            rhoProperties(mm,nn,3) = calcAngle(edgeDir(nn,:),rhoCPls(nn,:) );
            rhoProperties(mm,nn,4) = calcAngle(edgeDir(nn,:),rhoCMns(nn,:) );
            %----
%             %1st  symmetry
%             if (edgeCentreProperties(mm,nn,2) < 0)
%                 edgeCentreProperties(mm,nn,2) = edgeCentreProperties(mm,nn,2) + pi;
%             end
%             %avoid divide by zero
%             if (edgeCentreProperties(mm,nn,1) ~= 0)
%                 edgeCentreProperties(mm,nn,3) = calcAngle(edgeDir(mm,:),disp );
%                 % 1st symmetry
%                 if (edgeCentreProperties(mm,nn,3)< 0)
%                     edgeCentreProperties(mm,nn,3) = edgeCentreProperties(mm,nn,3) + pi;
%                 end
%                 % 2nd symmetry
%                 if (edgeCentreProperties(mm,nn,3) > pi/2)
%                     edgeCentreProperties(mm,nn,2) = pi - edgeCentreProperties(mm,nn,2) ;
%                     edgeCentreProperties(mm,nn,3) = pi - edgeCentreProperties(mm,nn,3);
%                 end
%             end
        end
    end
    % prop: 1 dist, 2 orientation, 3 radial
    
   
end