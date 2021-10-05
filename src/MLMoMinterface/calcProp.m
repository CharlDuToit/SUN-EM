%function prop = calcProp(Solver_setup, mmDir,nnDir,mmEdgeCenter, nnEdgeCenter )
function prop = calcProp(Solver_setup, mm, nn )
    % prop: 1 dist, 2 orientation, 3 radial
    mmDir = Solver_setup.nodes_xyz(Solver_setup.rwg_basis_functions_shared_edge_nodes(mm,1), :) - ...
        Solver_setup.nodes_xyz(Solver_setup.rwg_basis_functions_shared_edge_nodes(mm,2),: );
    mmDir = mmDir./norm(mmDir);
    
    nnDir = Solver_setup.nodes_xyz(Solver_setup.rwg_basis_functions_shared_edge_nodes(nn,1), :) - ...
        Solver_setup.nodes_xyz(Solver_setup.rwg_basis_functions_shared_edge_nodes(nn,2),: );
    nnDir = nnDir./norm(nnDir);
    
    mmEdgeCenter = Solver_setup.rwg_basis_functions_shared_edge_centre(mm, :);
    nnEdgeCenter = Solver_setup.rwg_basis_functions_shared_edge_centre(nn, :); 
    
    %Solver_setup.nodes_xyz
    prop = zeros(1,1,3);
    disp = nnEdgeCenter - mmEdgeCenter;
    prop(1,1,1) =  norm(disp);
    prop(1,1,2) = calcAngle(mmDir,nnDir );
    
    %1st  symmetry
    if (prop(1,1,2) < 0)
        prop(1,1,2) = prop(1,1,2) + pi;
    end
    %avoid divide by zero
    if (prop(1,1,1) ~= 0)
        prop(1,1,3) = calcAngle(mmDir,disp );
        % 1st symmetry
        if (prop(1,1,3)< 0)
            prop(1,1,3) = prop(1,1,3) + pi;
        end
        % 2nd symmetry
        if (prop(1,1,3) > pi/2)
            prop(1,1,2) = pi - prop(1,1,2) ;
            prop(1,1,3) = pi - prop(1,1,3);
        end
    end
end