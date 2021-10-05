
% 6 nodes
nodesXYZ = zeros(6, 3);
nodesXYZ(1,:) = [-1 0 0];
nodesXYZ(2,:) = [1 0 0];
nodesXYZ(3,:) = [0 sqrt(3) 0];
nodesXYZ(4,:) = [0 -sqrt(3) 0];
nodesXYZ(5,:) = [2 sqrt(3) 0];
nodesXYZ(6,:) = [-2 sqrt(3) 0];

% 3 internal edges
intEdgeNodes = zeros(3,2);
intEdgeNodes(1,:) = [1 2];
intEdgeNodes(2,:) = [2 3];
intEdgeNodes(3,:) = [3 1];

intEdgeLengths = zeros(3,1);
intEdgeLengths(1) =  2;
intEdgeLengths(2) =  2;
intEdgeLengths(3) =  2;

intEdgeCentres = zeros(3,3);
intEdgeCentres(1, :) =  0.5*nodesXYZ(intEdgeNodes(1,1) , :) + 0.5*nodesXYZ(intEdgeNodes(1,2) , :);
intEdgeCentres(2, :) =  0.5*nodesXYZ(intEdgeNodes(2,1) , :) + 0.5*nodesXYZ(intEdgeNodes(2,2) , :);
intEdgeCentres(3, :) =  0.5*nodesXYZ(intEdgeNodes(3,1) , :) + 0.5*nodesXYZ(intEdgeNodes(3,2) , :);


% 4 triangles
triangleVertices = zeros(4,3);
triangleVertices(1, :) = [1 2 3];
triangleVertices(2, :) = [1 2 4];
triangleVertices(3, :) = [2 3 5];
triangleVertices(4, :) = [3 1 6];

%---
trianglePlus = zeros(3, 1);
trianglePlus(1) = 1; % 1
trianglePlus(2) = 1;
trianglePlus(3) = 1;

trianglePlusFreeVertex = zeros(3, 1);
trianglePlusFreeVertex(1) = 3; %3
trianglePlusFreeVertex(2) = 1;
trianglePlusFreeVertex(3) = 2;

%---

triangleMinus = zeros(3, 1);
triangleMinus(1) = 2; % 2
triangleMinus(2) = 3;
triangleMinus(3) = 4;

triangleMinusFreeVertex = zeros(3, 1);
triangleMinusFreeVertex(1) = 4; %4
triangleMinusFreeVertex(2) = 5;
triangleMinusFreeVertex(3) = 6;


%---------------------

Solver_setup = [];
Solver_setup.frequencies = [];
Solver_setup.num_metallic_edges = 3;
Solver_setup.num_metallic_triangles = 4;

Solver_setup.nodes_xyz = nodesXYZ;
Solver_setup.triangle_vertices = triangleVertices;
Solver_setup.rwg_basis_functions_shared_edge_nodes = intEdgeNodes;
Solver_setup.rwg_basis_functions_length_m = intEdgeLengths;
Solver_setup.rwg_basis_functions_shared_edge_centre = intEdgeCentres;
Solver_setup.rwg_basis_functions_trianglePlus =trianglePlus;
Solver_setup.rwg_basis_functions_triangleMinus = triangleMinus;
Solver_setup.rwg_basis_functions_trianglePlusFreeVertex = trianglePlusFreeVertex;
Solver_setup.rwg_basis_functions_triangleMinusFreeVertex = triangleMinusFreeVertex;

%plot(nodesXYZ(:,1), nodesXYZ(:,2), '.', 'markerSize', 20);
[new_solver_setup] = addTriangles(Solver_setup);
%[new_solver_setup] = addTriangles(new_solver_setup);
%plot(new_solver_setup.nodes_xyz(:,1), new_solver_setup.nodes_xyz(:,2), '.', 'markerSize', 20);


