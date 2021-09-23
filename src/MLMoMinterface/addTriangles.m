function [new_solver_setup] = addTriangles(Solver_setup)
     
    % numbers
    numOldIntEdges = Solver_setup.num_metallic_edges;  
    [numOldNodes, ~] = size(Solver_setup.nodes_xyz);
    numOldTriangles = Solver_setup.num_metallic_triangles;
    nanInd = find(isnan(Solver_setup.nodes_xyz(:,1)));
    if (numel(nanInd) > 0)
        numOldNodes = nanInd(1) -1 ;
    end
    
    % append old nodesXYZ
    nodesXYZ = zeros(numOldNodes *3, 3);
    nodesXYZ(1:numOldNodes, :) = Solver_setup.nodes_xyz(1:numOldNodes,:);
    
    %old 
    oldTriangleVertices = Solver_setup.triangle_vertices;
    oldIntEdgeNodes = Solver_setup.rwg_basis_functions_shared_edge_nodes;
    oldIntEdgeLengths = Solver_setup.rwg_basis_functions_length_m;
    oldIntEdgeCentres = Solver_setup.rwg_basis_functions_shared_edge_centre;
    oldTrianglePlus = Solver_setup.rwg_basis_functions_trianglePlus;
    oldTriangleMinus = Solver_setup.rwg_basis_functions_triangleMinus;
    oldTrianglePlusFreeVertex = Solver_setup.rwg_basis_functions_trianglePlusFreeVertex;
    oldTriangleMinusFreeVertex = Solver_setup.rwg_basis_functions_triangleMinusFreeVertex;

    % new
    newTriangleVertices = zeros(numOldTriangles * 4, 3); % node labels
    newIntEdgeNodes = zeros(numOldIntEdges * 6, 2); % node labels
    newIntEdgeLengths = zeros(numOldIntEdges * 6, 1); % meter
    newIntEdgeCentres = zeros(numOldIntEdges * 6, 3); % xyz
    newTrianglePlus = zeros(numOldIntEdges * 6, 1); % triangle labels
    newTriangleMinus = zeros(numOldIntEdges * 6, 1); % triangle labels
    newTrianglePlusFreeVertex = zeros(numOldIntEdges * 6, 1); % node labels
    newTriangleMinusFreeVertex = zeros(numOldIntEdges * 6, 1); % node labels
    
    %links to old data
    % newEdgeLinkOldEdge(: , 1): old edge label
    % newEdgeLinkOldEdge(: , 2): link type : 
    % -1 : new original edge in Traingle Minus
    %  0 : new edge ontop of half of old edge
    %  1 : new original edge in Triangle Plus
    %  2 : new edge does not lie on half of old edge, and old edge is
    %       external
    newEdgeLinkOldEdge = zeros(numOldIntEdges * 4, 2); % old edge label, link type
    
    %links to old data
    %oldEdgeLinkFreeVertexTri(1) : old edge
    %oldEdgeLinkFreeVertexTri(2) : side : 1 = oldIntEdgeNodes(1), 2 =oldIntEdgeNodes(2)
    %oldEdgeLinkFreeVertexTri(3) : Sign : 1 = minus, 2 = plus
    %oldEdgeLinkFreeVertexTri(4) : Info : 1 = Free Vertex, 2 = Triangle label
   % oldEdgeLinkFreeVertexTri = zeros(numOldIntEdges * 4, 2, 2 , 2);
    
    oldEdgeLinkOverlapNewEdge = zeros(numOldIntEdges, 2); % new edge on index 1, new edge on index 2
    
    %links to old data
    oldEdgeLinkNewNode = zeros(numOldIntEdges , 1);
    
    % tri plus and tri minus for each edge has been checked
    % oldEdgeChecked(:, 1) = 0 : Triangle Minus has NOT been checked
    % oldEdgeChecked(:, 1) = 1 : Triangle Minus HAS checked
    % oldEdgeChecked(:, 2) = 0 : Triangle Plus has NOT been checked
    % oldEdgeChecked(:, 2) = 1 : Triangle Plus HAS been checked
    oldEdgeChecked = zeros(numOldIntEdges, 2);
    
   [oldTriangleEdges, oldExtEdgeNodes, oldExtEdgeCentres, oldExtEdgeLengths, oldExtEdgeOppositeNodes] =...
       verticesToEdges(oldTriangleVertices, oldIntEdgeNodes, nodesXYZ);
   
   newIntEdgeCount = 0;
   newTriangleCount = 0;
   nodeCount = numOldNodes;
   for oldEdgeInd = 1:numOldIntEdges
       
       for sign = 1:2
           
           %has negative/postive triangle of old edge been checked?
           if (oldEdgeChecked(oldEdgeInd, sign) == 0)
               
               if (sign == 1)
                   oldTri = oldTriangleMinus(oldEdgeInd);
               else
                   oldTri = oldTrianglePlus(oldEdgeInd);
               end
               %oldTri = oldTriangleMinus(oldEdgeInd);
               
               % 3 edge labels, one of these edges = oldEdgeInd
               if (oldTri == 0)
                   abc = 1; % error from previous iteration. will crash
               end
               oldTriEdges = oldTriangleEdges(oldTri, :);
               
               %Create or lookup the 3 nodes associated with old edge tri
               newNodesLabels = zeros(3,1);
               for k = 1:3
                   e = oldTriEdges(k);
                   if (e > numOldIntEdges) % external edge
                       nodeCount = nodeCount +1;
                       nodesXYZ(nodeCount, :) = oldExtEdgeCentres(e - numOldIntEdges, :);
                       newNodesLabels(k) = nodeCount;
                   elseif ( oldEdgeLinkNewNode(e) == 0) % edge has not been checked
                       nodeCount = nodeCount +1;
                       nodesXYZ(nodeCount, :) = oldIntEdgeCentres(e, :);
                       oldEdgeLinkNewNode(e) = nodeCount;
                       newNodesLabels(k) = nodeCount;
                   else % old edge already has a new node
                       newNodesLabels(k) = oldEdgeLinkNewNode(e);
                   end
               end % for k
               
               % guaranteed to Create 3 new triangles and edges
               fourthTriangleIndex = newTriangleCount + 4;
               for k = 1:3
                   %e1 = oldTriEdges(k);
                   e2 = oldTriEdges(mod(k, 3) +1 ) ;
                   e3 = oldTriEdges(mod(k +1, 3) +1 ) ;                   
                   
                   newIntEdgeCount = newIntEdgeCount + 1;
                   newIntEdgeNodes(newIntEdgeCount, :) = [newNodesLabels(k) newNodesLabels(mod(k ,3) + 1) ];
                   newIntEdgeCentres(newIntEdgeCount, :) = 0.5*nodesXYZ(newNodesLabels(k),:) + 0.5*nodesXYZ(newNodesLabels(mod(k ,3) + 1),:);
                   newIntEdgeLengths(newIntEdgeCount, 1) = norm(nodesXYZ(newNodesLabels(k),:) - nodesXYZ(newNodesLabels(mod(k ,3) + 1),:));
                   
                   %third vertex must still be determined
                   newTriangleCount = newTriangleCount +1;
                   newTriangleVertices(newTriangleCount, 1:2) = [newNodesLabels(k) newNodesLabels(mod(k ,3) + 1) ];
                   
                   % old edge parallel to new edge
                   newEdgeLinkOldEdge(newIntEdgeCount,1) = e3 ;
                    
                   if (newIntEdgeCount == 47)
                       abc = 5;
                   end
                   
                   if(e3 > numOldIntEdges) % old external edge parallel to new edge
                       newEdgeLinkOldEdge(newIntEdgeCount,2) = 2 ;
                       newTriangleVertices(newTriangleCount, 3) = oldExtEdgeOppositeNodes(e3 - numOldIntEdges);
                       % always assign triangle to triangle plus if parallel
                       % old edge is external (arbitrary)
                       newTrianglePlus(newIntEdgeCount) = newTriangleCount;
                       newTriangleMinus(newIntEdgeCount) = fourthTriangleIndex;
                       
                       newTrianglePlusFreeVertex(newIntEdgeCount) = newTriangleVertices(newTriangleCount, 3);
                       newTriangleMinusFreeVertex(newIntEdgeCount) = newNodesLabels(mod(k +1, 3) +1 ) ;
                   else % old internal edge parallel to new edge
                       %oldTriMinus (from oldEdgeInd loop) must be pos or neg
                       %triangle of other 2 old edges that make up oldTriMinus
                       if (oldTriangleMinus(e3) == oldTri) % Minus
                           newEdgeLinkOldEdge(newIntEdgeCount,2) = -1 ;
                           oldEdgeChecked(e3, 1) = 1;
                           newTriangleVertices(newTriangleCount, 3) = oldTriangleMinusFreeVertex(e3);
                           
                           newTriangleMinusFreeVertex(newIntEdgeCount) = oldTriangleMinusFreeVertex(e3);
                           newTrianglePlusFreeVertex(newIntEdgeCount) = newNodesLabels(mod(k +1 ,3) + 1);
                           
                           newTrianglePlus(newIntEdgeCount) = fourthTriangleIndex;
                           newTriangleMinus(newIntEdgeCount) = newTriangleCount;

                           
                       else % plus
                           if (oldTrianglePlus(e3) ~= oldTri)
                               abc = 5;
                           end
                           newEdgeLinkOldEdge(newIntEdgeCount,2) = 1 ;
                           oldEdgeChecked(e3, 2) = 1;
                           newTriangleVertices(newTriangleCount, 3) = oldTrianglePlusFreeVertex(e3);
                           
                           newTrianglePlusFreeVertex(newIntEdgeCount) = oldTrianglePlusFreeVertex(e3);
                           newTriangleMinusFreeVertex(newIntEdgeCount) = newNodesLabels(mod(k +1 ,3) + 1);
                           
                           newTriangleMinus(newIntEdgeCount) = fourthTriangleIndex;
                           newTrianglePlus(newIntEdgeCount) = newTriangleCount;
                           
                       end % if (oldTriangleMinus(e3) == oldTri) % Minus
                   end %if(e3 > numOldIntEdges)
                   
                   if (e2 <= numOldIntEdges )
                       nextTriangle = 1;
                       if (k ==3)
                           nextTriangle = -2;
                       end
                       
                       if (oldTriangleMinus(e2) == oldTri) % Minus
                           e2Sign =1;
                       else
                           if (oldTrianglePlus(e2) ~= oldTri)
                               abc = 5;
                           end
                           e2Sign =2;
                       end
                       
                       %if (sign == 1)
                       if (oldEdgeLinkOverlapNewEdge(e2, 1) == 0) % edge has not been created
                           
                           oldEdgeCentreNode = oldEdgeLinkNewNode(e2);
                           
                           newIntEdgeCount = newIntEdgeCount + 1;
                           newEdge1 = newIntEdgeCount;
                           newIntEdgeNodes(newIntEdgeCount, :) = [oldIntEdgeNodes(e2, 1) oldEdgeCentreNode ];
                           newIntEdgeCentres(newIntEdgeCount, :) = 0.5*nodesXYZ(oldIntEdgeNodes(e2, 1),:) + 0.5*nodesXYZ(oldEdgeCentreNode,:);
                           newIntEdgeLengths(newIntEdgeCount, 1) = norm(nodesXYZ(oldIntEdgeNodes(e2, 1),:) - nodesXYZ(oldEdgeCentreNode,:));
                           
                           if (newIntEdgeCount == 47)
                               abc = 5;
                           end
                           
                           newIntEdgeCount = newIntEdgeCount + 1;
                           newEdge2 = newIntEdgeCount;
                           newIntEdgeNodes(newIntEdgeCount, :) = [oldIntEdgeNodes(e2, 2) oldEdgeCentreNode ];
                           newIntEdgeCentres(newIntEdgeCount, :) = 0.5*nodesXYZ(oldIntEdgeNodes(e2, 2),:) + 0.5*nodesXYZ(oldEdgeCentreNode,:);
                           newIntEdgeLengths(newIntEdgeCount, 1) = norm(nodesXYZ(oldIntEdgeNodes(e2, 2),:) - nodesXYZ(oldEdgeCentreNode,:));
                           
                           if (newIntEdgeCount == 47)
                               abc = 5;
                           end
                           
                           if (newTriangleVertices(newTriangleCount, 3) == oldIntEdgeNodes(e2,1)) % index 1
                               oldEdgeLinkOverlapNewEdge(e2, 1) = newEdge1;
                               oldEdgeLinkOverlapNewEdge(e2, 2) = newEdge2;
                               
                           else % index 2
                               oldEdgeLinkOverlapNewEdge(e2, 2) = newEdge1;
                               oldEdgeLinkOverlapNewEdge(e2, 1) = newEdge2;
                           end
                           
                           newEdgeLinkOldEdge(newEdge1,:) = [e2  0];
                           newEdgeLinkOldEdge(newEdge2,:) = [e2  0];

                           %newTriangleMinus(newIntEdgeCount) = newTriangleCount + nextTriangle;
                           %newTriangleMinusFreeVertex(newIntEdgeCount) = newNodesLabels(mod(k+1,3) +1);
                           
                       end % if (oldEdgeLinkOverlapNewEdge(e2, 1) == 0) % edge has not been created
                       
                       
                       if (newTriangleVertices(newTriangleCount, 3) == oldIntEdgeNodes(e2,1)) % index 1
                           newEdge1 = oldEdgeLinkOverlapNewEdge(e2, 1);
                           newEdge2 = oldEdgeLinkOverlapNewEdge(e2, 2);
                           
                       else % index 2
                           newEdge1 = oldEdgeLinkOverlapNewEdge(e2, 2);
                           newEdge2 = oldEdgeLinkOverlapNewEdge(e2, 1);
                       end
                       
                       if (newEdge2 == 47 || newEdge1 == 47)
                           abc = 5;
                       end
                       
                       if (e2Sign == 1) % neg
                           newTriangleMinus(newEdge1) = newTriangleCount;
                           newTriangleMinusFreeVertex(newEdge1) = newNodesLabels(k);
                           
                           newTriangleMinus(newEdge2) = newTriangleCount + nextTriangle;
                           newTriangleMinusFreeVertex(newEdge2) = newNodesLabels(mod(k+1,3) +1);
                           if (newTriangleMinus(newEdge2) == 0)
                               abc = 1;
                           end
                       else
                           newTrianglePlus(newEdge1) = newTriangleCount;
                           newTrianglePlusFreeVertex(newEdge1) = newNodesLabels(k);
                           
                           newTrianglePlus(newEdge2) = newTriangleCount + nextTriangle;
                           newTrianglePlusFreeVertex(newEdge2) = newNodesLabels(mod(k+1,3) +1);
                       end

                   end %if e2 <= numOldIntEdges
               end % for k = 1:3
               
               %Create centre triangle. Already has Plus/ Minus assigned
               newTriangleCount = newTriangleCount +1;
               newTriangleVertices(newTriangleCount, :) = [newNodesLabels(1) newNodesLabels(2) newNodesLabels(3) ];
               
           end % triangle pls/minus checked?
           
       end% for sign
       
       %Create 2 edges on half of oldEdgeInd
%        oldEdgeCentreNode = oldEdgeLinkNewNode(oldEdgeInd);
%        newIntEdgeCount = newIntEdgeCount + 1;
%        newIntEdgeNodes(newIntEdgeCount, :) = [oldIntEdgeNodes(oldEdgeInd, 1) oldEdgeCentreNode ];
%        newIntEdgeCentres(newIntEdgeCount, :) = 0.5*nodesXYZ(oldIntEdgeNodes(oldEdgeInd, 1),:) + 0.5*nodesXYZ(oldEdgeCentreNode,:);
%        newIntEdgeLengths(newIntEdgeCount, 1) = norm(nodesXYZ(oldIntEdgeNodes(oldEdgeInd, 1),:) - nodesXYZ(oldEdgeCentreNode,:));
%        
%        newTriangleMinusFreeVertex(newIntEdgeCount) = oldEdgeLinkFreeVertexTri(oldEdgeInd,1 ,1,1);
%        newTrianglePlusFreeVertex(newIntEdgeCount) = oldEdgeLinkFreeVertexTri(oldEdgeInd,1 ,2,1);
%        newTriangleMinus(newIntEdgeCount) = oldEdgeLinkFreeVertexTri(oldEdgeInd,1 ,1,2 );
%        newTrianglePlus(newIntEdgeCount) = oldEdgeLinkFreeVertexTri(oldEdgeInd,1 ,2,2 );
% 
%        % Second
%        newIntEdgeCount = newIntEdgeCount + 1;
%        newIntEdgeNodes(newIntEdgeCount, :) = [oldIntEdgeNodes(oldEdgeInd, 2) oldEdgeCentreNode ];
%        newIntEdgeCentres(newIntEdgeCount, :) = 0.5*nodesXYZ(oldIntEdgeNodes(oldEdgeInd, 2),:) + 0.5*nodesXYZ(oldEdgeCentreNode,:);
%        newIntEdgeLengths(newIntEdgeCount, 1) = norm(nodesXYZ(oldIntEdgeNodes(oldEdgeInd, 2),:) - nodesXYZ(oldEdgeCentreNode,:));
%        
%        newTriangleMinusFreeVertex(newIntEdgeCount) = oldEdgeLinkFreeVertexTri(oldEdgeInd,2 ,1,1);
%        newTrianglePlusFreeVertex(newIntEdgeCount) = oldEdgeLinkFreeVertexTri(oldEdgeInd,2 ,2,1);
%        newTriangleMinus(newIntEdgeCount) = oldEdgeLinkFreeVertexTri(oldEdgeInd,2 ,1,2 );
%        newTrianglePlus(newIntEdgeCount) = oldEdgeLinkFreeVertexTri(oldEdgeInd,2 ,2,2 );
       
   end % for oldEdgeInd
   nodesXYZ = nodesXYZ(1:nodeCount, :);
   newTriangleVertices = newTriangleVertices(1:newTriangleCount, :);
   newIntEdgeNodes = newIntEdgeNodes(1:newIntEdgeCount, :);
   newIntEdgeLengths = newIntEdgeLengths(1:newIntEdgeCount, :);
   newIntEdgeCentres = newIntEdgeCentres(1:newIntEdgeCount, :);
   newTrianglePlus = newTrianglePlus(1:newIntEdgeCount, :);
   newTriangleMinus = newTriangleMinus(1:newIntEdgeCount, :);
   newTrianglePlusFreeVertex = newTrianglePlusFreeVertex(1:newIntEdgeCount, :);
   newTriangleMinusFreeVertex = newTriangleMinusFreeVertex(1:newIntEdgeCount, :);

    new_solver_setup = [];
    new_solver_setup.num_metallic_edges = newIntEdgeCount;
    new_solver_setup.num_metallic_triangles = newTriangleCount; 
    new_solver_setup.nodes_xyz = nodesXYZ;
    new_solver_setup.triangle_vertices = newTriangleVertices;
    new_solver_setup.rwg_basis_functions_shared_edge_nodes = newIntEdgeNodes;
    new_solver_setup.rwg_basis_functions_length_m = newIntEdgeLengths;
    new_solver_setup.rwg_basis_functions_shared_edge_centre = newIntEdgeCentres;
    new_solver_setup.rwg_basis_functions_trianglePlus = newTrianglePlus;
    new_solver_setup.rwg_basis_functions_triangleMinus = newTriangleMinus;
    new_solver_setup.rwg_basis_functions_trianglePlusFreeVertex = newTrianglePlusFreeVertex;
    new_solver_setup.rwg_basis_functions_triangleMinusFreeVertex = newTriangleMinusFreeVertex;

end


function [triangleEdges, extEdgeNodes, extEdgeCentres, extEdgeLengths, extEdgeOppositeNodes] = verticesToEdges(triangleVertices, intEdgeNodes, nodesXYZ)
    %input
    % triangleVertices: labels [numTri, 3]
    % intEdgeNodes: labels [numIntEdge, 2]
    % nodesXYZ : xyz [numNodes, 3]

    % output
    % triangleEdges: labels [numTri, 3]
    % extEdgeNodes: labels [numExtEdge, 2]
    % extEdgeCentres: xyz [numExtEdge, 3]
    % extEdgeLengths: m [numExtEdge, 1]
    % extEdgeOppositeNodes : labels [numExtEdge, 1]
    
    numTri = numel(triangleVertices(:,1));
    numIntEdge = numel(intEdgeNodes(:,1));
    triangleEdges = zeros(numTri, 3);
    extEdgeNodes = zeros(numIntEdge *4, 2);
    extEdgeCentres = zeros(numIntEdge *4, 3);
    extEdgeLengths = zeros(numIntEdge *4, 1);
    extEdgeOppositeNodes = zeros(numIntEdge *4, 1);
    
    extEdgeCount = 0;
    for triCount = 1:numTri
        for node = 1:3
            a = node;
            b = mod(node,3) +1;
            
            edgeLabel = find( ( intEdgeNodes(:,1) == triangleVertices(triCount,a) & intEdgeNodes(:,2) == triangleVertices(triCount,b) ) |...
                (intEdgeNodes(:,1) == triangleVertices(triCount,b) & intEdgeNodes(:,2) == triangleVertices(triCount,a)) );
            if (numel(edgeLabel) > 0) % internal edge
                triangleEdges(triCount, node) = edgeLabel(1);
            else % external edge
                c = mod(node+1,3) + 1;
                extEdgeCount = extEdgeCount + 1;
                triangleEdges(triCount, node) = numIntEdge + extEdgeCount;
                extEdgeNodes(extEdgeCount, :) = [triangleVertices(triCount,a)  triangleVertices(triCount,b)];
                extEdgeCentres(extEdgeCount, :) = 0.5*nodesXYZ(extEdgeNodes(extEdgeCount, 1) , :) + 0.5*nodesXYZ(extEdgeNodes(extEdgeCount, 2) , :);
                extEdgeLengths(extEdgeCount, 1) = norm(nodesXYZ(extEdgeNodes(extEdgeCount, 1) , :) - nodesXYZ(extEdgeNodes(extEdgeCount, 2) , :));
                extEdgeOppositeNodes(extEdgeCount, 1) = triangleVertices(triCount,c);
            end
        end % for node
    end % for triCount
    extEdgeNodes = extEdgeNodes(1:extEdgeCount, :);
    extEdgeCentres = extEdgeCentres(1:extEdgeCount, :);
    extEdgeLengths = extEdgeLengths(1:extEdgeCount, 1);
    extEdgeOppositeNodes = extEdgeOppositeNodes(1:extEdgeCount, 1);
end

