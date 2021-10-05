%clear;
Const = sunem_initialise('two_triangles',false);
Const.FEKOmatfilename          = 'two_triangles.mat'; 
Const.FEKOstrfilename          = 'two_triangles.str';
Const.FEKOrhsfilename          = 'two_triangles.rhs'; % ?
Const.FEKOoutfilename          = 'two_triangles.out'; % 
Const.FEKOefefilename          = 'two_triangles.efe'; % ?
Const.FEKOffefilename          = 'two_triangles.ffe'; % ?
Const.runMLMoMsolver              = true;

[Const, zMatrices, yVectors, xVectors] = extractFEKOMoMmatrixEq(Const);
[Const, Solver_setup] = parseFEKOoutfile(Const, yVectors);
Const.QUAD_PTS = 12;
%Const.MLMoMClusterSizeScale = 1;
%Const.MLMoMMinPercentImprov = 0;
%Const.MLMoMIncludeRealCalc = 0;
%[Solution] = runEMsolvers(Const, Solver_setup, zMatrices, yVectors, xVectors);

%mlmom_oldSolverSetup = Solution.mlmom;

%plot(Solver_setup.nodes_xyz(:,1), Solver_setup.nodes_xyz(:,2), '.', 'markerSize', 20);
[new_solver_setup, newEdgeLinkOldEdge, newEdgeParallelExternalEdgeLinkOldInternalEdge] = addTriangles(Solver_setup);
numNewEdges = new_solver_setup.num_mom_basis_functions;
%plot(new_solver_setup.nodes_xyz(:,1), new_solver_setup.nodes_xyz(:,2), '.', 'markerSize', 20);
[terms, singIndices, properties] = fillZmnTermsByEdge(Const,Solver_setup) ;

[newAllTerms, newAllProperties, newSingInd, calculatedNonSingInd, calculatedTriInd, calculatedSelfInd] = ...
    projectOldSolverSetup(new_solver_setup, Solver_setup,newEdgeLinkOldEdge, newEdgeParallelExternalEdgeLinkOldInternalEdge,  terms, properties, singIndices );

[zMatrices_newSolverSetup_fillByEdge]=FillZMatrixByEdge(Const,new_solver_setup) ;
[terms, singIndices, properties] = fillZmnTermsByEdge(Const,new_solver_setup) ;

zMN_projected = zeros(numNewEdges,numNewEdges);

numCalcNonSing = numel(find(calculatedNonSingInd == 1));
numCalcTri = numel(find(calculatedTriInd == 1));
numCalcSelf = numel(find(calculatedSelfInd == 1));

zMN_projected_calcNonSing = zeros(numCalcNonSing, 1);
zMN_projected_calcTri = zeros(numCalcTri,1);
zMN_projected_calcSelf = zeros(numCalcSelf,1);

zMN_fillbyEdge_calcNonSing = zeros(numCalcNonSing, 1);
zMN_fillbyEdge_calcTri = zeros(numCalcTri, 1);
zMN_fillbyEdge_calcSelf = zeros(numCalcSelf, 1);

calcSingCount = 0;
calcTriCount = 0;
calcSelfCount = 0;
for mm = 1:numNewEdges
    for nn = 1:numNewEdges
        zMN_projected(mm,nn) = sum(newAllTerms(mm,nn,:));
        if (calculatedNonSingInd(mm,nn))
            calcSingCount = calcSingCount +1;
            zMN_projected_calcNonSing(calcSingCount) = zMN_projected(mm,nn);
            zMN_fillbyEdge_calcNonSing(calcSingCount) = zMatrices_newSolverSetup_fillByEdge.values(mm,nn);
        end
        if (calculatedTriInd(mm,nn))
            calcTriCount = calcTriCount + 1;
            zMN_projected_calcTri(calcTriCount) = zMN_projected(mm,nn);
            zMN_fillbyEdge_calcTri(calcTriCount) = zMatrices_newSolverSetup_fillByEdge.values(mm,nn);
        end
        if (calculatedSelfInd(mm,nn))
            calcSelfCount = calcSelfCount + 1;
            zMN_projected_calcSelf(calcSelfCount) = zMN_projected(mm,nn);
            zMN_fillbyEdge_calcSelf(calcSelfCount) = zMatrices_newSolverSetup_fillByEdge.values(mm,nn);
        end
        
    end
end

%diag_zMN_projected = diag(zMN_projected);
%diag_zMatrices_newSolverSetup_fillByEdge = diag(zMatrices_newSolverSetup_fillByEdge.values);