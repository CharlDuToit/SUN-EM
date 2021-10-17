
mmInd = mlmomAddTriangles.groupIndices.twoUnique_pos_indices(:,1);
nnInd =  mlmomAddTriangles.groupIndices.twoUnique_pos_indices(:,2);
ref = zeros(numel(mmInd),1);
proj = zeros(numel(mmInd),1);
for k = 1:numel(mmInd)
    ref(k,1) = mlmomAddTriangles.refZmn(mmInd(k), nnInd(k), 1);
    proj(k,1) = mlmomAddTriangles.projZmn(mmInd(k), nnInd(k), 1);
end
ratio = zeros(numel(mmInd),2);
ratio(:,1) = real(proj) ./ real(ref);
ratio(:,2) = imag(proj) ./ imag(ref);
% a = [];
% if (eq(a,[]))
%     b = 5;
% else
%     b= 3;
% end

%activePlots = [0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0];
%activePlots = [0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0];
%activePlots = [0 0 0 0 0 0 1];
%plotMLMOM(mlmom, activePlots, 3, 1, 673);
 %[numClusters, ~]=size(mlmom.clusterMeans);
 %clusterCounts = mlmom.clusterCounts;
 %[~ , ind] = sort(clusterCounts);
% for k = 1:numClusters
     %plotMLMOM(mlmom, activePlots, 500, 1, ind(numClusters - k + 1));
 %end


% clear;
% Const = sunem_initialise('square_plate',false);
% Const.FEKOmatfilename          = 'square_plate.mat'; 
% Const.FEKOstrfilename          = 'square_plate.str';
% Const.FEKOrhsfilename          = 'square_plate.rhs'; % ?
% Const.FEKOoutfilename          = 'square_plate.out'; % 
% Const.FEKOefefilename          = 'square_plate.efe'; % ?
% Const.FEKOffefilename          = 'square_plate.ffe'; % ?
% Const.runMLMoMsolver              = true;
% 
% [Const, zMatrices, yVectors, xVectors] = extractFEKOMoMmatrixEq(Const);
% [Const, Solver_setup] = parseFEKOoutfile(Const, yVectors);
% 
%  tic;
%  Const.QUAD_PTS = 1;
%  [zMatrices1point]=FillZMatrixByEdge(Const,Solver_setup) ;
%  t = toc;
 %tic;
% %Const.QUAD_PTS = 6;
% %[zMatrices3point]=FillZMatrixByEdge(Const,Solver_setup) ;
% %b = toc;
% %[Solution] = runEMsolvers(Const, Solver_setup, zMatrices, yVectors, xVectors);