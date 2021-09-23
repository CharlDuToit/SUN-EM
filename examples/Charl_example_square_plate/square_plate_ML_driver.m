clear;
Const = sunem_initialise('square_plate',false);
Const.FEKOmatfilename          = 'square_plate.mat'; 
Const.FEKOstrfilename          = 'square_plate.str';
Const.FEKOrhsfilename          = 'square_plate.rhs'; % ?
Const.FEKOoutfilename          = 'square_plate.out'; % 
Const.FEKOefefilename          = 'square_plate.efe'; % ?
Const.FEKOffefilename          = 'square_plate.ffe'; % ?
Const.runMLMoMsolver              = true;

[Const, zMatrices, yVectors, xVectors] = extractFEKOMoMmatrixEq(Const);
[Const, Solver_setup] = parseFEKOoutfile(Const, yVectors);

Const.MLMoMClusterSizeScale = 1;
Const.QUAD_PTS = 1;
[Solution] = runEMsolvers(Const, Solver_setup, zMatrices, yVectors, xVectors);

mlmom = Solution.mlmom;


