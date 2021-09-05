function [mlmom] = runMLMoMsolver(Const, Solver_setup, zMatrices)
    %runMoMsolver
    %   Usage:
    %       [mlmom] = runMoMsolver(Const, Solver_setup, zMatrices)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       Solver_setup
    %           Solver specific struct, e.g. frequency range, basis function details, geometry details    
    %       zMatrices
    %           The reference Z-matrix data. This can be from FEKO (extracted from the *.mat file.
    %
    %   Output Arguments:
    %       mlmom
    %           Structs containing ML-MoM trained neural network and timing data
    %
    %   Description:
    %       Runs the ML-MoM solution based on the Z that was read / parsed.
    %       from the FEKO *.out, *.mat, *.str and *.rhs files
    %       Z is the training input data, Y is the training output data
    %
    %   Temporary development comments:
    %       Train model based on one frequency
    %
    %   =======================
    %   Written by Charl du Toit on August 20, 2021.
    %   Stellenbosch University
    %   Email: 21708886@sun.ac.za

    %narginchk(5,5);

    message_fc(Const,' ');
    message_fc(Const,'------------------------------------------------------------------------------------');
    message_fc(Const,sprintf('Running ML-MoM solver'));

    % Initialise the return values
    mlmom  = [];
    mlmom.name = 'mlmom';
    %Nmom = Solver_setup.num_mom_basis_functions;   % Total number of basis functions for whole problem
    %numSols = xVectors.numSols;                   % The number of solutions configurations
    %mlmom.numSols = 1; %numSols;                     % For now, set to 1. (TO-DO: Update)
    numFreq = Solver_setup.frequencies.freq_num;   % The number of frequency points to process
    %numRHSperFreq = mlmom.numSols / numFreq;         % The number of solutions per frequency point.
                                                   % For now, should be 1 (TO-DO: Update)
                                                   
    % Some info about the solution configurations
    %message_fc(Const,sprintf('  numSols : %d', mlmom.numSols));
    message_fc(Const,sprintf('  numFreq : %d', numFreq));
    %message_fc(Const,sprintf('  numRHSperFreq : %d', numRHSperFreq));

    % Calculate the solution vector (all frequency points, all RHSes)
    %mlmom.Isol = complex(zeros(Nmom,mlmom.numSols));
    
    mlmom.nonSingZmnTerms = zeros(0);
    mlmom.refNonSingZmn= zeros(0);
    mlmom.nonSingZmnUnityWeight = zeros(0);

    % The timing calculations also need to take into account that there is a
    % frequency loop
    %mlmom.setupTime = zeros(1,numFreq);
    %mlmom.factorisationTime = zeros(1,numFreq);
    % Zero also the total times (for all frequency iterations)
    mlmom.totsetupTime = 0.0;
    %mlmom.totfactorisationTime = 0.0;
    mlmom.totsolTime = 0.0;
    mlmom.excludeSelfTerms = 1;
    
    % 8 terms for each observation
    [mlmom.nonSingZmnTerms , singIndices] = extractNonSingZmnTerms(Const, Solver_setup);
    % reference Zmn values
    mlmom.refNonSingZmn = extractNonSingZmn(zMatrices.values, singIndices );
    % Same estimation as internal solver
    mlmom.nonSingZmnUnityWeight  = sum(mlmom.nonSingZmnTerms, 2);
    
    %training
    mlmom.weights = ( (mlmom.nonSingZmnTerms' * mlmom.nonSingZmnTerms) \ mlmom.nonSingZmnTerms')* mlmom.refNonSingZmn ;
    
    %predict
    mlmom.predNonSingZmn = mlmom.nonSingZmnTerms * mlmom.weights;
    diff = mlmom.refNonSingZmn - mlmom.predNonSingZmn;
    mlmom.predError = diff'*diff;
    
    %Error
    diff = mlmom.refNonSingZmn - mlmom.nonSingZmnUnityWeight;
    mlmom.unityWeightError = diff'*diff;
  
    
end
    

    