function [mlmom] = runMLMoMsolver(Const, Solver_setup, zMatrices)
    %runMLMoMsolver
    %   Usage:
    %       [mlmom] = runMoMsolver(Const, Solver_setup, zMatrices)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       Solver_setup
    %           Solver specific struct, e.g. frequency range, basis function details, geometry details    
    %       zMatrices
    %           The reference Z-matrix data. This MUST be from FEKO (extracted from the *.mat file)
    %           or some other accurate reference simulation data.
    %
    %   Output Arguments:
    %       mlmom
    %           Structs containing trained ML-MoM 
    %
    %   Description:
    %
    %       Runs the ML-MoM solution based on the Z that was read / parsed.
    %       Geometric assumptions:
    %       1. Planar body
    %       2. Consistent mesh size
    %       Algorithm:
    %       For each frequency, split into real and imaginary MLR models
    %       1. Performs quadrature integration and stores 8 terms
    %       2. Group data by 2 , 3  and 4 unique triangles
    %       3. For 4 triangles, performs clustering
    %       4. Use MLR to find optimal term weights of each group and cluster
    %
    %       Trained model can be used to predict a different Solver_setup,
    %       but with the same frequency (and mesh size? to be tested)
    %
    %   =======================
    %   Written by Charl du Toit on August 20, 2021.
    %   Stellenbosch University
    %   Email: 21708886@sun.ac.za

    %narginchk(5,5);

    message_fc(Const,' ');
    message_fc(Const,'------------------------------------------------------------------------------------');
    message_fc(Const,sprintf('Running ML-MoM solver'));
    
    % ------------ INIT MLMOM ------------
    
    % Variable convention:
    % dataset + singularism + ?
    %
    % dataset : ref, pred, unity
    % singularism : nonSing, sing ( 2 or 3 triangles), self, tri
    % ? : Context should be clear, first 2 parts could be missing as well
    
    
    
    
    % ------------ EXTRACT  ------------ 
    
    mlmom = [];
    [selfZmnTerms, triZmnTerms, nonSingZmnTerms, nonSingZmnProp, singInd] = extractZmnInfo(Const, Solver_setup);
    [refSelfZmn, refTriZmn, refNonSingZmn] = extractRefZmn(zMatrices.values, singInd );
    edgeLengths = Solver_setup.rwg_basis_functions_length_m;
    
    numFreq = zMatrices.numFreq;
    weightModels = cell(numFreq, 2); % real and imaginary models
    
    % ------------ REAL AND IMAG MODEL FOR EACH FREQ  ------------ 
    
    for i = 1:numFreq
        weightModels{i, 1} = calcWeightsStruct(refSelfZmn(:, i), refTriZmn(:, i), refNonSingZmn(:, i),...
    selfZmnTerms(:,:, i), triZmnTerms(:,:, i), nonSingZmnTerms(:,:, i), nonSingZmnProp,edgeLengths, 1);

        weightModels{i, 2} = calcWeightsStruct(refSelfZmn(:, i), refTriZmn(:, i), refNonSingZmn(:, i),...
    selfZmnTerms(:,:, i), triZmnTerms(:,:, i), nonSingZmnTerms(:,:, i), nonSingZmnProp,edgeLengths, 0);

    end
    
    % ------------ UPDATE STRUCT  ------------ 
    
    mlmom.weightModels =weightModels;
    
end
    

    