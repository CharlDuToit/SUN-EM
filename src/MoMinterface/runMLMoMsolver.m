function [mlmom] = runMLMoMsolver(Const, Solver_setup, zMatrices, yVectors, refIsol)
    %runMoMsolver
    %   Usage:
    %       [mom] = runMoMsolver(Const)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       Solver_setup
    %           Solver specific struct, e.g. frequency range, basis function details, geometry details    
    %       zMatrices
    %           The Z-matrices data. This can be from FEKO (extracted from the *.mat file, or internally
    %           calculated).
    %       yVectors
    %           The Yrhs-vector data
    %       refIsol
    %           The reference solution-vector data (e.g. MoM solution of FEKO or SUN-EM)
    %
    %   Output Arguments:
    %       mlmom
    %           Structs containing ML-MoM solution and timing data
    %
    %   Description:
    %       Runs the ML-MoM solution based on the Z and Y data that was read / parsed.
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

    narginchk(5,5);

    message_fc(Const,' ');
    message_fc(Const,'------------------------------------------------------------------------------------');
    message_fc(Const,sprintf('Running ML-MoM solver'));

    % Initialise the return values
    mlmom  = [];
    mlmom.name = 'mlmom';
    Nmom = Solver_setup.num_mom_basis_functions;   % Total number of basis functions for whole problem
    %numSols = xVectors.numSols;                   % The number of solutions configurations
    mlmom.numSols = 1; %numSols;                     % For now, set to 1. (TO-DO: Update)
    numFreq = Solver_setup.frequencies.freq_num;   % The number of frequency points to process
    numRHSperFreq = mlmom.numSols / numFreq;         % The number of solutions per frequency point.
                                                   % For now, should be 1 (TO-DO: Update)
                                                   
    % Some info about the solution configurations
    message_fc(Const,sprintf('  numSols : %d', mlmom.numSols));
    message_fc(Const,sprintf('  numFreq : %d', numFreq));
    message_fc(Const,sprintf('  numRHSperFreq : %d', numRHSperFreq));

    % Calculate the solution vector (all frequency points, all RHSes)
    mlmom.Isol = complex(zeros(Nmom,mlmom.numSols));

    % The timing calculations also need to take into account that there is a
    % frequency loop
    mlmom.setupTime = zeros(1,numFreq);
    mlmom.factorisationTime = zeros(1,numFreq);
    % Zero also the total times (for all frequency iterations)
    mlmom.totsetupTime = 0.0;
    mlmom.totfactorisationTime = 0.0;
    mlmom.totsolTime = 0.0;
    
% ===== MODULE FUNCTIONS ========
function [x,t] = extract_training_data(zMatrices, yVectors)
    
end

function network = construct_network(hiddenSizes,trainFcn)
    
end

function network = train_network(network, x, t)
    % might have to implement own back propagation algorithm
end
    
% Map Minimum and Maximum Input Processing Function
function y = mapminmax_apply(x,settings)
y = bsxfun(@minus,x,settings.xoffset);
y = bsxfun(@times,y,settings.gain);
y = bsxfun(@plus,y,settings.ymin);
end

% Sigmoid Symmetric Transfer Function
function a = tansig_apply(n,~)
a = 2 ./ (1 + exp(-2*n)) - 1;
end

% Map Minimum and Maximum Output Reverse-Processing Function
function x = mapminmax_reverse(y,settings)
x = bsxfun(@minus,y,settings.ymin);
x = bsxfun(@rdivide,x,settings.gain);
x = bsxfun(@plus,x,settings.xoffset);
end


    