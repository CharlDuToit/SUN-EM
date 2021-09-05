% Author: Danie Ludick (dludick@sun.ac.za)
% Project: Square plate
%
% Note: Each project directory / example directory needs to have a sunem_initialise.m
% script, that is used to setup the correct environment for that example.
%
% Refer to the /doc folder for more information

% --------------------------------------------------------------------------------------------------
% Initialise the environment
% --------------------------------------------------------------------------------------------------
% Project output directory: './square_plate/'
% Debug: True/False
Const = sunem_initialise('square_plate',false);

% --------------------------------------------------------------------------------------------------
% Program flow settings
% --------------------------------------------------------------------------------------------------

% Choose the solvers that will be executed
Const.runMoMsolver              = true;

% --------------------------------------------------------------------------------------------------
% Define input files for extracting FEKO data
% --------------------------------------------------------------------------------------------------
Const.FEKOmatfilename          = 'square_plate.mat'; 
Const.FEKOstrfilename          = 'square_plate.str';
Const.FEKOrhsfilename          = 'square_plate.rhs'; % ?
Const.FEKOoutfilename          = 'square_plate.out'; % 
Const.FEKOefefilename          = 'square_plate.efe'; % ?
Const.FEKOffefilename          = 'square_plate.ffe'; % ?

% The Following file is used to port solutions to FEKO 
% (for post-processing in POSTFEKO).
% TO-DO: [DL] Add this.
% Const.output_strfilename    = '';
% Const.writeFEKOstrfile = [0 0 0 0];

% --------------------------------------------------------------------------------------------------
% Read the MoM matrix equation from the file
% --------------------------------------------------------------------------------------------------
[Const, zMatrices, yVectors, xVectors] = extractFEKOMoMmatrixEq(Const);

% -- Read Z
%[zMatrices] = readFEKOZMatrixFromFile(Const, Const.FEKOmatfilename);

% -- Read X (FEKO MoM solution)
%[xVectors] = readFEKOXvectorFromFile(Const, Const.FEKOstrfilename);

% --------------------------------------------------------------------------------------------------
% Parse the setup files to extract the frequency sweep, the geometry and basis function setup 
% --------------------------------------------------------------------------------------------------
% TO-DO: At a later stage we can also add other meshing / geometry
% preprocessxing, e.g. Gmsh or GiD. For now the solver setup is read from FEKO.
[Const, Solver_setup] = parseFEKOoutfile(Const, yVectors);

% --------------------------------------------------------------------------------------------------
% Run the EM solver
% --------------------------------------------------------------------------------------------------
[Solution] = runEMsolvers(Const, Solver_setup, zMatrices, yVectors, xVectors);

% --------------------------------------------------------------------------------------------------
% Postprocess the results, e.g. calculate the Electric field
% --------------------------------------------------------------------------------------------------

figure;
hold on;
grid on;
samples = [xVectors.Isol(32);xVectors.Isol(39);xVectors.Isol(35);xVectors.Isol(3);
    xVectors.Isol(109);xVectors.Isol(88);xVectors.Isol(96);xVectors.Isol(80)];
samples  = abs(samples);
num_samples = 8;
x_val = [1; 2;3;4;5;6;7;8];
x_val = x_val - 0.5;
x_val = x_val/53.33;
%plot(1:num_samples,samples,'LineWidth',3);
plot(x_val,samples,'LineWidth',3);
%legend('SUN-EM','FEKO (*.efe file)', 'FEKO (*.ffe file)');
set(get(gca, 'XLabel'), 'String', ('y [m]'));
set(get(gca, 'YLabel'), 'String', ('|Current| [A]'))

%plot(1:FEKO_total_farfield_samples,FEKO_farfield_magnitude./max_FEKO_farfield_magnitude,'o','LineWidth',3);

