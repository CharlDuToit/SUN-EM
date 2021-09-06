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
    
    % Variable convention
    % dataset + singularism + ?
    %
    % dataset : ref, pred, unity
    % singularism : nonSing, sing ( 2 or 3 triangles), self, tri
    % ? : Context should be clear, first 2 parts could be missing as well
    % 
    mlmom  = [];
    mlmom.name = 'mlmom';      
    mlmom.nonSingZmnTerms = [];
    mlmom.refNonSingZmn= [];
    mlmom.unityNonSingZmn = [];
    mlmom.singInd = [];
    mlmom.clusterInd = [];
    mlmom.nonSingWeights = [];
    mlmom.predNonSingZmn = [];     
    mlmom.nonSingZmnProp = [];
    mlmom.dirDotDirInterval = [];
    mlmom.dirDotDispInterval =[];
    mlmom.distInterval = [];
    
    mlmom.totsetupTime = 0.0;
    mlmom.totsolTime = 0.0; 
    mlmom.numNonSing = 0;
    mlmom.numTerms = 0;
    mlmom.predNonSingMSE = 0;
    mlmom.unityNonSingMSE = 0;
    mlmom.predNonSingNormMSE = 0;
    mlmom.unityNonSingNormMSE = 0;    
    mlmom.threshDist = 0;
    mlmom.numNoClass = 0;
    
    % ------------ EXTRACT  ------------
    
    [nonSingZmnTerms ,nonSingZmnProp, singInd] = extractNonSingZmnInfo(Const, Solver_setup);
    [numNonSing , numTerms] = size(nonSingZmnTerms);
    
    % ------------ NON SINGULAR CLUSTERING  ------------
    
    refNonSingZmn = extractZmnRows(zMatrices.values, singInd );  
    unityNonSingZmn  = sum(nonSingZmnTerms, 2); % Same estimation as internal solver
    
% Threshold regression
%                 in = [0.1 0.005; 0.1 0.01; 0.1 0.015 ; 0.2 0.01 ; 0.2 0.02; 0.2 0.03;...
%                     0.4 0.02; 0.4 0.04; 0.4 0.06; 0.8 0.04; 0.8 0.08; 0.8 0.12; 1.6 0.08; 1.6 0.16; 1.6 0.24];
%                 out = [0.015; 0.025; 0.03; 0.04; 0.058; 0.068; 0.06; 0.1;0.131;0.12;0.19;0.23;0.2;0.38;0.5];
%                 w = ( (in' *in) \ in')* out ;
%                 w = [0.049472140762463  1.769794721407624];
    %square plate
    % 1.9586e-07                          10 15 5
    % 1.8848e-07   9.824429047596036e+02  12  8 10
    % 3.273e-07                           20  20 5
    
    %Initialise cluster sizes
    numDotClusters = 12;   
    numPreThreshDistClusters = 8;
    numPostThreshDistClusters = 10;
    
    % Threshold distance
    maxDist = max(nonSingZmnProp(:,1));
    avgEdgeLength = sum(Solver_setup.rwg_basis_functions_length_m)/Solver_setup.num_mom_basis_functions;
    threshDist = 0.035849377364104*maxDist + 1.769794721407624*avgEdgeLength;
    
    % intervals for properties
    distInterval = linspace(0, threshDist,numPreThreshDistClusters );
    numDistClusters = numPostThreshDistClusters + numPreThreshDistClusters -1;
    distInterval(numPreThreshDistClusters:numDistClusters) = linspace(threshDist, maxDist, numPostThreshDistClusters);
    dirDotDirInterval = linspace( -1 , 1 , numDotClusters);
    dirDotDispInterval = linspace( -1 , 1 , numDotClusters);  
    
    % cell clusters
    clusterInd = cell(numDistClusters- 1, numDotClusters -1 ,numDotClusters -1);
    nonSingWeights = cell(numDistClusters- 1, numDotClusters -1 ,numDotClusters -1);
    
    %initialise
    predNonSingZmn = zeros(numNonSing, 1);
    prop = nonSingZmnProp;
    numNoClass = 0;
    includeBias = 1;
    
    for i = 1:(numDistClusters-1)
        low_dist = distInterval(i);
        high_dist = distInterval(i+1);
        for j = 1:(numDotClusters-1)
            low_dir_dot_dir = dirDotDirInterval(j);
            high_dir_dot_dir = dirDotDirInterval(j+1);
            for k = 1:(numDotClusters-1)
                low_dir_dot_disp = dirDotDispInterval(k);
                high_dir_dot_disp = dirDotDispInterval(k+1);
                
                ind = find(prop(:,1) >= low_dist & prop(:,1) <= high_dist & ...
                    prop(:,2) >= low_dir_dot_dir & prop(:,2) <= high_dir_dot_dir &...
                    prop(:,3) >= low_dir_dot_disp & prop(:,3) <= high_dir_dot_disp  );
                
                clusterInd{i,j,k} = ind;              
                clusterTerms =  nonSingZmnTerms(ind, :);
                clusterRefZmn = refNonSingZmn(ind, :);
                [N, ~]= size(clusterTerms);
                
                if (N >= numTerms + includeBias)
                    if (includeBias)
                        clusterTerms(:,numTerms+1) = ones(N,1);
                    end
                    nonSingWeights{i,j,k} = ( (clusterTerms' *clusterTerms) \ clusterTerms')* clusterRefZmn ;
                    predNonSingZmn(ind) = clusterTerms * nonSingWeights{i,j,k};
                    
                else
                    %not enough terms, use unity
                    numNoClass = numNoClass + N;
                    predNonSingZmn(ind) = unityNonSingZmn(ind);
                    w = ones(numTerms,1);
                    if  (includeBias)
                        w(numTerms +1 ,1) = 0;
                    end
                    nonSingWeights{i,j,k} = w;
                end
                
            end % for k
        end % for j
    end % for i
    
    % ------------ ERROR ------------
    
    predDiff = predNonSingZmn -refNonSingZmn;
    predNonSingNormSquareDiff = predDiff.^2 ./ refNonSingZmn.^2;
    predNonSingNormMSE = sum(predNonSingNormSquareDiff) /numNonSing;
    predNonSingMSE = predDiff'*predDiff /numNonSing;
    
    unityNonSingDiff =  unityNonSingZmn - refNonSingZmn;
    unityNonSingNormSquareDiff = unityNonSingDiff.^2 ./ refNonSingZmn.^2;
    unityNonSingNormMSE = sum(unityNonSingNormSquareDiff)./numNonSing;
    unityNonSingMSE = unityNonSingDiff'*unityNonSingDiff ./numNonSing;
    
    % ------------ PLOT DATA ------------
    % for instant debugging
    
    gridSize = 500;
    
    %plotTitle = 'Unity weight error';
    %unityWeightDiff = log(abs(mlmom.nonSingUnityWeightZmn ./ mlmom.refNonSingZmn));
    %plotError(prop, unityNonSingDiff, gridSize, plotTitle);
    %plotError(prop, log(unityNonSingNormSquareDiff), gridSize, plotTitle);
    
    plotTitle = 'Predicted data error';
    %predDiff = log(abs(mlmom.predNonSingZmn ./ mlmom.refNonSingZmn));
    plotError(prop, predDiff, gridSize, plotTitle);
    %plotError(prop, log(predNonSingNormSquareDiff), gridSize, plotTitle);
    
    % ------------ UPDATE MLMOM ------------
    
    mlmom.nonSingZmnTerms = nonSingZmnTerms;  
    mlmom.refNonSingZmn = refNonSingZmn;
    mlmom.unityNonSingZmn = unityNonSingZmn;
    mlmom.singInd = singInd;
    mlmom.clusterInd = clusterInd;
    mlmom.nonSingWeights = nonSingWeights;
    mlmom.predNonSingZmn = predNonSingZmn; 
    mlmom.nonSingZmnProp = nonSingZmnProp;   
    mlmom.dirDotDirInterval = dirDotDirInterval;
    mlmom.dirDotDispInterval =dirDotDispInterval;
    mlmom.distInterval = distInterval;
    
    %mlmom.totsetupTime = 0.0;
    %mlmom.totsolTime = 0.0;
    mlmom.numNonSing = numNonSing;
    mlmom.numTerms = numTerms;
    mlmom.predNonSingMSE = predNonSingMSE;
    mlmom.unityNonSingMSE = unityNonSingMSE;
    mlmom.predNonSingNormMSE = predNonSingNormMSE;
    mlmom.unityNonSingNormMSE = unityNonSingNormMSE;       
    mlmom.threshDist = threshDist;
    mlmom.numNoClass = numNoClass;

end
    

    