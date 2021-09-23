activePlots =   [0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0];
%activePlots =   [0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0];
%plotMLMOM(mlmom, activePlots, 5, 1, 5);
 [numClusters, ~]=size(mlmom.clusterMeans);
 clusterCounts = mlmom.clusterCounts;
 [~ , ind] = sort(clusterCounts);
 for k = 1:numClusters
     plotMLMOM(mlmom, activePlots, 500, 1, ind(numClusters - k + 1));
 end
 

