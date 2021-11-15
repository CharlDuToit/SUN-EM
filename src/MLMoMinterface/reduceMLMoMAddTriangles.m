function reducedMLMoMAddTriangles = reduceMLMoMAddTriangles(mlmomAddTriangles)
    %only keeps mlmom member needed for prediction
    
    reducedMLMoMAddTriangles = [];
    
    reducedMLMoMAddTriangles.weightModels = mlmomAddTriangles.weightModels;
    reducedMLMoMAddTriangles.groupMeans = mlmomAddTriangles.groupMeans;
    reducedMLMoMAddTriangles.errorSummary = mlmomAddTriangles.errorSummary;
    reducedMLMoMAddTriangles.predRelNormPercentError = mlmomAddTriangles.predRelNormPercentError;
    reducedMLMoMAddTriangles.numFreq = mlmomAddTriangles.numFreq;
    reducedMLMoMAddTriangles.quadPts = mlmomAddTriangles.quadPts;
    reducedMLMoMAddTriangles.singDataThresh = mlmomAddTriangles.singDataThresh;
    reducedMLMoMAddTriangles.oldZmnTime = mlmomAddTriangles.oldZmnTime;
    reducedMLMoMAddTriangles.projectTime = mlmomAddTriangles.projectTime;
    reducedMLMoMAddTriangles.refZmnTime = mlmomAddTriangles.refZmnTime;
    reducedMLMoMAddTriangles.assignTime = mlmomAddTriangles.assignTime;
    reducedMLMoMAddTriangles.multiplyTime = mlmomAddTriangles.multiplyTime;
    reducedMLMoMAddTriangles.predictTime = mlmomAddTriangles.predictTime;
    reducedMLMoMAddTriangles.trainingTime = mlmomAddTriangles.trainingTime;
    
end