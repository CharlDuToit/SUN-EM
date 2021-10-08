function centreDistances = calcCentreDistance(Solver_setup )
    
    numEdges = Solver_setup.num_mom_basis_functions;
    
    centreDistances = zeros(numEdges,numEdges,4);
    for mm = 1:numEdges
        %for nn = 1:numEdges
        for nn = mm:numEdges
            
            mmPlusCentre = Solver_setup.triangle_centre_point(Solver_setup.rwg_basis_functions_trianglePlus(mm),:);
            mmMinusCentre = Solver_setup.triangle_centre_point(Solver_setup.rwg_basis_functions_triangleMinus(mm),:);
            nnPlusCentre =  Solver_setup.triangle_centre_point(Solver_setup.rwg_basis_functions_trianglePlus(nn),:);
            nnMinusCentre = Solver_setup.triangle_centre_point(Solver_setup.rwg_basis_functions_triangleMinus(nn),:);
            
            mmPlus_nnPlus = norm(mmPlusCentre-nnPlusCentre);
            mmPlus_nnMinus = norm(mmPlusCentre-nnMinusCentre);
            mmMinus_nnPlus= norm(mmMinusCentre-nnPlusCentre);
            mmMinus_nnMinus = norm(mmMinusCentre-nnMinusCentre);
            
            %switch 
            centreDistances(mm,nn, :) = [mmPlus_nnPlus, mmPlus_nnMinus, mmMinus_nnPlus, mmMinus_nnMinus];
            centreDistances(nn,mm, :) = [mmPlus_nnPlus, mmMinus_nnPlus, mmPlus_nnMinus, mmMinus_nnMinus];
        end
    end
    
end