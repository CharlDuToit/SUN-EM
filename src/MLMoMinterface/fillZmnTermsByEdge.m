% function [terms, singIndices,...
%     dist,edge_mm_dir_dot_edge_nn_dir, edge_mm_dir_dot_edge_nn_disp ] = fillZmnTermsByEdge(Const,Solver_setup)
%function [terms, singIndices, properties] = fillZmnTermsByEdge(Const,Solver_setup)   
function [terms, singIndices] = fillZmnTermsByEdge(Const,Solver_setup, centreDistances ) 
%function [terms, singIndices] = fillZmnTermsByEdge(Const,Solver_setup ) 
    %fillZmnTermsByEdge
    %   Date: 2021.08.07
    %   Usage:
    %       [Z] = fillZmnTermsByEdge(Const,Solver_setup)
    %
    %   Input Arguments:
    %       Const: 
    %           A global struct containing:
    %       Solver_setup
    %           A global struct containing solver setup details, e.g. frequency range,
    %           geometry data, basis function setup, etc.
    %   Output Arguments:
    %       Z
    %           The Z-matrices data calculated internally
    %
    %   Description:
    %       Fills the impedance matrix at each of the frequency samples.
    %
    %   =======================
    %   Written by Danie Ludick on 2018.05.253
    %   Modified by Charl du Toit on 2021.09.10
    %   Stellenbosch University
    %   Email: dludick.sun.ac.za 
    %   Credit to: Prof. David B. Davidson for the FillZMatrixByEdge routine that is based on his
    %   MATLAB implementation as detailed in [1]
    %   
    %   References: 
    %   [1] David B. Davidson, Computational Electromagnetics for RF and Microwave Engineering, 
    %       Second Edition, (see Chapter 6) - [DBD2011]
    %   [2] Xinlei Chen, Changqing Gu, Zhenyi Niu, and Zhuo L, "Fast Dipole Method for Electromagnetic Scattering From Perfect 
    %       Electric Conducting Targets", IEEE TRANSACTIONS ON ANTENNAS AND PROPAGATION, VOL. 60, NO. 2, FEBRUARY 2012

    message_fc(Const, sprintf('  Calculating input Z-matrix and its geometric properties '));
    
    % Populate now the variables as used by [1] - esp. replace the use of
    % global variables.
    num_dofs = Solver_setup.num_metallic_edges;       % Replacing global NUM_DOFS
    elements = Solver_setup.triangle_vertices;        % Replacing global ELEMENTS
    node_coord = Solver_setup.nodes_xyz;              % Replacing global NODE_COORD
    ell = Solver_setup.rwg_basis_functions_length_m;  % Replacing global ELL
    

    sing     = Const.SING;
    eps_0    = Const.EPS_0;
    mu_0     = Const.MU_0;
    
    % Set some general parameters
    number_of_frequencies = Solver_setup.frequencies.freq_num; % Number of frequencies
    Z.numFreq = number_of_frequencies;
    Z.mBasis  = num_dofs; % number of rows (testing/field functions)
    Z.nBasis  = num_dofs; % number of cols (basis/source  functions)
    
    % == start == Charl modification
    non_sing_quad_pts = Const.QUAD_PTS;
    %sing_quad_pts = 12;
    sing_quad_pts = Const.QUAD_PTS; 
    
    if (sing_quad_pts < 2)
        sing_quad_pts = 3;
    end
    
    %properties = calcProperties(Solver_setup);
    
    terms = complex(zeros(num_dofs,num_dofs, 8, number_of_frequencies));
    %terms(:,:,1) = A_m_pls_n_pls
    %terms(:,:,2) = Phi_m_pls_n_pls
    %terms(:,:,3) = A_m_pls_n_mns
    %terms(:,:,4) = Phi_m_pls_n_mns
    %terms(:,:,5) = A_m_mns_n_pls
    %terms(:,:,6) = Phi_m_mns_n_pls
    %terms(:,:,7) = A_m_mns_n_mns
    %terms(:,:,8) = Phi_m_mns_n_mns
    
    singIndices = zeros(num_dofs ,num_dofs);
    singIndicesCalculated = false;
    
    %  == end == Charl modification
    
    % Extract the triangle midpoints
    r_c = Solver_setup.triangle_centre_point;
    rho_c_pls = Solver_setup.rho_c_pls;
    rho_c_mns = Solver_setup.rho_c_mns;

    % Allocate some space for our impedance matrix
    %Z.values = complex(zeros(num_dofs,num_dofs,number_of_frequencies)); % Generalized impedance matrix. 

    % We will be calculating the Z matrix at each of the various frequencies:
    for freq_index = 1:number_of_frequencies

        % Extract the particular frequency value:
        freq = Solver_setup.frequencies.samples(freq_index);
        % Calculate some frequency dependent parameters required below
        omega = 2*pi*freq;       % Radial frequency
        lambda = Const.C0/freq;  % Wavelength in m
        k  = 2*pi/lambda;        % Wavenumber in rad/m
        %jk = 1i*k;
        
        message_fc(Const, sprintf('    Processing frequency %d of %d (%.2f Hz) ',freq_index,number_of_frequencies,freq))

        % Assemble by edges - not optimally fast computationally, but easy.
        % These are eqns. 32 and 33.
        num_std_zmn_entries_calculated = 0;
        for mm = 1:num_dofs

            %pp_pls = EDGECONXELEMS(mm,1);            
            pp_pls = Solver_setup.rwg_basis_functions_trianglePlus(mm);
            %pp_mns = EDGECONXELEMS(mm,2);
            pp_mns = Solver_setup.rwg_basis_functions_triangleMinus(mm);

           
            for nn = 1:num_dofs
                
                %Only 2 unique triangles
                if (~singIndicesCalculated)
                    if (mm == nn )
                        singIndices(mm, nn) = 1;
                    end
                    %continue
                end
                
                % 2018.06.13: We can use the Equivalent Dipole Method (EDM), if active and only if the distance 
                % between the two RWG dipole moments is sufficent (otherwise we still use normal matrix fill below)
                Zmn_calculated = false;
                    
                % Perform the standard Matrix calculation (using Gausian quadrature integration, if the EDM is not used above
                if (~Zmn_calculated)
                    num_std_zmn_entries_calculated = num_std_zmn_entries_calculated + 1;

                    %fprintf('Processing element %d,%d\n',mm,nn);
                    % There are four terms to add here;
                    % From integrals over source faces n+ & n- evaluated at centres of field faces m+ and m-.
                    % The n integrals associate with q, the m, with p.
                    % In datastructures EDGECONXELEMS & DOFLOCALNUM, second index 1
                    % associates with +, 2 with -. 

                    %qq_pls = EDGECONXELEMS(nn,1);
                    qq_pls = Solver_setup.rwg_basis_functions_trianglePlus(nn);
                    %qq_mns = EDGECONXELEMS(nn,2);
                    qq_mns = Solver_setup.rwg_basis_functions_triangleMinus(nn);
                    
                    %Only 3 unique triangles
                    if (~singIndicesCalculated)                      
                        if (pp_pls == qq_pls || pp_pls == qq_mns || pp_mns == qq_pls || pp_mns == qq_mns )
                            singIndices(mm, nn) = 1;
                        end
                        %continue
                    end
                
                    triangle_tn_plus_free_vertex = Solver_setup.rwg_basis_functions_trianglePlusFreeVertex(nn);
                    triangle_tn_minus_free_vertex = Solver_setup.rwg_basis_functions_triangleMinusFreeVertex(nn);
                    
                    % First, find contribution from n+ and n- faces evaluated at m+
                    
                    % Change quad_pts based on type of triangle pair
                    quad_pts = non_sing_quad_pts;
                    if (singIndices(mm,nn))
                        quad_pts = sing_quad_pts;
                    end
                                    
                    % --------------------------------------------------------------                
                    % Look at Tm+
                    % --------------------------------------------------------------
                    % -- Contribution from Tn+
                    [MagVecPot,ScalPot] = Potentials(elements,node_coord,ell,pp_pls,qq_pls,mm,nn,triangle_tn_plus_free_vertex,...
                        k,r_c,quad_pts,sing,eps_0,mu_0,omega,centreDistances(pp_pls,qq_pls));
                    %[MagVecPot,ScalPot] = Potentials(elements,node_coord,ell,pp_pls,qq_pls,mm,nn,triangle_tn_plus_free_vertex,...
                       % k,r_c,quad_pts,sing,eps_0,mu_0,omega);
                    
                    % [DBD2011] implementation below. I think there is a sign issue here.
                    %Amn_pls_source_pls = MagVecPot;
                    %Phi_mn_pls_source_pls = -ScalPot;
                    % [DL - 2018] implementation: Swop the signs around
                    %Amn_pls_source_pls = -MagVecPot;
                    %Phi_mn_pls_source_pls = +ScalPot;
                    
                    %terms(:,:,1) = A_m_pls_n_pls
                    %terms(:,:,2) = Phi_m_pls_n_pls
                    terms(mm,nn,1,freq_index) = ell(mm) * 0.5i*omega*dot(-MagVecPot',rho_c_pls(mm,:));
                    terms(mm,nn,2,freq_index) = - ell(mm) * ScalPot;

                    
                    % -- Contribution from Tn-
                    
                    
                    %[MagVecPot,ScalPot] = Potentials(elements,node_coord,ell,pp_pls,qq_mns,mm,nn,triangle_tn_minus_free_vertex,...
                        %k,r_c,quad_pts,sing,eps_0,mu_0,omega);
                    [MagVecPot,ScalPot] = Potentials(elements,node_coord,ell,pp_pls,qq_mns,mm,nn,triangle_tn_minus_free_vertex,...
                        k,r_c,quad_pts,sing,eps_0,mu_0,omega,centreDistances(pp_pls,qq_mns));
                    % [DBD2011] implementation below. I think there is a sign issue here.
                    %Amn_pls_source_mns = - MagVecPot;
                    %Phi_mn_pls_source_mns = +ScalPot;
                    % [DL - 2018] implementation: Swop the signs around
                    %Amn_pls_source_mns = + MagVecPot;
                    %Phi_mn_pls_source_mns = -ScalPot;
                    
                    %terms(:,:,3) = A_m_pls_n_mns
                    %terms(:,:,4) = Phi_m_pls_n_mns
                    terms(mm,nn,3,freq_index) = ell(mm) * 0.5i*omega*dot(MagVecPot',rho_c_pls(mm,:));
                    terms(mm,nn,4,freq_index) =  ell(mm) * ScalPot;
                                    
                    %Amn_pls = Amn_pls_source_pls + Amn_pls_source_mns;
                    %Phi_mn_pls = Phi_mn_pls_source_pls + Phi_mn_pls_source_mns;
                    
                    % --------------------------------------------------------------                
                    % Look at Tm-
                    % --------------------------------------------------------------
                    % -- Contribution from Tn+
                    [MagVecPot,ScalPot] = Potentials(elements,node_coord,ell,pp_mns,qq_pls,mm,nn,triangle_tn_plus_free_vertex,...
                        k,r_c,quad_pts,sing,eps_0,mu_0,omega,centreDistances(pp_mns,qq_pls));
                    %[MagVecPot,ScalPot] = Potentials(elements,node_coord,ell,pp_mns,qq_pls,mm,nn,triangle_tn_plus_free_vertex,...
                       %k,r_c,quad_pts,sing,eps_0,mu_0,omega);
                    % [DBD2011] implementation below. I think there is a sign issue here.
                    %Amn_mns_source_pls = MagVecPot;
                    %Phi_mn_mns_source_pls = -ScalPot;
                    % [DL - 2018] implementation: Swop the signs around
                    %Amn_mns_source_pls = -MagVecPot;
                    %Phi_mn_mns_source_pls = +ScalPot;
                    
                    %terms(:,:,5) = A_m_mns_n_pls
                    %terms(:,:,6) = Phi_m_mns_n_pls
                    terms(mm,nn,5,freq_index) = ell(mm) * 0.5i*omega*dot(-MagVecPot',rho_c_mns(mm,:));
                    terms(mm,nn,6,freq_index) =  ell(mm) * ScalPot;
                    
                    % -- Contribution from Tn-
                    %[MagVecPot,ScalPot] = Potentials(elements,node_coord,ell,pp_mns,qq_mns,mm,nn,triangle_tn_minus_free_vertex,...
                        %k,r_c,quad_pts,sing,eps_0,mu_0,omega);
                    [MagVecPot,ScalPot] = Potentials(elements,node_coord,ell,pp_mns,qq_mns,mm,nn,triangle_tn_minus_free_vertex,...
                        k,r_c,quad_pts,sing,eps_0,mu_0,omega,centreDistances(pp_mns,qq_mns));
                    % [DBD2011] implementation below. I think there is a sign issue here.
                    %Amn_mns_source_mns = - MagVecPot;
                    %Phi_mn_mns_source_mns = +ScalPot;
                    % [DL - 2018] implementation: Swop the signs around
                    %Amn_mns_source_mns = +MagVecPot;
                    %Phi_mn_mns_source_mns = -ScalPot;
                    
                    %terms(:,:,7) = A_m_mns_n_mns
                    %terms(:,:,8) = Phi_m_mns_n_mns
                    terms(mm,nn,7,freq_index) = ell(mm) * 0.5i*omega*dot(MagVecPot',rho_c_mns(mm,:));
                    terms(mm,nn,8,freq_index) =  - ell(mm) * ScalPot;
                    
                    %Amn_mns = Amn_mns_source_pls + Amn_mns_source_mns;
                    %Phi_mn_mns = Phi_mn_mns_source_pls + Phi_mn_mns_source_mns;
                           
                    % Assemble with eq. 17 in [RWG82]
                    %Phi_mn_pls = Phi_mn_pls_source_pls + Phi_mn_pls_source_mns;
                    %Phi_mn_mns = Phi_mn_mns_source_pls + Phi_mn_mns_source_mns;
                    %Amn_pls = Amn_pls_source_pls + Amn_pls_source_mns;
                    %Amn_mns = Amn_mns_source_pls + Amn_mns_source_mns;
                    

                    %Z.values(mm,nn) = 1i*omega*...
                        %(dot(Amn_pls',rho_c_pls(mm,:))/2 + dot(Amn_mns',rho_c_mns(mm,:))/2) + Phi_mn_mns - Phi_mn_pls;
                    
                    %mm
                    %nn
                    %Z.values(mm,nn) = ell(mm)* Z.values(mm,nn);                    
                end % if (~Zmn_calculated)
                
            end % for nn = 1:NUM_DOFS
        end %for mm = 1:NUM_DOFS
        % Singular indices is not dependant on frequency
        singIndicesCalculated = true;

    end %for freq_index = 1:Solver_setup.frequencies.freq_num 
   
end %function FillZMatrixByEdge

% =================================================================================
% Local function for evaluating the potential integrals, as done in [1]
% =================================================================================
function [MagVecPot,ScalPot] = Potentials(elements,node_coord,ell, field_pt,source_pt,field_edge, source_edge, ii, k,r_c,quad_pts,sing,eps_0,mu_0,omega,dist)
%function [MagVecPot,ScalPot] = Potentials(elements,node_coord,ell, field_pt,source_pt,field_edge, source_edge, ii, k,r_c,quad_pts,sing,eps_0,mu_0,omega)
    % This subfunction computes the magnetic vector and scalar potentials for 
    % field point face field_pt and source point face source_pt. 

    % Code for debugging singularity scheme below:
    %[Ipq,Ipq_xi,Ipq_eta,Ipq_zeta] = Int_pq_debug(field_pt,source_pt,r_c(field_pt,:),k,quad_pts,sing);
    [Ipq,Ipq_xi,Ipq_eta,Ipq_zeta] = Int_pq(elements,node_coord, field_pt,source_pt,r_c(field_pt,:),k,quad_pts,sing,dist); %r_c trinagle center
    %[Ipq,Ipq_xi,Ipq_eta,Ipq_zeta] = Int_pq(elements,node_coord, field_pt,source_pt,r_c(field_pt,:),k,quad_pts,sing);
    
    % Extract the nodes of the source triangle (q)
    qnodes = elements(source_pt,:);    
    r = zeros(3,3); % Store all three vertices of triangle q here in r()
    r(1,:) = node_coord(qnodes(1),:);
    r(2,:) = node_coord(qnodes(2),:);
    r(3,:) = node_coord(qnodes(3),:);
    
    % Extract the position of the ith free vertex (i.e. associated with the
    % ith RWG)
    %ii_nodes = elements(ii,:);
    rii = zeros(1,3);
    rii(1,1) = node_coord(ii,1);
    rii(1,2) = node_coord(ii,2);
    rii(1,3) = node_coord(ii,3);
        
    
    %ii = DOFLOCALNUM(source_edge,source_tri); % This is the free vertex
    %associated with the source_edge - which is now passed here as an
    %argument by the calling routine.    

    % [RWG82, Eq. (32) - without sign
    MagVecPot = mu_0*ell(source_edge)/(4*pi)*...
        ( r(1,:)*Ipq_xi + r(2,:)*Ipq_eta + r(3,:)*Ipq_zeta - rii(1,:)*Ipq);

    % [RWG82, Eq. (33) - without sign
    ScalPot = ell(source_edge)/(1i*2*pi*omega*eps_0) * Ipq;
        
end % function Potentials

%function [Ipq,Ipq_xi,Ipq_eta,Ipq_zeta] = Int_pq(elements,node_coord, p,q,r_cp,k,quad_pts,sing)
function [Ipq,Ipq_xi,Ipq_eta,Ipq_zeta] = Int_pq(elements,node_coord, p,q,r_cp,k,quad_pts,sing,dist) %r_cp is triangle p center
    % INT_PQ returns the four integrals for faces (=elements) p and q,
    % corresponding to [eqs.34a-34d,RWG82]. Note that these integrals are
    % evaluated in normalized coordinates, and are not scaled by the area.
    % The overall matrix element is assembled from these elsewhere.
    % Note also that p is the field point, q the source point. The integrals are
    % performed over triangle q.
    % Note further that xi, eta and zeta are equivalent to
    % lambda_1, lambda_2, and lambda_3 respectively.
    % If flag sing is set, then the singular terms are evaluated using a
    % special integration rule.

    % Author: D B Davidson, Dec 2009.
    % Corrections for singular intergral evaluation: 1 June 2010 DBD.

    %global ELEMENTS NODE_COORD

    qnodes = elements(q,:);
    n1 = node_coord(qnodes(1),:);
    n2 = node_coord(qnodes(2),:);
    n3 = node_coord(qnodes(3),:);
    area = tri_area3D(n1,n2,n3);

    Ipq=0;
    Ipq_xi=0;
    Ipq_eta=0;
    Ipq_zeta=0;

    if p==q && sing
        [Ipq,Ipq_xi,Ipq_eta,Ipq_zeta] = intg_sing_SGF(k,r_cp,n1,n2,n3,3,4);
        Ipq = Ipq/(2*area); 
        % The factor of 2 here is required to get good results, but it is not
        % clear why this is correct. If removed, then the w/2 below should also be
        % changed to w, but the results are then out by 2. 
        Ipq_xi = Ipq_xi/(2*area);
        Ipq_eta = Ipq_eta/(2*area);
        Ipq_zeta = Ipq_zeta/(2*area);
        % 2A factor above required to give result in simplex coordinates - see
        % [eq 31. RWG82]. (Function intg_sing_SGF returns the LHS thereof).

    else
        [w,lambda] = tri_quad(quad_pts);
        w=w/2; % Rule must be correctly normalized.

        %r_cp
        for nn=1:quad_pts
            if (quad_pts == 1)
                R_p = dist;
                %R_p = dist(p,q);
            else
                r_prime = lambda(nn,1)*n1 + lambda(nn,2)*n2 + lambda(nn,3)*n3; % [eq.30,RWG82]
                R_p = norm(r_cp-r_prime); % [eq.27,RWG82]
            end          

            % DJdbg --> plot r_cp and r_prime
            %plot(r_prime(1,1),r_prime(1,2),'or','MarkerFaceColor','g','MarkerSize',10);        
            %plot(r_cp(1,1),r_cp(1,2),'or','MarkerFaceColor','k','MarkerSize',10);   

            GF  = exp(-1i*k*R_p)/R_p; % Green's function
            %w(nn)
            Ipq     = Ipq+w(nn)*GF;
            %lambda(nn,1)
            Ipq_xi  = Ipq_xi+w(nn)*lambda(nn,1)*GF;
            %lambda(nn,2)
            Ipq_eta = Ipq_eta+w(nn)*lambda(nn,2)*GF;
            %Ipq_zeta = Ipq_zeta+w(nn)*lambda(nn,3)*GF; % Code gives same
            %answers as below.
        end
        Ipq_zeta = Ipq - Ipq_xi - Ipq_eta;
    end

end





