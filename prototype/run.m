%% Group 11 - April 13th 2020
% Finite element method for two-dimensional reaction-diffusion system.
%
% Either run simulation by name as
%     >> run( name )
% where 'name' is the string 'orchard', 'shelf life', 'refrigerator', 
% 'precooling', 'disorder inducing' or 'optimal CA',
%
% or call with specific parameters as
%    >> run( T_cel, n_u, n_v )
% where T_cel is the temperature in degrees Celsius, and n_u and n_v are the 
% oxygen and carbon dioxide concentrations as a fraction between 0 and 1,
%
% or run default 'refrigerator' simulation ( T_cel = 7, n_u = 0.208, n_v = 0 )
%    >> run
%
% Inspired on https://www.math.hu-berlin.de/~cc/cc_homepage/download/1999-AJ_CC_FS-50_Lines_of_Matlab.pdf


function run( varargin )
    clc

    %% Read input
    if ( nargin < 2 )
        % run default simulation 'precooling'
        if ( nargin == 0 )
            name = 'refrigerator' ;
        
        % run simulation specified by name
        elseif ( nargin == 1 )
            name = varargin{1} ;
        end
        
        switch lower(name)
            case 'orchard'
                T_cel = 25 ;
                n_u   = 0.208 ;
                n_v   = 0.0004 ;
            
            case 'shelf life'
                T_cel = 20 ;
                n_u   = 0.208 ;
                n_v   = 0 ;
                
            case 'refrigerator'
                T_cel = 7 ;
                n_u   = 0.208 ;
                n_v   = 0 ;
                
            case 'precooling'
                T_cel = -1 ;
                n_u   = 0.208 ;
                n_v   = 0 ;
                
            case 'disorder inducing'
                T_cel = -1 ;
                n_u   = 0.02 ;
                n_v   = 0.05 ;
                
            case 'optimal ca'
                T_cel = -1 ;
                n_u   = 0.02 ;
                n_v   = 0.007 ;
                
            otherwise
                error( "Did not understand which simulation to run. Run 'help run' for more information." )
        end
        
    % run simulation specified by parameters    
    elseif ( nargin == 3 )
        T_cel = varargin{1} ;
        n_u   = varargin{2} ;
        n_v   = varargin{3} ;
        
        % store parameters of simulation
        name = join(['conditions : ', num2str(100*n_u), '% O_2, ', num2str(100*n_v), '% CO_2, and ', num2str(T_cel), '°C' ]) ;
    else
        error( "Did not understand which simulation to run. Run 'help run' for more information." )
    end
    %
    disp( "Run simulation with T_cel = " + num2str(T_cel) + ", n_u = " + num2str(n_u)+ ", n_v = " + num2str(n_v) )
    disp( "	" )
    
    
    %% Create worspace
    workspace ;
    load workspace.mat ;

    
    %% Load domain
    load mesh/HCT_Mesh_Data.mat
    %
    coordinates = Nodes(:, 2:3)/50 ;
    % coordinates = mesh.Nodes'/50 ;
    elements3   = Elements( : , 2:end ) ;
    % number of vertices
    M           = size(coordinates, 1) ;
    % coordinates of mesh vertices
    r           = coordinates(:, 1) ;
    z           = coordinates(:, 2) ;
    % edge information
    G1_edges    = InnerBEdges( :, 2:end ) ;
    G2_edges    = OuterBEdges( :, 2:end ) ;

    
    %% Solve FEM model
    % parameters of homotopy continuation
    dt = 0.05 ;
    maxit = 50 ;

    % intialize concentrations
    C = zeros(2*M, 1) ;
    % K = [ K_u , 0 ; 0 , K_v ]
    K = assemble_K( coordinates, elements3, G2_edges, s_ur, s_vr, s_uz, s_vz, r_u, r_v ) ;
    % f = [ f_u ; f_v ]
    f = assemble_f( coordinates, G2_edges, r_u, r_v, C_u_amb, C_v_amb ) ;

    % keep track of progress
    fprintf( '        %12s \n', 'converged at' ) ;
    fprintf( '%3s     %10s      %8s \n','t',  'iteration', 'residual' ) ;
    fprintf( ' ----   ------------    -------- \n') ;
          
    % perform hopotopy continuation
    for t = 0:dt:1

        fprintf( ' %3.2f',t ) ;
        
        % prediction step
        % nonlinearity H = [ H_u(C) ; H_v(C) ]
        H       = assemble_H( coordinates, elements3, C, R_u, R_v ) ;
        % Jacobian J = dH/dC
        J       = assemble_J( coordinates, elements3, C, dR_u_u, dR_u_v, dR_v_u, dR_v_v ) ;
        % update concentrations with forward euler
        dC_dt   = ( K + t*J ) \ ( -H ) ;
        C       = C + dt * dC_dt ;

        % correction step with Newton method
        for n=1:maxit

            % no need to recompute K and f

            % nonlinearity H = [ H_u(C) ; H_v(C) ]
            H = assemble_H( coordinates, elements3, C, R_u, R_v ) ;
            % Jacobian J = dH/dC
            J = assemble_J( coordinates, elements3, C, dR_u_u, dR_u_v, dR_v_u, dR_v_v ) ;

            % Variational
            G = K*C - f + t*H ;

            % solving one Newton step (J_G)^-1 * G
            P = ( K + t*J ) \ G ;

            % check for convergence
            if norm(P) < 10^(-12)
                fprintf( '       %2d          %6.2e\n', n, norm(P) ) ;
                break
            end

            % backtracking
            b = 1 ;
            for k = 1:50
                % store proposed new concentrations
                temp = C - b*P ;

                % recompute Variational
                H = assemble_H( coordinates, elements3, temp, R_u, R_v ) ;
                res = K*temp - f + t*H ;

                % check if new concentrations indeed reduce the Variational
                if ( norm(res) > norm(G) )
                    b = b/2 ;
                else
                    break
                end
            end

            % update estimate for concentrations
            C = C - b*P;
            
            % check for convergence
            if n == maxit
                fprintf( '    %8s       %6.2e\n', 'maximum', norm(P) ) ;
                break
            end
        end
    end


    %% Show oxygen and carbon dioxide solutions
    addpath( '../matlab' )
    %
    figure('position', [300 100 450 400])
    subplot(1, 2, 1)
    show(elements3,[],coordinates,full( C(1:M) ), [0, 10], C_u_amb);
    title('Oxygen [mol/m³]', 'FontSize', 10)
    %
    subplot(1, 2, 2)
    show(elements3,[],coordinates,full( C(M+1:end) ), [0, 5], C_v_amb);
    title('Carbon dioxide [mol/m³]', 'FontSize', 10)
    %
    sgtitle( join(['Simulated ', name]), 'FontSize', 12 )
    
end