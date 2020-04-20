%% Group 11 - April 13th 2020
%
% Solve reaction-diffusion system with finite element method in Matlab.
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
    addpath('../util/')
    [T_cel, n_u, n_v, name] = read_input( varargin{:} ) ;
    
    
    %% Create worspace
    workspace ;
    load workspace.mat ;

    
    %% Load domain
    addpath('../data/meshes/')
    load pear.mat
    
    coordinates = Nodes(:, 2:3) ;
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
    
    tic
    % perform hopotopy continuation
%     for t = 0:dt:1
    t  = 0 ;
    dt = 1 ;
    maxit = 10 ;
    while t < 1
        % update homotopy continuation
        t = t + dt ;
        % store current concentrations if need to restart this step
        C_hold = C ;
            
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
            if norm(P) < 10^(-10)
                fprintf( ' %3.2f    %2d          %6.2e\n', t, n, norm(P) ) ;
                n = 0 ;
                break
            end

            % backtracking
            b = 1 ;
            for k = 1:10
                % store proposed new concentrations
                temp = C - b*P ;

                % recompute variational
                H = assemble_H( coordinates, elements3, temp, R_u, R_v ) ;
                res = K*temp - f + t*H ;

                % check if new concentrations indeed reduce the variational
                if ( norm(res) > norm(G) )
                    b = b/2 ;
                else
                    break
                end
            end
        
            % update estimate for concentrations
            C = C - b*P;
        end
        
        % if solver did not converge
        if n == maxit
            % reset homotopy parameter
            t = t - dt ;
            % reset concentrations to beginning of this step
            C = C_hold ;
            % reduce step size
            dt = min(dt/2, 1-t) ;
            % skip remaining of while-loop
            continue
        end
        
        % reset stepsize to large value
        dt = 2*dt ;
    end
    toc

    %% Compute residuals of nonlinear system
    res = K*C - f + H ;

    
    %% Show oxygen and carbon dioxide solutions
    addpath( '../util' )
    
    % show solutions
    show( C, name, elements3, coordinates, T_cel, n_u, n_v )
        
    % show residuals
    % show( res, 'residuals', elements3, coordinates )
end