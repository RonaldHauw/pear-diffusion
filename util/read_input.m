%% Group 11 - April 13th 2020
%
% Read inputs for running prototype and software.

function [T_cel, n_u, n_v, name] = read_input( varargin )

    if ( nargin < 2 )
        % run default simulation 'refrigerator'
        if ( nargin == 0 )
            name = 'refrigerator' ;
        % run simulation specified by name
        elseif ( nargin == 1 )
            name = varargin{1} ;
        end
        
        switch lower(name)
            case {'orchard', 'or'}
                % parameters of simulation
                T_cel = 25 ;
                n_u   = 0.208 ;
                n_v   = 0.0004 ;
                name  = "Orchard conditions" ; 
            
            case {'shelf life', 'sl', 'shelflife'}
                % parameters of simulation
                T_cel = 20 ;
                n_u   = 0.208 ;
                n_v   = 0 ;
                name  = "Shelf life conditions" ;
                
            case {'refrigerator', 'r'}
                % parameters of simulation
                T_cel = 7 ;
                n_u   = 0.208 ;
                n_v   = 0 ;
                name  = "Refrigerator conditions" ;
                
            case {'precooling', 'p'}
                % parameters of simulation
                T_cel = -1 ;
                n_u   = 0.208 ;
                n_v   = 0 ;
                name  = "Precooling conditions" ;
                
            case {'disorder inducing', 'di', 'disorderinducing'}
                % parameters of simulation
                T_cel = -1 ;
                n_u   = 0.02 ;
                n_v   = 0.05 ;
                name  = "Disorder inducing conditions" ;
                
            case {'optimal ca', 'oca', 'optimalca'}
                % parameters of simulation
                T_cel = -1 ;
                n_u   = 0.02 ;
                n_v   = 0.007 ;
                name  = "Optimal CA conditions" ;
                
            otherwise
                error( "Did not understand which simulation to run. Run 'help run' for more information." )
        end
        
    % run simulation specified by parameters    
    elseif ( nargin == 3 )
        T_cel = varargin{1} ;
        n_u   = varargin{2} ;
        n_v   = varargin{3} ;
        
        % store parameters of simulation
        name = join(['Conditions : ', num2str(100*n_u), '% O_2, ', num2str(100*n_v), '% CO_2, and ', num2str(T_cel), '°C' ]) ;
    else
        error( "Did not understand which simulation to run. Run 'help run' for more information." )
    end
    %
    disp( "Simulating T_cel = " + num2str(T_cel) + ", n_u = " + num2str(n_u)+ ", n_v = " + num2str(n_v) )
    disp( " " )
end