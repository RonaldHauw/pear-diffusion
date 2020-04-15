%% Group 11 - April 13th 2020
%
% Test all simulations, boudary conditions and convergence of jacobian.
%
% Call as
%     >> test_all( )

function test_all()
    clc
    clear all
    addpath( '../util' )
    
    %% Run each simulation
    run( 'orchard' )
    run( 'shelf life' )
    run( 'refrigerator' )
    run( 'precooling' )
    run( 'disorder inducing' )
    run( 'optimal ca' )
    
    %% Test finite difference approximation of jacobians
    
    
end