%% Group 11 - March 4th 2020
% compute derivative of reaction functions

clear all
clc

% declare variables
syms K_mv V_mu c_u K_mu c_v r_q K_mfu V_mfv

% declare function
R_u = K_mv .* V_mu .* c_u ./ ( (K_mu+c_u) .* (K_mv+c_v) ) ;
R_v = r_q.* (K_mv .* V_mu .* c_u ./ ( (K_mu+c_u) .* (K_mv+c_v) )) +  ...
      K_mfu .* V_mfv ./ (K_mfu+c_u) ;

% vizualize functions
disp('R_u = '); pretty(R_u); disp(' ')
disp('R_v = '); pretty(R_v); disp(' ')

% differentiate functions
dR_u_u = simplify( diff( R_u, c_u ) ) ;
dR_u_v = simplify( diff( R_u, c_v ) ) ;
dR_v_u = simplify( diff( R_v, c_u ) ) ;
dR_v_v = simplify( diff( R_v, c_v ) ) ;

% visualize derivatives
disp('dR_u_u = '); pretty(dR_u_u); disp(' ')
disp('dR_u_v = '); pretty(dR_u_v); disp(' ')
disp('dR_v_u = '); pretty(dR_v_u); disp(' ')
disp('dR_v_v = '); pretty(dR_v_v); disp(' ')
