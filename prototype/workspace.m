%% Group 11 - March 9th 2020
% variables of diffusion-reaction system

% DECLARE PARAMETERS and functions
% variable
s_ur            = 2.8e-10 ;
s_uz            = 1.1e-9 ;
s_vr            = 2.32e-9 ;
s_vz            = 6.97e-9 ;
%
T_ref           = 293.15 ;
R_g             = 8.314 ;
V_mu_ref        = 2.39e-4 ;
E_a_vmu_ref     = 80200 ;
T               = T_cel + 273.15 ;
V_mu            = V_mu_ref * exp( E_a_vmu_ref/R_g * ( 1/T_ref - 1/T ) ) ;
%
V_mfv_ref       = 1.61e-4 ;
E_a_vmfv_ref    = 56700 ;
V_mfv           = V_mfv_ref * exp( E_a_vmfv_ref/R_g * ( 1/T_ref-1/T ) ) ;
%
K_mu            = 0.4103 ;
K_mv            = 27.2438 ;
K_mfu           = 0.1149 ;
%
r_q             = 0.97 ;
r_u             = 7e-7 ; 
r_v             = 7.5e-7 ;
%
p_atm           = 101300 ;
C_u_amb         = p_atm*n_u/R_g/T ;
C_v_amb         = p_atm*n_v/R_g/T ;
% respiration
R_u     		= @(c_u, c_v) K_mv * V_mu * c_u ./ ( (K_mu+c_u) .* (K_mv+c_v) ) ;
R_v     		= @(c_u, c_v) r_q * R_u(c_u, c_v) +  K_mfu * V_mfv ./ (K_mfu+c_u) ;
% derivative of respiration
dR_u_u  		= @(c_u, c_v) K_mu * K_mv * V_mu ./ ( (K_mu+c_u).^2 .* (K_mv+c_v) ) ;
dR_u_v  		= @(c_u, c_v) - K_mv * V_mu .* c_u ./ ( (K_mu+c_u) .* (K_mv+c_v).^2 ) ;
dR_v_u  		= @(c_u, c_v) r_q .* dR_u_u(c_u, c_v) - K_mfu .* V_mfv ./ (K_mfu+c_u).^2 ;
dR_v_v  		= @(c_u, c_v) r_q .* dR_u_v(c_u, c_v) ;

save('workspace')