\ENCODING=ISO-8859-1
\Problem name: 

Minimize
 obj1: p_thermal_1_0 + p_thermal_2_0 + 0 p_flow_network_3_0
       + 0 p_flow_network_3_1 + 0 p_flow_network_3_2 + 0 p_flow_network_3_3
       + 0 W_voltage_real_3_3 + 0 W_voltage_real_3_4 + 0 W_voltage_real_3_6
       + 0 W_voltage_real_3_8 + 0 W_voltage_real_3_9 + 0 W_voltage_real_3_12
       + 0 W_voltage_real_3_13 + 0 W_voltage_real_3_14 + 0 W_voltage_imag_3_3
       + 0 W_voltage_imag_3_4 + 0 W_voltage_imag_3_6 + 0 W_voltage_imag_3_8
       + 0 W_voltage_imag_3_9 + 0 W_voltage_imag_3_12 + 0 W_voltage_imag_3_13
       + 0 W_voltage_imag_3_14 + 0 W_voltage_imag_3_15
Subject To
 node_injection_c_0_0:           p_thermal_2_0 - s_network_3_0  = 0
 node_injection_c_0_1:           - s_network_3_1  = 0
 node_injection_c_0_2:           - s_network_3_2  = 0
 node_injection_c_0_3:           p_thermal_1_0 - s_network_3_3  = 0
 StartUp_ShutDown_Variables_Const_Thermal_1_0: u_thermal_1_0 - v_thermal_1_0
                                 + w_thermal_1_0  = 1
 StartUp_Commitment_Const_Thermal_1_0: u_thermal_1_0 - v_thermal_1_0 >= 0
 ShutDown_Commitment_Const_Thermal_1_0: u_thermal_1_0 + w_thermal_1_0 <= 1
 MinPower_Const_Thermal_1_0:     p_thermal_1_0 >= 0
 MaxPower_Const_Thermal_1_0:     318 u_thermal_1_0 - p_thermal_1_0 >= 0
 MaxPower_Const_Thermal_1_1:     318 u_thermal_1_0 - p_thermal_1_0
                                 - 318 v_thermal_1_0 >= 0
 MaxPower_Const_Thermal_1_2:     318 u_thermal_1_0 - p_thermal_1_0 >= 0
 StartUp_ShutDown_Variables_Const_Thermal_2_0: u_thermal_2_0 - v_thermal_2_0
                                 + w_thermal_2_0  = 1
 StartUp_Commitment_Const_Thermal_2_0: u_thermal_2_0 - v_thermal_2_0 >= 0
 ShutDown_Commitment_Const_Thermal_2_0: u_thermal_2_0 + w_thermal_2_0 <= 1
 MinPower_Const_Thermal_2_0:     p_thermal_2_0 >= 0
 MaxPower_Const_Thermal_2_0:     - p_thermal_2_0 >= 0
 MaxPower_Const_Thermal_2_1:     - p_thermal_2_0 >= 0
 MaxPower_Const_Thermal_2_2:     - p_thermal_2_0 >= 0
 AC_angle_bounds_limit_3_0:      - 1.56655753389007e-14 W_voltage_real_3_1
                                 + W_voltage_imag_3_1 >= 0
 AC_angle_bounds_limit_3_1:      - 1.56655753389007e-14 W_voltage_real_3_2
                                 + W_voltage_imag_3_2 >= 0
 AC_angle_bounds_limit_3_2:      - 1.56655753389007e-14 W_voltage_real_3_7
                                 + W_voltage_imag_3_7 >= 0
 AC_angle_bounds_limit_3_3:      - 1.56655753389007e-14 W_voltage_real_3_11
                                 + W_voltage_imag_3_11 >= 0
 AC_angle_bounds_limit_3_4:      1.56655753389007e-14 W_voltage_real_3_1
                                 + W_voltage_imag_3_1 <= 0
 AC_angle_bounds_limit_3_5:      1.56655753389007e-14 W_voltage_real_3_2
                                 + W_voltage_imag_3_2 <= 0
 AC_angle_bounds_limit_3_6:      1.56655753389007e-14 W_voltage_real_3_7
                                 + W_voltage_imag_3_7 <= 0
 AC_angle_bounds_limit_3_7:      1.56655753389007e-14 W_voltage_real_3_11
                                 + W_voltage_imag_3_11 <= 0
 AC_voltage_defintion_const_3_0: 3.81562881562882 W_voltage_real_3_0
                                 - 3.81562881562882 W_voltage_real_3_1
                                 - 19.0268940781441 W_voltage_imag_3_0
                                 + 19.0781440781441 W_voltage_imag_3_1
                                 - S_power_flow_3_0 - S_power_flow_3_8  = 0
 AC_voltage_defintion_const_3_1: 5.16956162117453 W_voltage_real_3_0
                                 - 5.16956162117453 W_voltage_real_3_2
                                 - 25.8090581058726 W_voltage_imag_3_0
                                 + 25.8478081058726 W_voltage_imag_3_2
                                 - S_power_flow_3_1 - S_power_flow_3_9  = 0
 AC_voltage_defintion_const_3_2: 5.16956162117453 W_voltage_real_3_5
                                 - 5.16956162117453 W_voltage_real_3_7
                                 - 25.8090581058726 W_voltage_imag_3_5
                                 + 25.8478081058726 W_voltage_imag_3_7
                                 - S_power_flow_3_2 - S_power_flow_3_10  = 0
 AC_voltage_defintion_const_3_3: 3.02370585389453 W_voltage_real_3_10
                                 - 3.02370585389453 W_voltage_real_3_11
                                 - 15.0547792694727 W_voltage_imag_3_10
                                 + 15.1185292694727 W_voltage_imag_3_11
                                 - S_power_flow_3_3 - S_power_flow_3_11  = 0
 AC_voltage_defintion_const_3_4: - S_power_flow_3_4 - S_power_flow_3_12  = 0
 AC_voltage_defintion_const_3_5: - S_power_flow_3_5 - S_power_flow_3_13  = 0
 AC_voltage_defintion_const_3_6: - S_power_flow_3_6 - S_power_flow_3_14  = 0
 AC_voltage_defintion_const_3_7: - S_power_flow_3_7 - S_power_flow_3_15  = 0
 AC_voltage_defintion_const_3_8: 19.0268940781441 W_voltage_real_3_0
                                 - 19.0781440781441 W_voltage_real_3_1
                                 + 3.81562881562882 W_voltage_imag_3_0
                                 - 3.81562881562882 W_voltage_imag_3_1  = 0
 AC_voltage_defintion_const_3_9: 25.8090581058726 W_voltage_real_3_0
                                 - 25.8478081058726 W_voltage_real_3_2
                                 + 5.16956162117453 W_voltage_imag_3_0
                                 - 5.16956162117453 W_voltage_imag_3_2  = 0
 AC_voltage_defintion_const_3_10: 25.8090581058726 W_voltage_real_3_5
                                 - 25.8478081058726 W_voltage_real_3_7
                                 + 5.16956162117453 W_voltage_imag_3_5
                                 - 5.16956162117453 W_voltage_imag_3_7  = 0
 AC_voltage_defintion_const_3_11: 15.0547792694727 W_voltage_real_3_10
                                 - 15.1185292694727 W_voltage_real_3_11
                                 + 3.02370585389453 W_voltage_imag_3_10
                                 - 3.02370585389453 W_voltage_imag_3_11  = 0
 AC_voltage_defintion_const_3_12:  = 0
 AC_voltage_defintion_const_3_13:  = 0
 AC_voltage_defintion_const_3_14:  = 0
 AC_voltage_defintion_const_3_15:  = 0
 AC_thermal_limit_const_3_0:     S_power_flow_3_0 + S_power_flow_3_8 <= 62500
 AC_thermal_limit_const_3_1:     S_power_flow_3_1 + S_power_flow_3_9 <= 62500
 AC_thermal_limit_const_3_2:     S_power_flow_3_2 + S_power_flow_3_10 <= 62500
 AC_thermal_limit_const_3_3:     S_power_flow_3_3 + S_power_flow_3_11 <= 62500
 AC_thermal_limit_const_3_4:     S_power_flow_3_4 + S_power_flow_3_12 <= 62500
 AC_thermal_limit_const_3_5:     S_power_flow_3_5 + S_power_flow_3_13 <= 62500
 AC_thermal_limit_const_3_6:     S_power_flow_3_6 + S_power_flow_3_14 <= 62500
 AC_thermal_limit_const_3_7:     S_power_flow_3_7 + S_power_flow_3_15 <= 62500
 AC_socp_const_3_0:              0.5 W_voltage_real_3_0
                                 + 0.25 W_voltage_real_3_5 + quad_aux_var_0
                                 <= 0
 AC_socp_const_3_1:              0.5 W_voltage_real_3_0
                                 + 0.25 W_voltage_real_3_10 + quad_aux_var_1
                                 <= 0
 AC_socp_const_3_2:              0.5 W_voltage_real_3_5
                                 + 0.25 W_voltage_real_3_15 + quad_aux_var_2
                                 <= 0
 AC_socp_const_3_3:              0.5 W_voltage_real_3_10
                                 + 0.25 W_voltage_real_3_15 + quad_aux_var_3
                                 <= 0
 quad_aux_con_0: quad_aux_var_0 + [ - 2 W_voltage_real_3_0 ^2
                 - W_voltage_real_3_1 ^2 - 4 W_voltage_real_3_5 ^2 ] >= 0
 quad_aux_con_1: quad_aux_var_1 + [ - 2 W_voltage_real_3_0 ^2
                 - W_voltage_real_3_2 ^2 - 4 W_voltage_real_3_10 ^2 ] >= 0
 quad_aux_con_2: quad_aux_var_2 + [ - 2 W_voltage_real_3_5 ^2
                 - W_voltage_real_3_7 ^2 - 4 W_voltage_real_3_15 ^2 ] >= 0
 quad_aux_con_3: quad_aux_var_3 + [ - 2 W_voltage_real_3_10 ^2
                 - W_voltage_real_3_11 ^2 - 4 W_voltage_real_3_15 ^2 ] >= 0
Bounds
 0 <= u_thermal_1_0 <= 1
      v_thermal_1_0 = 0
 0 <= w_thermal_1_0 <= 1
 0 <= u_thermal_2_0 <= 1
      v_thermal_2_0 = 0
 0 <= w_thermal_2_0 <= 1
      s_network_3_0 Free
      s_network_3_1 Free
      s_network_3_2 Free
      s_network_3_3 Free
      p_flow_network_3_0 Free
      p_flow_network_3_1 Free
      p_flow_network_3_2 Free
      p_flow_network_3_3 Free
 0.81 <= W_voltage_real_3_0 <= 1.21
      W_voltage_real_3_1 Free
      W_voltage_real_3_2 Free
      W_voltage_real_3_3 Free
      W_voltage_real_3_4 Free
 0.81 <= W_voltage_real_3_5 <= 1.21
      W_voltage_real_3_6 Free
      W_voltage_real_3_7 Free
      W_voltage_real_3_8 Free
      W_voltage_real_3_9 Free
 0.81 <= W_voltage_real_3_10 <= 1.21
      W_voltage_real_3_11 Free
      W_voltage_real_3_12 Free
      W_voltage_real_3_13 Free
      W_voltage_real_3_14 Free
 0.81 <= W_voltage_real_3_15 <= 1.21
      W_voltage_imag_3_0 Free
      W_voltage_imag_3_1 Free
      W_voltage_imag_3_2 Free
      W_voltage_imag_3_3 Free
      W_voltage_imag_3_4 Free
      W_voltage_imag_3_5 Free
      W_voltage_imag_3_6 Free
      W_voltage_imag_3_7 Free
      W_voltage_imag_3_8 Free
      W_voltage_imag_3_9 Free
      W_voltage_imag_3_10 Free
      W_voltage_imag_3_11 Free
      W_voltage_imag_3_12 Free
      W_voltage_imag_3_13 Free
      W_voltage_imag_3_14 Free
      W_voltage_imag_3_15 Free
      S_power_flow_3_0 Free
      S_power_flow_3_1 Free
      S_power_flow_3_2 Free
      S_power_flow_3_3 Free
      S_power_flow_3_4 Free
      S_power_flow_3_5 Free
      S_power_flow_3_6 Free
      S_power_flow_3_7 Free
      S_power_flow_3_8 Free
      S_power_flow_3_9 Free
      S_power_flow_3_10 Free
      S_power_flow_3_11 Free
      S_power_flow_3_12 Free
      S_power_flow_3_13 Free
      S_power_flow_3_14 Free
      S_power_flow_3_15 Free
      quad_aux_var_0 Free
      quad_aux_var_1 Free
      quad_aux_var_2 Free
      quad_aux_var_3 Free
Binaries
 u_thermal_1_0  v_thermal_1_0  w_thermal_1_0  u_thermal_2_0  v_thermal_2_0 
 w_thermal_2_0 
End
