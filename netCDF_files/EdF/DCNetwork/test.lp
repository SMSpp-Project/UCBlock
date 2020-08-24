\ENCODING=ISO-8859-1
\Problem name: testCPX

Minimize
 obj: 4 p_1_0 + 4 p_1_1 + 4 p_1_2 + 4 p_1_3 + 0 pr_1_0 + 0 pr_1_1 + 0 pr_1_2
      + 0 pr_1_3 + 0 sr_1_0 + 0 sr_1_1 + 0 sr_1_2 + 0 sr_1_3 + 10 p_2_0
      + 10 p_2_1 + 10 p_2_2 + 10 p_2_3 + 0 pr_2_0 + 0 pr_2_1 + 0 pr_2_2
      + 0 pr_2_3 + 0 sr_2_0 + 0 sr_2_1 + 0 sr_2_2 + 0 sr_2_3 + 20 p_3_0
      + 20 p_3_1 + 20 p_3_2 + 20 p_3_3 + 0 pr_3_0 + 0 pr_3_1 + 0 pr_3_2
      + 0 pr_3_3 + 0 sr_3_0 + 0 sr_3_1 + 0 sr_3_2 + 0 sr_3_3 + 10 p_4_0
      + 10 p_4_1 + 10 p_4_2 + 10 p_4_3 + 0 pr_4_0 + 0 pr_4_1 + 0 pr_4_2
      + 0 pr_4_3 + 0 sr_4_0 + 0 sr_4_1 + 0 sr_4_2 + 0 sr_4_3 + 4 p_5_0
      + 4 p_5_1 + 4 p_5_2 + 4 p_5_3 + 0 pr_5_0 + 0 pr_5_1 + 0 pr_5_2 + 0 pr_5_3
      + 0 sr_5_0 + 0 sr_5_1 + 0 sr_5_2 + 0 sr_5_3
Subject To
 node_injection_c_0_0:          - S_6_0  = -50
 node_injection_c_0_1:          - S_6_1  = 0
 node_injection_c_0_2:          - S_6_2  = 0
 node_injection_c_0_3:          - S_6_3  = 0
 node_injection_c_0_4:          - S_6_4  = 0
 node_injection_c_0_5:          - S_7_0  = -50
 node_injection_c_0_6:          - S_7_1  = 0
 node_injection_c_0_7:          - S_7_2  = 0
 node_injection_c_0_8:          - S_7_3  = 0
 node_injection_c_0_9:          - S_7_4  = 0
 node_injection_c_0_10:         - S_8_0  = -50
 node_injection_c_0_11:         - S_8_1  = 0
 node_injection_c_0_12:         - S_8_2  = 0
 node_injection_c_0_13:         - S_8_3  = 0
 node_injection_c_0_14:         - S_8_4  = 0
 node_injection_c_0_15:         - S_9_0  = -50
 node_injection_c_0_16:         - S_9_1  = 0
 node_injection_c_0_17:         - S_9_2  = 0
 node_injection_c_0_18:         - S_9_3  = 0
 node_injection_c_0_19:         - S_9_4  = 0
 StartUp_ShutDown_Vars_1_0:     u_1_0 - v_1_0 + w_1_0  = 1
 StartUp_ShutDown_Vars_1_1:     - u_1_0 + u_1_1 - v_1_1 + w_1_1  = 0
 StartUp_ShutDown_Vars_1_2:     - u_1_1 + u_1_2 - v_1_2 + w_1_2  = 0
 StartUp_ShutDown_Vars_1_3:     - u_1_2 + u_1_3 - v_1_3 + w_1_3  = 0
 RampUp_1_0:                    p_1_0 - 200 v_1_0 - RgRampUp_1_0  = 0
 RampUp_1_1:                    900 u_1_0 + p_1_0 - p_1_1 + 200 v_1_1 >= 0
 RampUp_1_2:                    900 u_1_1 + p_1_1 - p_1_2 + 200 v_1_2 >= 0
 RampUp_1_3:                    900 u_1_2 + p_1_2 - p_1_3 + 200 v_1_3 >= 0
 RampDown_1_0:                  900 u_1_0 + p_1_0 + 200 w_1_0 >= 500
 RampDown_1_1:                  900 u_1_1 - p_1_0 + p_1_1 + 200 w_1_1 >= 0
 RampDown_1_2:                  900 u_1_2 - p_1_1 + p_1_2 + 200 w_1_2 >= 0
 RampDown_1_3:                  900 u_1_3 - p_1_2 + p_1_3 + 200 w_1_3 >= 0
 MinPower_1_0:                  - 200 u_1_0 + p_1_0 >= 0
 MinPower_1_1:                  - 200 u_1_1 + p_1_1 >= 0
 MinPower_1_2:                  - 200 u_1_2 + p_1_2 >= 0
 MinPower_1_3:                  - 200 u_1_3 + p_1_3 >= 0
 MaxPower_1_0:                  925 u_1_0 - p_1_0 >= 0
 MaxPower_1_1:                  925 u_1_1 - p_1_1 >= 0
 MaxPower_1_2:                  925 u_1_2 - p_1_2 >= 0
 MaxPower_1_3:                  925 u_1_3 - p_1_3 >= 0
 StartUp_ShutDown_Vars_2_0:     u_2_0 - v_2_0 + w_2_0  = 1
 StartUp_ShutDown_Vars_2_1:     - u_2_0 + u_2_1 - v_2_1 + w_2_1  = 0
 StartUp_ShutDown_Vars_2_2:     - u_2_1 + u_2_2 - v_2_2 + w_2_2  = 0
 StartUp_ShutDown_Vars_2_3:     - u_2_2 + u_2_3 - v_2_3 + w_2_3  = 0
 StartUp_2_0:                   u_2_1 - v_2_1 >= 0
 StartUp_2_1:                   u_2_2 - v_2_2 >= 0
 StartUp_2_2:                   u_2_3 - v_2_3 >= 0
 ShutDown_2_0:                  u_2_1 + w_2_1 - RgShutDown_2_0  = 0
 ShutDown_2_1:                  u_2_2 + w_2_2 - RgShutDown_2_1  = 0
 ShutDown_2_2:                  u_2_3 + w_2_3 - RgShutDown_2_2  = 0
 RampUp_2_0:                    p_2_0 - 100 v_2_0 - RgRampUp_2_0  = 0
 RampUp_2_1:                    500 u_2_0 + p_2_0 - p_2_1 + 100 v_2_1 >= 0
 RampUp_2_2:                    500 u_2_1 + p_2_1 - p_2_2 + 100 v_2_2 >= 0
 RampUp_2_3:                    500 u_2_2 + p_2_2 - p_2_3 + 100 v_2_3 >= 0
 RampDown_2_0:                  500 u_2_0 + p_2_0 + 100 w_2_0 >= 500
 RampDown_2_1:                  500 u_2_1 - p_2_0 + p_2_1 + 100 w_2_1 >= 0
 RampDown_2_2:                  500 u_2_2 - p_2_1 + p_2_2 + 100 w_2_2 >= 0
 RampDown_2_3:                  500 u_2_3 - p_2_2 + p_2_3 + 100 w_2_3 >= 0
 MinPower_2_0:                  - 100 u_2_0 + p_2_0 >= 0
 MinPower_2_1:                  - 100 u_2_1 + p_2_1 >= 0
 MinPower_2_2:                  - 100 u_2_2 + p_2_2 >= 0
 MinPower_2_3:                  - 100 u_2_3 + p_2_3 >= 0
 MaxPower_2_0:                  1000 u_2_0 - p_2_0 >= 0
 MaxPower_2_1:                  1000 u_2_1 - p_2_1 >= 0
 MaxPower_2_2:                  1000 u_2_2 - p_2_2 >= 0
 MaxPower_2_3:                  1000 u_2_3 - p_2_3 >= 0
 StartUp_ShutDown_Vars_3_0:     u_3_0 - v_3_0 + w_3_0  = 0
 StartUp_ShutDown_Vars_3_1:     - u_3_0 + u_3_1 - v_3_1 + w_3_1  = 0
 StartUp_ShutDown_Vars_3_2:     - u_3_1 + u_3_2 - v_3_2 + w_3_2  = 0
 StartUp_ShutDown_Vars_3_3:     - u_3_2 + u_3_3 - v_3_3 + w_3_3  = 0
 StartUp_3_0:                   u_3_1 - v_3_1 >= 0
 StartUp_3_1:                   u_3_2 - v_3_2 >= 0
 StartUp_3_2:                   u_3_3 - v_3_3 >= 0
 ShutDown_3_0:                  u_3_1 + w_3_1 - RgShutDown_3_0  = 0
 ShutDown_3_1:                  u_3_2 + w_3_2 - RgShutDown_3_1  = 0
 ShutDown_3_2:                  u_3_3 + w_3_3 - RgShutDown_3_2  = 0
 RampUp_3_0:                    p_3_0 - 100 v_3_0  = 0
 RampUp_3_1:                    500 u_3_0 + p_3_0 - p_3_1 + 100 v_3_1 >= 0
 RampUp_3_2:                    500 u_3_1 + p_3_1 - p_3_2 + 100 v_3_2 >= 0
 RampUp_3_3:                    500 u_3_2 + p_3_2 - p_3_3 + 100 v_3_3 >= 0
 RampDown_3_0:                  500 u_3_0 + p_3_0 + 100 w_3_0 >= 0
 RampDown_3_1:                  500 u_3_1 - p_3_0 + p_3_1 + 100 w_3_1 >= 0
 RampDown_3_2:                  500 u_3_2 - p_3_1 + p_3_2 + 100 w_3_2 >= 0
 RampDown_3_3:                  500 u_3_3 - p_3_2 + p_3_3 + 100 w_3_3 >= 0
 MinPower_3_0:                  - 100 u_3_0 + p_3_0 >= 0
 MinPower_3_1:                  - 100 u_3_1 + p_3_1 >= 0
 MinPower_3_2:                  - 100 u_3_2 + p_3_2 >= 0
 MinPower_3_3:                  - 100 u_3_3 + p_3_3 >= 0
 MaxPower_3_0:                  500 u_3_0 - p_3_0 >= 0
 MaxPower_3_1:                  500 u_3_1 - p_3_1 >= 0
 MaxPower_3_2:                  500 u_3_2 - p_3_2 >= 0
 MaxPower_3_3:                  500 u_3_3 - p_3_3 >= 0
 StartUp_ShutDown_Vars_4_0:     u_4_0 - v_4_0 + w_4_0  = 1
 StartUp_ShutDown_Vars_4_1:     - u_4_0 + u_4_1 - v_4_1 + w_4_1  = 0
 StartUp_ShutDown_Vars_4_2:     - u_4_1 + u_4_2 - v_4_2 + w_4_2  = 0
 StartUp_ShutDown_Vars_4_3:     - u_4_2 + u_4_3 - v_4_3 + w_4_3  = 0
 StartUp_4_0:                   u_4_1 - v_4_1 >= 0
 StartUp_4_1:                   u_4_2 - v_4_2 >= 0
 StartUp_4_2:                   u_4_3 - v_4_3 >= 0
 ShutDown_4_0:                  u_4_1 + w_4_1 - RgShutDown_4_0  = 0
 ShutDown_4_1:                  u_4_2 + w_4_2 - RgShutDown_4_1  = 0
 ShutDown_4_2:                  u_4_3 + w_4_3 - RgShutDown_4_2  = 0
 RampUp_4_0:                    p_4_0 - 100 v_4_0 - RgRampUp_4_0  = 0
 RampUp_4_1:                    500 u_4_0 + p_4_0 - p_4_1 + 100 v_4_1 >= 0
 RampUp_4_2:                    500 u_4_1 + p_4_1 - p_4_2 + 100 v_4_2 >= 0
 RampUp_4_3:                    500 u_4_2 + p_4_2 - p_4_3 + 100 v_4_3 >= 0
 RampDown_4_0:                  500 u_4_0 + p_4_0 + 100 w_4_0 >= 500
 RampDown_4_1:                  500 u_4_1 - p_4_0 + p_4_1 + 100 w_4_1 >= 0
 RampDown_4_2:                  500 u_4_2 - p_4_1 + p_4_2 + 100 w_4_2 >= 0
 RampDown_4_3:                  500 u_4_3 - p_4_2 + p_4_3 + 100 w_4_3 >= 0
 MinPower_4_0:                  - 100 u_4_0 + p_4_0 >= 0
 MinPower_4_1:                  - 100 u_4_1 + p_4_1 >= 0
 MinPower_4_2:                  - 100 u_4_2 + p_4_2 >= 0
 MinPower_4_3:                  - 100 u_4_3 + p_4_3 >= 0
 MaxPower_4_0:                  1000 u_4_0 - p_4_0 >= 0
 MaxPower_4_1:                  1000 u_4_1 - p_4_1 >= 0
 MaxPower_4_2:                  1000 u_4_2 - p_4_2 >= 0
 MaxPower_4_3:                  1000 u_4_3 - p_4_3 >= 0
 StartUp_ShutDown_Vars_5_0:     u_5_0 - v_5_0 + w_5_0  = 1
 StartUp_ShutDown_Vars_5_1:     - u_5_0 + u_5_1 - v_5_1 + w_5_1  = 0
 StartUp_ShutDown_Vars_5_2:     - u_5_1 + u_5_2 - v_5_2 + w_5_2  = 0
 StartUp_ShutDown_Vars_5_3:     - u_5_2 + u_5_3 - v_5_3 + w_5_3  = 0
 RampUp_5_0:                    p_5_0 - 200 v_5_0 - RgRampUp_5_0  = 0
 RampUp_5_1:                    900 u_5_0 + p_5_0 - p_5_1 + 200 v_5_1 >= 0
 RampUp_5_2:                    900 u_5_1 + p_5_1 - p_5_2 + 200 v_5_2 >= 0
 RampUp_5_3:                    900 u_5_2 + p_5_2 - p_5_3 + 200 v_5_3 >= 0
 RampDown_5_0:                  900 u_5_0 + p_5_0 + 200 w_5_0 >= 500
 RampDown_5_1:                  900 u_5_1 - p_5_0 + p_5_1 + 200 w_5_1 >= 0
 RampDown_5_2:                  900 u_5_2 - p_5_1 + p_5_2 + 200 w_5_2 >= 0
 RampDown_5_3:                  900 u_5_3 - p_5_2 + p_5_3 + 200 w_5_3 >= 0
 MinPower_5_0:                  - 200 u_5_0 + p_5_0 >= 0
 MinPower_5_1:                  - 200 u_5_1 + p_5_1 >= 0
 MinPower_5_2:                  - 200 u_5_2 + p_5_2 >= 0
 MinPower_5_3:                  - 200 u_5_3 + p_5_3 >= 0
 MaxPower_5_0:                  925 u_5_0 - p_5_0 >= 0
 MaxPower_5_1:                  925 u_5_1 - p_5_1 >= 0
 MaxPower_5_2:                  925 u_5_2 - p_5_2 >= 0
 MaxPower_5_3:                  925 u_5_3 - p_5_3 >= 0
 HVDC_power_flow_limit_6_0:     f_6_0 - RgHVDC_power_flow_limit_6_0  = -100
 HVDC_power_flow_limit_6_1:     f_6_1 - RgHVDC_power_flow_limit_6_1  = -500
 HVDC_power_flow_limit_6_2:     f_6_2 - RgHVDC_power_flow_limit_6_2  = -250
 HVDC_power_flow_limit_6_3:     f_6_3 - RgHVDC_power_flow_limit_6_3  = -500
 HVDC_power_flow_limit_6_4:     f_6_4 - RgHVDC_power_flow_limit_6_4  = -100
 HVDC_power_flow_injection_6_0: - S_6_0 + f_6_0  = 0
 HVDC_power_flow_injection_6_1: - S_6_1 - f_6_0 + f_6_1 + f_6_2  = 0
 HVDC_power_flow_injection_6_2: - S_6_2 - f_6_1 + f_6_3  = 0
 HVDC_power_flow_injection_6_3: - S_6_3 - f_6_2 - f_6_3 + f_6_4  = 0
 HVDC_power_flow_injection_6_4: - S_6_4 - f_6_4  = 0
 HVDC_power_flow_limit_7_0:     f_7_0 - RgHVDC_power_flow_limit_7_0  = -100
 HVDC_power_flow_limit_7_1:     f_7_1 - RgHVDC_power_flow_limit_7_1  = -500
 HVDC_power_flow_limit_7_2:     f_7_2 - RgHVDC_power_flow_limit_7_2  = -250
 HVDC_power_flow_limit_7_3:     f_7_3 - RgHVDC_power_flow_limit_7_3  = -500
 HVDC_power_flow_limit_7_4:     f_7_4 - RgHVDC_power_flow_limit_7_4  = -100
 HVDC_power_flow_injection_7_0: - S_7_0 + f_7_0  = 0
 HVDC_power_flow_injection_7_1: - S_7_1 - f_7_0 + f_7_1 + f_7_2  = 0
 HVDC_power_flow_injection_7_2: - S_7_2 - f_7_1 + f_7_3  = 0
 HVDC_power_flow_injection_7_3: - S_7_3 - f_7_2 - f_7_3 + f_7_4  = 0
 HVDC_power_flow_injection_7_4: - S_7_4 - f_7_4  = 0
 HVDC_power_flow_limit_8_0:     f_8_0 - RgHVDC_power_flow_limit_8_0  = -100
 HVDC_power_flow_limit_8_1:     f_8_1 - RgHVDC_power_flow_limit_8_1  = -500
 HVDC_power_flow_limit_8_2:     f_8_2 - RgHVDC_power_flow_limit_8_2  = -250
 HVDC_power_flow_limit_8_3:     f_8_3 - RgHVDC_power_flow_limit_8_3  = -500
 HVDC_power_flow_limit_8_4:     f_8_4 - RgHVDC_power_flow_limit_8_4  = -100
 HVDC_power_flow_injection_8_0: - S_8_0 + f_8_0  = 0
 HVDC_power_flow_injection_8_1: - S_8_1 - f_8_0 + f_8_1 + f_8_2  = 0
 HVDC_power_flow_injection_8_2: - S_8_2 - f_8_1 + f_8_3  = 0
 HVDC_power_flow_injection_8_3: - S_8_3 - f_8_2 - f_8_3 + f_8_4  = 0
 HVDC_power_flow_injection_8_4: - S_8_4 - f_8_4  = 0
 HVDC_power_flow_limit_9_0:     f_9_0 - RgHVDC_power_flow_limit_9_0  = -100
 HVDC_power_flow_limit_9_1:     f_9_1 - RgHVDC_power_flow_limit_9_1  = -500
 HVDC_power_flow_limit_9_2:     f_9_2 - RgHVDC_power_flow_limit_9_2  = -250
 HVDC_power_flow_limit_9_3:     f_9_3 - RgHVDC_power_flow_limit_9_3  = -500
 HVDC_power_flow_limit_9_4:     f_9_4 - RgHVDC_power_flow_limit_9_4  = -100
 HVDC_power_flow_injection_9_0: - S_9_0 + f_9_0  = 0
 HVDC_power_flow_injection_9_1: - S_9_1 - f_9_0 + f_9_1 + f_9_2  = 0
 HVDC_power_flow_injection_9_2: - S_9_2 - f_9_1 + f_9_3  = 0
 HVDC_power_flow_injection_9_3: - S_9_3 - f_9_2 - f_9_3 + f_9_4  = 0
 HVDC_power_flow_injection_9_4: - S_9_4 - f_9_4  = 0
Bounds
 0 <= u_1_0 <= 1
 0 <= u_1_1 <= 1
 0 <= u_1_2 <= 1
 0 <= u_1_3 <= 1
      v_1_0 = 0
      v_1_1 = 0
      v_1_2 = 0
      v_1_3 = 0
 0 <= w_1_0 <= 1
 0 <= w_1_1 <= 1
 0 <= w_1_2 <= 1
 0 <= w_1_3 <= 1
 0 <= u_2_0 <= 1
 0 <= u_2_1 <= 1
 0 <= u_2_2 <= 1
 0 <= u_2_3 <= 1
      v_2_0 = 0
 0 <= v_2_1 <= 1
 0 <= v_2_2 <= 1
 0 <= v_2_3 <= 1
 0 <= w_2_0 <= 1
 0 <= w_2_1 <= 1
 0 <= w_2_2 <= 1
 0 <= w_2_3 <= 1
 0 <= u_3_0 <= 1
 0 <= u_3_1 <= 1
 0 <= u_3_2 <= 1
 0 <= u_3_3 <= 1
 0 <= v_3_0 <= 1
 0 <= v_3_1 <= 1
 0 <= v_3_2 <= 1
 0 <= v_3_3 <= 1
      w_3_0 = 0
 0 <= w_3_1 <= 1
 0 <= w_3_2 <= 1
 0 <= w_3_3 <= 1
 0 <= u_4_0 <= 1
 0 <= u_4_1 <= 1
 0 <= u_4_2 <= 1
 0 <= u_4_3 <= 1
      v_4_0 = 0
 0 <= v_4_1 <= 1
 0 <= v_4_2 <= 1
 0 <= v_4_3 <= 1
 0 <= w_4_0 <= 1
 0 <= w_4_1 <= 1
 0 <= w_4_2 <= 1
 0 <= w_4_3 <= 1
 0 <= u_5_0 <= 1
 0 <= u_5_1 <= 1
 0 <= u_5_2 <= 1
 0 <= u_5_3 <= 1
      v_5_0 = 0
      v_5_1 = 0
      v_5_2 = 0
      v_5_3 = 0
 0 <= w_5_0 <= 1
 0 <= w_5_1 <= 1
 0 <= w_5_2 <= 1
 0 <= w_5_3 <= 1
      S_6_0 Free
      S_6_1 Free
      S_6_2 Free
      S_6_3 Free
      S_6_4 Free
      f_6_0 Free
      f_6_1 Free
      f_6_2 Free
      f_6_3 Free
      f_6_4 Free
      S_7_0 Free
      S_7_1 Free
      S_7_2 Free
      S_7_3 Free
      S_7_4 Free
      f_7_0 Free
      f_7_1 Free
      f_7_2 Free
      f_7_3 Free
      f_7_4 Free
      S_8_0 Free
      S_8_1 Free
      S_8_2 Free
      S_8_3 Free
      S_8_4 Free
      f_8_0 Free
      f_8_1 Free
      f_8_2 Free
      f_8_3 Free
      f_8_4 Free
      S_9_0 Free
      S_9_1 Free
      S_9_2 Free
      S_9_3 Free
      S_9_4 Free
      f_9_0 Free
      f_9_1 Free
      f_9_2 Free
      f_9_3 Free
      f_9_4 Free
 0 <= RgRampUp_1_0 <= 1400
 0 <= RgShutDown_2_0 <= 1
 0 <= RgShutDown_2_1 <= 1
 0 <= RgShutDown_2_2 <= 1
 0 <= RgRampUp_2_0 <= 1000
 0 <= RgShutDown_3_0 <= 1
 0 <= RgShutDown_3_1 <= 1
 0 <= RgShutDown_3_2 <= 1
 0 <= RgShutDown_4_0 <= 1
 0 <= RgShutDown_4_1 <= 1
 0 <= RgShutDown_4_2 <= 1
 0 <= RgRampUp_4_0 <= 1000
 0 <= RgRampUp_5_0 <= 1400
 0 <= RgHVDC_power_flow_limit_6_0 <= 200
 0 <= RgHVDC_power_flow_limit_6_1 <= 1000
 0 <= RgHVDC_power_flow_limit_6_2 <= 500
 0 <= RgHVDC_power_flow_limit_6_3 <= 1000
 0 <= RgHVDC_power_flow_limit_6_4 <= 200
 0 <= RgHVDC_power_flow_limit_7_0 <= 200
 0 <= RgHVDC_power_flow_limit_7_1 <= 1000
 0 <= RgHVDC_power_flow_limit_7_2 <= 500
 0 <= RgHVDC_power_flow_limit_7_3 <= 1000
 0 <= RgHVDC_power_flow_limit_7_4 <= 200
 0 <= RgHVDC_power_flow_limit_8_0 <= 200
 0 <= RgHVDC_power_flow_limit_8_1 <= 1000
 0 <= RgHVDC_power_flow_limit_8_2 <= 500
 0 <= RgHVDC_power_flow_limit_8_3 <= 1000
 0 <= RgHVDC_power_flow_limit_8_4 <= 200
 0 <= RgHVDC_power_flow_limit_9_0 <= 200
 0 <= RgHVDC_power_flow_limit_9_1 <= 1000
 0 <= RgHVDC_power_flow_limit_9_2 <= 500
 0 <= RgHVDC_power_flow_limit_9_3 <= 1000
 0 <= RgHVDC_power_flow_limit_9_4 <= 200
Binaries
 u_1_0  u_1_1  u_1_2  u_1_3  v_1_0  v_1_1  v_1_2  v_1_3  w_1_0  w_1_1  w_1_2 
 w_1_3  u_2_0  u_2_1  u_2_2  u_2_3  v_2_0  v_2_1  v_2_2  v_2_3  w_2_0  w_2_1 
 w_2_2  w_2_3  u_3_0  u_3_1  u_3_2  u_3_3  v_3_0  v_3_1  v_3_2  v_3_3  w_3_0 
 w_3_1  w_3_2  w_3_3  u_4_0  u_4_1  u_4_2  u_4_3  v_4_0  v_4_1  v_4_2  v_4_3 
 w_4_0  w_4_1  w_4_2  w_4_3  u_5_0  u_5_1  u_5_2  u_5_3  v_5_0  v_5_1  v_5_2 
 v_5_3  w_5_0  w_5_1  w_5_2  w_5_3 
End
