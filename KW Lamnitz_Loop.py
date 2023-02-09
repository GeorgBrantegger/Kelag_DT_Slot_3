# code cell 0
import os
import sys
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np

current = os.path.dirname(os.path.realpath('Main_Programm.ipynb'))
parent = os.path.dirname(current)
sys.path.append(parent)
from Ausgleichsbecken.Ausgleichsbecken_class_file import Ausgleichsbecken_class
from Druckrohrleitung.Druckrohrleitung_class_file import Druckrohrleitung_class
from functions.pressure_conversion import pressure_conversion
from Kraftwerk.Kraftwerk_class_file import Kraftwerk_class
from Regler.Regler_class_file import PI_controller_class
from Turbinen.Turbinen_class_file import Turbine

# code cell 1
# for loop creation

Area_list = np.round(np.arange(20.,30.,5.),1)
Kp_list   = np.round(np.arange(0.7,1.3,0.2),1)
Ti_list   = np.round(np.arange(200.,220.,25.),1)

# # if one wants to use the loop to save 1 specific configuration:
# desired_area    = 60
# desired_KP      = 0.7
# desired_ti      = 200.

# Area_list   = np.round(np.arange(desired_area,desired_area+1.,1.),1)
# Kp_list     = np.round(np.arange(desired_KP,desired_KP+1.,1),1)
# Ti_list     = np.round(np.arange(desired_ti,desired_ti+1.,1.),1)

for i in range(np.size(Area_list)):
    for j in range(np.size(Kp_list)):
        for k in range(np.size(Ti_list)):
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            print("Current Time =", current_time)

            print('i = ',i, '/ ', str(np.size(Area_list)-1))
            print('j = ',j, '/ ', str(np.size(Kp_list)-1))
            print('k = ',k, '/ ', str(np.size(Ti_list)-1))
            print('area = ',Area_list[i])
            print('K_p = ',Kp_list[j])
            print('T_i = ',Ti_list[k])

            with open('log.txt','a') as f:
                f.write("Current Time =" + current_time + '\n')
                f.write('i = '+str(i)+ '/ '+ str(np.size(Area_list)-1)+ '\n')
                f.write('j = '+str(j)+ '/ '+ str(np.size(Kp_list)-1)+ '\n')
                f.write('k = '+str(k)+ '/ '+ str(np.size(Ti_list)-1)+ '\n')
                f.write('area = '+str(Area_list[i])+ '\n')
                f.write('K_p = '+str(Kp_list[j])+ '\n')
                f.write('T_i = '+str(Ti_list[k])+ '\n')


            # code cell 2
            # define constants

                # for physics
            g                   = 9.81                                          # [m/s²]    gravitational acceleration 
            rho                 = 0.9982067*1e3                                 # [kg/m³]   density of water 
            pUnit_calc          = 'Pa'                                          # [string]  DO NOT CHANGE! for pressure conversion in print statements and plot labels 
            pUnit_conv          = 'mWS'                                         # [string]  for pressure conversion in print statements and plot labels

                # for KW OL 
            OL_T1_Q_nenn        = 1.0                                           # [m³/s]    nominal flux of turbine 
            OL_T1_p_nenn        = pressure_conversion(1,'bar',pUnit_calc)       # [Pa]      nominal pressure of turbine  ## p_nenn wird konstant gehalten, Wert ist also fiktiv
            OL_p_pseudo         = 1.1*OL_T1_p_nenn                              # ficticious pressure applied to OL turbines to avoid LA>1 error caused by unfortunate rounding
            OL_T1_closingTime   = 600.                                           # [s]       closing time of turbine

                # for KW UL
            UL_T1_Q_nenn        = 1.1                                           # [m³/s]    nominal flux of turbine 
            UL_T1_p_nenn        = pressure_conversion(120.,'mWS',pUnit_calc)     # [Pa]      nominal pressure of turbine 
            UL_T1_closingTime   = 60.                                           # [s]       closing time of turbine

                # for PI controller
            Con_targetLevel     = 1.25                                          # [m]       target level of the PI controller
            Con_K_p             = Kp_list[j]                                    # [-]       proportionality constant of PI controller
            Con_T_i             = Ti_list[k]                                          # [s]       timespan in which a steady state error is corrected by the intergal term
            Con_deadbandRange   = 0.00                                          # [m]       Deadband range around targetLevel for which the controller does NOT intervene

                # for pipeline
            Pip_length          = 2000.                                         # [m]       length of pipeline
            Pip_dia             = 0.9                                           # [m]       diameter of pipeline
            Pip_area            = Pip_dia**2/4*np.pi                            # [m²]      crossectional area of pipeline
            Pip_head            = 130.                                           # [m]       hydraulic head of pipeline without reservoir
            Pip_angle           = np.arcsin(Pip_head/Pip_length)                # [rad]     elevation angle of pipeline 
            Pip_n_seg           = 50                                            # [-]       number of pipe segments in discretization
            Pip_f_D             = 0.015                                         # [-]       Darcy friction factor
            Pip_pw_vel          = 600.                                          # [m/s]     propagation velocity of the pressure wave (pw) in the given pipeline
                # derivatives of the pipeline constants
            Pip_dx              = Pip_length/Pip_n_seg                          # [m]       length of each pipe segment
            Pip_dt              = Pip_dx/Pip_pw_vel                             # [s]       timestep according to method of characteristics
            Pip_nn              = Pip_n_seg+1                                   # [1]       number of nodes
            Pip_x_vec           = np.arange(0,Pip_nn,1)*Pip_dx                  # [m]       vector holding the distance of each node from the upstream reservoir along the pipeline
            Pip_h_vec           = np.arange(0,Pip_nn,1)*Pip_head/Pip_n_seg      # [m]       vector holding the vertical distance of each node from the upstream reservoir

                # for reservoir
            Res_area_base       = Area_list[i]                                  # [m²]      total base are of the cuboid reservoir   
            Res_area_out        = Pip_area                                      # [m²]      outflux area of the reservoir, given by pipeline area
            Res_level_crit_lo   = Con_targetLevel-0.5                           # [m]       for yet-to-be-implemented warnings
            Res_level_crit_hi   = np.inf                                        # [m]       for yet-to-be-implemented warnings
            Res_dt_approx       = 1e-3                                          # [s]       approx. timestep of reservoir time evolution to ensure numerical stability (see Res_nt why approx.)
            Res_nt              = max(1,int(Pip_dt//Res_dt_approx))             # [1]       number of timesteps of the reservoir time evolution within one timestep of the pipeline
            Res_dt              = Pip_dt/Res_nt                                 # [s]       harmonised timestep of reservoir time evolution

                # for general simulation
            flux_init           = OL_T1_Q_nenn                                  # [m³/s]    initial flux through whole system for steady state initialization  
            #OL_LAs_init         = [1.,0.3]                                     # [vec]     initial guide vane openings of OL-KW
            level_init          = Con_targetLevel                               # [m]       initial water level in upstream reservoir for steady state initialization
            simTime_target      = 1200.                                         # [s]       target for total simulation time (will vary slightly to fit with Pip_dt)
            nt                  = int(simTime_target//Pip_dt)                   # [1]       Number of timesteps of the whole system
            t_vec               = np.arange(0,nt+1,1)*Pip_dt                    # [s]       time vector. At each step of t_vec the system parameters are stored


            # code cell 3
            # create objects

            # influx setting turbines
            OL_T1 = Turbine(OL_T1_Q_nenn,OL_T1_p_nenn,OL_T1_closingTime,Pip_dt,pUnit_conv)

            KW_OL = Kraftwerk_class()
            KW_OL.add_turbine(OL_T1)

            KW_OL.set_steady_state_by_flux(flux_init,OL_p_pseudo)

            # KW_OL.set_steady_state_by_LA(OL_LAs_init,OL_p_pseudo)
            # flux_init = KW_OL.get_current_Q()

            # Upstream reservoir
            reservoir = Ausgleichsbecken_class(Res_area_base,Res_area_out,Res_dt,pUnit_conv,Res_level_crit_lo,Res_level_crit_hi,rho)
            reservoir.set_steady_state(flux_init,level_init)

            # pipeline
            pipe = Druckrohrleitung_class(Pip_length,Pip_dia,Pip_head,Pip_n_seg,Pip_f_D,Pip_pw_vel,Pip_dt,pUnit_conv,rho)
            pipe.set_steady_state(flux_init,reservoir.get_current_pressure())

            # downstream turbines
            UL_T1 = Turbine(UL_T1_Q_nenn,UL_T1_p_nenn,UL_T1_closingTime,Pip_dt,pUnit_conv)

            KW_UL = Kraftwerk_class()
            KW_UL.add_turbine(UL_T1)

            KW_UL.set_steady_state_by_flux(flux_init,pipe.get_current_pressure_distribution()[-1])

            # level controller
            level_control = PI_controller_class(Con_targetLevel,Con_deadbandRange,Con_K_p,Con_T_i,Pip_dt)
            level_control.set_control_variable(UL_T1.get_current_LA(),display_warning=False)


            # code cell 5
            # initialization for Timeloop

            # OL KW
                # manual input to modulate influx
            OL_T1_LA_soll_vec = np.full_like(t_vec,OL_T1.get_current_LA())          # storing the target value for the guide van opening
            OL_T1_LA_soll_vec[np.argmin(np.abs(t_vec-100)):] = 0.                   # changing the target value for the guide vane opening at t = 100 s
            OL_T1_LA_soll_vec[np.argmin(np.abs(t_vec-600)):] = OL_T1_LA_soll_vec[0] # changing the target value for the guide vane opening at t = 600 s        

            # creating a bunch of vectors that are used to store usefull information - either for analysis or for the following step in the timeloop

            # reservoir
            Q_in_vec = np.zeros_like(t_vec)                                 # for storing the influx to the reservoir
            Q_in_vec[0] = flux_init                                         # storing the initial influx to the reservoir
            # Outflux from reservoir is stored in Q_boundary_res
            level_vec  = np.zeros_like(t_vec)                               # for storing the level in the reservoir at the end of each pipeline timestep
            level_vec[0] = level_init                                       # storing the initial level in the reservoir
            volume_vec = np.zeros_like(t_vec)                               # for storing the volume in the reservoir at the end of each pipeline timestep
            volume_vec[0] = reservoir.get_current_volume()                  # storing the initial volume in the reservoir

            # pipeline
            v_old = pipe.get_current_velocity_distribution()                # for storing the velocity from the last timestep
            v_min = pipe.get_lowest_velocity_per_node()                     # for storing minimal flux velocity at each node
            v_max = pipe.get_highest_velocity_per_node()                    # for storing maximal flux velocity at each node
            Q_old = pipe.get_current_flux_distribution()                    # for storing the flux from the last timestep
            Q_min = pipe.get_lowest_flux_per_node()                         # for storing minimal flux at each node
            Q_max = pipe.get_highest_flux_per_node()                        # for storing maximal flux at each node
            p_old = pipe.get_current_pressure_distribution()                # for storing the pressure from the last timestep
            p_min = pipe.get_lowest_pressure_per_node()                     # for storing minimal pressure at each node
            p_max = pipe.get_highest_pressure_per_node()                    # for storing maximal pressure at each node
            p_0   = pipe.get_initial_pressure_distribution()                # storing initial pressure at each node

            v_boundary_res  = np.zeros_like(t_vec)                          # for storing the boundary velocity at the reservoir
            v_boundary_tur  = np.zeros_like(t_vec)                          # for storing the boundary velocity at the turbine
            Q_boundary_res  = np.zeros_like(t_vec)                          # for storing the boundary flux at the reservoir
            Q_boundary_tur  = np.zeros_like(t_vec)                          # for storing the boundary flux at the turbine
            p_boundary_res  = np.zeros_like(t_vec)                          # for storing the boundary pressure at the reservoir
            p_boundary_tur  = np.zeros_like(t_vec)                          # for storing the boundary pressure at the turbine

            v_boundary_res[0] = v_old[0]                                    # storing the initial value for the boundary velocity at the reservoir
            v_boundary_tur[0] = v_old[-1]                                   # storing the initial value for the boundary velocity at the turbine
            Q_boundary_res[0] = Q_old[0]                                    # storing the initial value for the boundary flux at the reservoir
            Q_boundary_tur[0] = Q_old[-1]                                   # storing the initial value for the boundary flux at the turbine
            p_boundary_res[0] = p_old[0]                                    # storing the initial value for the boundary pressure at the reservoir
            p_boundary_tur[0] = p_old[-1]                                   # storing the initial value for the boundary pressure at the turbine

            # OL KW
            OL_T1_LA_ist_vec = np.zeros_like(t_vec)                         # for storing the actual value of the guide vane opening
            OL_T1_LA_ist_vec[0] = OL_T1.get_current_LA()                    # storing the initial value of the guide vane opening

            # UL KW
            UL_T1_LA_soll_vec = np.zeros_like(t_vec)                        # for storing the target value of the guide vane opening
            UL_T1_LA_soll_vec[0] = UL_T1.get_current_LA()                   # storing the initial value of the guide vane opening

            UL_T1_LA_ist_vec = np.zeros_like(t_vec)                         # for storing the actual value of the guide vane opening
            UL_T1_LA_ist_vec[0] = UL_T1.get_current_LA()                    # storing the initial value of the guide vane opening


            # code cell 8
            # time loop
            # needed for turbine convergence
            convergence_parameters = [p_old[-2],v_old[-2],Pip_dia,Pip_area,Pip_angle,Pip_f_D,Pip_pw_vel,rho,Pip_dt,p_old[-1]]

            # loop through time steps of the pipeline
            for it_pipe in range(1,nt+1):

                # update OL_KW and the influx into the reservoir
                KW_OL.update_LAs([OL_T1_LA_soll_vec[it_pipe]])
                KW_OL.set_pressure(OL_p_pseudo)
                Q_in_vec[it_pipe] = KW_OL.get_current_Q()
                reservoir.set_influx(Q_in_vec[it_pipe])

            # for each pipeline timestep, execute Res_nt timesteps of the reservoir code
                # set initial condition for the reservoir time evolution calculted with the timestep_reservoir_evolution() method
                reservoir.set_pressure(p_old[0],display_warning=False)
                reservoir.set_outflux(Q_old[0],display_warning=False)
                # calculate the time evolution of the reservoir level within each pipeline timestep to avoid runaway numerical error
                for it_res in range(Res_nt):
                    reservoir.timestep_reservoir_evolution()   
                # save the level and the volume in the reservoir                                                              
                level_vec[it_pipe]  = reservoir.get_current_level()                                                 
                volume_vec[it_pipe] = reservoir.get_current_volume()         

                # update target value for UL_KW from the level controller
                level_control.update_control_variable(level_vec[it_pipe])
                UL_T1_LA_soll_vec[it_pipe] = level_control.get_current_control_variable()                                            
                
                # change the guide vane opening based on the target value and closing time limitation
                KW_UL.update_LAs([UL_T1_LA_soll_vec[it_pipe]])
                    # save the actual guide vane openings
                OL_T1_LA_ist_vec[it_pipe] = KW_OL.get_current_LAs()
                UL_T1_LA_ist_vec[it_pipe] = KW_UL.get_current_LAs()

                # set boundary condition for the next timestep of the characteristic method
                convergence_parameters[0] = p_old[-2]
                convergence_parameters[1] = v_old[-2]
                convergence_parameters[9] = p_old[-1]
                KW_UL.set_pressure(p_old[-1])
                    # use the convergence method to avoid numerical errors
                KW_UL.converge(convergence_parameters)
                    # save the first set of boundary conditions
                p_boundary_res[it_pipe] = reservoir.get_current_pressure()
                v_boundary_tur[it_pipe] = 1/Pip_area*KW_UL.get_current_Q()
                Q_boundary_tur[it_pipe] = KW_UL.get_current_Q()

                # set the the boundary condition in the pipe and thereby calculate boundary pressure at turbine
                pipe.set_boundary_conditions_next_timestep(p_boundary_res[it_pipe],v_boundary_tur[it_pipe])
                    # save the second set of boundary conditions
                p_boundary_tur[it_pipe] = pipe.get_current_pressure_distribution()[-1]
                v_boundary_res[it_pipe] = pipe.get_current_velocity_distribution()[0]
                Q_boundary_res[it_pipe] = pipe.get_current_flux_distribution()[0]

                # perform the next timestep via the characteristic method
                    # use vectorized method for performance
                pipe.timestep_characteristic_method_vectorized()

                # prepare for next loop
                p_old = pipe.get_current_pressure_distribution()
                v_old = pipe.get_current_velocity_distribution()
                Q_old = pipe.get_current_flux_distribution()


            # code cell 10
            # code for plotting and safing the figures generated in the loop

            level_plot_min = 0
            level_plot_max = 3
            volume_plot_min = level_plot_min*Res_area_base
            volume_plot_max = level_plot_max*Res_area_base

            fig3,axs3 = plt.subplots(2,2,figsize=(16,9))
            fig3.suptitle('Fläche = '+str(Res_area_base)+'\n'+'Kp = '+str(Con_K_p)+'  Ti = '+str(Con_T_i))
            axs3[0,0].set_title('Level and Volume reservoir')
            axs3[0,0].plot(t_vec,level_vec,label='level')
            axs3[0,0].plot(t_vec,np.full_like(t_vec,Res_level_crit_lo),label='level_limit',c='r')
            axs3[0,0].set_xlabel(r'$t$ [$\mathrm{s}$]')
            axs3[0,0].set_ylabel(r'$h$ [m]')
            axs3[0,0].set_ylim(level_plot_min,level_plot_max)
            x_twin_00 = axs3[0,0].twinx()
            x_twin_00.set_ylabel(r'$V$ [$\mathrm{m}^3$]')
            x_twin_00.plot(t_vec,volume_vec)
            x_twin_00.set_ylim(volume_plot_min,volume_plot_max)
            axs3[0,0].legend()

            axs3[0,1].set_title('LA')
            axs3[0,1].plot(t_vec,100*OL_T1_LA_soll_vec,label='OL_T1 Target',c='b')
            axs3[0,1].scatter(t_vec[::200],100*OL_T1_LA_ist_vec[::200],label='OL_T1 Actual',c='b',marker='+')
            axs3[0,1].plot(t_vec,100*UL_T1_LA_soll_vec,label='UL_T1 Target',c='r')
            axs3[0,1].scatter(t_vec[::200],100*UL_T1_LA_ist_vec[::200],label='UL_T1 Actual',c='r',marker='+')
            axs3[0,1].set_xlabel(r'$t$ [$\mathrm{s}$]')
            axs3[0,1].set_ylabel(r'$LA$ [%]')
            axs3[0,1].legend()

            axs3[1,0].set_title('Fluxes')
            axs3[1,0].plot(t_vec,Q_in_vec,label='Influx')
            axs3[1,0].plot(t_vec,Q_boundary_res,label='Outflux')
            axs3[1,0].scatter(t_vec[::200],Q_boundary_tur[::200],label='Flux Turbine',c='g',marker='+')
            axs3[1,0].set_xlabel(r'$t$ [$\mathrm{s}$]')
            axs3[1,0].set_ylabel(r'$Q$ [$\mathrm{m}^3/\mathrm{s}$]')
            axs3[1,0].legend()

            axs3[1,1].set_title('Pressure change vs t=0 at reservoir and turbine')
            axs3[1,1].plot(t_vec,pressure_conversion(p_boundary_res-p_boundary_res[0],pUnit_calc, pUnit_conv),label='Reservoir')
            axs3[1,1].plot(t_vec,pressure_conversion(p_boundary_tur-p_boundary_tur[0],pUnit_calc, pUnit_conv),label='Turbine')
            axs3[1,1].set_xlabel(r'$t$ [$\mathrm{s}$]')
            axs3[1,1].set_ylabel(r'$p$ ['+pUnit_conv+']')
            axs3[1,1].legend()

            fig3.tight_layout()
            plt.close()

            figname = 'Simulation Lamnitz\KW_Lamnitz_Fläche_'+str(Res_area_base)+'_Kp_'+str(round(Con_K_p,1))+'_Ti_'+str(Con_T_i)+'.png'
            fig3.savefig(figname)


