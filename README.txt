
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MATLAB Code for simulating the COVID-19 University epidemic evolution with or without the optimal control scenarios
in the work

Ranking the Effectiveness of Non-pharmaceutical Interventions to Counter COVID-19 in Universities with Vaccinated Population
by Zirui Niu and Giordano Scarciotti

Zirui Niu, October 25, 2021
Contact: z.niu@outlook.com
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


This folder has two subfolders which contain MATLAB scripts for reproducing the results in our article. More specifically:


Subfoler 1: UncontrolledModel_Simulation
This folder contains MATLAB scripts that can be used to reproduce the results in the baseline scenario when no epidemic mitigation measures are conducted in the University. It 
includes several scripts:

	*MAIN_UncontrolledModelSimulation.m
	The MATLAB script is used to simulate the epidemic evolution in the University with vaccinated population when no countermeasures are conducted. 
	The model prameter values are set based on various research into COVID-19 spread and vaccine efficacy. Number of students and staff are set according to the 
	EEE department of Imperial College London. These parameter values and initial conditions can all be tuned to adapt to any other actual case.
	
	*Other MATLAB scripts
	The other scripts are supporting functions for constructing the system dynamics and plotting the un-controlled epidemic evolution as well as R0 variations.


Subfoler 2: OptimalControl_FourScenarios
This folder contains MATLAB scripts that can be used to reproduce the results of the four designated epidemic mitigation scenarios in the University. Five intervention measures
can be implemented by the University after re-opening: mask wearing, social distancing, environmental disinfection, quarantine on infected students, and quarantine of infected staff.
Their effectiveness are considered and included as input variables. Considering the efforts made by different measures, four scenarios are constructed with four different objectives of 
optimal control. Constraints are designed to require the epidemic must end within 120 days and thatat least 95% of students and staff are not infected. Results can be obtained by 
solving the optimal control problem with the help of the library function "OptimTraj" designed by [1]. For this part, it includes several scripts:

	*MAIN_MinimumCase.m
	The MATLAB script is used to solve the optimal control problem when the University spares no effort to extinguish the COVID-19 pandemic. Optimal trajectories of 
	both the resulting epidemic evolution and the expected effectiveness of five countermeasures are plotted by this script. 

	*MAIN_MinimumControlNorm.m
	The MATLAB script is used to solve the optimal control problem when the University aims to minimise the useof all interventions (including quaranting). Optimal
	trajectories of both the resulting epidemic evolution and the expected effectiveness of five countermeasures are plotted by this script. 
	
	*MAIN_MinimumNonQuarantine.m
	The MATLAB script is used to solve the optimal control problem when we want to minimize  the  use  of  masks,  social  distancing  and  environmental  disinfection.
	Optimal trajectories of both the resulting epidemic evolution and the expected effectiveness of five countermeasures are plotted by this script. 
	
	*MAIN_MinimumQuarantine.m
	The MATLAB script is used to solve the optimal control problem when the University plans to minimise the use of quarantine but allow a strong use of mask wearing, 
	social distancing andenvironmental  disinfection. Optimal trajectories of both the resulting epidemic evolution and the expected effectiveness of five countermeasures 
	are plotted by this script. 
	
	*Other MATLAB scripts
	The other scripts are supporting functions for constructing the system dynamics, defining ojective function of optimal control problem, and plotting the un-controlled
	epidemic evolution as well as R0 variations.

	***The "OptimTraj-master" folder is the MATLAB library for constructing and solving optimal control problems [1].


******Reference*****
[1] M. Kelly, OptimTraj - Trajectory Optimization for Matlab, (2017), GitHub repository, https://github.com/MatthewPeterKelly/OptimTraj