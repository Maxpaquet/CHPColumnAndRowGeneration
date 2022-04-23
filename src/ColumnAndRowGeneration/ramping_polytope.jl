using CSV;
using DataFrames;
using JuMP;
using Gurobi;
using Plots;
using GLPK;
using TickTock;
using PyPlot;
using LaTeXStrings
using JSON;
using StatsPlots;

include("src_test_function.jl");
include("src.jl");

# const GUROBI_ENV = Gurobi.Env()  # this is to avoid having msgs from Gurobi all the time

#########################################################################################################
#########################################################################################################
############################################ RUNNING METHODS ############################################
#########################################################################################################
#########################################################################################################

function test_posterior_data()
	nb_gen = 1;
	MinRunCapacity 		= [6 10 10 4 4 2];
	MaxRunCapacity 		= [16 16 16 7 7 6];
	RampUp				= [5 5 5 1 1 2];
	RampDown 			= [5 5 5 1 1 2];;
	UT 					= [1 1 1 1 1 1];
	DT		 			= [1 1 1 1 1 1];
	NoLoadConsumption 	= [10 0 0 0 0 0];
	F 					= [53 53 53 30 30 0];
	C 					= [30 3 3 2 2 7];
	SU = MinRunCapacity;
	SD = MinRunCapacity;
	# Demand/Load
	VOLL = 3000;
	# data_demand = [6 11 16 11 13 15]; # 6 10 16 9 13 15 -- 6 10 16 11 13 15 -- 6 11 16 11 13 15
	# data_demand = [6 11 16 11];
	data_demand = [20 20 20];
	T_max = length(data_demand);
	
	### Optimization programs
	u_prior 		= zeros(nb_gen);
	v_prior 		= zeros(nb_gen);
	w_prior 		= zeros(nb_gen);
	p_prior 		= zeros(nb_gen);
	pbar_prior 		= zeros(nb_gen);

	u_posterior 	= zeros(nb_gen);
	v_posterior 	= zeros(nb_gen);
	w_posterior 	= zeros(nb_gen);
	p_posterior 	= zeros(nb_gen);
	pbar_posterior 	= zeros(nb_gen);
	#=u_posterior 	= [];
	v_posterior 	= [];
	w_posterior 	= [];
	p_posterior 	= [];
	pbar_posterior 	= [];=#

	(p_EF, pbar_EF,gamma_EF, p_time_EF, pbar_time_EF, u_EF, v_EF, w_EF, price_EF, l_EF, obj_EF, Solving_time_EF) = Extended_Formulation(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
	println("p_time_EF $(p_time_EF)");
	println("u_EF : $(u_EF)");
	println("v_EF : $(v_EF)");
	println("w_EF : $(w_EF)");
	
	Printer = true;
	G_c = Int64[g for g=1:nb_gen]
	# gather_gen, indices_gen = Gather_Generators(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, G_c);
	# gather_gen, indices_gen = Gather_Generators(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, MinRunCapacity, MinRunCapacity, G_c);
	gather_gen = Int64[i for i=1:nb_gen]; indices_gen = Int64[i for i=1:nb_gen];
	delta_criterion = 0;
	(iter_max, price_BD, obj_vec_BD, p_BD, pbar_BD, u_BD, v_BD, w_BD, l_BD, obj_real_BD, Solving_time, nb_cuts_tot, nb_salves_solve, prices_vect, p_vect) = Bender_Decomposition(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, gather_gen, indices_gen, G_c, delta_criterion, Printer);
	(p_time_restricted, pbar_time_restricted, u_restricted, v_restricted, w_restricted, price_restricted, l_restricted, obj_restricted, obj_pricing, Solving_time_CG, obj_restricted_vector, obj_pricing_vector, prices_vect_cr) = Column_And_Row_Generation(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer);
	println();println();
	println("price_BD : $(price_BD)");
	println("price_restricted : $(price_restricted)");
	println();
	println("obj_real_BD    : $(obj_real_BD)");
	println("obj_restricted : $(obj_restricted)");
	println("obj_pricing    : $(obj_pricing)");
end





function Check_Column_And_Row_Generation()
	###########################################################################################
	###########################################################################################
	#=MinRunCapacity = [10 10];
	MaxRunCapacity = [20 20];
	RampUp = [5 1]; 				# [MW/time period]
	RampDown = [6 4];					# [MW/time period] RD = [20 20];
	UT = [1 1];					# [time period]
	DT = [1 1];					# [time period]
	SU = [12 10]; # Doit être au moins aussi grand que MinRunCapacity
	SD = [12 15]; # Doit être au moins aussi grand que MinRunCapacity
	StartupCost = [5 5];				# [EUR]
	MargCost = [10 10];				# [EUR/MW]
	NoLoadConsumption = [1 1];
	# data_demand = [15 15 15 20];
	data_demand = [5 16 89 21 13];
	VOLL = 3000;
	T_max = length(data_demand);
	nb_gen = 2;=#
	
	###########################################################################################
	###########################################################################################
	#=MinRunCapacity = [10 10];
	MaxRunCapacity = [20 20];
	RampUp = [20 1]; 				# [MW/time period]
	RampDown = [20 4];					# [MW/time period] RD = [20 20];
	UT = [1 1];					# [time period]
	DT = [1 1];					# [time period]
	SU = [20 10]; # Doit être au moins aussi grand que MinRunCapacity
	SD = [20 15]; # Doit être au moins aussi grand que MinRunCapacity
	StartupCost = [5 5];				# [EUR]
	MargCost = [10 10];				# [EUR/MW]
	NoLoadConsumption = [1 1];
	data_demand = [15 15 15 20];
	# data_demand = [15]
	VOLL = 3000;
	T_max = length(data_demand);
	nb_gen = 1;=#

	###########################################################################################
	###########################################################################################
	#=
	df1 = CSV.read("D:/Documents/Mémoire Papavasiliou/Unit commitment Toy/UCDataBE/Demand.csv", DataFrame);
	Demand = convert(Matrix, df1);
	df2 = CSV.read("D:/Documents/Mémoire Papavasiliou/Unit commitment Toy/UCDataBE/Generators.csv", DataFrame);
	Generators = convert(Matrix, df2);
	######## Demand ########
	Data_types = Demand[:,1];
	(a,b) = size(Demand)
	Demand = convert(Matrix{Float64}, Demand[:,2:b])
	AutumnWE = Demand[5,:];
	WinterWE = Demand[8,:];
	SpringWE = Demand[6,:];
	SummerWE = Demand[7,:];
	AutumnWD = Demand[1,:];
	WinterWD = Demand[4,:];
	SpringWD = Demand[2,:];
	SummerWD = Demand[3,:];

	######## Generators ########
	nb_gen = 68; # There are normally 68 generators

	(a,b) = size(Generators);
	dt = 1; # [hour]
	#dt = 0.25;
	Generators = convert(Matrix{Float64}, Generators[1:nb_gen,2:b]);
	MinRunCapacity 		= Generators[:,1];								# [MW]
	MaxRunCapacity 		= Generators[:,2];								# [MW]
	RampUp				= (60*dt)*Generators[:,3];						# [MW/min] must multiply by 15 to have [MW/15-min]
	RampDown 			= (60*dt)*Generators[:,4];						# [MW/min] must multiply by 15
	UT_1 				= convert(Array{Int64}, (1/dt)*Generators[:,5]);# [15 min]
	DT_1		 		= convert(Array{Int64}, (1/dt)*Generators[:,6]);# [15 min]
	NoLoadConsumption 	= Generators[:,7];								# [MW], NoLoadCost = NoLoadConsumption * MargCost
	StartupCost 		= Generators[:,8];								# [€]
	MargCost 			= dt*Generators[:,9];							# [€/MWh]
	SU = MinRunCapacity;
	SD = MinRunCapacity;
	
	VOLL = 3000;

	# T_max = 35; # There normally 96 time-steps  25
	# Demand in hours
	data_demand_vec = WinterWE; # WinterWE --> not working with Bender_Decomposition T_max = 56 and timestep of 15 minuts
	data_demand = zeros(24);
	c = 1;
	for i=1:24
		data_demand[i] = data_demand_vec[c];
		c+=4;
	end
	T_max = 24;

	# T_max = 35;
	# T_max = 80;
	#=T_max = 56;
	data_demand_vec = AutumnWE;
	data_demand = data_demand_vec[1:T_max];=#

	UT = check_UT_DT(UT_1, T_max);
	DT = check_UT_DT(DT_1, T_max);
	=#
	###########################################################################################
	###########################################################################################
	#=
	# Generators
	nb_gen = 6;
	MinRunCapacity 		= [10 10 10 4 4 2];
	MaxRunCapacity 		= [16 16 16 7 7 6];
	RampUp				= [5 5 5 1 1 2];
	RampDown 			= [5 5 5 1 1 2];;
	UT 					= [1 1 1 1 1 1];
	DT		 			= [1 1 1 1 1 1];
	NoLoadConsumption 	= [0 0 0 0 0 0];
	StartupCost 		= [52 51 50 30 30 0];
	MargCost 			= [3 3 3 2 2 7];
	SU = MinRunCapacity;
	SD = MinRunCapacity;
	# Demand/Load
	VOLL = 3000;
	data_demand = [30 20 33 40];
	# data_demand = [30 40];
	T_max = length(data_demand);
	=#
	###########################################################################################
	###########################################################################################
	nb_gen = 1;
	MinRunCapacity 		= [6];
	MaxRunCapacity 		= [16];
	RampUp				= [5];
	RampDown 			= [5];
	UT 					= [1];
	DT			 		= [1];
	NoLoadConsumption 	= [10];
	StartupCost 		= [53];
	MargCost 			= [30];
	SU = [6];
	SD = [6];
	data_demand = [6 11 16 11];
	T_max = length(data_demand);
	################################################################################################
	################################### PRIOR AND POSTERIOR DATA ###################################
	################################################################################################
	VOLL = 3000;
	u_prior 		= zeros(nb_gen);
	v_prior 		= zeros(nb_gen);
	w_prior 		= zeros(nb_gen);
	p_prior 		= zeros(nb_gen);
	pbar_prior 		= zeros(nb_gen);

	u_posterior 	= zeros(nb_gen);
	v_posterior 	= zeros(nb_gen);
	w_posterior 	= zeros(nb_gen);
	p_posterior 	= zeros(nb_gen);
	pbar_posterior 	= zeros(nb_gen);

	#(p_matching, pbar_matching, u_matching, v_matching, w_matching, l_matching, obj_matching, Solving_time_matching) = Matching(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, MinRunCapacity, MinRunCapacity, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
	#println("obj_matching : $(obj_matching)");

	Printer = false;

	#=G_c = Int64[i for i=1:nb_gen];
	# gather_gen, indices_gen = Gather_Generators(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, MinRunCapacity, MinRunCapacity, G_c);
	gather_gen = Int64[i for i=1:nb_gen]; indices_gen = Int64[i for i=1:nb_gen];
	delta_criterion = 0;
	started_time = time();
	(iter_max, price_BD, obj_vec_BD, p_BD, pbar_BD, u_BD, v_BD, w_BD, l_BD, obj_real_BD, Solving_time_BD, nb_cuts_tot, nb_salves_solve) = Bender_Decomposition(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, gather_gen, indices_gen, G_c, delta_criterion, Printer);
	elapsed = time()-started_time;
	println("Bender Decomposition time elapsed : $(elapsed)");
	println("obj_real_BD : ",obj_real_BD);

	started_time = time();
	(iter_max, price_opt, obj_vec, p_time_opt, pbar_time_opt, u_opt, v_opt, w_opt, l_BD, obj_real_BD, Solving_time_BD_1) = Bender_Decomposition_1(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, delta_criterion, false);
	elapsed = time()-started_time;
	println("Bender Decomposition_1 elapsed : $(elapsed)");
	println("obj_real_BD : ",obj_real_BD);
	# println("Solving time BD : ",Solving_time_BD);println();

	println();
	started_time = time();
	(p_relax, pbar_relax, u_relax, v_relax, w_relax, price_UB, l_relax, obj_real_relax, Solving_time_LP) = LP_Relaxation(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
	elapsed = time()-started_time;
	println("Linear Relaxation time elapsed : $(elapsed)");
	println("obj_real_relax : ",obj_real_relax);
	
	println();
	started_time = time();
	(p_time_restricted_1, pbar_time_restricted_1, u_restricted_1, v_restricted_1, w_restricted_1, price_restricted_1, l_restricted_1, obj_restricted_1, obj_pricing_1, Solving_time_CG_1, obj_restricted_vec_1, obj_pricing_vector_1) = Column_And_Row_Generation_1(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer);
	elapsed = time()-started_time;
	println("Column-and-row generation algorithm (regular one)  : $(elapsed)");
	println("Solving Solving time (regular one) : ",Solving_time_CG_1);
	println("obj_restricted_1 : ",obj_restricted_1);
	println("obj_pricing_1    : ",obj_pricing_1);
	println();=#

	println();
	started_time = time();
	(p_time_restricted, pbar_time_restricted, u_restricted, v_restricted, w_restricted, price_restricted, l_restricted, obj_restricted, obj_pricing, Solving_time_CG, obj_restricted_vec, obj_pricing_vector, price_vect, iter_stop, Solving_time_master, Solving_time_slave, nb_var_vec, nb_cons_vec, Solving_time_master_vec, Solving_time_slave_vec) = Column_And_Row_Generation(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer);
	elapsed = time()-started_time;
	println("Column-and-row generation algorithm (hot-start model) elapsed : $(elapsed)");
	println("Solving Solving time (hot-start model) : ",Solving_time_CG);
	println("iter_stop : $(iter_stop)");
	println("obj_restricted : ",obj_restricted);
	println("obj_pricing    : ",obj_pricing);
	println("nb_var_vec = $(nb_var_vec);")
	println("nb_cons_vec = $(nb_cons_vec);");
	println("Solving_time_master_vec = $(Solving_time_master_vec);");
	println("Solving_time_slave_vec = $(Solving_time_slave_vec);");

	println();
	delta_criterion = -10^(-5);
	Started_time = time();
	(obj_column_generation, price_column_generation, iter_stop_column_generation, Solving_time, Solving_time_master, Solving_time_slaves, obj_vect, prices_iterates, nb_var_vec ,nb_cons_vec, Solving_time_master_vec, Solving_time_slave_vec) = Column_Generation_DW(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, delta_criterion,Printer)
	elapsed = time() - Started_time;
	println("Column generation algorithm elapsed : $(elapsed)");
	println("Solving_time : $(Solving_time)");
	println("iter_stop_column_generation : $(iter_stop_column_generation)");
	println("nb_var_vec = $(nb_var_vec);");
	println("nb_cons_vec = $(nb_cons_vec);");
	println("Solving_time_master_vec = $(Solving_time_master_vec);");
	println("Solving_time_slave_vec = $(Solving_time_slave_vec);");

	#=println();
	Started_time = time();
	(p_EF, pbar_EF, gamma_EF, p_time_EF, pbar_time_EF, u_EF, v_EF, w_EF, price_EF, l_EF, obj_real_EF, Solving_time_EF) = Extended_Formulation(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
	elasped = time() - Started_time;
	println("Extended Formulation time elasped : $(elasped)");
	println("obj_real_EF : ",obj_real_EF);
	println();=#

	#=
	#phi_criteria_vec = 0:10:100;
	phi_criteria_vec = 0:0.000001:0.00001;
	(p_matching, pbar_matching, u_matching, v_matching, w_matching, l_matching, obj_real_matching, Solving_time_matching) = Matching(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, false);
	# phi_criteria_vec = 0:10:20;
	elapsed_vec = zeros(length(phi_criteria_vec));
	Solving_time_CG_vec = zeros(length(phi_criteria_vec));
	iter_stop_vec = zeros(length(phi_criteria_vec));
	obj_restricted_vec = zeros(length(phi_criteria_vec));
	uplifts_vec = zeros(length(phi_criteria_vec));
	
	for (i,phi_criteria) in enumerate(phi_criteria_vec)
		started_time = time();
		(p_time_restricted, pbar_time_restricted, u_restricted, v_restricted, w_restricted, price_restricted, l_restricted, obj_restricted, obj_pricing, Solving_time_CG, obj_restricted_vector, obj_pricing_vector, price_vect, iter_stop, Solving_time_master, Solving_time_slave) = Column_And_Row_Generation(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer,phi_criteria);
		elapsed = time()-started_time;
		elapsed_vec[i] = elapsed;
		Solving_time_CG_vec[i] = Solving_time_CG;
		iter_stop_vec[i] = iter_stop;
		obj_restricted_vec[i] = obj_restricted;
		uplifts_vec[i] = Compute_Uplifts_Illustrative_Example(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price_restricted, p_matching, u_matching, v_matching, l_matching);
	end
	println("delta_vec = $(phi_criteria_vec);");
	println("elapsed_vec = $(elapsed_vec);");
	println("Solving_time_vec = $(Solving_time_CG_vec);");
	println("iteration_vec = $(iter_stop_vec);");
	println("obj_master_vec = $(obj_restricted_vec);");
	println("uplifts_vec = $(uplifts_vec);")
	=#
end

#Check_Column_And_Row_Generation();


function Run_Bender_Decomposition()
	df1 = CSV.read("D:/Documents/Mémoire Papavasiliou/Unit commitment Toy/UCDataBE/Demand.csv", DataFrame);
	Demand = convert(Matrix, df1);
	df2 = CSV.read("D:/Documents/Mémoire Papavasiliou/Unit commitment Toy/UCDataBE/Generators.csv", DataFrame);
	Generators = convert(Matrix, df2);
	######## Demand ########
	Data_types = Demand[:,1];
	(a,b) = size(Demand)
	Demand = convert(Matrix{Float64}, Demand[:,2:b])
	AutumnWE = Demand[5,:];
	WinterWE = Demand[8,:];
	SpringWE = Demand[6,:];
	SummerWE = Demand[7,:];
	AutumnWD = Demand[1,:];
	WinterWD = Demand[4,:];
	SpringWD = Demand[2,:];
	SummerWD = Demand[3,:];

	######## Generators ########
	nb_gen = 68; # There are normally 68 generators

	(a,b) = size(Generators);
	dt = 1; # [hour]
	# dt = 0.25;
	Generators = convert(Matrix{Float64}, Generators[1:nb_gen,2:b]);
	MinRunCapacity 		= Generators[:,1];								# [MW]
	MaxRunCapacity 		= Generators[:,2];								# [MW]
	RampUp				= (60*dt)*Generators[:,3];						# [MW/min] must multiply by 15 to have [MW/15-min]
	RampDown 			= (60*dt)*Generators[:,4];						# [MW/min] must multiply by 15
	UT_1 				= convert(Array{Int64}, (1/dt)*Generators[:,5]);# [15 min]
	DT_1		 		= convert(Array{Int64}, (1/dt)*Generators[:,6]);# [15 min]
	NoLoadConsumption 	= Generators[:,7];								# [MW], NoLoadCost = NoLoadConsumption * MargCost
	StartupCost 		= Generators[:,8];								# [€]
	MargCost 			= dt*Generators[:,9];							# [€/MWh]
	SU = MinRunCapacity;
	SD = MinRunCapacity;
	
	VOLL = 3000;

	# T_max = 35; # There normally 96 time-steps  25
	# Demand in hours
	data_demand_vec = WinterWE; # WinterWE --> not working with Bender_Decomposition T_max = 56 and timestep of 15 minuts
	data_demand = zeros(24);
	c = 1;
	for i=1:24
		data_demand[i] = data_demand_vec[c];
		c+=4;
	end
	T_max = 24;

	UT = check_UT_DT(UT_1, T_max);
	DT = check_UT_DT(DT_1, T_max);

	u_prior 		= zeros(nb_gen);
	v_prior 		= zeros(nb_gen);
	w_prior 		= zeros(nb_gen);
	p_prior 		= zeros(nb_gen);
	pbar_prior 		= zeros(nb_gen);
	u_posterior 	= zeros(nb_gen);
	v_posterior 	= zeros(nb_gen);
	w_posterior 	= zeros(nb_gen);
	p_posterior 	= zeros(nb_gen);
	pbar_posterior 	= zeros(nb_gen);

	G_c = Int64[i for i=1:nb_gen];gather_gen = Int64[i for i=1:nb_gen]; indices_gen = Int64[i for i=1:nb_gen];
	
	#=delta_criterion = 0;
	started_time = time();
	(iter_max, price_BD, obj_vec_BD, p_BD, pbar_BD, u_BD, v_BD, w_BD, l_BD, obj_BD, Solving_time_BD, nb_cuts_tot, nb_salves_solve) = Bender_Decomposition(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, gather_gen, indices_gen, G_c, delta_criterion, false);
	elapsed = time() - started_time;
	println("Bender Decomposition time - Total running time : $(elapsed)");
	println("Solving_time_BD : $(Solving_time_BD)");
	println("obj_BD : $(obj_BD)");
	println("price_BD : $(price_BD)");
	println();println();
	started_time = time();
	(iter_max, price_opt, obj_vec, p_time_opt, pbar_time_opt, u_opt, v_opt, w_opt, l_BD, obj_real_BD, Solving_time_BD_1) = Bender_Decomposition_1(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
	elapsed = time() - started_time;
	println("Bender Decomposition time - Total running time : $(elapsed)");
	println("Solving_time_BD : $(Solving_time_BD_1)");
	println("obj_real_BD : $(obj_real_BD)");=#
	
	println();println();
	(p_matching, pbar_matching, u_matching, v_matching, w_matching, l_matching, obj_real_matching, Solving_time_matching) = Matching(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, true);
	println("obj_real_matching : $(obj_real_matching)");
	# println("uplifts (BD) = $(obj_real_matching - obj_real_BD)");
	
	delta_criterion_vec = 0:0.1:1;
	objective_vec = zeros(length(delta_criterion_vec));
	uplifts_vec = zeros(length(delta_criterion_vec));
	iteration_vec = zeros(length(delta_criterion_vec));
	elapsed_vec = zeros(length(delta_criterion_vec));
	Solving_time_vec = zeros(length(delta_criterion_vec));
	for (i,delta_criterion) in enumerate(delta_criterion_vec)
		println();
		started_time = time();
		(iter_max, price_BD, obj_vec_BD, p_BD, pbar_BD, u_BD, v_BD, w_BD, l_BD, obj_real_BD, Solving_time_BD, nb_cuts_tot, nb_salves_solve) = Bender_Decomposition(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, gather_gen, indices_gen, G_c, delta_criterion, false);
		elapsed = time() - started_time;
		println("Bender Decomposition time - Total running time : $(elapsed)")
		# println("price_BD : $(price_BD)");
		objective_vec[i] = obj_real_BD;
		uplifts_vec[i] = obj_real_matching - obj_real_BD;
		iteration_vec[i] = iter_max;
		elapsed_vec[i] = elapsed;
		Solving_time_vec[i] = Solving_time_BD;
	end
	println("delta_criterion_vec = $(delta_criterion_vec);");
	println("objective_vec = $(objective_vec);");
	println("uplifts_vec = $(uplifts_vec);");
	println("elapsed_vec = $(elapsed_vec);");
	println("Solving_time_vec = $(Solving_time_vec);");
	println("iteration_vec = $(iteration_vec);");
	
end

# Run_Bender_Decomposition()


function Run_Extended_Formulation()
	#=df1 = CSV.read("D:/Documents/Mémoire Papavasiliou/Unit commitment Toy/UCDataBE/Demand.csv", DataFrame);
	Demand = convert(Matrix, df1);
	df2 = CSV.read("D:/Documents/Mémoire Papavasiliou/Unit commitment Toy/UCDataBE/Generators.csv", DataFrame);
	Generators = convert(Matrix, df2);
	######## Demand ########
	Data_types = Demand[:,1];
	(a,b) = size(Demand)
	Demand = convert(Matrix{Float64}, Demand[:,2:b])
	AutumnWD = Demand[1,:];
	SpringWD = Demand[2,:];
	SummerWD = Demand[3,:];
	WinterWD = Demand[4,:];
	AutumnWE = Demand[5,:];
	SpringWE = Demand[6,:];
	SummerWE = Demand[7,:];
	WinterWE = Demand[8,:];

	######## Generators ########
	nb_gen = 68; # There are normally 68 generators

	(a,b) = size(Generators);
	# dt = 1; # [hour]
	dt = 0.25;
	Generators = convert(Matrix{Float64}, Generators[1:nb_gen,2:b]);
	MinRunCapacity 		= Generators[:,1];								# [MW]
	MaxRunCapacity 		= 4*Generators[:,2];							# [MW]  10*
	RampUp				= (60*dt)*Generators[:,3];						# [MW/min] must multiply by 15 to have [MW/15-min]
	RampDown 			= (60*dt)*Generators[:,4];						# [MW/min] must multiply by 15
	UT_1 				= convert(Array{Int64}, (1/dt)*Generators[:,5]);# [15 min]
	DT_1		 		= convert(Array{Int64}, (1/dt)*Generators[:,6]);# [15 min]
	NoLoadConsumption 	= Generators[:,7];								# [MW], NoLoadCost = NoLoadConsumption * MargCost
	StartupCost 		= Generators[:,8];								# [€]
	MargCost 			= dt*Generators[:,9];							# [€/MWh]
	SU = MinRunCapacity;
	SD = MinRunCapacity;
	VOLL = 3000;
	T_max = 35; # There normally 96 time-steps  25
	data_demand_vec = WinterWE;
	data_demand = data_demand_vec[1:T_max];
	# Demand in hours
	#=data_demand_vec = WinterWE;
	data_demand = zeros(24);
	c = 1;
	for i=1:24
		data_demand[i] = data_demand_vec[c];
		c+=4;
	end
	VOLL = 3000;
	T_max = 24;=#
	=#
	########################################################################################
	########################################################################################
	#=MinRunCapacity = [10 12 20];
	MaxRunCapacity = [20 40 30];
	RampUp = [5 10 10]; 				# [MW/time period]
	RampDown = [5 15 10];					# [MW/time period] RD = [20 20];
	UT_1 = [1 1 2];					# [time period]
	DT_1 = [1 1 4];					# [time period]
	SU = [10 12 20]; # Doit être au moins aussi grand que MinRunCapacity
	SD = [10 12 20]; # Doit être au moins aussi grand que MinRunCapacity
	F = [70 75 5];				# [EUR]
	C = [50 56 2];				# [EUR/MW]
	NoLoadConsumption = [1 2 6];

	# data_demand = 2*[3.21029 13.3892 17.5783 18.5449 5.62669 16.5953 0.869666 11.8104 3.28692 17.9711]; ## Example 5
	# data_demand =  2*[3.21029 13.3892 17.5783 18.5449];
	# data_demand =  2*[3.21029 13.3892 17.5783 18.5449 5.62669 16.5953];
	data_demand = [10]
	VOLL = 3000;
	T_max = length(data_demand);
	nb_gen = 2; # 2=#
	########################################################################################
	########################################################################################
	nb_gen = 1;
	MinRunCapacity 		= [6];
	MaxRunCapacity 		= [16];
	RampUp				= [5];
	RampDown 			= [5];
	UT_1 				= [1];
	DT_1		 		= [1];
	NoLoadConsumption 	= [10];
	StartupCost 		= [53];
	MargCost 			= [30];
	SU = [6];
	SD = [6];
	data_demand = [6 11 16 11];
	T_max = length(data_demand);
	########################################################################################
	########################################################################################

	UT = check_UT_DT(UT_1, T_max);
	DT = check_UT_DT(DT_1, T_max);
	VOLL = 3000;

	u_prior 		= zeros(nb_gen);
	v_prior 		= zeros(nb_gen);
	w_prior 		= zeros(nb_gen);
	p_prior 		= zeros(nb_gen);
	pbar_prior 		= zeros(nb_gen);
	u_posterior 	= zeros(nb_gen);
	v_posterior 	= zeros(nb_gen);
	w_posterior 	= zeros(nb_gen);
	p_posterior 	= zeros(nb_gen);
	pbar_posterior 	= zeros(nb_gen);
	started_time = time();
	(_, _, gamma_EF, p_EF, pbar_EF, u_EF, v_EF, w_EF, price, l_val, obj_EF, Solving_time) = Extended_Formulation(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, true);
	elapsed = time()-started_time;
	println("obj_EF : $(obj_EF)");
	println("elapsed : $(elapsed)");
	println("price : $(price)");

	println();
	(p_relax, pbar_relax, u_relax, v_relax, w_relax, price_UB, l_relax, obj_real_relax, Solving_time_LP) = LP_Relaxation(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, false);
	println("obj_real_relax : $(obj_real_relax)");
	println("price_UB : $(price_UB)");
	println();
	println();
	G_c = Int64[i for i=1:nb_gen];
	gather_gen = Int64[i for i=1:nb_gen];
	indices_gen = Int64[i for i=1:nb_gen];
	delta_criterion = 0;
	started_time = time();
	(iter_max, price_BD, obj_vec_BD, p_BD, pbar_BD, u_BD, v_BD, w_BD, l_BD, obj_real_BD, Solving_time_BD, nb_cuts_tot, nb_salves_solve) = Bender_Decomposition(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, gather_gen, indices_gen, G_c, delta_criterion, false);
	elapsed = time() - started_time;
	println("obj_real_BD : $(obj_real_BD)");
	println("Bender Decomposition time - Total running time : $(elapsed)");
	println("Solving_time_BD : $(Solving_time_BD)");
	println("price_BD : $(price_BD)");

	(p_matching, pbar_matching, u_matching, v_matching, w_matching, l_matching, obj_real_matching, Solving_time_matching) = Matching(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, true);
	println("obj_real_matching : $(obj_real_matching)");
	println("uplifts (EF) = obj_real_matching - obj_EF         = $(obj_real_matching-obj_EF)");
	println("uplifts (LP) = obj_real_matching - obj_real_relax = $(obj_real_matching-obj_real_relax)");
end

# Run_Extended_Formulation()


function ComputeUplift(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, delta_criterion, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior)
	### STEP 2
	# Bender Decomposition
	# G_c = Compute_G_c(MinRunCapacity, MaxRunCapacity, RU, RD, UT);
	G_c = Int64[i for i=1:nb_gen];
	# gather_gen, indices_gen = Gather_Generators(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, G_c);
	gather_gen = Int64[i for i=1:nb_gen]; indices_gen = Int64[i for i=1:nb_gen];
	println();
	started_time = time();
	(iter_max, price_BD, obj_vec_BD, p_BD, pbar_BD, u_BD, v_BD, w_BD, l_BD, obj_real_BD, Solving_time_BD, nb_cuts_tot, nb_salves_solve) = Bender_Decomposition(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, gather_gen, indices_gen, G_c, delta_criterion, false);
	elapsed = time() - started_time;
	println("Bender Decomposition time - Total running time : $(elapsed)")
	println("price_BD : $(price_BD)");
	# obj_real_BD = 0;obj_vec_BD = [0]
	#println();println("p_BD : ",p_BD);println();

	Printer = false;
	println();
	started_time = time();
	(p_time_restricted, pbar_time_restricted, u_restricted, v_restricted, w_restricted, price_restricted, l_restricted, obj_restricted, obj_pricing, Solving_time_CG, obj_restricted_vec, obj_pricing_vector) = Column_And_Row_Generation(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer);
	elapsed = time() - started_time;
	println("Column-and-row generation algorithm - Total running time : $(elapsed)");
	println("Solving Solving time (hot-start model) : ",Solving_time_CG);println();
	println("price_restricted : $(price_restricted)");

	#=println();
	println("Column-and-row generation algorithm (regular one)");
	@time (p_time_restricted_1, pbar_time_restricted_1, u_restricted_1, v_restricted_1, w_restricted_1, price_restricted_1, l_restricted_1, obj_restricted_1, obj_pricing_1, Solving_time_CG_1, obj_restricted_vec_1, obj_pricing_vector_1) = Column_And_Row_Generation_1(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer);

	println("Solving Solving time (regular one) : ",Solving_time_CG_1);println();
	#println();println("p_time_restricted : ",p_time_restricted);println();=#

	println();
	started_time = time();
	(obj_column_generation, price_column_generation, iter_stop_column_generation) = Column_Generation_DW(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer)
	elapsed = time() - started_time;
	println("Column generation algorithm - Total running time : $(elapsed)");
	println("iter_stop_column_generation : $(iter_stop_column_generation)");
	println("price_column_generation : $(price_column_generation)");

	# Compute CHP solving the Extended Formulation
	started_time = time();
	(p_EF, pbar_EF, gamma_EF, p_time_EF, pbar_time_EF, u_EF, v_EF, w_EF, price_EF, l_EF, obj_real_EF, Solving_time_EF) = Extended_Formulation(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, true);
	elapsed = time() - started_time;
	println("Extended Formulation time - Total running time : $(elapsed)");
	println("Solving_time_EF : $(Solving_time_EF)");
	println("price_EF : $(price_EF)")
	println("obj_real_EF : $(obj_real_EF)");
	println();


	
	#obj_real_EF = 0

	### STEP 1
	# Solve matching problem : unit commitment problem
	println();
	started_time = time();
	(p_matching, pbar_matching, u_matching, v_matching, w_matching, l_matching, obj_real_matching, Solving_time_matching) = Matching(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, true);
	elapsed = time() - started_time;
	println("Matching time - Total running time : $(elapsed)");
	println("Solving_time_matching : $(Solving_time_matching)");

	println();
	started_time = time();
	(p_relax, pbar_relax, u_relax, v_relax, w_relax, price_UB, l_relax, obj_real_relax, Solving_time_LP) = LP_Relaxation(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, true);
	elapsed = time() - started_time;
	println("Linear Relaxation time - Total running time : $(elapsed)");
	println("Solving_time_LP : $(Solving_time_LP)");
	println("price_UB : $(price_UB)");

	println();
	println("------------------------------------------------------------------------------");
	println("obj_vec_BD : ");
	println(obj_vec_BD);
	println("------------------------------------------------------------------------------");
	println("obj_restricted_vec : ");
	println(obj_restricted_vec);
	println("------------------------------------------------------------------------------");
	println("obj_pricing_vector");
	println(obj_pricing_vector);
	println("------------------------------------------------------------------------------");
	
	### STEP 3
	### Step 3.1 : Compute maximum profit with CHP from EF
	received_profit_producer_EF = zeros(nb_gen); max_profit_producer_EF = zeros(nb_gen); uplift_producer_EF = zeros(nb_gen);
	max_profit_consumer_EF = 0;received_profit_consumer_EF = 0;uplift_consumer_EF = 0;
	for g=1:nb_gen
		(max_profit_EF, p_max_profit, pbar_max_profit, u_max_profit, v_max_profit, w_max_profit) = MaxProfit_Producer(MinRunCapacity[g], MaxRunCapacity[g], RU[g], RD[g], UT[g], DT[g], SU[g], SD[g], T_max, F[g], C[g], NoLoadConsumption[g], price_EF, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, g);
		max_profit_producer_EF[g] = max_profit_EF;

		received_profit_producer_EF[g] = sum(price_EF[t]*p_matching[g,t+1] - (NoLoadConsumption[g]*C[g]*u_matching[g,t+1] + F[g]*v_matching[g,t+1] + C[g]*p_matching[g,t+1]) for t=1:T_max);
		uplift_producer_EF[g] = max_profit_producer_EF[g] - received_profit_producer_EF[g];
	end
	# Computer maximum profit for Consumer with CHP from EF
	(max_profit_consumer, l_max_profit_consumer) = MaxProfit_Consumer(data_demand, T_max, price_EF, VOLL);
	max_profit_consumer_EF = max_profit_consumer;
	received_profit_consumer_EF = sum( (VOLL - price_EF[t])*l_matching[t] for t=1:T_max);
	uplift_consumer_EF = max_profit_consumer_EF - received_profit_consumer_EF;

	### Step 3.2 : Compute maximum profit with CHP from BD
	received_profit_producer_BD = zeros(nb_gen); max_profit_producer_BD = zeros(nb_gen); uplift_producer_BD = zeros(nb_gen);
	max_profit_consumer_BD = 0; received_profit_consumer_BD = 0; uplift_consumer_BD = 0;
	for g=1:nb_gen
		(max_profit_BD, p_max_profit, pbar_max_profit, u_max_profit, v_max_profit, w_max_profit) = MaxProfit_Producer(MinRunCapacity[g], MaxRunCapacity[g], RU[g], RD[g], UT[g], DT[g], SU[g], SD[g], T_max, F[g], C[g], NoLoadConsumption[g], price_BD, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, g);
		max_profit_producer_BD[g] = max_profit_BD;

		received_profit_producer_BD[g] = sum(price_BD[t]*p_matching[g,t+1] -  (NoLoadConsumption[g]*C[g]*u_matching[g,t+1] + F[g]*v_matching[g,t+1] + C[g]*p_matching[g,t+1]) for t=1:T_max);
		uplift_producer_BD[g] = max_profit_producer_BD[g] - received_profit_producer_BD[g];
	end
	# Computer maximum profit for Consumer with CHP from BD
	(max_profit_consumer, l_max_profit_consumer) = MaxProfit_Consumer(data_demand, T_max, price_BD, VOLL);
	max_profit_consumer_BD = max_profit_consumer;
	received_profit_consumer_BD = sum( (VOLL - price_BD[t])*l_matching[t] for t=1:T_max);
	uplift_consumer_BD = (max_profit_consumer_BD - received_profit_consumer_BD);

	### Step 4 : Computing upper bound on uplifts using price from LP relaxation
	received_profit_producer_LP = zeros(nb_gen); max_profit_producer_LP = zeros(nb_gen); uplift_producer_LP = zeros(nb_gen);max_profit_consumer_LP = 0; received_profit_consumer_LP = 0; uplift_consumer_LP = 0;
	for g=1:nb_gen
		(max_profit_LP, p_max_profit, pbar_max_profit, u_max_profit, v_max_profit, w_max_profit) = MaxProfit_Producer(MinRunCapacity[g], MaxRunCapacity[g], RU[g], RD[g], UT[g], DT[g], SU[g], SD[g], T_max, F[g], C[g], NoLoadConsumption[g], price_UB, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, g);
		max_profit_producer_LP[g] = max_profit_LP;

		received_profit_producer_LP[g] = sum(price_UB[t]*p_matching[g,t+1] -  (NoLoadConsumption[g]*C[g]*u_matching[g,t+1] + F[g]*v_matching[g,t+1] + C[g]*p_matching[g,t+1]) for t=1:T_max);
		uplift_producer_LP[g] = max_profit_producer_LP[g] - received_profit_producer_LP[g];
	end
	# Computer maximum profit for Consumer with CHP from LP_Relaxation
	(max_profit_consumer, l_max_profit_consumer) = MaxProfit_Consumer(data_demand, T_max, price_UB, VOLL);
	max_profit_consumer_LP = max_profit_consumer;
	received_profit_consumer_LP = sum( (VOLL - price_UB[t])*l_matching[t] for t=1:T_max);
	uplift_consumer_LP = (max_profit_consumer_LP - received_profit_consumer_LP);


	### Step 5 : Compute maximum profit with CHP from column-and-row gen
	received_profit_producer_CG = zeros(nb_gen); max_profit_producer_CG = zeros(nb_gen); uplift_producer_CG = zeros(nb_gen);max_profit_consumer_CG = 0; received_profit_consumer_CG = 0; uplift_consumer_CG = 0;
	for g=1:nb_gen
		(max_profit_CG, p_max_profit_CG, pbar_max_profit_CG, u_max_profit_CG, v_max_profit_CG, w_max_profit_CG) = MaxProfit_Producer(MinRunCapacity[g], MaxRunCapacity[g], RU[g], RD[g], UT[g], DT[g], SU[g], SD[g], T_max, F[g], C[g], NoLoadConsumption[g], price_restricted, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, g);
		max_profit_producer_CG[g] = max_profit_CG;

		received_profit_producer_CG[g] = sum(price_restricted[t]*p_matching[g,t+1] -  (NoLoadConsumption[g]*C[g]*u_matching[g,t+1] + F[g]*v_matching[g,t+1] + C[g]*p_matching[g,t+1]) for t=1:T_max);
		uplift_producer_CG[g] = max_profit_producer_CG[g] - received_profit_producer_CG[g];
	end
	# Computer maximum profit for Consumer with CHP from LP_Relaxation
	(max_profit_consumer, l_max_profit_consumer) = MaxProfit_Consumer(data_demand, T_max, price_restricted, VOLL);
	max_profit_consumer_CG = max_profit_consumer;
	received_profit_consumer_CG = sum( (VOLL - price_restricted[t])*l_matching[t] for t=1:T_max);
	uplift_consumer_CG = (max_profit_consumer_CG - received_profit_consumer_CG);


	### Step 5 : Compute maximum profit with CHP from column generation
	received_profit_producer_ColG = zeros(nb_gen); max_profit_producer_ColG = zeros(nb_gen); uplift_producer_ColG = zeros(nb_gen);max_profit_consumer_ColG = 0; received_profit_consumer_ColG = 0; uplift_consumer_ColG = 0;
	for g=1:nb_gen
		(max_profit_ColG, _, _, _, _, _) = MaxProfit_Producer(MinRunCapacity[g], MaxRunCapacity[g], RU[g], RD[g], UT[g], DT[g], SU[g], SD[g], T_max, F[g], C[g], NoLoadConsumption[g], price_column_generation, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, g);
		max_profit_producer_ColG[g] = max_profit_ColG;

		received_profit_producer_ColG[g] = sum(price_restricted[t]*p_matching[g,t+1] -  (NoLoadConsumption[g]*C[g]*u_matching[g,t+1] + F[g]*v_matching[g,t+1] + C[g]*p_matching[g,t+1]) for t=1:T_max);
		uplift_producer_ColG[g] = max_profit_producer_ColG[g] - received_profit_producer_ColG[g];
	end
	# Computer maximum profit for Consumer with CHP from LP_Relaxation
	(max_profit_consumer, l_max_profit_consumer) = MaxProfit_Consumer(data_demand, T_max, price_column_generation, VOLL);
	max_profit_consumer_ColG = max_profit_consumer;
	received_profit_consumer_ColG = sum( (VOLL - price_restricted[t])*l_matching[t] for t=1:T_max);
	uplift_consumer_ColG = (max_profit_consumer_ColG - received_profit_consumer_ColG);

	### Sanity checks
	#=obj_sanity_check_MaxProfit_LP = sanity_check_MaxProfit(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price_UB);
	obj_sanity_check_MaxProfit_EF = sanity_check_MaxProfit(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price_EF);
	println();
	println("Test on dual multipliers for Max Profit Producer obj_sanity_check_MaxProfit using price EF : ",obj_sanity_check_MaxProfit_EF);
	println("Test on dual multipliers for Max Profit Producer obj_sanity_check_MaxProfit using price LP : ",obj_sanity_check_MaxProfit_LP);
	println();
	obj_sanity_Extended_Formulation = Sanity_Extended_Formulation(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price_EF);
	obj_sanity_LP_Relaxation = Sanity_LP_Relaxation(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price_UB);
	println();
	println("Test for dual multipliers obj_sanity_Extended_Formulation (should be equal to obj EF) : ",obj_sanity_Extended_Formulation," (",obj_real_EF,")");	
	println("Test for dual multipliers obj_sanity_LP_Relaxation (should be equal to obj LP Relaxation) : ",obj_sanity_LP_Relaxation," (",obj_real_relax,")");
	println();=#

	return (received_profit_producer_EF, max_profit_producer_EF, uplift_producer_EF, received_profit_consumer_EF, max_profit_consumer_EF, uplift_consumer_EF,
		received_profit_producer_BD, max_profit_producer_BD, uplift_producer_BD, received_profit_consumer_BD, max_profit_consumer_BD, uplift_consumer_BD,
		received_profit_producer_LP, max_profit_producer_LP, uplift_producer_LP, received_profit_consumer_LP, max_profit_consumer_LP, uplift_consumer_LP,
		received_profit_producer_CG, max_profit_producer_CG, uplift_producer_CG, received_profit_consumer_CG, max_profit_consumer_CG, uplift_consumer_CG,
		received_profit_producer_ColG, max_profit_producer_ColG, uplift_producer_ColG, received_profit_consumer_ColG, max_profit_consumer_ColG, uplift_consumer_ColG,
		obj_real_matching, obj_real_relax, obj_real_EF, obj_real_BD, obj_restricted, Solving_time_EF, Solving_time_BD, Solving_time_LP, Solving_time_matching,
		nb_cuts_tot, nb_salves_solve, Solving_time_CG, obj_column_generation);
end


### GENERATOR POUR LEQUEL IL Y A UNE FEASIBLE SOLUTION DANS LE DISPATCH POLYTOPE
### i.e. il y a une solution au problème décrit dans la fonction "Extended_Formulation"


function Check_Algorithms()
	#######################################################################
	#######################################################################
	#=MinRunCapacity = [10 12 20];
	MaxRunCapacity = [20 40 30];
	RU = [10 10 10]; 				# [MW/time period]
	RD = [5 15 10];					# [MW/time period] RD = [20 20];
	UT = [1 1 2];					# [time period]
	DT = [1 1 4];					# [time period]
	SU = [10 12 20]; # Doit être au moins aussi grand que MinRunCapacity
	SD = [15 12 20]; # Doit être au moins aussi grand que MinRunCapacity
	F = [70 75 5];				# [EUR]
	C = [50 56 2];				# [EUR/MW]
	NoLoadConsumption = [1 2 6];

	# data_demand = rand(1,24)*sum(MaxRunCapacity);
	# data_demand = [15 10 20 15 10 15]; ## Example 1
	## data_demand = [19.9827  2.64094  11.6176  12.9832  17.7714  5.62079  8.10247  19.7018  2.68533  6.51353]; ## Example 3
	data_demand = [18.1187  2.29986  1.49833  0.00790886  15.3403  15.8105  11.7689  13.755  3.24349  6.00786]; ## Example 4
	## data_demand = [3.21029 13.3892 17.5783 18.5449 5.62669 16.5953 0.869666 11.8104 3.28692 17.9711]; ## Example 5
	VOLL = 3000;
	T_max = length(data_demand);
	nb_gen = 1; # 2

	UT = check_UT_DT(UT, T_max)
	DT = check_UT_DT(DT, T_max)=#
	#######################################################################
	#######################################################################
	#=MinRunCapacity = [10 10];
	MaxRunCapacity = [20 20];
	RU = [5 1]; 				# [MW/time period]
	RD = [6 4];					# [MW/time period] RD = [20 20];
	UT = [1 1];					# [time period]
	DT = [1 1];					# [time period]
	SU = [12 10]; # Doit être au moins aussi grand que MinRunCapacity
	SD = [12 15]; # Doit être au moins aussi grand que MinRunCapacity
	F = [5 5];				# [EUR]
	C = [10 10];				# [EUR/MW]
	NoLoadConsumption = [1 1];
	data_demand = [15 15 15 20];
	# data_demand = [15]
	VOLL = 3000;
	T_max = length(data_demand);
	nb_gen = 1;=#
	#######################################################################
	#######################################################################
	nb_gen = 1;
	MinRunCapacity = [6];
	MaxRunCapacity = [16];
	RU = [5];
	RD = [5];
	UT = [1];
	DT = [1];
	SU = [6];
	SD = [6];
	NoLoadConsumption 	= [10];
	F 		= [53];
	C 			= [30];
	data_demand = [6 11 16 11];
	VOLL = 3000
	T_max = length(data_demand);

	#######################################################################
	#######################################################################

	u_prior = zeros(nb_gen);
	v_prior = zeros(nb_gen);
	w_prior = zeros(nb_gen);
	p_prior = zeros(nb_gen);
	pbar_prior = zeros(nb_gen);

	u_posterior = zeros(nb_gen);
	v_posterior = zeros(nb_gen);
	w_posterior = zeros(nb_gen);
	p_posterior = zeros(nb_gen);
	pbar_posterior = zeros(nb_gen);

	delta_criterion = 0;
	# (received_profit_producer_EF, max_profit_producer_EF, uplift_producer_EF, received_profit_consumer_EF, max_profit_consumer_EF, uplift_consumer_EF, received_profit_producer_BD, max_profit_producer_BD, uplift_producer_BD, received_profit_consumer_BD, max_profit_consumer_BD, uplift_consumer_BD, received_profit_producer_LP, max_profit_producer_LP, uplift_producer_LP, received_profit_consumer_LP, max_profit_consumer_LP, uplift_consumer_LP, obj_real_matching, obj_real_relax, obj_real_EF, obj_real_BD, Solving_time_EF, Solving_time_BD, Solving_time_LP, Solving_time_matching, nb_cuts_tot, nb_salves_solve, Solving_time_CG) = ComputeUplift(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, delta_criterion, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
	(received_profit_producer_EF, max_profit_producer_EF, uplift_producer_EF, received_profit_consumer_EF, max_profit_consumer_EF, uplift_consumer_EF, received_profit_producer_BD, max_profit_producer_BD, uplift_producer_BD, received_profit_consumer_BD, max_profit_consumer_BD, uplift_consumer_BD, received_profit_producer_LP, max_profit_producer_LP, uplift_producer_LP, received_profit_consumer_LP, max_profit_consumer_LP, uplift_consumer_LP, received_profit_producer_CG, max_profit_producer_CG, uplift_producer_CG, received_profit_consumer_CG, max_profit_consumer_CG, uplift_consumer_CG, received_profit_producer_ColG, max_profit_producer_ColG, uplift_producer_ColG, received_profit_consumer_ColG, max_profit_consumer_ColG, uplift_consumer_ColG, obj_real_matching, obj_real_relax, obj_real_EF, obj_real_BD, obj_restricted, Solving_time_EF, Solving_time_BD, Solving_time_LP, Solving_time_matching, nb_cuts_tot, nb_salves_solve, Solving_time_CG) = ComputeUplift(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, delta_criterion, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
	println();
	println();
	println("data_demand : ",data_demand);
	println("obj matching : ",obj_real_matching);
	println("obj EF       : ",obj_real_EF);
	println("obj BD       : ",obj_real_BD);
	println("obj column   : ",obj_restricted);
	println("obj LP       : ",obj_real_relax);
	println();
	println();
	println("Duality GAP (EF)  (obj_matching - obj_EF)  : ", obj_real_matching-obj_real_EF);
	println("Duality GAP (BD)  (obj_matching - obj_BD)  : ", obj_real_matching-obj_real_BD);
	println("Duality GAP (CRG) (obj_matching - obj_CRG) : ", obj_real_matching-obj_restricted);
	println("Duality GAP (LP)  (obj_matching - obj_LP)  : ", obj_real_matching-obj_real_relax); #obj_real_matching-obj_real_relax
	println();
	println();
	#=println("max_profit_producer_EF : ",max_profit_producer_EF);
	println("received_profit_producer_EF : ",received_profit_producer_EF);
	println("uplift_producer_EF : ", uplift_producer_EF);
	println();
	println("max_profit_consumer_EF : ",max_profit_consumer_EF);
	println("received_profit_consumer_EF : ",received_profit_consumer_EF);
	println("uplift_consumer_EF : ", uplift_consumer_EF);
	println();=#
	println("sum of the maximum profit (should equal obj EF) ",sum(max_profit_producer_EF) + sum(max_profit_consumer_EF)," (",obj_real_EF,")");
	println("sum of the received profit (should equal obj matching) ",sum(received_profit_producer_EF) + sum(received_profit_consumer_EF)," (",obj_real_matching,")");
	println("sum uplifts using CHP EF (shoud equal dual GAP EF) ",sum(uplift_producer_EF) + sum(uplift_consumer_EF)," (",obj_real_matching-obj_real_EF,")");
	println("Solving time EF : ",Solving_time_EF);
	println();
	println();
	#=println("max_profit_producer_BD : ", max_profit_producer_BD);
	println("received_profit_producer_BD : ", received_profit_producer_BD);
	println("uplift_producer_BD : ", uplift_producer_BD);
	println();
	println("max_profit_consumer_BD : ", max_profit_consumer_BD);
	println("received_profit_consumer_BD : ", received_profit_consumer_BD);
	println("uplift_consumer_BD : ", uplift_consumer_BD);
	println();=#
	println("sum of the maximum profit (should equal obj BD) ",sum(max_profit_producer_BD) + sum(max_profit_consumer_BD)," (",obj_real_BD,")");
	println("sum of the received profit (should equal obj matching) ",sum(received_profit_producer_BD) + sum(received_profit_consumer_BD)," (",obj_real_matching,")");
	println("sum uplifts using CHP BD (should equal dual GAP BD) ",sum(uplift_producer_BD) + sum(uplift_consumer_BD)," (",obj_real_matching-obj_real_BD,")");
	println("Solving time BD : ",Solving_time_BD);
	println("slaves violated : ",nb_salves_solve);
	println("cuts added : ",nb_cuts_tot);
	println();
	println();
	#=println("max_profit_producer_CG : ", max_profit_producer_CG);
	println("received_profit_producer_CG : ", received_profit_producer_CG);
	println("uplift_producer_CG : ", uplift_producer_CG);
	println();
	println("max_profit_consumer_CG : ", max_profit_consumer_CG);
	println("received_profit_consumer_CG : ", received_profit_consumer_CG);
	println("uplift_consumer_CG : ", uplift_consumer_CG);
	println();=#
	println("sum of the maximum profit (should equal obj CG) ",sum(max_profit_producer_CG) + sum(max_profit_consumer_CG)," (",obj_restricted,")");
	println("sum of the received profit (should equal obj matching) ",sum(received_profit_producer_CG) + sum(received_profit_consumer_CG)," (",obj_real_matching,")");
	println("sum uplifts using CHP CG (should equal dual GAP CG) ",sum(uplift_producer_CG) + sum(uplift_consumer_CG)," (",obj_real_matching-obj_restricted,")");
	println("Solving time CG : ",Solving_time_CG);
	println();
	println();
	#=println("max_profit_producer_LP : ",max_profit_producer_LP);
	println("received_profit_producer_LP : ",received_profit_producer_LP);
	println("uplift_producer_LP : ",uplift_producer_LP);
	println();
	println("max_profit_consumer_LP : ",max_profit_consumer_LP);
	println("received_profit_consumer_LP : ",received_profit_consumer_LP);
	println("uplift_consumer_LP : ",uplift_consumer_LP);
	println();=#
	println("sum of the maximum profit (should equal obj LP) ",sum(max_profit_producer_LP) + sum(max_profit_consumer_LP)," (",obj_real_relax,")");
	println("sum of the received profit (should equal obj matching) ",sum(received_profit_producer_LP) + sum(received_profit_consumer_LP)," (",obj_real_matching,")");
	println("sum uplifts using CHP LP (should equal dual GAP LP) ",sum(uplift_producer_LP) + sum(uplift_consumer_LP)," (",obj_real_matching-obj_real_relax,")");
	println("Solving time LP : ",Solving_time_LP);
end

# Check_Algorithms()


function Belgian_System_Network()
	df1 = CSV.read("D:/Documents/Mémoire Papavasiliou/Unit commitment Toy/UCDataBE/Demand.csv", DataFrame);
	Demand = convert(Matrix, df1);
	df2 = CSV.read("D:/Documents/Mémoire Papavasiliou/Unit commitment Toy/UCDataBE/Generators.csv", DataFrame);
	Generators = convert(Matrix, df2);
	######## Demand ########
	Data_types = Demand[:,1];
	(a,b) = size(Demand)
	Demand = convert(Matrix{Float64}, Demand[:,2:b]);
	AutumnWE = Demand[5,:];
	WinterWE = Demand[8,:];
	SpringWE = Demand[6,:];
	SummerWE = Demand[7,:];
	AutumnWD = Demand[1,:];
	WinterWD = Demand[4,:];
	SpringWD = Demand[2,:];
	SummerWD = Demand[3,:];

	######## Generators ########
	nb_gen = 68; # There are normally 68 generators

	(a,b) = size(Generators);
	dt = 1; # [hour]
	Generators = convert(Matrix{Float64}, Generators[1:nb_gen,2:b]);
	MinRunCapacity 		= Generators[:,1];								# [MW]
	MaxRunCapacity 		= Generators[:,2];
	RampUp				= (60*dt)*Generators[:,3];						# [MW/min] must multiply by 15 to have [MW/15-min]
	RampDown 			= (60*dt)*Generators[:,4];						# [MW/min] must multiply by 15
	UT_1 				= convert(Array{Int64}, (1/dt)*Generators[:,5]);# [15 min]
	DT_1		 		= convert(Array{Int64}, (1/dt)*Generators[:,6]);# [15 min]
	NoLoadConsumption 	= Generators[:,7];								# [MW], NoLoadCost = NoLoadConsumption * MargCost
	StartupCost 		= Generators[:,8];								# [€]
	MargCost 			= dt*Generators[:,9];							# [€/MWh]
	VOLL = 3000;
	SU = MinRunCapacity;
	SD = MinRunCapacity;

	# T_max = 35; # There normally 96 time-steps  25

	# Demand in hours
	data_demand_vec = WinterWE;
	data_demand = zeros(24);
	c = 1;
	for i=1:24
		data_demand[i] = data_demand_vec[c];
		c+=4;
	end
	T_max = 24;
	
	
	# T_max = 35;
	#=T_max = 35;
	data_demand_vec = AutumnWE; # WinterWE
	data_demand = data_demand_vec[1:T_max];=#

	UT = check_UT_DT(UT_1, T_max);
	DT = check_UT_DT(DT_1, T_max);
	# delta_criterion = 10^(-2);
	delta_criterion = 0;

	u_prior 		= zeros(nb_gen);
	v_prior 		= zeros(nb_gen);
	w_prior 		= zeros(nb_gen);
	p_prior 		= zeros(nb_gen);
	pbar_prior 		= zeros(nb_gen);
	u_posterior 	= zeros(nb_gen);
	v_posterior 	= zeros(nb_gen);
	w_posterior 	= zeros(nb_gen);
	p_posterior 	= zeros(nb_gen);
	pbar_posterior 	= zeros(nb_gen);

	#=(received_profit_producer_EF, max_profit_producer_EF, uplift_producer_EF, received_profit_consumer_EF, max_profit_consumer_EF, uplift_consumer_EF, received_profit_producer_BD, max_profit_producer_BD, uplift_producer_BD, received_profit_consumer_BD, max_profit_consumer_BD, uplift_consumer_BD, received_profit_producer_LP, max_profit_producer_LP, uplift_producer_LP, received_profit_consumer_LP, max_profit_consumer_LP, uplift_consumer_LP, received_profit_producer_CG, max_profit_producer_CG, uplift_producer_CG, received_profit_consumer_CG, max_profit_consumer_CG, uplift_consumer_CG, received_profit_producer_ColG, max_profit_producer_ColG, uplift_producer_ColG, received_profit_consumer_ColG, max_profit_consumer_ColG, uplift_consumer_ColG, obj_real_matching, obj_real_relax, obj_real_EF, obj_real_BD, obj_restricted, Solving_time_EF, Solving_time_BD, Solving_time_LP, Solving_time_matching, nb_cuts_tot, nb_salves_solve, Solving_time_CG, obj_column_generation) = ComputeUplift(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, delta_criterion, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
	println();println();println();println();
	println("received_profit_producer_ColG = $(received_profit_producer_ColG)");
	println("max_profit_producer_ColG = $(max_profit_producer_ColG)");
	println("uplift_producer_ColG = $(uplift_producer_ColG)");
	println("received_profit_consumer_ColG = $(received_profit_consumer_ColG)");
	println("max_profit_consumer_ColG = $(max_profit_consumer_ColG)");
	println("uplift_consumer_ColG = $(uplift_consumer_ColG)");
	println();
	println("received_profit_producer_BD : $(received_profit_producer_BD)");
	println("max_profit_producer_BD : $(max_profit_producer_BD)");
	println("uplift_producer_BD : $(uplift_producer_BD)");
	println("received_profit_consumer_BD : $(received_profit_consumer_BD)");
	println("max_profit_consumer_BD : $(max_profit_consumer_BD)");
	println("uplift_consumer_BD : $(uplift_consumer_BD)");
	println();
	println("received_profit_producer_CG = $(received_profit_producer_CG)")
	println("max_profit_producer_CG = $(max_profit_producer_CG)")
	println("uplift_producer_CG = $(uplift_producer_CG)")
	println("received_profit_consumer_CG = $(received_profit_consumer_CG)")
	println("max_profit_consumer_CG = $(max_profit_consumer_CG)")
	println("uplift_consumer_CG = $(uplift_consumer_CG)")

	println();
	println();
	println("data_demand    : ",data_demand);
	println("obj matching   : ",obj_real_matching);
	println("obj EF         : ",obj_real_EF);
	println("obj BD         : ",obj_real_BD);
	println("obj_restricted : $(obj_restricted)");
	println("obj col gen    : $(obj_column_generation)");
	println("obj LP         : ",obj_real_relax);
	println();
	println();
	println("Extended Formuation");
	println("sum uplifts using CHP EF (shoud equal dual GAP EF) ",sum(uplift_producer_EF) + sum(uplift_consumer_EF)," (",obj_real_matching-obj_real_EF,")");
	println("Solving time EF : ",Solving_time_EF);
	println();
	println();
	println("Bender Decomposition");
	println("sum uplifts using CHP BD (should equal dual GAP BD) ",sum(uplift_producer_BD) + sum(uplift_consumer_BD)," (",obj_real_matching-obj_real_BD,")");
	println("Solving time BD : ",Solving_time_BD);
	#println("nb_cuts_tot : ",nb_cuts_tot);
	#println("nb_salves_solve : ",nb_salves_solve);
	println();
	println();
	println("Column-and-Row Generation");
	println("sum uplifts using CHP CRG (should equal dual GAP CRG) ",sum(uplift_producer_CG) + sum(uplift_consumer_CG)," (",obj_real_matching-obj_restricted,")");
	println("Solving_time_CG : ",Solving_time_CG);
	println();
	println();
	println("Column Generation");
	println("sum uplifts using CHP RG (should equal dual GAP CG) ",sum(uplift_producer_ColG) + sum(uplift_consumer_ColG)," (",obj_real_matching-obj_column_generation,")");
	#println("Solving_time_CG : ",Solving_time_CG);
	println();
	println();
	println("Linear Relaxation");
	println("sum uplifts using CHP LP (should equal dual GAP LP) ",sum(uplift_producer_LP) + sum(uplift_consumer_LP)," (",obj_real_matching-obj_real_relax,")");
	println("Solving time LP : ",Solving_time_LP);
	println();
	println();
	println("Solving time matching : ",Solving_time_matching);
	println();println();println();=#

	max_profit_producer_EF = [-0.0, 1.04519e6, 1.02607e6, 1.15747e7, 1.37827e6, 1.04617e6, 1.74427e6, 1.52849e6, 1.51673e6, 9.3159e6, 424923.0, -0.0, 1.26147e7, 1.02826e7, 5.20847e6, 5.18638e6, 5.19765e6, 5.19683e6, 5.20892e6, 5.21173e6, 1.80381e6, 1.15371e6, -0.0, 6.57562e6, 1.87849e6, 7.62674e5, 2.83593e6, 1.09087e7, 3.81097e6, 3.84805e6, 2.17181e6, 1.61051e6, 1.39341e6, 5.21506e6, 1.28164e6, 1.26706e6, 1.28532e6, 6.79236e5, 1.15097e7, 1.06647e7, 9.2194e6, 3.23243e6, 4.91025e5, 2.61628e5, 5.21338e6, 1.44336e6, 1.20358e6, 1.33255e7, 9.81202e6, 5.98251e6, 6.47633e6, 4.83243e6, 3.00081e6, 2.81055e6, 1.15332e6, 4.64487e6, 2.48086e6, 1.08743e6, 5.39478e5, 1.32359e7, 7.93892e5, -0.0, -0.0, -0.0, -0.0, 1.27212e6, 5.84001e6, 5.44382e6]
	received_profit_producer_EF = [0.0, 1.04519e6, 1.02607e6, 1.15747e7, 1.37827e6, 1.04617e6, 1.73828e6, 1.52849e6, 1.51625e6, 9.3159e6, 424923.0, 0.0, 1.26147e7, 1.02826e7, 5.20847e6, 5.18638e6, 5.19765e6, 5.19683e6, 5.20892e6, 5.21173e6, 1.79474e6, 1.15371e6, 0.0, 6.57562e6, 1.87849e6, 7.62674e5, 2.83593e6, 1.09087e7, 3.81097e6, 3.84805e6, 2.17181e6, 1.61051e6, 1.39341e6, 5.21506e6, 1.28164e6, 1.26706e6, 1.28532e6, 6.79236e5, 1.15097e7, 1.06647e7, 9.2194e6, 3.2309e6, 4.91025e5, 2.61628e5, 5.21338e6, 1.44336e6, 1.20358e6, 1.33255e7, 9.81202e6, 5.98251e6, 6.47633e6, 4.83243e6, 3.00081e6, 2.81055e6, 1.14922e6, 4.64487e6, 2.48086e6, 1.08743e6, 5.39478e5, 1.32359e7, 7.93892e5, 0.0, 0.0, 0.0, 0.0, 1.27212e6, 5.84001e6, 5.44382e6]
	uplift_producer_EF = [-0.0, 3.49246e-10, 0.0, -1.86265e-9, 4.65661e-10, -1.16415e-10, 5985.18, 4.65661e-10, 478.819, 1.86265e-9, -5.82077e-11, -0.0, -1.86265e-9, 1.86265e-9, -9.31323e-10, 0.0, 9.31323e-10, 0.0, -9.31323e-10, 0.0, 9070.23, 2.32831e-10, -0.0, 1.86265e-9, -4.65661e-10, -1.04774e-9, -4.65661e-10, 3.72529e-9, -1.39698e-9, 1.39698e-9, 0.0, 0.0, -2.32831e-10, 0.0, -2.32831e-10, 4.65661e-10, 0.0, 0.0, -1.86265e-9, 1.86265e-9, -1.86265e-9, 1531.96, 5.82077e-11, 0.0, 0.0, -2.32831e-10, 2.32831e-10, 1.86265e-9, -1.86265e-9, 0.0, 0.0, 9.31323e-10, -4.65661e-10, 0.0, 4105.73, 0.0, 0.0, 0.0, -1.04774e-9, 0.0, 3.49246e-10, -0.0, -0.0, -0.0, -0.0, -2.32831e-10, 2.79397e-9, 1.86265e-9]
	max_profit_consumer_EF = 3.480844521331947e8;
	received_profit_consumer_EF = 3.480844521331947e8;
	uplift_consumer_EF = 0.0;
	
	max_profit_producer_LP = [-0.0, 1.04553e6, 1.02641e6, 1.15855e7, 1.37984e6, 1.04729e6, 1.74315e6, 1.52899e6, 1.51691e6, 9.32448e6, 4.25279e5, -0.0, 1.2626e7, 1.0292e7, 5.21345e6, 5.19136e6, 5.20263e6, 5.20181e6, 5.2139e6, 5.21671e6, 1.80224e6, 1.15481e6, -0.0, 6.58283e6, 1.88021e6, 7.61367e5, 2.83824e6, 1.09187e7, 3.81441e6, 3.85136e6, 2.1739e6, 1.61193e6, 1.39459e6, 5.22004e6, 1.2827e6, 1.26812e6, 1.28637e6, 6.79795e5, 1.152e7, 1.06739e7, 9.22689e6, 3.23057e6, 4.91182e5, 2.61886e5, 5.21836e6, 1.44458e6, 1.20456e6, 1.33374e7, 9.82061e6, 5.98972e6, 6.48211e6, 4.83689e6, 3.004e6, 2.81354e6, 1.153e6, 4.64865e6, 2.48188e6, 1.08849e6, 5.38574e5, 1.32477e7, 7.94628e5, -0.0, -0.0, -0.0, -0.0, 1.27317e6, 5.84634e6, 5.45015e6]
	received_profit_producer_LP = [0.0, 1.04553e6, 1.02641e6, 1.15855e7, 1.37984e6, 1.04729e6, 1.73616e6, 1.52899e6, 1.51337e6, 9.32448e6, 4.25279e5, 0.0, 1.2626e7, 1.0292e7, 5.21345e6, 5.19136e6, 5.20263e6, 5.20181e6, 5.2139e6, 5.21671e6, 1.78947e6, 1.15481e6, 0.0, 6.5828e6, 1.87844e6, 7.61367e5, 2.83824e6, 1.09187e7, 3.81441e6, 3.85136e6, 2.1739e6, 1.61193e6, 1.39459e6, 5.22004e6, 1.2827e6, 1.26812e6, 1.28637e6, 6.79795e5, 1.152e7, 1.06739e7, 9.22689e6, 3.23057e6, 4.91182e5, 2.61886e5, 5.21836e6, 1.44458e6, 1.20456e6, 1.33374e7, 9.82061e6, 5.98972e6, 6.47838e6, 4.83687e6, 3.004e6, 2.81354e6, 1.14649e6, 4.64865e6, 2.48168e6, 1.08849e6, 5.38574e5, 1.32477e7, 7.94628e5, 0.0, 0.0, 0.0, 0.0, 1.27317e6, 5.84634e6, 5.45015e6]
	uplift_producer_LP = [-0.0, 2.32831e-10, 1.16415e-10, 0.0, 2.32831e-10, 0.0, 6983.48, 2.32831e-10, 3533.68, 1.86265e-9, -1.16415e-10, -0.0, 1.86265e-9, 0.0, 0.0, 0.0, 0.0, 9.31323e-10, 0.0, 9.31323e-10, 12769.8, 2.32831e-10, -0.0, 30.3147, 1767.08, -1.16415e-9, 0.0, 0.0, -4.65661e-10, 9.31323e-10, 4.65661e-10, 0.0, -6.98492e-10, -9.31323e-10, -4.65661e-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.82077e-11, 2.91038e-11, -9.31323e-10, 0.0, 0.0, 3.72529e-9, -1.86265e-9, -9.31323e-10, 3723.06, 22.6394, -4.65661e-10, 4.65661e-10, 6515.5, 0.0, 205.177, 0.0, -1.04774e-9, 0.0, 2.32831e-10, -0.0, -0.0, -0.0, -0.0, 2.32831e-10, 1.86265e-9, 1.86265e-9]
	max_profit_consumer_LP = 3.4787824750056434e8;
	received_profit_consumer_LP = 3.4787824750056434e8;
	uplift_consumer_LP = 0.0;

	received_profit_producer_ColG = [0.0, 1.04519e6, 1.02607e6, 1.15747e7, 1.37827e6, 1.04617e6, 1.73828e6, 1.52849e6, 1.51625e6, 9.3159e6, 424923.0, 0.0, 1.26147e7, 1.02826e7, 5.20847e6, 5.18638e6, 5.19765e6, 5.19683e6, 5.20892e6, 5.21173e6, 1.79474e6, 1.15371e6, 0.0, 6.57562e6, 1.87849e6, 7.62674e5, 2.83593e6, 1.09087e7, 3.81097e6, 3.84805e6, 2.17181e6, 1.61051e6, 1.39341e6, 5.21506e6, 1.28164e6, 1.26706e6, 1.28532e6, 6.79236e5, 1.15097e7, 1.06647e7, 9.2194e6, 3.2309e6, 4.91025e5, 2.61628e5, 5.21338e6, 1.44336e6, 1.20358e6, 1.33255e7, 9.81202e6, 5.98251e6, 6.47633e6, 4.83243e6, 3.00081e6, 2.81055e6, 1.14922e6, 4.64487e6, 2.48086e6, 1.08743e6, 5.39478e5, 1.32359e7, 7.93892e5, 0.0, 0.0, 0.0, 0.0, 1.27212e6, 5.84001e6, 5.44382e6]
	max_profit_producer_ColG = [-0.0, 1.04519e6, 1.02607e6, 1.15747e7, 1.37827e6, 1.04617e6, 1.74427e6, 1.52849e6, 1.51673e6, 9.3159e6, 424923.0, -0.0, 1.26147e7, 1.02826e7, 5.20847e6, 5.18638e6, 5.19765e6, 5.19683e6, 5.20892e6, 5.21173e6, 1.80381e6, 1.15371e6, -0.0, 6.57562e6, 1.87849e6, 7.62674e5, 2.83593e6, 1.09087e7, 3.81097e6, 3.84805e6, 2.17181e6, 1.61051e6, 1.39341e6, 5.21506e6, 1.28164e6, 1.26706e6, 1.28532e6, 6.79236e5, 1.15097e7, 1.06647e7, 9.2194e6, 3.23243e6, 4.91025e5, 2.61628e5, 5.21338e6, 1.44336e6, 1.20358e6, 1.33255e7, 9.81202e6, 5.98251e6, 6.47633e6, 4.83243e6, 3.00081e6, 2.81055e6, 1.15332e6, 4.64487e6, 2.48086e6, 1.08743e6, 5.39478e5, 1.32359e7, 7.93892e5, -0.0, -0.0, -0.0, -0.0, 1.27212e6, 5.84001e6, 5.44382e6]
	uplift_producer_ColG = [-0.0, 6.98492e-10, 1.16415e-10, 0.0, 6.98492e-10, 2.32831e-10, 5985.18, 6.98492e-10, 478.819, 5.58794e-9, -5.82077e-11, -0.0, 1.86265e-9, 1.86265e-9, 1.86265e-9, 2.79397e-9, 2.79397e-9, 9.31323e-10, 9.31323e-10, 1.86265e-9, 9070.23, 4.65661e-10, -0.0, 3.72529e-9, 0.0, -6.98492e-10, 0.0, 3.72529e-9, -4.65661e-10, 2.79397e-9, 9.31323e-10, 4.65661e-10, 0.0, 9.31323e-10, 2.32831e-10, 6.98492e-10, 2.32831e-10, 2.32831e-10, 0.0, 3.72529e-9, 1.86265e-9, 1531.96, 1.16415e-10, 5.82077e-11, 1.86265e-9, 4.65661e-10, 4.65661e-10, 3.72529e-9, 1.86265e-9, 1.86265e-9, 0.0, 2.79397e-9, 9.31323e-10, 9.31323e-10, 4105.73, 9.31323e-10, 9.31323e-10, 2.32831e-10, -8.14907e-10, 0.0, 3.49246e-10, -0.0, -0.0, -0.0, -0.0, 4.65661e-10, 4.65661e-9, 3.72529e-9]
	received_profit_consumer_ColG = 3.480844521331946e8
	max_profit_consumer_ColG = 3.4808445213319457e8
	uplift_consumer_ColG = -5.960464477539063e-8

	received_profit_producer_BD = [0.0, 1.04519e6, 1.02607e6, 1.15747e7, 1.37827e6, 1.04617e6, 1.73828e6, 1.52849e6, 1.51625e6, 9.3159e6, 424923.0, 0.0, 1.26147e7, 1.02826e7, 5.20847e6, 5.18638e6, 5.19765e6, 5.19683e6, 5.20892e6, 5.21173e6, 1.79474e6, 1.15371e6, 0.0, 6.57562e6, 1.87849e6, 7.62674e5, 2.83593e6, 1.09087e7, 3.81097e6, 3.84805e6, 2.17181e6, 1.61051e6, 1.39341e6, 5.21506e6, 1.28164e6, 1.26706e6, 1.28532e6, 6.79236e5, 1.15097e7, 1.06647e7, 9.2194e6, 3.2309e6, 4.91025e5, 2.61628e5, 5.21338e6, 1.44336e6, 1.20358e6, 1.33255e7, 9.81202e6, 5.98251e6, 6.47633e6, 4.83243e6, 3.00081e6, 2.81055e6, 1.14922e6, 4.64487e6, 2.48086e6, 1.08743e6, 5.39478e5, 1.32359e7, 7.93892e5, 0.0, 0.0, 0.0, 0.0, 1.27212e6, 5.84001e6, 5.44382e6]
	max_profit_producer_BD = [-0.0, 1.04519e6, 1.02607e6, 1.15747e7, 1.37827e6, 1.04617e6, 1.74427e6, 1.52849e6, 1.51673e6, 9.3159e6, 424923.0, -0.0, 1.26147e7, 1.02826e7, 5.20847e6, 5.18638e6, 5.19765e6, 5.19683e6, 5.20892e6, 5.21173e6, 1.80381e6, 1.15371e6, -0.0, 6.57562e6, 1.87849e6, 7.62674e5, 2.83593e6, 1.09087e7, 3.81097e6, 3.84805e6, 2.17181e6, 1.61051e6, 1.39341e6, 5.21506e6, 1.28164e6, 1.26706e6, 1.28532e6, 6.79236e5, 1.15097e7, 1.06647e7, 9.2194e6, 3.23243e6, 4.91025e5, 2.61628e5, 5.21338e6, 1.44336e6, 1.20358e6, 1.33255e7, 9.81202e6, 5.98251e6, 6.47633e6, 4.83243e6, 3.00081e6, 2.81055e6, 1.15332e6, 4.64487e6, 2.48086e6, 1.08743e6, 5.39478e5, 1.32359e7, 7.93892e5, -0.0, -0.0, -0.0, -0.0, 1.27212e6, 5.84001e6, 5.44382e6]
	uplift_producer_BD = [-0.0, 3.49246e-10, -1.16415e-10, 0.0, 2.32831e-10, 1.16415e-10, 5985.18, 2.32831e-10, 478.819, 1.86265e-9, -1.16415e-10, -0.0, 0.0, -1.86265e-9, 0.0, -9.31323e-10, 9.31323e-10, 0.0, -9.31323e-10, 9.31323e-10, 9070.23, 0.0, -0.0, 9.31323e-10, -6.98492e-10, -1.04774e-9, -4.65661e-10, 1.86265e-9, -1.39698e-9, 1.39698e-9, 4.65661e-10, 4.65661e-10, -4.65661e-10, 0.0, 0.0, 2.32831e-10, 0.0, 0.0, -1.86265e-9, 0.0, 0.0, 1531.96, 5.82077e-11, 0.0, -9.31323e-10, 2.32831e-10, 4.65661e-10, 0.0, -1.86265e-9, -9.31323e-10, 0.0, 0.0, 0.0, 0.0, 4105.73, 9.31323e-10, 4.65661e-10, 0.0, -1.16415e-9, -1.86265e-9, 1.16415e-10, -0.0, -0.0, -0.0, -0.0, 0.0, 1.86265e-9, 1.86265e-9]
	received_profit_consumer_BD = 3.480844521331946e8
	max_profit_consumer_BD = 3.480844521331946e8
	uplift_consumer_BD = 0.0

	received_profit_producer_CG = [0.0, 1.04519e6, 1.02607e6, 1.15747e7, 1.37827e6, 1.04617e6, 1.73828e6, 1.52849e6, 1.51625e6, 9.3159e6, 424923.0, 0.0, 1.26147e7, 1.02826e7, 5.20847e6, 5.18638e6, 5.19765e6, 5.19683e6, 5.20892e6, 5.21173e6, 1.79474e6, 1.15371e6, 0.0, 6.57562e6, 1.87849e6, 7.62674e5, 2.83593e6, 1.09087e7, 3.81097e6, 3.84805e6, 2.17181e6, 1.61051e6, 1.39341e6, 5.21506e6, 1.28164e6, 1.26706e6, 1.28532e6, 6.79236e5, 1.15097e7, 1.06647e7, 9.2194e6, 3.2309e6, 4.91025e5, 2.61628e5, 5.21338e6, 1.44336e6, 1.20358e6, 1.33255e7, 9.81202e6, 5.98251e6, 6.47633e6, 4.83243e6, 3.00081e6, 2.81055e6, 1.14922e6, 4.64487e6, 2.48086e6, 1.08743e6, 5.39478e5, 1.32359e7, 7.93892e5, 0.0, 0.0, 0.0, 0.0, 1.27212e6, 5.84001e6, 5.44382e6]
	max_profit_producer_CG = [-0.0, 1.04519e6, 1.02607e6, 1.15747e7, 1.37827e6, 1.04617e6, 1.74427e6, 1.52849e6, 1.51673e6, 9.3159e6, 424923.0, -0.0, 1.26147e7, 1.02826e7, 5.20847e6, 5.18638e6, 5.19765e6, 5.19683e6, 5.20892e6, 5.21173e6, 1.80381e6, 1.15371e6, -0.0, 6.57562e6, 1.87849e6, 7.62674e5, 2.83593e6, 1.09087e7, 3.81097e6, 3.84805e6, 2.17181e6, 1.61051e6, 1.39341e6, 5.21506e6, 1.28164e6, 1.26706e6, 1.28532e6, 6.79236e5, 1.15097e7, 1.06647e7, 9.2194e6, 3.23243e6, 4.91025e5, 2.61628e5, 5.21338e6, 1.44336e6, 1.20358e6, 1.33255e7, 9.81202e6, 5.98251e6, 6.47633e6, 4.83243e6, 3.00081e6, 2.81055e6, 1.15332e6, 4.64487e6, 2.48086e6, 1.08743e6, 5.39478e5, 1.32359e7, 7.93892e5, -0.0, -0.0, -0.0, -0.0, 1.27212e6, 5.84001e6, 5.44382e6]
	uplift_producer_CG = [-0.0, 4.65661e-10, 0.0, -1.86265e-9, 2.32831e-10, 0.0, 5985.18, 4.65661e-10, 478.819, 1.86265e-9, -1.16415e-10, -0.0, -1.86265e-9, 0.0, 0.0, 9.31323e-10, 1.86265e-9, 9.31323e-10, -9.31323e-10, 0.0, 9070.23, 2.32831e-10, -0.0, 1.86265e-9, -6.98492e-10, -1.16415e-9, 0.0, 1.86265e-9, -9.31323e-10, 9.31323e-10, 4.65661e-10, 0.0, -2.32831e-10, 0.0, -2.32831e-10, 4.65661e-10, 0.0, 0.0, -1.86265e-9, 1.86265e-9, 1.86265e-9, 1531.96, 0.0, 0.0, 9.31323e-10, 2.32831e-10, 4.65661e-10, 0.0, 0.0, -9.31323e-10, -9.31323e-10, 0.0, -9.31323e-10, -4.65661e-10, 4105.73, 9.31323e-10, 0.0, 0.0, -1.04774e-9, 0.0, 2.32831e-10, -0.0, -0.0, -0.0, -0.0, -2.32831e-10, 1.86265e-9, 2.79397e-9]
	received_profit_consumer_CG = 3.480844521331946e8
	max_profit_consumer_CG = 3.480844521331946e8
	uplift_consumer_CG = 0.0

	#######
	### Graph for uplifts
	obj_bar = zeros(nb_gen,3); # nb_gen+1
	ctg = repeat(["Max Profit", "Received Profit", "Uplifts"], inner=nb_gen); # nb_gen+1
	for g=1:nb_gen
		obj_bar[g,1] = max_profit_producer_EF[g];
		obj_bar[g,2] = received_profit_producer_EF[g];
		obj_bar[g,3] = uplift_producer_EF[g];
	end
	#obj_bar[nb_gen+1,1] = max_profit_consumer_EF;
	#obj_bar[nb_gen+1,2] = received_profit_consumer_EF;
	#obj_bar[nb_gen+1,3] = uplift_consumer_EF;
	
	uplifts_vec = zeros(nb_gen,2);
	ctg = repeat(["Uplifts from EF", "Uplifts from LP"], inner=nb_gen); # nb_gen+1
	for g = 1:nb_gen
		uplifts_vec[g,1] = uplift_producer_EF[g]
		uplifts_vec[g,2] = uplift_producer_LP[g]
	end

	uplifts_vec2 = zeros(nb_gen,3);
	ctg2 = repeat(["Uplifts from EF", "Uplifts from LP", "Uplifts from RG"], inner=nb_gen); # nb_gen+1
	for g = 1:nb_gen
		uplifts_vec2[g,1] = uplift_producer_EF[g]
		uplifts_vec2[g,2] = uplift_producer_LP[g]
		uplifts_vec2[g,3] = uplift_producer_BD[g]
	end

	uplifts_vec3 = zeros(nb_gen,4);
	ctg3 = repeat(["Uplifts from EF", "Uplifts from LP", "Uplifts from RG","Uplifts from CG"], inner=nb_gen); # nb_gen+1
	for g = 1:nb_gen
		uplifts_vec3[g,1] = uplift_producer_EF[g]
		uplifts_vec3[g,2] = uplift_producer_LP[g]
		uplifts_vec3[g,3] = uplift_producer_BD[g]
		uplifts_vec3[g,4] = uplift_producer_ColG[g]
	end

	uplifts_vec4 = zeros(nb_gen,5);
	ctg4 = repeat(["Uplifts from EF", "Uplifts from LP", "Uplifts from RG","Uplifts from CG","Uplifts from CRG"], inner=nb_gen); # nb_gen+1
	for g = 1:nb_gen
		uplifts_vec4[g,1] = uplift_producer_EF[g]
		uplifts_vec4[g,2] = uplift_producer_LP[g]
		uplifts_vec4[g,3] = uplift_producer_BD[g]
		uplifts_vec4[g,4] = uplift_producer_ColG[g]
		uplifts_vec4[g,5] = uplift_producer_CG[g]
	end

	# x = repeat(1:nb_gen,outer=3);
	#nam = repeat("G" .* string.(1:nb_gen+1), outer = 3)
	pyplot();
	#figure();
	#p1 = StatsPlots.groupedbar(nam, obj_bar,bar_position=:dodge,bar_width=0.7,xlabel="Generator",group=ctg,xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=20,ylabel="[€]"); # ,color=[:blue :red]
	p1 = StatsPlots.groupedbar(obj_bar,bar_position=:dodge,bar_width=0.7,xlabel="Generator",group=ctg,xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=20,ylabel="[€]"); # ,color=[:blue :red]
	display(p1);
	#sfigure();
	p2 = Plots.bar(uplift_producer_EF,bar_width=0.7,xaxis=false,xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=20,ylabel="Uplift [€]", legend=false); # ,color=[:blue :red]
	display(p2);

	#nam = repeat("G" .* string.(1:nb_gen), outer = 2)
	#nam = repeat(string.(1:nb_gen), outer = 2)
	#figure();
	p3 = StatsPlots.groupedbar(uplifts_vec,bar_position=:dodge,bar_width=1.2,group=ctg,xtickfontsize=20,ytickfontsize=20,xguidefontsize=40,yguidefontsize=40,legendfontsize=30,ylabel="Uplifts [€]",xlabel="Generators");
	display(p3);

	pyplot();
	figure();
	p4 = StatsPlots.groupedbar(uplifts_vec2,bar_position=:dodge,bar_width=1.2,group=ctg2,xtickfontsize=20,ytickfontsize=20,xguidefontsize=40,yguidefontsize=40,legendfontsize=30,ylabel="Uplifts [€]",xlabel="Generators");
	display(p4);

	pyplot();
	figure();
	p5 = StatsPlots.groupedbar(uplifts_vec3,bar_position=:dodge,bar_width=1.2,group=ctg3,xtickfontsize=20,ytickfontsize=20,xguidefontsize=40,yguidefontsize=40,legendfontsize=30,ylabel="Uplifts [€]",xlabel="Generators");
	display(p5);

	pyplot();
	figure();
	p6 = StatsPlots.groupedbar(uplifts_vec4,bar_position=:dodge,bar_width=1.2,group=ctg4,xtickfontsize=20,ytickfontsize=20,xguidefontsize=40,yguidefontsize=40,legendfontsize=30,ylabel="Uplifts [€]",xlabel="Generators");
	display(p6);

end

# Belgian_System_Network()



function Process_Example_GitHub_not_saving(filename, Printer=false)
	#data_file = "GitHub_Example/rts_gmlc/2020-02-09.json"
	# data_file = "D:/Documents/Mémoire Papavasiliou/CHP_code/GitHub_Example/rts_gmlc/2020-01-27.json"
	# data_file = "D:/Documents/Mémoire Papavasiliou/CHP_code/GitHub_Example/rts_gmlc/$(filename)"
	data = JSON.parsefile(filename);
	thermal_gens = keys(data["thermal_generators"]);
	# renewable_gens = keys(data["renewable_generators"]);
	time_periods = 1:data["time_periods"]

	nb_gen = length(thermal_gens);
	T_max = data["time_periods"];
	MinRunCapacity = zeros(nb_gen);
	MaxRunCapacity = zeros(nb_gen);
	RampUp = zeros(nb_gen);
	RampDown = zeros(nb_gen);
	UT = zeros(nb_gen);
	DT = zeros(nb_gen);
	SU = zeros(nb_gen);
	SD = zeros(nb_gen);
	StartupCost = zeros(nb_gen);
	MarginalCost = zeros(nb_gen);
	NoLoadConsumption = zeros(nb_gen);
	L = zeros(T_max);
	VOLL = 3000;
	for t in time_periods
		L[t] = data["demand"][t];
	end

	for (i,(g,gen)) in enumerate(data["thermal_generators"])
		# println("gen : $(gen)");
		MinRunCapacity[i] = gen["power_output_minimum"];
		MaxRunCapacity[i] = gen["power_output_maximum"];
		RampUp[i] = gen["ramp_up_limit"];
		RampDown[i] = gen["ramp_down_limit"];
		UT[i] = gen["time_up_minimum"];
		DT[i] = gen["time_down_minimum"];
		SU[i] = gen["ramp_startup_limit"];
		SD[i] = gen["ramp_shutdown_limit"];
		#println("gen[startup] : $(values(gen["startup"]))");
		count = 0.0;
		index = 0.0;
		# piecewise_production['cost']
		for (k,gen_startup) in enumerate(gen["startup"])
			count+=gen_startup["cost"];
			index+=1.0;
		end
		count = count/index;
		StartupCost[i] = count;
		count = 0.0;
		index = 0.0;
		for (k,gen_piecewise_production) in enumerate(gen["piecewise_production"])
			count+=gen_piecewise_production["cost"];
			index+=1.0;
		end
		count = count/index;
		MarginalCost[i] = count;
		#println("gen[unit_on_t0] : $(gen["unit_on_t0"])");
		#println("gen[power_output_t0] : $(gen["power_output_t0"])");
		# push!(df_Generators, ("$(i)", "SLOW", "Belgium",MinRunCapacity[i],MaxRunCapacity[i],RampUp[i],RampDown[i],UT[i],DT[i],NoLoadConsumption[i],StartupCost[i],MarginalCost[i],gen["unit_on_t0"],gen["power_output_t0"]));
	end
	if Printer
		println("Example Guthub ($(nb_gen) generators and $(T_max) hours)");
		println("MinRunCapacity : $(MinRunCapacity)");
		println();
		println("MaxRunCapacity : $(MaxRunCapacity)");
		println();
		println("RampUp : $(RampUp)");
		println();
		println("RampDown : $(RampDown)");
		println();
		println("UT : $(UT)");
		println();
		println("DT : $(DT)");
		println();
		println("SU : $(SU)");
		println();
		println("SD : $(SD)");
		println();
		println("StartupCost : $(StartupCost)");
		println();
		println("MarginalCost : $(MarginalCost)");
		println();
		println("L : $(L)");
		println();
		println("df_Generators : $(df_Generators)");
		println("df_Demand : $(df_Demand)");
	end
	return (nb_gen, T_max, MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, StartupCost, MarginalCost, NoLoadConsumption, L, VOLL);
end

function Process_Example_GitHub(filename, save_path, Printer=false, Build_Dataframe=false)
	#data_file = "GitHub_Example/rts_gmlc/2020-02-09.json"
	# data_file = "D:/Documents/Mémoire Papavasiliou/CHP_code/GitHub_Example/rts_gmlc/2020-01-27.json"
	# data_file = "D:/Documents/Mémoire Papavasiliou/CHP_code/GitHub_Example/rts_gmlc/$(filename)"
	data = JSON.parsefile(filename);
	thermal_gens = keys(data["thermal_generators"]);
	# renewable_gens = keys(data["renewable_generators"]);
	time_periods = 1:data["time_periods"]

	nb_gen = length(thermal_gens);
	T_max = data["time_periods"];
	MinRunCapacity = zeros(nb_gen);
	MaxRunCapacity = zeros(nb_gen);
	RampUp = zeros(nb_gen);
	RampDown = zeros(nb_gen);
	UT = zeros(nb_gen);
	DT = zeros(nb_gen);
	SU = zeros(nb_gen);
	SD = zeros(nb_gen);
	StartupCost = zeros(nb_gen);
	MarginalCost = zeros(nb_gen);
	NoLoadConsumption = zeros(nb_gen);
	L = zeros(T_max);
	VOLL = 3000;
	if Build_Dataframe
		df_Generators = DataFrame(A = String[], B = String[], C = String[], D = Float64[],E = Float64[],F = Float64[],G = Float64[],H = Float64[],I = Float64[],J = Float64[],K = Float64[],L = Float64[],M = Float64[],N = Float64[]);
		df_Demand = DataFrame();
	end
	for t in time_periods
		L[t] = data["demand"][t];
		if Build_Dataframe
			colname="hour$(t)"
			df_Demand[!,colname] = [L[t]];
		end
	end

	for (i,(g,gen)) in enumerate(data["thermal_generators"])
		# println("gen : $(gen)");
		MinRunCapacity[i] = gen["power_output_minimum"];
		MaxRunCapacity[i] = gen["power_output_maximum"];
		RampUp[i] = gen["ramp_up_limit"];
		RampDown[i] = gen["ramp_down_limit"];
		UT[i] = gen["time_up_minimum"];
		DT[i] = gen["time_down_minimum"];
		SU[i] = gen["ramp_startup_limit"];
		SD[i] = gen["ramp_shutdown_limit"];
		#println("gen[startup] : $(values(gen["startup"]))");
		count = 0.0;
		index = 0.0;
		# piecewise_production['cost']
		for (k,gen_startup) in enumerate(gen["startup"])
			count+=gen_startup["cost"];
			index+=1.0;
		end
		count = count/index;
		StartupCost[i] = count;
		count = 0.0;
		index = 0.0;
		for (k,gen_piecewise_production) in enumerate(gen["piecewise_production"])
			count+=gen_piecewise_production["cost"];
			index+=1.0;
		end
		count = count/index;
		MarginalCost[i] = count;
		if Build_Dataframe
			push!(df_Generators, ("$(i)", "SLOW", "Belgium",MinRunCapacity[i],MaxRunCapacity[i],RampUp[i],RampDown[i],UT[i],DT[i],NoLoadConsumption[i],StartupCost[i],MarginalCost[i],gen["unit_on_t0"],gen["power_output_t0"]));
		end
	end
	if Printer
		println("Example Guthub ($(nb_gen) generators and $(T_max) hours)");
		println("MinRunCapacity : $(MinRunCapacity)");
		println();
		println("MaxRunCapacity : $(MaxRunCapacity)");
		println();
		println("RampUp : $(RampUp)");
		println();
		println("RampDown : $(RampDown)");
		println();
		println("UT : $(UT)");
		println();
		println("DT : $(DT)");
		println();
		println("SU : $(SU)");
		println();
		println("SD : $(SD)");
		println();
		println("StartupCost : $(StartupCost)");
		println();
		println("MarginalCost : $(MarginalCost)");
		println();
		println("L : $(L)");
		println();
		println("df_Generators : $(df_Generators)");
		println("df_Demand : $(df_Demand)");
	end
	if Build_Dataframe
		save_path_generators = "$(save_path)/Generators.csv";
		save_path_demand = "$(save_path)/Demand.csv";
		touch(save_path_generators);
		touch(save_path_demand);
		CSV.write(save_path_generators,  df_Generators, header=true);
		CSV.write(save_path_demand,  df_Demand, header=true);
	end
	return (nb_gen, T_max, MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, StartupCost, MarginalCost, NoLoadConsumption, L, VOLL);
end

function pre_process_data()
	filename_vec = ["2020-01-27.json", "2020-02-09.json", "2020-03-05.json", "2020-04-03.json", "2020-05-05.json", "2020-06-09.json", "2020-07-06.json", "2020-08-12.json", "2020-09-20.json", "2020-10-27.json", "2020-11-25.json", "2020-12-23.json"];
	#filename = "2020-01-27.json";
	for filename in filename_vec
		println("filename : $(filename)")
		file_n = filename[1:10];
		f = "D:/Documents/Mémoire Papavasiliou/CHP_code/GitHub_Example/rts_gmlc/$(filename)";
		save_path = "D:/Documents/Mémoire Papavasiliou/Code - Nicolas Stevens/chp_dw_maxime/data/rts_gmlc_$(file_n)";
		Process_Example_GitHub(f, save_path, false);
	end
end
# pre_process_data()

function Run_Example_GitHub_rts_gmlc(save_results=0)
	println("=====================================================================");
	println("RTS_CMLC");println();
	save_path = "D:/Documents/Mémoire Papavasiliou/CHP_code/output/rts_gmlc/"

	filename_vec = ["2020-01-27.json", "2020-02-09.json", "2020-03-05.json", "2020-04-03.json", "2020-05-05.json", "2020-06-09.json", "2020-07-06.json", "2020-08-12.json", "2020-09-20.json", "2020-10-27.json", "2020-11-25.json", "2020-12-23.json"]
	#filename_vec = ["2020-01-27.json"];


	ST_master_mean_CG = 0;
	ST_slave_mean_CG = 0;
	ST_master_mean_CRG = 0;
	ST_slave_mean_CRG = 0;
	res = Dict();
	for (ind,filename) in enumerate(filename_vec)
		f = "D:/Documents/Mémoire Papavasiliou/CHP_code/GitHub_Example/rts_gmlc/$(filename)"
		println("filename : $(filename)");
		(nb_gen, T_max, MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT_1, DT_1, SU, SD, StartupCost, MarginalCost, NoLoadConsumption, data_demand, VOLL) = Process_Example_GitHub_not_saving(f, false);
		data_demand = data_demand/2;
		UT = check_UT_DT(UT_1, T_max);
		DT = check_UT_DT(DT_1, T_max);

		#nb_gen = 5;

		u_prior = zeros(nb_gen);
		v_prior = zeros(nb_gen);
		w_prior = zeros(nb_gen);
		p_prior = zeros(nb_gen);
		pbar_prior = zeros(nb_gen);

		u_posterior = zeros(nb_gen);
		v_posterior = zeros(nb_gen);
		w_posterior = zeros(nb_gen);
		p_posterior = zeros(nb_gen);
		pbar_posterior = zeros(nb_gen);

		Printer = false;
		G_c = Int64[g for g=1:nb_gen]
		# gather_gen, indices_gen = Gather_Generators(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, G_c);
		gather_gen = Int64[i for i=1:nb_gen]; indices_gen = Int64[i for i=1:nb_gen];
		delta_criterion = 0;

		#=
		# (iter_max, price_BD, obj_vec_BD, p_BD, pbar_BD, u_BD, v_BD, w_BD, l_BD, obj_real_BD, Solving_time, nb_cuts_tot, nb_salves_solve) = Bender_Decomposition(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MarginalCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, gather_gen, indices_gen, G_c, delta_criterion, Printer);
		println("Bender_Decomposition_1 total time");
		started_time = time();
		(iter_max, price_opt, obj_vec, p_time_opt, pbar_time_opt, u_opt, v_opt, w_opt, l_BD, obj_real_BD, Solving_time_BD_1) = Bender_Decomposition_1(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MarginalCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
		elasped = time() - started_time;
		println("Total running time : $(elasped)");
		println("Solving_time_BD_1 : $(Solving_time_BD_1)");
		println("iter_max : $(iter_max)");
		println("price_opt : $(price_opt)");
		=#
		
		println("Column-and-row generation algorithm (hot-start model)");
		started_time = time();
		(_, _, _, _, _, price_restricted, _, obj_restricted, obj_pricing, Solving_time_CRG, obj_restricted_vector, obj_pricing_vector, price_vect, iter_stop, Solving_time_master, Solving_time_slave, nb_var_vec, nb_cons_vec, Solving_time_master_vec, Solving_time_slave_vec) = Column_And_Row_Generation(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MarginalCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer);
		elasped = time() - started_time;
		println("Total running time : $(elasped)");
		println("Solving Solving time (hot-start model) : ",Solving_time_CRG);
		println("iter_stop : $(iter_stop)");
		println("obj_restricted_vec = $(obj_restricted_vector);")
		#println("price_restricted : $(price_restricted)");
		println();
		ST_master_mean_CRG+=Solving_time_master;
		ST_slave_mean_CRG+=Solving_time_slave;


		delta_criterion = -10^(-5);
		started_time = time();
		(obj_master, price, iter_stop, Solving_time, Solving_time_master, Solving_time_slaves, obj_vect_CG, prices_iterates, nb_var_vec ,nb_cons_vec, Solving_time_master_vec, Solving_time_slave_vec) = Column_Generation_DW(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MarginalCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, delta_criterion,Printer);
		elapsed = time() - started_time;
		println("Total running time : $(elapsed)");
		println("Solving_time : $(Solving_time)");
		println("iter_stop : $(iter_stop)");
		println("obj_vect_CG = $(obj_vect_CG);");
		ST_master_mean_CG+=Solving_time_master;
		ST_slave_mean_CG+=Solving_time_slaves;

	end
	println("ST_master_mean_CRG = $(ST_master_mean_CRG/length(filename_vec))");
	println("ST_slave_mean_CRG = $(ST_slave_mean_CRG/length(filename_vec))");
	println("ST_master_mean_CG = $(ST_master_mean_CG/length(filename_vec))");
	println("ST_slave_mean_CG = $(ST_slave_mean_CG/length(filename_vec))");
	if save_results==1
		CSV.write("$(save_path)results_table.csv",res);
	end
end

# Run_Example_GitHub_rts_gmlc()


function Run_Example_GitHub_ferc()
	println("=====================================================================");
	println("FERC");println();
	# filename_vec = ["2015-01-01_hw.json", "2015-01-01_lw.json", "2015-02-01_hw.json", "2015-02-01_lw.json", "2015-03-01_hw.json", "2015-03-01_lw.json", "2015-04-01_hw.json", "2015-04-01_lw.json", "2015-05-01_hw.json", "2015-05-01_lw.json", "2015-06-01_hw.json", "2015-06-01_lw.json", "2015-07-01_hw.json", "2015-07-01_lw.json", "2015-08-01_hw.json", "2015-08-01_lw.json", "2015-09-01_hw.json", "2015-09-01_lw.json", "2015-10-01_hw.json", "2015-10-01_lw.json", "2015-11-02_hw.json", "2015-11-02_lw.json", "2015-12-01_hw.json", "2015-12-01_lw.json"];
	# filename_vec = ["2015-03-01_hw.json", "2015-03-01_lw.json", "2015-04-01_hw.json", "2015-04-01_lw.json", "2015-05-01_hw.json", "2015-05-01_lw.json", "2015-06-01_hw.json", "2015-06-01_lw.json", "2015-07-01_hw.json", "2015-07-01_lw.json", "2015-08-01_hw.json", "2015-08-01_lw.json"];
	filename_vec = ["2015-01-01_hw.json"];
	
	# filename_vec = ["2015-03-01_hw.json", "2015-03-01_lw.json", "2015-04-01_hw.json", "2015-04-01_lw.json", "2015-05-01_hw.json", "2015-05-01_lw.json"];
	# filename_vec = ["2015-06-01_hw.json", "2015-06-01_lw.json", "2015-07-01_hw.json", "2015-07-01_lw.json", "2015-08-01_hw.json", "2015-08-01_lw.json"];
	
	master_ST_mean_CG = 0;
	slave_ST_mean_CG = 0;
	master_ST_mean_CRG = 0;
	slave_ST_mean_CRG = 0;
	for (ind,filename) in enumerate(filename_vec)
		println("$(filename)");
		f = "D:/Documents/Mémoire Papavasiliou/CHP_code/GitHub_Example/ferc/$(filename)"
		(nb_gen, T_max, MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT_1, DT_1, SU, SD, StartupCost, MarginalCost, NoLoadConsumption, data_demand, VOLL) = Process_Example_GitHub_not_saving(f);

		UT = check_UT_DT(UT_1, T_max);
		DT = check_UT_DT(DT_1, T_max);
		data_demand = data_demand/12;
		
		#println(length(MaxRunCapacity));
		#println(length(data_demand));println();

		u_prior = zeros(nb_gen);
		v_prior = zeros(nb_gen);
		w_prior = zeros(nb_gen);
		p_prior = zeros(nb_gen);
		pbar_prior = zeros(nb_gen);

		u_posterior = zeros(nb_gen);
		v_posterior = zeros(nb_gen);
		w_posterior = zeros(nb_gen);
		p_posterior = zeros(nb_gen);
		pbar_posterior = zeros(nb_gen);

		Printer = false;
		G_c = Int64[g for g=1:nb_gen]
		# gather_gen, indices_gen = Gather_Generators(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, G_c);
		gather_gen = Int64[i for i=1:nb_gen]; indices_gen = Int64[i for i=1:nb_gen];
		delta_criterion = 0;

		#=
		# (iter_max, price_BD, obj_vec_BD, p_BD, pbar_BD, u_BD, v_BD, w_BD, l_BD, obj_real_BD, Solving_time, nb_cuts_tot, nb_salves_solve) = Bender_Decomposition(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MarginalCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, gather_gen, indices_gen, G_c, delta_criterion, Printer);
		println("");
		started_time = time();
		(iter_max, price_opt, obj_vec, p_time_opt, pbar_time_opt, u_opt, v_opt, w_opt, l_BD, obj_real_BD, Solving_time_BD_1) = Bender_Decomposition_1(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MarginalCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
		elasped = time() - started_time;
		println(" Bender_Decomposition_1 Total running time : $(elasped)");
		println("Solving_time_BD_1 : $(Solving_time_BD_1)");
		println("iter_max : $(iter_max)");
		println("price_opt : $(price_opt)");=#
		
		started_time = time();
		(_, _, _, _, _, price_restricted, _, obj_restricted, obj_pricing, Solving_time_CRG, obj_restricted_vec_CRG, obj_pricing_vector, prices_iterates, iter_stop_CRG, solving_time_master, solving_time_slave) = Column_And_Row_Generation(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MarginalCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer);
		elapsed_CRG = time() - started_time;
		master_ST_mean_CRG+=solving_time_master;
		slave_ST_mean_CRG+=solving_time_slave;

		delta_criterion = -10^(-5);
		started_time = time();
		(obj_master, price, iter_stop, Solving_time, solving_time_master, solving_time_slaves, obj_vect_CG) = Column_Generation_DW(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MarginalCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, delta_criterion, Printer);
		elapsed = time() - started_time;
		println("Total running time : $(elapsed)");
		println("Solving_time : $(Solving_time)");
		println("iter_stop : $(iter_stop)");
		println("obj_vect_CG = $(obj_vect_CG);");
		println();
		println("Total running time : $(elapsed_CRG)");
		println("Solving Solving time : ",Solving_time_CRG);
		println("iter_stop_CRG : $(iter_stop_CRG)");
		println("obj_restricted_vec_CRG = $(obj_restricted_vec_CRG);");
		println();
		println();
		master_ST_mean_CG+=solving_time_master;
		slave_ST_mean_CG+=solving_time_slaves;
	end
	println();
	println();
	println("master_ST_mean_CG : $(master_ST_mean_CG/length(filename_vec))");
	println("slave_ST_mean_CG : $(slave_ST_mean_CG/length(filename_vec))");
	println();
	println("master_ST_mean_CRG : $(master_ST_mean_CRG/length(filename_vec))");
	println("slave_ST_mean_CRG : $(slave_ST_mean_CRG/length(filename_vec))");
end

# Run_Example_GitHub_ferc()

function Export_Data_For_DW(save_data=1)
	df1 = CSV.read("D:/Documents/Mémoire Papavasiliou/Unit commitment Toy/UCDataBE/Demand.csv", DataFrame);
	Demand = convert(Matrix, df1);
	######## Demand ########
	Data_types = Demand[:,1];
	(a,b) = size(Demand)
	Demand = convert(Matrix{Float64}, Demand[:,2:b])
	AutumnWD = Demand[1,:];
	SpringWD = Demand[2,:];
	SummerWD = Demand[3,:];
	WinterWD = Demand[4,:];
	AutumnWE = Demand[5,:];
	SpringWE = Demand[6,:];
	SummerWE = Demand[7,:];
	WinterWE = Demand[8,:];

	T_max = 24;
	AutumnWD_24h_array = zeros(T_max);
	SpringWD_24h_array = zeros(T_max);
	SummerWD_24h_array = zeros(T_max);
	WinterWD_24h_array = zeros(T_max);
	AutumnWE_24h_array = zeros(T_max);
	SpringWE_24h_array = zeros(T_max);
	SummerWE_24h_array = zeros(T_max);
	WinterWE_24h_array = zeros(T_max);

	c = 1;
	for i=1:24
		AutumnWD_24h_array[i] = AutumnWD[c];
		SpringWD_24h_array[i] = SpringWD[c];
		SummerWD_24h_array[i] = SummerWD[c];
		WinterWD_24h_array[i] = WinterWD[c];
		AutumnWE_24h_array[i] = AutumnWE[c];
		SpringWE_24h_array[i] = SpringWE[c];
		SummerWE_24h_array[i] = SummerWE[c];
		WinterWE_24h_array[i] = WinterWE[c];
		c+=4;
	end
	demand_dict = Dict();
	demand_dict["AutumnWD_24h"] = AutumnWD_24h_array;
	demand_dict["SpringWD_24h"] = SpringWD_24h_array;
	demand_dict["SummerWD_24h"] = SummerWD_24h_array;
	demand_dict["WinterWD_24h"] = WinterWD_24h_array;
	demand_dict["AutumnWE_24h"] = AutumnWE_24h_array;
	demand_dict["SpringWE_24h"] = SpringWE_24h_array;
	demand_dict["SummerWE_24h"] = SummerWE_24h_array;
	demand_dict["WinterWE_24h"] = WinterWE_24h_array;
	println(demand_dict);
	if save_data==1
		# save_path = "D:/Documents/Mémoire Papavasiliou/CHP_code/Data_Nicolas_Stevens"
		save_path = "D:/Documents/Mémoire Papavasiliou/Code - Nicolas Stevens/chp_dw_maxime/data/UCDataBE_adapted_Maxime2"
		CSV.write("$save_path/Demand_to_do.csv", demand_dict);
	end
end

# Export_Data_For_DW();

function check_data()
	df1 = CSV.read("D:/Documents/Mémoire Papavasiliou/Unit commitment Toy/UCDataBE/Demand.csv", DataFrame);
	Demand = convert(Matrix, df1);
	df2 = CSV.read("D:/Documents/Mémoire Papavasiliou/Unit commitment Toy/UCDataBE/Generators.csv", DataFrame);
	Generators = convert(Matrix, df2);
	######## Demand ########
	Data_types = Demand[:,1];
	(a,b) = size(Demand)
	Demand = convert(Matrix{Float64}, Demand[:,2:b])
	AutumnWD = Demand[1,:];
	SpringWD = Demand[2,:];
	SummerWD = Demand[3,:];
	WinterWD = Demand[4,:];
	AutumnWE = Demand[5,:];
	SpringWE = Demand[6,:];
	SummerWE = Demand[7,:];
	WinterWE = Demand[8,:];

	######## Generators ########
	nb_gen = 68; # There are normally 68 generators

	(a,b) = size(Generators);
	dt = 1; # [hour]
	Generators = convert(Matrix{Float64}, Generators[1:nb_gen,2:b]);
	MinRunCapacity 		= Generators[:,1];								# [MW]
	MaxRunCapacity 		= Generators[:,2];							# [MW]  10*
	RampUp				= (60*dt)*Generators[:,3];						# [MW/min] must multiply by 15 to have [MW/15-min]
	RampDown 			= (60*dt)*Generators[:,4];						# [MW/min] must multiply by 15
	UT_1 				= convert(Array{Int64}, (1/dt)*Generators[:,5]);# [15 min]
	DT_1		 		= convert(Array{Int64}, (1/dt)*Generators[:,6]);# [15 min]
	NoLoadConsumption 	= Generators[:,7];								# [MW], NoLoadCost = NoLoadConsumption * MargCost
	StartupCost 		= Generators[:,8];								# [€]
	MargCost 			= dt*Generators[:,9];							# [€/MWh]
	VOLL = 3000;
	SU = MinRunCapacity; # MinRunCapacity MaxRunCapacity
	SD = MinRunCapacity; # MinRunCapacity MaxRunCapacity
	# T_max = 35; # There normally 96 time-steps  25
	println("RampUp : $(RampUp)");
	println("RampDown : $(RampDown)");
end
# check_data()


function Run_Column_and_Row_Generation()
	df1 = CSV.read("D:/Documents/Mémoire Papavasiliou/Unit commitment Toy/UCDataBE/Demand.csv", DataFrame);
	Demand = convert(Matrix, df1);
	df2 = CSV.read("D:/Documents/Mémoire Papavasiliou/Unit commitment Toy/UCDataBE/Generators.csv", DataFrame);
	Generators = convert(Matrix, df2);
	######## Demand ########
	Data_types = Demand[:,1];
	(a,b) = size(Demand)
	Demand = convert(Matrix{Float64}, Demand[:,2:b])
	AutumnWE = Demand[5,:];
	WinterWE = Demand[8,:];
	SpringWE = Demand[6,:];
	SummerWE = Demand[7,:];
	AutumnWD = Demand[1,:];
	WinterWD = Demand[4,:];
	SpringWD = Demand[2,:];
	SummerWD = Demand[3,:];

	######## Generators ########
	nb_gen = 68; # There are normally 68 generators

	(a,b) = size(Generators);
	#dt = 1; # [hour]
	dt = 0.25;
	Generators = convert(Matrix{Float64}, Generators[1:nb_gen,2:b]);
	MinRunCapacity 		= Generators[:,1];								# [MW]
	MaxRunCapacity 		= Generators[:,2];								# [MW]
	RampUp				= (60*dt)*Generators[:,3];						# [MW/min] must multiply by 15 to have [MW/15-min]
	RampDown 			= (60*dt)*Generators[:,4];						# [MW/min] must multiply by 15
	UT_1 				= convert(Array{Int64}, (1/dt)*Generators[:,5]);# [15 min]
	DT_1		 		= convert(Array{Int64}, (1/dt)*Generators[:,6]);# [15 min]
	NoLoadConsumption 	= Generators[:,7];								# [MW], NoLoadCost = NoLoadConsumption * MargCost
	StartupCost 		= Generators[:,8];								# [€]
	MargCost 			= dt*Generators[:,9];							# [€/MWh]
	SU = MinRunCapacity;
	SD = MinRunCapacity;
	


	VOLL = 3000;
	T_max = length(AutumnWE);

	UT = check_UT_DT(UT_1, T_max);
	DT = check_UT_DT(DT_1, T_max);
	################################################################################################
	################################### PRIOR AND POSTERIOR DATA ###################################
	################################################################################################
	u_prior 		= zeros(nb_gen);
	v_prior 		= zeros(nb_gen);
	w_prior 		= zeros(nb_gen);
	p_prior 		= zeros(nb_gen);
	pbar_prior 		= zeros(nb_gen);

	u_posterior 	= zeros(nb_gen);
	v_posterior 	= zeros(nb_gen);
	w_posterior 	= zeros(nb_gen);
	p_posterior 	= zeros(nb_gen);
	pbar_posterior 	= zeros(nb_gen);

	G_c = Int64[g for g=1:nb_gen];gather_gen = Int64[i for i=1:nb_gen]; indices_gen = Int64[i for i=1:nb_gen];
	#delta_criterion = 0;

	master_ST_mean_CG = 0;
	slave_ST_mean_CG = 0;
	master_ST_mean_CRG = 0;
	slave_ST_mean_CRG = 0;
	slave_ST_mean_RG = 0;
	master_ST_mean_RG = 0;
	Printer = false;
	demand_vec = [AutumnWE]; # [AutumnWD,SpringWD,SummerWD,WinterWD,AutumnWE,SpringWE,SummerWE,WinterWE]
	for (index,data_demand) in enumerate(demand_vec)
		#=demand = zeros(24);
		c = 1;
		for i=1:24
			demand[i] = data_demand[c];
			c+=4;
		end
		T_max = 24;=#

		demand = data_demand;
		
		
		println("index : $(index)");
		started_time = time();
		(_, _, _, _, _, price_restricted, _, obj_restricted, obj_pricing, Solving_time_CRG, obj_restricted_vector, obj_pricing_vector, price_vect, iter_stop, Solving_time_master, Solving_time_slave, nb_var_vec, nb_cons_vec, Solving_time_master_vec, Solving_time_slave_vec) = Column_And_Row_Generation(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer);
		elapsed = time() - started_time;
		println("Total running time = $(elapsed)");
		println("Solving Solving time = ",Solving_time_CRG);
		println("iter_stop = $(iter_stop)");
		println("obj_restricted_vector = $(obj_restricted_vector);");
		println("");println("");
		println("price_restricted : $(price_restricted)");
		println("");println("");
		master_ST_mean_CRG+=Solving_time_master;
		slave_ST_mean_CRG+=Solving_time_slave;

		#delta_criterion = -0.001;
		delta_criterion = -10^(-5);
		started_time = time();
		(obj_master, price, iter_stop, Solving_time, Solving_time_master, Solving_time_slaves, obj_vect_CG, prices_iterates, nb_var_vec ,nb_cons_vec, Solving_time_master_vec, Solving_time_slave_vec) = Column_Generation_DW(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, delta_criterion, Printer);
		elapsed = time() - started_time;
		println("Total running time = $(elapsed)");
		println("Solving_time = $(Solving_time)");
		println("iter_stop = $(iter_stop)");
		println("obj_vect_CG = $(obj_vect_CG);")
		println();
		master_ST_mean_CG+=Solving_time_master;
		slave_ST_mean_CG+=Solving_time_slaves;
		

		#=println("");
		started_time = time();
		(iter_max, price_opt, obj_vec, p_time_opt, pbar_time_opt, u_opt, v_opt, w_opt, l_BD, obj_real_BD, Solving_time_BD_1) = Bender_Decomposition_1(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
		elasped = time() - started_time;
		println(" Bender_Decomposition_1 Total running time : $(elasped)");
		println("Solving_time_BD_1 : $(Solving_time_BD_1)");
		println("iter_max : $(iter_max)");
		println("price_opt : $(price_opt)");=#

		#=started_time = time();
		(p_EF, pbar_EF, gamma_EF, p_time_EF, pbar_time_EF, u_EF, v_EF, w_EF, price_EF, l_EF, obj_real_EF, Solving_time_EF) = Extended_Formulation(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer);
		elapsed = time() - started_time;
		println("Extended Formulation time - Total running time : $(elapsed)");
		println("Solving_time_EF : $(Solving_time_EF)");=#
	end
	println("master_ST_mean_RG= $(master_ST_mean_RG/length(demand_vec))");
	println("slave_ST_mean_RG= $(slave_ST_mean_RG/length(demand_vec))");
	println();
	println("master_ST_mean_CG= $(master_ST_mean_CG/length(demand_vec))");
	println("slave_ST_mean_CG= $(slave_ST_mean_CG/length(demand_vec))");
	println();
	println("master_ST_mean_CRG= $(master_ST_mean_CRG/length(demand_vec))");
	println("slave_ST_mean_CRG= $(slave_ST_mean_CRG/length(demand_vec))");
end

# Run_Column_and_Row_Generation();


function Run_Column_and_Row_Generation_24()
	df1 = CSV.read("D:/Documents/Mémoire Papavasiliou/Unit commitment Toy/UCDataBE/Demand.csv", DataFrame);
	Demand = convert(Matrix, df1);
	df2 = CSV.read("D:/Documents/Mémoire Papavasiliou/Unit commitment Toy/UCDataBE/Generators.csv", DataFrame);
	Generators = convert(Matrix, df2);
	######## Demand ########
	Data_types = Demand[:,1];
	(a,b) = size(Demand)
	Demand = convert(Matrix{Float64}, Demand[:,2:b])
	AutumnWE = Demand[5,:];
	WinterWE = Demand[8,:];
	SpringWE = Demand[6,:];
	SummerWE = Demand[7,:];
	AutumnWD = Demand[1,:];
	WinterWD = Demand[4,:];
	SpringWD = Demand[2,:];
	SummerWD = Demand[3,:];

	######## Generators ########
	nb_gen = 68; # There are normally 68 generators

	(a,b) = size(Generators);
	dt = 1; # [hour]
	#dt = 0.25;
	Generators = convert(Matrix{Float64}, Generators[1:nb_gen,2:b]);
	MinRunCapacity 		= Generators[:,1];								# [MW]
	MaxRunCapacity 		= Generators[:,2];								# [MW]
	RampUp				= (60*dt)*Generators[:,3];						# [MW/min] must multiply by 15 to have [MW/15-min]
	RampDown 			= (60*dt)*Generators[:,4];						# [MW/min] must multiply by 15
	UT_1 				= convert(Array{Int64}, (1/dt)*Generators[:,5]);# [15 min]
	DT_1		 		= convert(Array{Int64}, (1/dt)*Generators[:,6]);# [15 min]
	NoLoadConsumption 	= Generators[:,7];								# [MW], NoLoadCost = NoLoadConsumption * MargCost
	StartupCost 		= Generators[:,8];								# [€]
	MargCost 			= dt*Generators[:,9];							# [€/MWh]
	SU = MinRunCapacity;
	SD = MinRunCapacity;
	


	VOLL = 3000;
	#T_max = length(AutumnWE);
	
	T_max = 24;
	UT = check_UT_DT(UT_1, T_max);
	DT = check_UT_DT(DT_1, T_max);
	################################################################################################
	################################### PRIOR AND POSTERIOR DATA ###################################
	################################################################################################
	u_prior 		= zeros(nb_gen);
	v_prior 		= zeros(nb_gen);
	w_prior 		= zeros(nb_gen);
	p_prior 		= zeros(nb_gen);
	pbar_prior 		= zeros(nb_gen);

	u_posterior 	= zeros(nb_gen);
	v_posterior 	= zeros(nb_gen);
	w_posterior 	= zeros(nb_gen);
	p_posterior 	= zeros(nb_gen);
	pbar_posterior 	= zeros(nb_gen);

	G_c = Int64[g for g=1:nb_gen];gather_gen = Int64[i for i=1:nb_gen]; indices_gen = Int64[i for i=1:nb_gen];
	#delta_criterion = 0;

	master_ST_mean_CG = 0;
	slave_ST_mean_CG = 0;
	master_ST_mean_CRG = 0;
	slave_ST_mean_CRG = 0;
	slave_ST_mean_RG = 0;
	master_ST_mean_RG = 0;
	Printer = false;

	nb_var_mean = 0;
	nb_con_mean = 0;

	demand_vec = [WinterWE]; # [AutumnWD,SpringWD,SummerWD,WinterWD,AutumnWE,SpringWE,SummerWE,WinterWE]
	for (index,data_demand) in enumerate(demand_vec)
		demand = zeros(24);
		c = 1;
		for i=1:24
			demand[i] = data_demand[c];
			c+=4;
		end
		T_max = 24;
		#demand = data_demand;
		
		#=println("index : $(index)");
		started_time = time();
		(_, _, _, _, _, price_restricted, _, obj_restricted, obj_pricing, Solving_time_CRG, obj_restricted_vector, obj_pricing_vector, price_vect, iter_stop, Solving_time_master, Solving_time_slave, nb_var_vec, nb_cons_vec, Solving_time_master_vec, Solving_time_slave_vec) = Column_And_Row_Generation(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer);
		elapsed = time() - started_time;
		println("Total running time = $(elapsed)");
		println("Solving Solving time = ",Solving_time_CRG);
		println("iter_stop = $(iter_stop)");
		println("obj_restricted_vector = $(obj_restricted_vector);");
		master_ST_mean_CRG+=Solving_time_master;
		slave_ST_mean_CRG+=Solving_time_slave;=#

		#=
		#delta_criterion = -0.001;
		delta_criterion = -10^(-5);
		started_time = time();
		(obj_master, price, iter_stop, Solving_time, Solving_time_master, Solving_time_slaves, obj_vect_CG, prices_iterates, nb_var_vec ,nb_cons_vec, Solving_time_master_vec, Solving_time_slave_vec) = Column_Generation_DW(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, delta_criterion, Printer);
		elapsed = time() - started_time;
		println("Total running time = $(elapsed)");
		println("Solving_time = $(Solving_time)");
		println("iter_stop = $(iter_stop)");
		println("obj_vect_CG = $(obj_vect_CG);")
		println();
		master_ST_mean_CG+=Solving_time_master;
		slave_ST_mean_CG+=Solving_time_slaves;=#
		
		println("");
		started_time = time();
		(iter_max, price_opt, obj_vec, p_time_opt, pbar_time_opt, u_opt, v_opt, w_opt, l_BD, obj_real_BD, Solving_time_BD_1) = Bender_Decomposition_1(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
		elasped = time() - started_time;
		println(" Bender_Decomposition_1 Total running time : $(elasped)");
		println("Solving_time_BD_1 : $(Solving_time_BD_1)");
		println("iter_max : $(iter_max)");
		println("price_opt : $(price_opt)");

		#=started_time = time();
		(p_EF, pbar_EF, gamma_EF, p_time_EF, pbar_time_EF, u_EF, v_EF, w_EF, price_EF, l_EF, obj_real_EF, Solving_time_EF,nb_var, nb_con) = Extended_Formulation(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer);
		elapsed = time() - started_time;
		println("Extended Formulation time - Total running time : $(elapsed)");
		println("Solving_time_EF : $(Solving_time_EF)");
		nb_var_mean+=nb_var;
		nb_con_mean+=nb_con;=#
	end
	println("master_ST_mean_RG= $(master_ST_mean_RG/length(demand_vec))");
	println("slave_ST_mean_RG= $(slave_ST_mean_RG/length(demand_vec))");
	println();
	println("master_ST_mean_CG= $(master_ST_mean_CG/length(demand_vec))");
	println("slave_ST_mean_CG= $(slave_ST_mean_CG/length(demand_vec))");
	println();
	println("master_ST_mean_CRG= $(master_ST_mean_CRG/length(demand_vec))");
	println("slave_ST_mean_CRG= $(slave_ST_mean_CRG/length(demand_vec))");
	println();
	println("nb_var_mean : $(nb_var_mean/length(demand_vec))");
	println("nb_con_mean : $(nb_con_mean/length(demand_vec))");
end

# Run_Column_and_Row_Generation_24();

function Column_And_Row_Generation_Different_Inititalization(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price_BD ,Printer=false)
	### Solve Matching and get "on-intervals" to initialize \bar{S}
	#(p_matching, pbar_matching, u_matching, v_matching, w_matching, l_matching, obj_real_matching, Solving_time) = Matching(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
	
	m_matching = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
	JuMP.set_optimizer_attribute(m_matching,"OutputFlag",0);
	@variable(m_matching, 0 <= p[g=1:nb_gen, t=0:T_max+1]);		# power output of each generator g at time period t
	@variable(m_matching, 0 <= pbar[g=1:nb_gen, t=0:T_max+1]);	# maximum power available of each generator g at time period t
	@variable(m_matching, u[g=1:nb_gen, t=0:T_max+1], Bin);		# commitment status of geenrator g at time period t
	@variable(m_matching, v[g=1:nb_gen, t=0:T_max+1], Bin); 	# startup status of generator g at time period t
	@variable(m_matching, w[g=1:nb_gen, t=0:T_max+1], Bin);		# shutdown status of generator g at time period t
	@variable(m_matching, 0 <= l[t=1:T_max]);
	@constraint(m_matching, [t=1:T_max], l[t] <= data_demand[t]);			# upper bound for elasticity of the demand
	### The following constraints described the feasible set for technical constraints, the relaxation of 3-bin space
	logical_constraint 	= @constraint(m_matching, [g=1:nb_gen, t=1:T_max+1], u[g,t] - u[g,t-1] == v[g,t] - w[g,t]);					# Constraint 2, logical constraint on the different status (commitment status, startup status, shutdown status)
	minimum_up_time 	= @constraint(m_matching, [g=1:nb_gen, t=UT[g]:T_max+1], sum(v[g,i]  for i=t-UT[g]+1:t) <= u[g,t]);			# Constraint 3, minimum up time constraints
	minimum_down_time 	= @constraint(m_matching, [g=1:nb_gen, t=DT[g]:T_max+1], sum(w[g,i]  for i=t-DT[g]+1:t) <= (1 - u[g,t])); 	# Constraint 4, minimum down time constraints
	generation_limits_1	= @constraint(m_matching, [g=1:nb_gen, t=0:T_max+1], MinRunCapacity[g]*u[g,t] <= p[g,t]);					# Constraint 5_1, generation limits constraints
	generation_limits_2 = @constraint(m_matching, [g=1:nb_gen, t=0:T_max+1], p[g,t] <= pbar[g,t]);									# Constraint 5_2, generation limits constraints
	generation_limits_3 = @constraint(m_matching, [g=1:nb_gen, t=0:T_max+1], pbar[g,t] <= MaxRunCapacity[g]*u[g,t]);				# Constraint 5_3, generation limits constraints
	ramp_up_constraint 	= @constraint(m_matching, [g=1:nb_gen, t=1:T_max+1], pbar[g,t] - p[g,t-1] <= RU[g]*u[g,t-1] + SU[g]*v[g,t]);# Constraint 6, ramp-up constraints and start-up mode
	ramp_down_consraint = @constraint(m_matching, [g=1:nb_gen, t=1:T_max+1], pbar[g,t-1] - p[g,t] <= RD[g]*u[g,t] + SD[g]*w[g,t]);	# Constraint 7, ramp-down constraints and shutdown mode
	#loads = @constraint(m_matching, [t=1:T_max], sum( p[g,t] for g=1:nb_gen) == l[t]);
	pricing_problem = PricingProblem(m_matching,p,pbar,u,v,w,l,logical_constraint,minimum_up_time,minimum_down_time,generation_limits_1,generation_limits_2,generation_limits_3,ramp_up_constraint,ramp_down_consraint);
	### Prior data
	if length(u_prior) > 0
		for g=1:nb_gen
			JuMP.fix(pricing_problem.u[g,0], u_prior[g]; force = true);
			JuMP.fix(pricing_problem.v[g,0], v_prior[g]; force = true);
			JuMP.fix(pricing_problem.w[g,0], w_prior[g]; force = true);
			JuMP.fix(pricing_problem.p[g,0], p_prior[g]; force = true);
			JuMP.fix(pricing_problem.pbar[g,0], pbar_prior[g]; force = true);
		end
	end
	### Posterior data
	if length(u_posterior) > 0
		for g=1:nb_gen
			JuMP.fix(pricing_problem.u[g,T_max+1],u_posterior[g]; force = true);
			JuMP.fix(pricing_problem.v[g,T_max+1],v_posterior[g]; force = true);
			JuMP.fix(pricing_problem.w[g,T_max+1],w_posterior[g]; force = true);
			JuMP.fix(pricing_problem.p[g,T_max+1],p_posterior[g]; force = true);
			JuMP.fix(pricing_problem.pbar[g,T_max+1],pbar_posterior[g]; force = true);
		end
	end

	#delete(pricing_problem.model, loads); # Only need the constraint for initialization. Don't need after because it will be dualize.
	### Objective function
	### Regular init
	#=
	@objective(pricing_problem.model, Min, sum( sum( NoLoadConsumption[g]*C[g]*pricing_problem.u[g,t] + F[g]*pricing_problem.v[g,t] + C[g]*pricing_problem.p[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum( pricing_problem.l[t] for t=1:T_max) );
	optimize!(pricing_problem.model);
	u_matching = value.(pricing_problem.u).data;
	delete(pricing_problem.model, loads); # Only need the constraint for initialization. Don't need after because it will be dualize.
	# println("u_matching : ",u_matching);
	if nb_gen>1
		u_matching = u_matching[:,2:T_max+1];
	else
		u_matching = u_matching[2:T_max+1];
	end	
	(A,B,nb_intervals_gen) = Compute_A_B(u_matching, nb_gen);
	if Printer
		println("[Column_And_Row_Generation] A : ",A);
		println("[Column_And_Row_Generation] B : ",B);
		println("[Column_And_Row_Generation] nb_intervals_gen : ",nb_intervals_gen);
	end
	=#
	
	### Price init
	price = copy(price_BD);
	@objective(pricing_problem.model, Min, sum( sum( NoLoadConsumption[g]*C[g]*pricing_problem.u[g,t] + F[g]*pricing_problem.v[g,t] + C[g]*pricing_problem.p[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum( pricing_problem.l[t] for t=1:T_max) - sum( price[t]*(sum( pricing_problem.p[g,t] for g=1:nb_gen) - pricing_problem.l[t]) for t=1:T_max) );
	optimize!(pricing_problem.model);
	u_pricing = value.(pricing_problem.u).data;
	if nb_gen>1
		u_pricing = u_pricing[:,2:T_max+1];
	else
		u_pricing = u_pricing[2:T_max+1];
	end	
	(A,B,nb_intervals_gen) = Compute_A_B(u_pricing, nb_gen);
	if Printer
		println("[Column_And_Row_Generation] A : ",A);
		println("[Column_And_Row_Generation] B : ",B);
		println("[Column_And_Row_Generation] nb_intervals_gen : ",nb_intervals_gen);
	end

	### Intervals init


	beta = 0;
	iter_max = 10; # 5000
	obj_restricted_vector = [];
	obj_pricing_vector = [];
	p_time_restricted = 0;
	pbar_time_restricted = 0;
	u_restricted = 0;
	v_restricted = 0;
	w_restricted = 0;
	price = 0;
	price_iterates = [];
	l_restricted = 0;
	obj_restricted = 0;
	iter_stop = 0;
	obj_pricing = 0;
	# phi = 1e-3;
	phi = 10^(-3);
	Solving_time = 0;
	for iter=1:iter_max
		println();
		println();
		if Printer
			println("[Column_And_Row_Generation] iteration ",iter);
		end
		iter_stop = iter;
		### Solve the restricted problem
		#=if iter==1
			price = copy(price_BD);
		else
			(p_restricted, pbar_restricted, gamma_restricted, p_time_restricted, pbar_time_restricted, u_restricted, v_restricted, w_restricted, price, l_restricted, obj_restricted, time) = Restricted_Extended_Formulation(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, A, B, nb_intervals_gen,u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
			push!(obj_restricted_vector, obj_restricted);
			push!(price_iterates, price);
			Solving_time += time;
		end=#
		(p_restricted, pbar_restricted, gamma_restricted, p_time_restricted, pbar_time_restricted, u_restricted, v_restricted, w_restricted, price, l_restricted, obj_restricted, time) = Restricted_Extended_Formulation(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, A, B, nb_intervals_gen,u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
		push!(obj_restricted_vector, obj_restricted);
		push!(price_iterates, price);
		Solving_time += time;

		### Solve the pricing problem
		@objective(pricing_problem.model, Min, sum( sum( NoLoadConsumption[g]*C[g]*pricing_problem.u[g,t] + F[g]*pricing_problem.v[g,t] + C[g]*pricing_problem.p[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum( pricing_problem.l[t] for t=1:T_max) - sum( price[t]*(sum( pricing_problem.p[g,t] for g=1:nb_gen) - pricing_problem.l[t]) for t=1:T_max) );
		optimize!(pricing_problem.model);
		Solving_time += MOI.get(pricing_problem.model, MOI.SolveTime());
		obj_pricing = objective_value(pricing_problem.model);
		push!(obj_pricing_vector, obj_pricing);
		u_pricing = value.(pricing_problem.u).data;

		if Printer
			println("[Column_And_Row_Generation] obj_restricted : ",obj_restricted);
			println("[Column_And_Row_Generation] obj_pricing    : ",obj_pricing);
			println("[Column_And_Row_Generation : while loop] price : $(price)");
			println("[Column_And_Row_Generation : while loop] p_opt : $(value.(pricing_problem.p))");
			#println(pricing_problem.model)
		end
		### Compute the Lagrangian dual bound
		if iter==1
			beta = obj_pricing;
		else
			beta = maximum([obj_pricing beta]);
		end
		#println("beta           : $(beta)");
		#println("obj_restricted : $(obj_restricted)");
		#println("obj_restricted <= beta+phi : $(obj_restricted) <= $(beta + phi) : $(obj_restricted <= beta + phi)");
		if obj_restricted <= beta + phi
			println();
			println();
			if Printer
				println("[Column_And_Row_Generation] OPTIMAL SOLUTION FOUND");
			end
			break; # STOP algorithm.
		end
		### Update current bundle
		if nb_gen>1
			u = u_pricing[:,2:T_max+1];
		else
			u = u_pricing[2:T_max+1];
		end
		
		#=if iter==1
			(A_new,B_new,nb_intervals_gen_new) = Compute_A_B(u, nb_gen);
		else
			A_new,B_new,nb_intervals_gen_new, A_added, B_added, nb_intervals_added = Updata_A_B(u, nb_gen, A, B, nb_intervals_gen);
		end=#
		A_new,B_new,nb_intervals_gen_new, A_added, B_added, nb_intervals_added = Updata_A_B(u, nb_gen, A, B, nb_intervals_gen);
		A = copy(A_new);
		B = copy(B_new);
		nb_intervals_gen = copy(nb_intervals_gen_new);
		if Printer
			println("[Column_And_Row_Generation : while loop] A : ",A);
			println("[Column_And_Row_Generation : while loop] B : ",B);
			println("[Column_And_Row_Generation : while loop] nb_intervals_gen : ",nb_intervals_gen);
			if iter>1
				println("A_added : ",A_added);
				println("B_added : ",B_added);
				println("nb_intervals_added : ",nb_intervals_added);
			end
		end
	end
	if Printer
		println("[Column_And_Row_Generation] Stop at iteration ",iter_stop);
	end
	return (p_time_restricted, pbar_time_restricted, u_restricted, v_restricted, w_restricted, price, l_restricted, obj_restricted, obj_pricing, Solving_time, obj_restricted_vector, obj_pricing_vector, price_iterates);
end


function Column_Generation_Different_Inititalization(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price_BD, Printer=false)
	p_schedules = Dict();
	u_schedules = Dict();
	v_schedules = Dict();
	w_schedules = Dict();
	cost_schedules = Dict();
	counter_schedules = Int64[0 for g=1:nb_gen];

	### Build Subproblems
	price = zeros(T_max);
	tab_subProblems = Array{SubProblem}(undef, nb_gen);
	count = 0;
	for g=1:nb_gen
		sub_prob = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
		JuMP.set_optimizer_attribute(sub_prob,"OutputFlag",0);
	
		@variable(sub_prob, p[t=0:T_max+1],lower_bound = 0);		# power output of each generator g at time period t
		@variable(sub_prob, pbar[t=0:T_max+1],lower_bound = 0);	# maximum power available of each generator g at time period t
		@variable(sub_prob, u[t=0:T_max+1], Bin);		# commitment status of geenrator g at time period t
		@variable(sub_prob, v[t=0:T_max+1], Bin); 	# startup status of generator g at time period t
		@variable(sub_prob, w[t=0:T_max+1], Bin);		# shutdown status of generator g at time period t
		@variable(sub_prob, cost);

	
		### The following constriants described the feasible set for technical constraints, the relaxation of 3-bin space
		@constraint(sub_prob, logical_constraint[t=1:T_max+1], u[t] - u[t-1] == v[t] - w[t]);					# Constraint 2, logical constraint on the different status (commitment status, startup status, shutdown status)
		@constraint(sub_prob, minimum_up_time[t=UT[g]:T_max+1], sum(v[i]  for i=t-UT[g]+1:t) <= u[t]);			# Constraint 3, minimum up time constraints
		@constraint(sub_prob, minimum_down_time[t=DT[g]:T_max+1], sum(w[i]  for i=t-DT[g]+1:t) <= (1 - u[t])); 	# Constraint 4, minimum down time constraints
		
		@constraint(sub_prob, generation_limits_1[t=0:T_max+1], MinRunCapacity[g]*u[t] <= p[t]);				# Constraint 5_1, generation limits constraints
		@constraint(sub_prob, generation_limits_2[t=0:T_max+1], p[t] <= pbar[t]);								# Constraint 5_2, generation limits constraints
		@constraint(sub_prob, generation_limits_3[t=0:T_max+1], pbar[t] <= MaxRunCapacity[g]*u[t]);				# Constraint 5_3, generation limits constraints
		@constraint(sub_prob, ramp_up_constraint[t=1:T_max+1], pbar[t] - p[t-1] <= RU[g]*u[t-1] + SU[g]*v[t]);	# Constraint 6, ramp-up constraints and start-up mode
		@constraint(sub_prob, ramp_down_consraint[t=1:T_max+1], pbar[t-1] - p[t] <= RD[g]*u[t] + SD[g]*w[t]);	# Constraint 7, ramp-down constraints and shutdown mode
		
		### Prior data
		if length(u_prior) > 0
			for g=1:nb_gen
				JuMP.fix(u[0], u_prior[g]; force = true);
				JuMP.fix(v[0], v_prior[g]; force = true);
				JuMP.fix(w[0], w_prior[g]; force = true);
				JuMP.fix(p[0], p_prior[g]; force = true);
				JuMP.fix(pbar[0], pbar_prior[g]; force = true);
			end
		end
		### Posterior data
		if length(u_posterior) > 0
			for g=1:nb_gen
				JuMP.fix(u[T_max+1],u_posterior[g]; force = true);
				JuMP.fix(v[T_max+1],v_posterior[g]; force = true);
				JuMP.fix(w[T_max+1],w_posterior[g]; force = true);
				JuMP.fix(p[T_max+1],p_posterior[g]; force = true);
				JuMP.fix(pbar[T_max+1],pbar_posterior[g]; force = true);
			end
		end
		@constraint(sub_prob, cost - sum( NoLoadConsumption[g]*C[g]*u[t] + F[g]*v[t] + C[g]*p[t] for t=1:T_max)==0);
		### Objective function
		@objective(sub_prob, Min, sum( NoLoadConsumption[g]*C[g]*u[t] + F[g]*v[t] + C[g]*p[t] for t=1:T_max) - sum( price[t]*p[t] for t=1:T_max));
		sub_problem = SubProblem(sub_prob, p, pbar, u, v, w, cost, logical_constraint, minimum_up_time, minimum_down_time, generation_limits_1, generation_limits_2, generation_limits_3, ramp_up_constraint, ramp_down_consraint);
		count+=1;
		tab_subProblems[count] = sub_problem;
	end
	
	### Initialisation (regular one)
	#=(p_matching, pbar_matching, u_matching, v_matching, w_matching, l_matching, obj_matching, solving_time_matching) = Matching(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
	for g=1:nb_gen
		counter_schedules[g]+=1;
		if nb_gen>1
			p_schedules[g,counter_schedules[g]] = p_matching[g,2:T_max+1];
			u_schedules[g,counter_schedules[g]] = u_matching[g,2:T_max+1];
			v_schedules[g,counter_schedules[g]] = v_matching[g,2:T_max+1];
			w_schedules[g,counter_schedules[g]] = w_matching[g,2:T_max+1];
			cost_schedules[g,counter_schedules[g]] = sum( NoLoadConsumption[g]*C[g]*u_matching[g,t] + F[g]*v_matching[g,t] + C[g]*p_matching[g,t] for t=2:T_max+1);
		else
			p_schedules[g,counter_schedules[g]] = p_matching[2:T_max+1];
			u_schedules[g,counter_schedules[g]] = u_matching[2:T_max+1];
			v_schedules[g,counter_schedules[g]] = v_matching[2:T_max+1];
			w_schedules[g,counter_schedules[g]] = w_matching[2:T_max+1];
			cost_schedules[g,counter_schedules[g]] = sum( NoLoadConsumption[g]*C[g]*u_matching[t] + F[g]*v_matching[t] + C[g]*p_matching[t] for t=2:T_max+1);
		end
	end=#
	
	### Initialization with prices
	price = copy(price_BD);
	for g in 1:nb_gen
		sub_problem = tab_subProblems[g];
		@objective(sub_problem.model, Min, sum( NoLoadConsumption[g]*C[g]*sub_problem.u[t] + F[g]*sub_problem.v[t] + C[g]*sub_problem.p[t] for t=1:T_max) - sum( price[t]*sub_problem.p[t] for t=1:T_max));
		optimize!(sub_problem.model);
		p_sub = value.(sub_problem.p).data;
		u_sub = value.(sub_problem.u).data;
		v_sub = value.(sub_problem.v).data;
		w_sub = value.(sub_problem.w).data;
		counter_schedules[g]+=1;
		p_schedules[g,counter_schedules[g]] = p_sub[2:T_max+1];
		u_schedules[g,counter_schedules[g]] = u_sub[2:T_max+1];
		v_schedules[g,counter_schedules[g]] = v_sub[2:T_max+1];
		w_schedules[g,counter_schedules[g]] = w_sub[2:T_max+1];
		cost_schedules[g,counter_schedules[g]] = sum( NoLoadConsumption[g]*C[g]*u_sub[t] + F[g]*v_sub[t] + C[g]*p_sub[t] for t=2:T_max+1);
	end

	if Printer
		println("Initial schedule considered.");
		println(p_schedules);
		println(u_schedules);
		println(v_schedules);
		println(w_schedules);
		println(cost_schedules);
	end

	### Iterations
	iter_max = 500; # 500
	eps = -0.001; # -10^(-5)
	stopping_criteria = 0;
	obj_master = 0;
	iter_stop = 0;
	Solving_time_master = 0;
	Solving_time_slaves = 0;
	obj_vec = [];

	for iter=1:iter_max
		iter_stop = iter;
		obj_master, price, pi_dual, sol_time_master = Restricted_Master_Program_Column_Generation(data_demand, counter_schedules, p_schedules, cost_schedules, nb_gen, T_max, VOLL);
		Solving_time_master += sol_time_master;
		push!(obj_vec, obj_master);
		if Printer
			println();println();
			println("iter : $(iter)");
			println("obj_master : $(obj_master)");
			println("price : $(price)");
		end
		stopping_criteria = 1;
		for g in 1:nb_gen
			sub_problem = tab_subProblems[g];
			@objective(sub_problem.model, Min, sum( NoLoadConsumption[g]*C[g]*sub_problem.u[t] + F[g]*sub_problem.v[t] + C[g]*sub_problem.p[t] for t=1:T_max) - sum( price[t]*sub_problem.p[t] for t=1:T_max));
			optimize!(sub_problem.model);
			Solving_time_slaves += MOI.get(sub_problem.model, MOI.SolveTime());
			if Printer
				println("g : $(g) -> obj = $(objective_value(sub_problem.model))");
			end
			#reduced_cost = ComputeReducedCost(objective_value(sub_problem.model), pi_dual[g]);
			reduced_cost = objective_value(sub_problem.model) - pi_dual[g];
			if reduced_cost < eps
				stopping_criteria = 0
				p_sub = value.(sub_problem.p).data;
				u_sub = value.(sub_problem.u).data;
				v_sub = value.(sub_problem.v).data;
				w_sub = value.(sub_problem.w).data;
				counter_schedules[g]+=1;
				p_schedules[g,counter_schedules[g]] = p_sub[2:T_max+1];
				u_schedules[g,counter_schedules[g]] = u_sub[2:T_max+1];
				v_schedules[g,counter_schedules[g]] = v_sub[2:T_max+1];
				w_schedules[g,counter_schedules[g]] = w_sub[2:T_max+1];
				cost_schedules[g,counter_schedules[g]] = sum( NoLoadConsumption[g]*C[g]*u_sub[t] + F[g]*v_sub[t] + C[g]*p_sub[t] for t=2:T_max+1);
			end
		end
		if Printer
			println(p_schedules);
			println(u_schedules);
			println(v_schedules);
			println(w_schedules);
			println(cost_schedules);
		end
		if stopping_criteria == 1 # no column were added.
			println("iter_stop : $(iter_stop)");
			break;
		end
	end
	return (obj_master, price, iter_stop, Solving_time_master+Solving_time_slaves, Solving_time_master, Solving_time_slaves, obj_vec);
end


function run_Column_And_Row_Generation_Different_Inititalization()
	df1 = CSV.read("D:/Documents/Mémoire Papavasiliou/Unit commitment Toy/UCDataBE/Demand.csv", DataFrame);
	Demand = convert(Matrix, df1);
	df2 = CSV.read("D:/Documents/Mémoire Papavasiliou/Unit commitment Toy/UCDataBE/Generators.csv", DataFrame);
	Generators = convert(Matrix, df2);
	######## Demand ########
	Data_types = Demand[:,1];
	(a,b) = size(Demand)
	Demand = convert(Matrix{Float64}, Demand[:,2:b])
	AutumnWD = Demand[1,:];
	SpringWD = Demand[2,:];
	SummerWD = Demand[3,:];
	WinterWD = Demand[4,:];
	AutumnWE = Demand[5,:];
	SpringWE = Demand[6,:];
	SummerWE = Demand[7,:];
	WinterWE = Demand[8,:];

	######## Generators ########
	nb_gen = 68; # There are normally 68 generators

	(a,b) = size(Generators);
	dt = 1; # [hour]
	Generators = convert(Matrix{Float64}, Generators[1:nb_gen,2:b]);
	MinRunCapacity 		= Generators[:,1];								# [MW]
	MaxRunCapacity 		= Generators[:,2];							# [MW]  10*
	RampUp				= (60*dt)*Generators[:,3];						# [MW/min] must multiply by 15 to have [MW/15-min]
	RampDown 			= (60*dt)*Generators[:,4];						# [MW/min] must multiply by 15
	UT_1 				= convert(Array{Int64}, (1/dt)*Generators[:,5]);# [15 min]
	DT_1		 		= convert(Array{Int64}, (1/dt)*Generators[:,6]);# [15 min]
	NoLoadConsumption 	= Generators[:,7];								# [MW], NoLoadCost = NoLoadConsumption * MargCost
	StartupCost 		= Generators[:,8];								# [€]
	MargCost 			= dt*Generators[:,9];							# [€/MWh]
	VOLL = 3000;
	SU = MinRunCapacity; # MinRunCapacity MaxRunCapacity
	SD = MinRunCapacity; # MinRunCapacity MaxRunCapacity
	# T_max = 35; # There normally 96 time-steps  25

	# Demand in hours
	data_demand_vec = WinterWE; # SpringWE
	#data_demand_vec = SpringWE;
	data_demand = zeros(24);
	c = 1;
	for i=1:24
		data_demand[i] = data_demand_vec[c];
		c+=4;
	end
	T_max = 24;

	# T_max = 35;
	#=T_max = 25;
	data_demand_vec = WinterWE;
	data_demand = data_demand_vec[1:T_max];=#
	u_prior 		= zeros(nb_gen);
	v_prior 		= zeros(nb_gen);
	w_prior 		= zeros(nb_gen);
	p_prior 		= zeros(nb_gen);
	pbar_prior 		= zeros(nb_gen);

	u_posterior 	= zeros(nb_gen);
	v_posterior 	= zeros(nb_gen);
	w_posterior 	= zeros(nb_gen);
	p_posterior 	= zeros(nb_gen);
	pbar_posterior 	= zeros(nb_gen);

	UT = check_UT_DT(UT_1, T_max);
	DT = check_UT_DT(DT_1, T_max);
	##########################################################################################
	##########################################################################################
	#=nb_gen = 1;
	MinRunCapacity = [6];
	MaxRunCapacity = [16];
	RampUp = [5];
	RampDown = [5];
	UT = [1];
	DT = [1];
	SU = [6];
	SD = [6];
	NoLoadConsumption 	= [10];
	StartupCost 		= [53];
	MargCost 			= [30];
	data_demand = [6 11 16 11];
	VOLL = 3000
	T_max = length(data_demand);=#
	##########################################################################################
	##########################################################################################
	#=nb_gen = 1;
	MinRunCapacity = [6];
	MaxRunCapacity = [16];
	RampUp = [5];
	RampDown = [5];
	UT = [1];
	DT = [1];
	SU = [6];
	SD = [6];
	NoLoadConsumption 	= [10];
	StartupCost 		= [53];
	MargCost 			= [30];
	data_demand = [6 11 16 11 15 24 10 12 4];
	VOLL = 3000
	T_max = length(data_demand);=#
	##########################################################################################
	##########################################################################################
	u_prior 		= zeros(nb_gen);
	v_prior 		= zeros(nb_gen);
	w_prior 		= zeros(nb_gen);
	p_prior 		= zeros(nb_gen);
	pbar_prior 		= zeros(nb_gen);

	u_posterior 	= zeros(nb_gen);
	v_posterior 	= zeros(nb_gen);
	w_posterior 	= zeros(nb_gen);
	p_posterior 	= zeros(nb_gen);
	pbar_posterior 	= zeros(nb_gen);

	delta_criterion = 0;
	Printer = false;
	G_c = Int64[g for g=1:nb_gen]
	# gather_gen, indices_gen = Gather_Generators(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, G_c);
	gather_gen = Int64[i for i=1:nb_gen]; indices_gen = Int64[i for i=1:nb_gen];
	delta_criterion = 0;

	#(iter_max, price_opt, obj_vec, _, _, _, _, _, _, obj_real_BD, Solving_time, _, _, prices_vect, _, _, _, _) = Bender_Decomposition(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, gather_gen, indices_gen, G_c, delta_criterion, Printer);
	#println("prices_vect : $(prices_vect)");
	#println("prices_vect : $(prices_vect[1,:])");
	
	#price_BD = [3000.0, 3000.0, 69.4517, 68.5876, 68.5876, 64.4349, 68.5876, 68.5876, 68.5876, 89.1048, 160.544, 3000.0, 3000.0, 3000.0, 154.245, 156.776, 88.765, 82.6683, 3000.0, 3000.0, 3000.0, 3000.0, 3000.0, 3000.0];
	#price_BD = 20*ones(24);
	price_BD = zeros(24);
	(p_time_restricted, pbar_time_restricted, u_restricted, v_restricted, w_restricted, price, l_restricted, obj_restricted_init, obj_pricing_init, Solving_time, obj_restricted_vector, obj_pricing_vector, price_iterates) = Column_And_Row_Generation_Different_Inititalization(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price_BD ,Printer)
	println("obj_restricted_vector : $(obj_restricted_vector)");
	println("obj_pricing_vector : $(obj_pricing_vector)");
	println();
	
	(obj_master, price, iter_stop, Solving_time, Solving_time_master, Solving_time_slaves, obj_vec) = Column_Generation_Different_Inititalization(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price_BD, Printer);
	println("obj_master : $(obj_master)");
	println("price : $(price)");
	println();

	(p_time_restricted, pbar_time_restricted, u_restricted, v_restricted, w_restricted, price, l_restricted, obj_restricted, obj_pricing, Solving_time, obj_restricted_vector, obj_pricing_vector, price_vect, iter_stop) = Column_And_Row_Generation(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, StartupCost, MargCost, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer);
	println("obj_restricted_init : $(obj_restricted_init)");
	println("obj_pricing_init    : $(obj_pricing_init)");
	println("obj_restricted      : $(obj_restricted)");
	println("obj_pricing         : $(obj_pricing)");
	println("obj_master          : $(obj_master)");
end

# run_Column_And_Row_Generation_Different_Inititalization()

function run_Column_Generation()
	#######################################################################
	#######################################################################
	#=nb_gen = 1;
	MinRunCapacity 		= [6];
	MaxRunCapacity 		= [16];
	RampUp				= [5];
	RampDown 			= [5];
	UT 					= [1];
	DT		 			= [1];
	NoLoadConsumption 	= [10];
	F 					= [53];
	C 					= [30];
	SU = MinRunCapacity;
	SD = MinRunCapacity;
	# Demand/Load
	VOLL = 3000;
	# data_demand = [6 11 16 11 13 15]; # 6 10 16 9 13 15 -- 6 10 16 11 13 15 -- 6 11 16 11 13 15
	data_demand = [6 11 16 11];
	T_max = length(data_demand);=#
	#######################################################################
	#######################################################################
	#=nb_gen = 1;
	MinRunCapacity 		= [10];
	MaxRunCapacity 		= [20];
	RampUp				= [10];
	RampDown 			= [5];
	UT 					= [1];
	DT		 			= [1];
	NoLoadConsumption 	= [1];
	F 					= [70];
	C 					= [50];
	SU = [10];
	SD = [15];
	# Demand/Load
	VOLL = 3000;
	# data_demand = [6 11 16 11 13 15]; # 6 10 16 9 13 15 -- 6 10 16 11 13 15 -- 6 11 16 11 13 15
	data_demand = [18.1187 2.29986 1.49833 0.00790886 15.3403 15.8105 11.7689 13.755 3.24349 6.00786];
	T_max = length(data_demand);=#
	#######################################################################
	#######################################################################
	#=nb_gen = 2;
	MinRunCapacity 		= [10 10];
	MaxRunCapacity 		= [20 20];
	RampUp				= [5 1];
	RampDown 			= [6 4];
	UT 					= [1 1];
	DT		 			= [1 1];
	NoLoadConsumption 	= [1 1];
	F 					= [5 5];
	C 					= [10 10];
	SU = [12 10];
	SD = [12 15];
	# Demand/Load
	VOLL = 3000;
	# data_demand = [6 11 16 11 13 15]; # 6 10 16 9 13 15 -- 6 10 16 11 13 15 -- 6 11 16 11 13 15
	data_demand = [22.5 22.5 22.5 30];
	T_max = length(data_demand);=#
	#######################################################################
	#######################################################################
	df1 = CSV.read("D:/Documents/Mémoire Papavasiliou/Unit commitment Toy/UCDataBE/Demand.csv", DataFrame);
	Demand = convert(Matrix, df1);
	df2 = CSV.read("D:/Documents/Mémoire Papavasiliou/Unit commitment Toy/UCDataBE/Generators.csv", DataFrame);
	Generators = convert(Matrix, df2);
	######## Demand ########
	Data_types = Demand[:,1];
	(a,b) = size(Demand)
	Demand = convert(Matrix{Float64}, Demand[:,2:b])
	AutumnWD = Demand[1,:];
	SpringWD = Demand[2,:];
	SummerWD = Demand[3,:];
	WinterWD = Demand[4,:];
	AutumnWE = Demand[5,:];
	SpringWE = Demand[6,:];
	SummerWE = Demand[7,:];
	WinterWE = Demand[8,:];

	######## Generators ########
	nb_gen = 68; # There are normally 68 generators

	(a,b) = size(Generators);
	dt = 1; # [hour]
	Generators = convert(Matrix{Float64}, Generators[1:nb_gen,2:b]);
	MinRunCapacity 		= Generators[:,1];								# [MW]
	MaxRunCapacity 		= Generators[:,2];							# [MW]  10*
	RampUp				= (60*dt)*Generators[:,3];						# [MW/min] must multiply by 15 to have [MW/15-min]
	RampDown 			= (60*dt)*Generators[:,4];						# [MW/min] must multiply by 15
	UT_1 				= convert(Array{Int64}, (1/dt)*Generators[:,5]);# [15 min]
	DT_1		 		= convert(Array{Int64}, (1/dt)*Generators[:,6]);# [15 min]
	NoLoadConsumption 	= Generators[:,7];								# [MW], NoLoadCost = NoLoadConsumption * MargCost
	F 					= Generators[:,8];								# [€]
	C		 			= dt*Generators[:,9];							# [€/MWh]
	VOLL = 3000;
	SU = MinRunCapacity; # MinRunCapacity MaxRunCapacity
	SD = MinRunCapacity; # MinRunCapacity MaxRunCapacity
	# T_max = 35; # There normally 96 time-steps  25

	# Demand in hours
	#data_demand = WinterWE; # SpringWE
	#T_max = length(data_demand);
	
	data_demand_vec = WinterWE; # SpringWE
	#data_demand_vec = SpringWE;
	data_demand = zeros(24);
	c = 1;
	for i=1:24
		data_demand[i] = data_demand_vec[c];
		c+=4;
	end
	T_max = 24;
	
	UT = check_UT_DT(UT_1, T_max);
	DT = check_UT_DT(DT_1, T_max);
	#######################################################################
	#######################################################################
	### Optimization programs
	u_prior 		= zeros(nb_gen);
	v_prior 		= zeros(nb_gen);
	w_prior 		= zeros(nb_gen);
	p_prior 		= zeros(nb_gen);
	pbar_prior 		= zeros(nb_gen);
	u_posterior 	= zeros(nb_gen);
	v_posterior 	= zeros(nb_gen);
	w_posterior 	= zeros(nb_gen);
	p_posterior 	= zeros(nb_gen);
	pbar_posterior 	= zeros(nb_gen);

	(p_matching, pbar_matching, u_matching, v_matching, w_matching, l_matching, obj_matching, Solving_time_matching) = Matching(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
	delta_criterion = -10^(-3);
	started_time = time();
	(obj_master, price, iter_stop, Solving_time_tot, Solving_time_master, Solving_time_slaves, obj_vect, prices_iterates) = Column_Generation_DW(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, delta_criterion,false);
	elapsed = time() - started_time;
	println("elapsed : $(elapsed)");
	println("Solving_time_tot : $(Solving_time_tot)")
	println("obj_master : $(obj_master)");
	println("price : $(price)");
	println("iter_stop : $(iter_stop)");
	println("uplifts : $(obj_matching - obj_master)");

	#=delta_criterion_vec = range(0;stop=10000,step=100);
	elapsed_vec = zeros(length(delta_criterion_vec));
	Solving_time_tot_vec = zeros(length(delta_criterion_vec));
	obj_master_vec = zeros(length(delta_criterion_vec));
	iter_stop_vec = zeros(length(delta_criterion_vec));
	uplifts_vec = zeros(length(delta_criterion_vec));

	for (i,d) in enumerate(delta_criterion_vec)
		delta_criterion = -d;
		started_time = time();
		(obj_master, price, iter_stop, Solving_time_tot, Solving_time_master, Solving_time_slaves, obj_vect, prices_iterates) = Column_Generation_DW(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, delta_criterion,false);
		elapsed = time() - started_time;
		#=println("elapsed : $(elapsed)");
		println("Solving_time_tot : $(Solving_time_tot)")
		println("obj_master : $(obj_master)");
		println("price : $(price)");
		println("iter_stop : $(iter_stop)");
		println("uplifts : $(obj_matching - obj_master)");=#
		elapsed_vec[i] = elapsed;
		Solving_time_tot_vec[i] = Solving_time_tot;
		obj_master_vec[i] = obj_master;
		iter_stop_vec[i] = iter_stop;
		#uplifts_vec[i] = obj_matching - obj_master;
		uplifts_vec[i] = Compute_Uplifts_Illustrative_Example(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price, p_matching, u_matching, v_matching, l_matching);
	end
	println("delta_criterion = $(delta_criterion_vec);");
	println("elapsed_vec = $(elapsed_vec);");
	println("Solving_time_tot_vec = $(Solving_time_tot_vec);");
	println("obj_master_vec = $(obj_master_vec);");
	println("iter_stop_vec = $(iter_stop_vec);");
	println("uplifts_vec = $(uplifts_vec);");=#
	
	#=(p_time_restricted, pbar_time_restricted, u_restricted, v_restricted, w_restricted, price, l_restricted, obj_restricted, obj_pricing, Solving_time, obj_restricted_vector, obj_pricing_vector, price_vect, iter_stop) = Column_And_Row_Generation(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, false);
	println("obj_restricted      : $(obj_restricted)");
	println("obj_pricing         : $(obj_pricing)");=#

end

# run_Column_Generation()

#######################################################################################################
#######################################################################################################
############################################ GRAPH METHODS ############################################
#######################################################################################################
#######################################################################################################

"""
Count the number of variables and constraints the extended formulation actually needs.
"""
function count_variables_constraints(UT, T_max, nb_gen)
	(A,B,nb_intervals_gen) = set_T(UT,T_max, nb_gen);

	nb_var = 2*sum(nb_intervals_gen)*T_max + sum(nb_intervals_gen) + 3*nb_gen*T_max; # number of p and pbar + gamma + (u,v,w)

	nb_constraints = 0;
	for g=1:nb_gen
		for i=1:nb_intervals_gen[g]
			nb_constraints += 2*(A[g,i]-1 + T_max - B[g,i]) + 5*(B[g,i] - A[g,i]) + 4*(B[g,i] - A[g,i] - 1);
		end
	end
	nb_constraints += 6*nb_gen*T_max + T_max;
	return nb_constraints,nb_var;
end


function count_nb_intervals(T_max)
	T = 1:T_max;
	UT = 1:2:10;
	y = zeros(length(UT),length(T));
	nb_var = zeros(length(UT),length(T));
	nb_con = zeros(length(UT),length(T));
	c_ut = 1;
	for ut in UT
		c = 1;
		for t in T
			if t<ut
				A = [0]; B = [0]; index_gen = [0];
			else
				(A,B,index_gen) = set_T(ut,t, 1);
			end
			# println("A : $(A)"); println("B : $(B)"); println("index_gen : $(index_gen)"); println("c : $(c) and c_ut : $(c_ut)");
			y[c_ut,c] = index_gen[1];
			nb_var[c_ut,c] = 3*index_gen[1]*t;
			nb_con[c_ut,c] = (2*t+3*(B[1]-A[1]+1)+4*(B[1]-A[1]))*index_gen[1];
			c+=1;
		end
		c_ut+=1;
	end
	# p1 = plot(T, y, title = "Number of intervals according to T",legend=:bottomleft,labels=["number of intervals"]);
	p1 = plot(T, transpose(y), title = "Number of intervals according to T",legend=:topleft,labels=["UT = $(UT[1])" "UT = $(UT[2])" "UT = $(UT[3])" "UT = $(UT[4])" "UT = $(UT[5])"],titlefontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=20,xtickfontsize=12,ytickfontsize=12);
	xlabel!(p1,"T");
	ylabel!(p1, "Number of intervals");
	display(p1);
	savefig(p1,"D:/Documents/Mémoire Papavasiliou/Master Thesis - Manuscript/Images/nb_intervals");

	p2 = plot(T, transpose(nb_var), title = "Number of variables according to T",legend=:topleft, labels=["UT = $(UT[1])" "UT = $(UT[2])" "UT = $(UT[3])" "UT = $(UT[4])" "UT = $(UT[5])"],titlefontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=20,xtickfontsize=12,ytickfontsize=12);
	xlabel!(p2,"T");
	ylabel!(p2, "Number of variables");
	display(p2);
	savefig(p2,"D:/Documents/Mémoire Papavasiliou/Master Thesis - Manuscript/Images/nb_variables");

	p3 = plot(T, transpose(nb_con), title = "Number of constraints according to T",legend=:topleft, labels=["UT = $(UT[1])" "UT = $(UT[2])" "UT = $(UT[3])" "UT = $(UT[4])" "UT = $(UT[5])"],titlefontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=20,xtickfontsize=12,ytickfontsize=12);
	xlabel!(p3,"T");
	ylabel!(p3, "Number of constraints");
	display(p3);
	savefig(p3,"D:/Documents/Mémoire Papavasiliou/Master Thesis - Manuscript/Images/nb_constraints");
end


"""
Plot graphs using the function count_variables_constraints.
"""
function graph_variables_constraints()
	nb_gen = 68;
	T_max = 200;
	nb_var = zeros(T_max);
	nb_cons = zeros(T_max);
	for t=1:T_max
		(cons,var) = count_variables_constraints(UT, t, nb_gen);
		nb_var[t] = var; nb_cons[t] = cons;
	end
	p1 = plot(1:length(nb_var), nb_var, title = "Number of variables according to T");
	xlabel!(p1,"T");
	ylabel!(p1, "Number of variables");
	display(p1)
	savefig(p1,"D:/Home/Desktop/Mémoire Papavasiliou/Ramping Polytope and Cut Generation for UCP/Image/Number variables and constraints/number_variables");

	p2 = plot(1:length(nb_cons), nb_cons, title = "Number of constraints according to T");
	xlabel!(p2,"T");
	ylabel!(p2, "Number of constraints");
	display(p2)
	savefig(p2,"D:/Home/Desktop/Mémoire Papavasiliou/Ramping Polytope and Cut Generation for UCP/Image/Number variables and constraints/number_constraints");
end


function Compute_Uplifts_Illustrative_Example(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price_to_use, p_matching, u_matching, v_matching, l_matching)
	received_profit_producer_EF = zeros(nb_gen); max_profit_producer_EF = zeros(nb_gen); uplift_producer_EF = zeros(nb_gen);
	max_profit_consumer_EF = 0;received_profit_consumer_EF = 0;uplift_consumer_EF = 0;
	for g=1:nb_gen
		(max_profit_EF, p_max_profit, pbar_max_profit, u_max_profit, v_max_profit, w_max_profit) = MaxProfit_Producer(MinRunCapacity[g], MaxRunCapacity[g], RampUp[g], RampDown[g], UT[g], DT[g], SU[g], SD[g], T_max, F[g], C[g], NoLoadConsumption[g], price_to_use, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, g);
		max_profit_producer_EF[g] = max_profit_EF;

		received_profit_producer_EF[g] = sum(price_to_use[t]*p_matching[g,t+1] - (NoLoadConsumption[g]*C[g]*u_matching[g,t+1] + F[g]*v_matching[g,t+1] + C[g]*p_matching[g,t+1]) for t=1:T_max);
		uplift_producer_EF[g] = max_profit_producer_EF[g] - received_profit_producer_EF[g];
	end
	# Computer maximum profit for Consumer with CHP from EF
	(max_profit_consumer, l_max_profit_consumer) = MaxProfit_Consumer(data_demand, T_max, price_to_use, VOLL);
	max_profit_consumer_EF = max_profit_consumer;
	received_profit_consumer_EF = sum( (VOLL - price_to_use[t])*l_matching[t] for t=1:T_max);
	uplift_consumer_EF = max_profit_consumer_EF - received_profit_consumer_EF;
	return sum(uplift_producer_EF) + uplift_consumer_EF;
end



function Small_Example_Manuscript(save_fig=0)
	### Data of the example
	# Generators
	nb_gen = 1;
	MinRunCapacity = [6];
	MaxRunCapacity = [16];
	RampUp = [5];
	RampDown = [5];
	UT = [1];
	DT = [1];
	SU = [6];
	SD = [6];
	C = [30];
	F = [53];
	NoLoadConsumption = [10];
	data_demand = [6 11 16 11];
	VOLL = 3000
	T_max = length(data_demand);
	# Demand/Load
	VOLL = 3000;
	# data_demand = [6 11 16 11 13 15]; # 6 10 16 9 13 15 -- 6 10 16 11 13 15 -- 6 11 16 11 13 15
	data_demand = [6 11 16 11];
	T_max = length(data_demand);
	
	### Optimization programs
	u_prior 		= zeros(nb_gen);
	v_prior 		= zeros(nb_gen);
	w_prior 		= zeros(nb_gen);
	p_prior 		= zeros(nb_gen);
	pbar_prior 		= zeros(nb_gen);
	u_posterior 	= zeros(nb_gen);
	v_posterior 	= zeros(nb_gen);
	w_posterior 	= zeros(nb_gen);
	p_posterior 	= zeros(nb_gen);
	pbar_posterior 	= zeros(nb_gen);

	Printer = false;
	G_c = Int64[g for g=1:nb_gen]
	# gather_gen, indices_gen = Gather_Generators(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, G_c);
	gather_gen = Int64[i for i=1:nb_gen]; indices_gen = Int64[i for i=1:nb_gen];
	delta_criterion = 0;

	#=(iter_max, price_BD, obj_vec_BD, p_BD, pbar_BD, u_BD, v_BD, w_BD, l_BD, obj_real_BD, Solving_time, nb_cuts_tot, nb_salves_solve, prices_vect, p_vect) = Bender_Decomposition(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, gather_gen, indices_gen, G_c, delta_criterion, Printer);
	(p_time_restricted, pbar_time_restricted, u_restricted, v_restricted, w_restricted, price_restricted, l_restricted, obj_restricted, obj_pricing, Solving_time_CG, obj_restricted_vector, obj_pricing_vector, prices_vect_cr) = Column_And_Row_Generation(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer);
	(p_EF, pbar_EF,gamma_EF, p_time_EF, pbar_time_EF, u_EF, v_EF, w_EF, price_EF, l_EF, obj_EF, Solving_time_EF) = Extended_Formulation(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
	price_LP = [88.8333 45.4375 3000 3000];
	# price_EF = [];
	println();println();
	println("price_BD         = $(price_BD)");
	println("price_restricted = $(price_restricted)");
	println("price_EF         = $(price_EF)");
	println("prices_vect = $(prices_vect)");
	println("prices_vect_cr = $(prices_vect_cr)");
	println();
	println("obj_real_BD    : $(obj_real_BD)");
	println("obj_restricted : $(obj_restricted)");
	println("obj_pricing    : $(obj_pricing)");
	println("obj_EF         : $(obj_EF)");=#

	#price_BD         = [98.7121, 52.4545, 3000.0, 3000.0]
	#price_restricted = [38.3333, 80.0, 3000.0, 3000.0]
	#price_EF         = [30.0, 90.0, 3000.0, 3000.0]

	#=
	# price = price_BD;
	price = [30 30 3000 3000];
	(p_dual_lag, pbar_dual_lag, u_dual_lag, v_dual_lag, w_dual_lag, l_dual_lag, obj_dual_lag, Solving_time_dual_lag) = Dual_Lagrangian(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price)
	println("obj_dual_lag hour 1 : $(obj_dual_lag)");
	
	price = [30 84.5455 3000 3000];
	(p_dual_lag, pbar_dual_lag, u_dual_lag, v_dual_lag, w_dual_lag, l_dual_lag, obj_dual_lag, Solving_time_dual_lag) = Dual_Lagrangian(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price)
	println("obj_dual_lag hour 2 : $(obj_dual_lag)");

	price = [38.3333 80 3000 3000];
	(p_dual_lag, pbar_dual_lag, u_dual_lag, v_dual_lag, w_dual_lag, l_dual_lag, obj_dual_lag, Solving_time_dual_lag) = Dual_Lagrangian(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price)
	println("obj_dual_lag hour 3 : $(obj_dual_lag)");=#
	
	#=
	path = "Examples/Illustrative_Example_2"
	### Objective value according to iterations
	p_a = Plots.plot(1:length(obj_vec_BD),obj_vec_BD,marker=:circle,label="Row Generation",color="blue",xlabel="Iterations",ylabel="[€] (objective value)",legend=:bottomright,xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=15, titlefontsize=20);
	#Plots.plot!(p_a,1:length(obj_vec_BD),obj_vec_BD,linestyle=:dot,label=false,color="black");
	if save_fig==1
		Plots.savefig(p_a,"$(path)/row_generation_iterations.png");
	end
	p_b = Plots.plot(1:length(obj_restricted_vector),obj_restricted_vector,marker=:circle,label="Restricted Problem",color="red",xlabel="Iterations",ylabel="[€] (objective value)",legend=:bottomright,xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=15, titlefontsize=20);
	#Plots.plot!(p_b,1:length(obj_restricted_vector),obj_restricted_vector,linestyle=:dot,label=false,color="black");
	Plots.plot!(p_b,1:length(obj_pricing_vector),obj_pricing_vector,marker=:circle,color="blue",label="Pricing Problem");
	#Plots.plot!(p_b,1:length(obj_pricing_vector),obj_pricing_vector,linestyle=:dot,label=false,color="red");
	if save_fig==1
		Plots.savefig(p_b,"$(path)/column_and_row_generation_iterations.png");
	end
	
	### According to price of the first hour
	# For Row-Generation
	obj_RG_1 = zeros(length(prices_vect));
	obj_Row_Generation = zeros(length(prices_vect));
	price_RG_1 = zeros(length(prices_vect));
	for (i,price_RG) in enumerate(prices_vect)
		price_RG_1[i] = price_RG[1];
		# obj_RG_1[i] = obj_vec_BD[i];
		price = copy(price_BD);
		price[1] = price_RG[1];
		(_, _, _, _, _, _, obj_dual_lag, _) = Dual_Lagrangian(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price);
		obj_RG_1[i] = obj_dual_lag;
		obj_Row_Generation[i] = obj_vec_BD[i];
	end

	# For Column-and-Row Generation
	obj_pricing_CRG_1 = zeros(length(prices_vect_cr));
	obj_restricted_CRG_1 = zeros(length(prices_vect_cr));
	obj_DL_CRG = zeros(length(prices_vect_cr));
	price_CRG_1 = zeros(length(prices_vect_cr));
	for (i,price_CRG) in enumerate(prices_vect_cr)
		price_CRG_1[i] = price_CRG[1];
		price = copy(price_restricted);
		price[1] = price_CRG[1];
		(_, _, _, _, _, _, obj_dual_lag, _) = Dual_Lagrangian(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price);
		obj_DL_CRG[i] = obj_dual_lag;
		obj_pricing_CRG_1[i] = obj_pricing_vector[i];
		obj_restricted_CRG_1[i] = obj_restricted_vector[i];
	end

	# price_hour_1 = range(70; stop=100, step=0.01);
	price_hour_1 = range(70; stop=100, step=0.1);
	obj_dual_lag_vec = zeros(length(price_hour_1));
	for (j,p_h_1) in enumerate(price_hour_1)
		#price = [p_h_1 price_BD[2] price_BD[3] price_BD[4]];
		price = copy(price_BD);
		price[1] = p_h_1;
		(_, _, _, _, _, _, obj_dual_lag, _) = Dual_Lagrangian(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price);
		obj_dual_lag_vec[j] = obj_dual_lag;
	end

	# price_RCG_hour_1 = range(20; stop=50, step=0.01);
	price_RCG_hour_1 = range(20; stop=50, step=0.1);
	obj_dual_lag_RCG_vec = zeros(length(price_RCG_hour_1));
	for (j,p_h_1) in enumerate(price_RCG_hour_1)
		#price = [p_h_1 price_restricted[2] price_restricted[3] price_restricted[4]];
		price = copy(price_restricted);
		price[1] = p_h_1;
		(_, _, _, _, _, _, obj_dual_lag, _) = Dual_Lagrangian(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price);
		obj_dual_lag_RCG_vec[j] = obj_dual_lag;
	end

	pyplot();
	figure();
	p1 = Plots.scatter(price_RG_1, obj_RG_1,markersize = 7,title="Analysis of Dual Lagrangian for price hour 1 from Row Generation",label=["Projection of Dual Lagrangian"],xlabel="p [€/MWh] (Price hour 1)",ylabel="[€] (objective value)",xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=15, titlefontsize=20); # markershape=:star6
	Plots.plot!(p1,price_hour_1, obj_dual_lag_vec,label=["Dual Lagragian with price (p,52.4545,3000,3000)"],linestyle=:dash);
	Plots.scatter!(p1, price_RG_1,obj_Row_Generation, markersize = 10,label=["Master Program"],markershape=:x); # markershape=:xcross
	Plots.plot!(p1,price_RG_1,obj_Row_Generation,color="green", linestyle=:dot,label=false);
	display(p1);
	if save_fig==1
		Plots.savefig("$(path)/dual_lagrangian_prices_row_generation");
	end
	
	pyplot();
	fig = figure();
	p2 = Plots.scatter(price_CRG_1, obj_DL_CRG, markersize = 7, title="Analysis of Dual Lagrangian for price hour 1 from\nColumn-and-Row Generation", label=["Projection of Dual Lagrangian"],xlabel="p [€/MWh] (Price hour 1)",ylabel="Dual Lagrangian",xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=15, titlefontsize=20,legend= :outerbottom); # legend= :bottomright
	Plots.plot!(p2,price_RCG_hour_1, obj_dual_lag_RCG_vec,label=["Dual Lagragian with price (p,80,3000,3000)"],linestyle=:dash);
	Plots.scatter!(p2,price_CRG_1, obj_pricing_CRG_1, markersize = 13,markershape=:+ ,label=["Pricing Problem"]); # ,seriestypce=:scatter, markershape=:star6
	Plots.scatter!(p2,price_CRG_1, obj_restricted_CRG_1, markersize = 10,label=["Restricted Problem"], markershape=:x); # ,seriestypce=:scattermarkershape=:xcross
	Plots.plot!(p2,price_CRG_1, obj_pricing_CRG_1,color="green", linestyle=:dot,label=false);
	Plots.plot!(p2,price_CRG_1, obj_restricted_CRG_1,color="purple", linestyle=:dot,label=false);
	display(p2);
	if save_fig==1
		Plots.savefig("$(path)/dual_lagrangian_prices_column_and_row_generation.png");
	end


	### According to price of the second hour
	# For Row-Generation
	obj_RG_2 = zeros(length(prices_vect));
	obj_Row_Generation = zeros(length(prices_vect));
	price_RG_2 = zeros(length(prices_vect));
	for (i,price_RG) in enumerate(prices_vect)
		price_RG_2[i] = price_RG[2];
		price = copy(price_BD);
		price[2] = price_RG[2];
		(_, _, _, _, _, _, obj_dual_lag, _) = Dual_Lagrangian(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price);
		obj_RG_2[i] = obj_dual_lag;
		obj_Row_Generation[i] = obj_vec_BD[i];
	end

	# For Column-and-Row Generation
	obj_pricing_CRG_2 = zeros(length(prices_vect_cr));
	obj_restricted_CRG_2 = zeros(length(prices_vect_cr));
	obj_DL_CRG = zeros(length(prices_vect_cr));
	price_CRG_2 = zeros(length(prices_vect_cr));
	for (i,price_CRG) in enumerate(prices_vect_cr)
		price_CRG_2[i] = price_CRG[2];
		price = copy(price_restricted);
		price[2] = price_CRG[2];
		(_, _, _, _, _, _, obj_dual_lag, _) = Dual_Lagrangian(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price);
		obj_DL_CRG[i] = obj_dual_lag;
		obj_pricing_CRG_2[i] = obj_pricing_vector[i];
		obj_restricted_CRG_2[i] = obj_restricted_vector[i];
	end

	# price_hour_2 = range(40; stop=60, step=0.01); # y = 25:1:90;
	price_hour_2 = range(40; stop=60, step=0.1);
	obj_dual_lag_vec = zeros(length(price_hour_2));
	for (j,p_h_2) in enumerate(price_hour_2)
		price = copy(price_BD);
		price[2] = p_h_2;
		(_, _, _, _, _, _, obj_dual_lag, _) = Dual_Lagrangian(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price);
		obj_dual_lag_vec[j] = obj_dual_lag;
	end

	# price_RCG_hour_2 = range(20; stop=90, step=0.01);
	price_RCG_hour_2 = range(20; stop=90, step=0.1);
	obj_dual_lag_RCG_vec = zeros(length(price_RCG_hour_2));
	for (j,p_h_2) in enumerate(price_RCG_hour_2)
		price = copy(price_restricted);
		price[2] = p_h_2;
		(_, _, _, _, _, _, obj_dual_lag, _) = Dual_Lagrangian(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price);
		obj_dual_lag_RCG_vec[j] = obj_dual_lag;
	end

	pyplot();
	figure();
	p1 = Plots.scatter(price_RG_2, obj_RG_2,markersize = 7,title="Analysis of Dual Lagrangian for price hour 2 from Row Generation",label=["Projection of Dual Lagrangian"],xlabel="p [€/MWh] (Price hour 2)",ylabel="[€] (objective value)",xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=15, titlefontsize=20); # markershape=:star6
	Plots.plot!(p1,price_hour_2, obj_dual_lag_vec,label=["Dual Lagragian with price (98.7121,p,3000,3000)"],linestyle=:dash);
	Plots.scatter!(p1, price_RG_2,obj_Row_Generation, markersize = 10,label=["Master Program"],markershape=:x); # markershape=:xcross
	Plots.plot!(p1,price_RG_2,obj_Row_Generation,color="green", linestyle=:dot,label=false);
	display(p1);
	if save_fig==1
		Plots.savefig("$(path)/dual_lagrangian_prices_row_generation_hour2.png");
	end
	
	pyplot();
	fig = figure();
	p2 = Plots.scatter(price_CRG_2, obj_DL_CRG, markersize = 7, title="Analysis of Dual Lagrangian for price hour 2 from\nColumn-and-Row Generation", label=["Projection of Dual Lagrangian"],xlabel="p [€/MWh] (Price hour 2)",ylabel="Dual Lagrangian",xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=15, titlefontsize=20,legend= :outerbottom); # legend= :bottomright
	Plots.plot!(p2,price_RCG_hour_2, obj_dual_lag_RCG_vec,label=["Dual Lagragian with price (38.33,p,3000,3000)"],linestyle=:dash);
	Plots.scatter!(p2,price_CRG_2, obj_pricing_CRG_2, markersize = 13,markershape=:+ ,label=["Pricing Problem"]); # ,seriestypce=:scatter, markershape=:star6
	Plots.scatter!(p2,price_CRG_2, obj_restricted_CRG_2, markersize = 10,label=["Restricted Problem"], markershape=:x); # ,seriestypce=:scattermarkershape=:xcross
	Plots.plot!(p2,price_CRG_2, obj_pricing_CRG_2,color="green", linestyle=:dot,label=false);
	Plots.plot!(p2,price_CRG_2, obj_restricted_CRG_2,color="purple", linestyle=:dot,label=false);
	display(p2);
	if save_fig==1
		Plots.savefig("$(path)/dual_lagrangian_prices_column_and_row_generation_hour2.png");
	end
	=#
	
	### Contour Row-Generation and Column-and-Row Generation for Dual Lagrangian
	# Row Generation : 
	# Iteration 1 : [88.8333, 45.4375, 3000.0, 3000.0]
	# Iteration 2 : [98.7121, 52.4545, 3000.0, 3000.0]
	# Column-and-Row Generation
	# Iteration 1 : [30.0, 30.0, 3000.0, 3000.0]
	# Iteration 2 : [30.0, 84.5455, 3000.0, 3000.0]
	# Iteration 3 : [38.3333, 80.0, 3000.0, 3000.0]
	# Extended Formulation
	# [53 3000 3000 3000]
	# Continuous Linear Relaxation
	# [146.333 84.25 3000 3000]

	price_BD = [98.7121, 52.4545, 3000.0, 3000.0];
	prices_vect = [[88.8333, 45.4375, 3000.0, 3000.0],[98.7121, 52.4545, 3000.0, 3000.0]];
	price_EF = [30 90 3000 3000];
	price_LP = [88.8333 45.4375 3000 3000];
	prices_vect_cr = [[3000 -3123.36 3000 3000], [-2822.33 52.4545 3000 3000], [88.8333 52.4545 3000 3000]];
	prices_vect_CRG = [[30 30 3000 3000],[30 84.5455 3000 3000],[38.3333 80 3000 3000]];
	price_restricted = [88.8333 52.4545 3000 3000];

	f(x,y) = begin
		price = copy(price_BD);
		price[1] = x;
		price[2] = y;
		(_, _, _, _, _, _, obj_dual_lag, _) = Dual_Lagrangian(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price);
		obj_dual_lag;
	end

	price_RG_hour1 = zeros(length(prices_vect));
	price_RG_hour2 = zeros(length(prices_vect));
	obj_DL_RG = zeros(length(prices_vect));
	for (i,price_RG) in enumerate(prices_vect)
		price_RG_hour1[i] = price_RG[1];
		price_RG_hour2[i] = price_RG[2];
		price = copy(price_BD);
		price[1] = price_RG[1];
		price[2] = price_RG[2];
		#(_, _, _, _, _, _, obj_dual_lag, _) = Dual_Lagrangian(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price);
		#obj_DL_RG[i] = obj_dual_lag;
	end

	price_CRG_hour1 = zeros(length(prices_vect_cr));
	price_CRG_hour2 = zeros(length(prices_vect_cr));
	obj_DL_CRG = zeros(length(prices_vect_cr));
	for (i,price_CRG) in enumerate(prices_vect_cr)
		price_CRG_hour1[i] = price_CRG[1];
		price_CRG_hour2[i] = price_CRG[2];
		price = copy(price_restricted);
		price[1] = price_CRG[1];
		price[2] = price_CRG[2];
		#(_, _, _, _, _, _, obj_dual_lag, _) = Dual_Lagrangian(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price);
		#obj_DL_CRG[i] = obj_dual_lag;
	end

	x = 25:1:110; # price of hour 1
	y = 25:1:100; # price of hour 2
	#x = 60:1:90; # price of hour 1
	#y = 60:1:90; # price of hour 2
	X = repeat(reshape(x, 1, :), length(y), 1);
    Y = repeat(y, 1, length(x));
    Z = map(f, X, Y);
	
	### Contour Row-Generation and Column-and-Row Generation for Uplifts
	(p_matching, pbar_matching, u_matching, v_matching, w_matching, l_matching, obj_real_matching, Solving_time_matching) = Matching(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
	g(x,y) = begin
		price = copy(price_BD);
		price[1] = x;
		price[2] = y;
		uplift = Compute_Uplifts_Illustrative_Example(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price, p_matching, u_matching, v_matching, l_matching)
		uplift;
	end

	Z2 = map(g, X, Y);

	### Price EF and LP
	pyplot();
	figure();
	p4 = Plots.contour(x, y, Z,label="Dual Lagragian with prices from Row-Generation", title="Dual Lagragian according to Price of hour 1 and hour 2" ,fill=true, levels=-68050:10:-67350,xlabel="Price hour 1 [€/MWh]",ylabel="Price hour 2 [€/MWh]",xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,zguidefontsize=20,legendfontsize=15,titlefontsize=20);
	Plots.scatter!(p4, [price_EF[1]], [price_EF[2]], markershape=:circle,markersize=10,label="EF prices", color="blue");
	Plots.scatter!(p4, [price_LP[1]], [price_LP[2]], markershape=:rect,markersize=10,label="LP prices"); # , color="blue"
	display(p4);

	pyplot();
	figure();
	p5 = Plots.contour(x, y, Z2, title="Uplifts according to Price of hour 1 and hour 2" ,fill=true,xlabel="Price hour 1 [€/MWh]",ylabel="Price hour 2 [€/MWh]",xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,zguidefontsize=20,legendfontsize=15,titlefontsize=20,levels = 0:1:700);
	Plots.scatter!(p5, [price_EF[1]], [price_EF[2]], markershape=:circle,markersize=10,label="EF prices", color="blue");
	Plots.scatter!(p5, [price_LP[1]], [price_LP[2]], markershape=:rect,markersize=10,label="LP prices"); # , color="blue"
	display(p5);
	
	### prices LP, EF and RG
	pyplot();
	figure();
	#p3 = heatmap(x, y, f, alpha=0.3)
	p6 = Plots.contour(x, y, Z,label="Dual Lagragian with prices from Row-Generation", title="Dual Lagragian according to price of hour 1 and hour 2" ,fill=true, levels=-68050:10:-67350,xlabel="Price hour 1 [€/MWh]",ylabel="Price hour 2 [€/MWh]",xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,zguidefontsize=20,legendfontsize=15,titlefontsize=20); # , fill = true , levels = 100:1:101
	Plots.scatter!(p6, [price_EF[1]], [price_EF[2]], markershape=:circle,markersize=10,label="EF prices", color="blue");
	Plots.scatter!(p6, [price_LP[1]], [price_LP[2]], markershape=:rect,markersize=10,label="LP prices"); # , color="blue"
	Plots.scatter!(p6, price_RG_hour1, price_RG_hour2, markershape=:x,markersize=10,label="Row Generation iterations", color="black");
	#Plots.scatter!(p6, price_CRG_hour1, price_CRG_hour2, markershape=:+,markersize=10,label="Column-and-Row Generation iterations", color="red");
	Plots.plot!(p6, price_RG_hour1, price_RG_hour2,color="black",label=false,linestyle=:dot);
	#Plots.plot!(p3, price_CRG_hour1, price_CRG_hour2,color="red",label=false,linestyle=:dot);
	display(p6);
	if save_fig==1
		Plots.savefig("$(path)/dual_lagrangian_hour_1_2.png");
	end

	# Z2 = map(g, X, Y);
	pyplot();
	figure();
	p7 = Plots.contour(x, y, Z2, title="Uplifts according to price of hour 1 and hour 2" ,fill=true,xlabel="Price hour 1 [€/MWh]",ylabel="Price hour 2 [€/MWh]",xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,zguidefontsize=20,legendfontsize=15,titlefontsize=20,levels = 0:1:700); # , levels=-68050:10:-67350
	Plots.scatter!(p7, [price_EF[1]], [price_EF[2]], markershape=:circle,markersize=10,label="EF prices", color="blue");
	Plots.scatter!(p7, [price_LP[1]], [price_LP[2]], markershape=:rect,markersize=10,label="LP prices"); # , color="blue"
	Plots.scatter!(p7, price_RG_hour1, price_RG_hour2, markershape=:x,markersize=10,label="Row Generation iterations", color="white");
	#Plots.scatter!(p7, price_CRG_hour1, price_CRG_hour2, markershape=:+,markersize=10,label="Column-and-Row Generation iterations", color="red");
	Plots.plot!(p7, price_RG_hour1, price_RG_hour2,color="white",label=false,linestyle=:dot);
	#Plots.plot!(p4, price_CRG_hour1, price_CRG_hour2,color="red",label=false,linestyle=:dot);
	display(p7);
	if save_fig==1
		Plots.savefig("$(path)/uplifts_hour_1_2.png");
	end

	### prices LP, EF, RG and CG
	prices_CG_hour1 = Float64[0 0 0];
	prices_CG_hour2 = Float64[0 0 0];
	for (i,p) in enumerate(prices_vect_cr)
		prices_CG_hour1[i] = p[1];
		prices_CG_hour2[i] = p[2];
	end
	pyplot();
	figure();
	#p3 = heatmap(x, y, f, alpha=0.3)
	p8 = Plots.contour(x, y, Z,label="Dual Lagragian with prices from Row-Generation", title="Dual Lagragian according to price of hour 1 and hour 2" ,fill=true, levels=-68050:10:-67350,xlabel="Price hour 1 [€/MWh]",ylabel="Price hour 2 [€/MWh]",xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,zguidefontsize=20,legendfontsize=15,titlefontsize=20); # , fill = true , levels = 100:1:101
	Plots.scatter!(p8, [price_EF[1]], [price_EF[2]], markershape=:circle,markersize=10,label="EF prices", color="blue");
	Plots.scatter!(p8, [price_LP[1]], [price_LP[2]], markershape=:rect,markersize=10,label="LP prices"); # , color="blue"
	Plots.scatter!(p8, price_RG_hour1, price_RG_hour2, markershape=:x,markersize=10,label="Row Generation iterations", color="black");
	Plots.scatter!(p8, [prices_CG_hour1[3]], [prices_CG_hour2[3]], markershape=:star5,markersize=18,label="Column Generation iterations", color="orange");
	Plots.plot!(p8, price_RG_hour1, price_RG_hour2,color="black",label=false,linestyle=:dot);
	display(p8);
	if save_fig==1
		Plots.savefig("$(path)/dual_lagrangian_hour_1_2.png");
	end

	# Z2 = map(g, X, Y);
	pyplot();
	figure();
	p9 = Plots.contour(x, y, Z2, title="Uplifts according to price of hour 1 and hour 2" ,fill=true,xlabel="Price hour 1 [€/MWh]",ylabel="Price hour 2 [€/MWh]",xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,zguidefontsize=20,legendfontsize=15,titlefontsize=20,levels = 0:1:700); # , levels=-68050:10:-67350
	Plots.scatter!(p9, [price_EF[1]], [price_EF[2]], markershape=:circle,markersize=10,label="EF prices", color="blue");
	Plots.scatter!(p9, [price_LP[1]], [price_LP[2]], markershape=:rect,markersize=10,label="LP prices"); # , color="blue"
	Plots.scatter!(p9, price_RG_hour1, price_RG_hour2, markershape=:x,markersize=10,label="Row Generation iterations", color="white");
	Plots.scatter!(p9, [prices_CG_hour1[3]], [prices_CG_hour2[3]], markershape=:star5,markersize=18,label="Column Generation iterations", color="orange");
	Plots.plot!(p9, price_RG_hour1, price_RG_hour2,color="white",label=false,linestyle=:dot);
	display(p9);
	if save_fig==1
		Plots.savefig("$(path)/uplifts_hour_1_2.png");
	end

	### prices LP, EF, RG, CG and CRG
	# prices_CG = [[3000 -3123.36 3000 3000], [-2822.33 52.4545 3000 3000], [88.8333 52.4545 3000 3000]];
	prices_vect_CRG = [[30,30,3000,3000],[30,84.5455,3000,3000],[38.3333,80,3000,3000]];
	prices_CRG_hour1 = zeros(length(prices_vect_CRG));
	prices_CRG_hour2 = zeros(length(prices_vect_CRG));
	for (i,p) in enumerate(prices_vect_CRG)
		prices_CRG_hour1[i] = p[1];
		prices_CRG_hour2[i] = p[2];
	end
	println("prices_CRG_hour1 : $(prices_CRG_hour1)");println("prices_CRG_hour2 : $(prices_CRG_hour2)");
	pyplot();
	figure();
	#p3 = heatmap(x, y, f, alpha=0.3)
	p8 = Plots.contour(x, y, Z,label="Dual Lagragian with prices from Row-Generation", title="Dual Lagragian according to price of hour 1 and hour 2" ,fill=true, levels=-68050:10:-67350,xlabel="Price hour 1 [€/MWh]",ylabel="Price hour 2 [€/MWh]",xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,zguidefontsize=20,legendfontsize=15,titlefontsize=20); # , fill = true , levels = 100:1:101
	Plots.scatter!(p8, [price_EF[1]], [price_EF[2]], markershape=:circle,markersize=10,label="EF prices", color="blue");
	Plots.scatter!(p8, [price_LP[1]], [price_LP[2]], markershape=:rect,markersize=10,label="LP prices"); # , color="blue"
	Plots.scatter!(p8, price_RG_hour1, price_RG_hour2, markershape=:x,markersize=10,label="Row Generation iterations", color="black");
	Plots.scatter!(p8, [prices_CG_hour1[3]], [prices_CG_hour2[3]], markershape=:star5,markersize=18,label="Column Generation iterations", color="orange");
	Plots.scatter!(p8, prices_CRG_hour1, prices_CRG_hour2, markershape=:+,markersize=10,label="Column-and-Row Generation iterations", color="red");
	Plots.plot!(p8, price_RG_hour1, price_RG_hour2,color="black",label=false,linestyle=:dot);
	Plots.plot!(p8, prices_CRG_hour1, prices_CRG_hour2,color="red",label=false,linestyle=:dot);
	display(p8);
	if save_fig==1
		Plots.savefig("$(path)/dual_lagrangian_hour_1_2.png");
	end

	#=x = 25:1:110; # price of hour 1
	y = 25:1:100; # price of hour 2
	#x = 60:1:90; # price of hour 1
	#y = 60:1:90; # price of hour 2
	X = repeat(reshape(x, 1, :), length(y), 1);
    Y = repeat(y, 1, length(x));
	Z2 = map(g, X, Y);=#

	pyplot();
	figure();
	p10 = Plots.contour(x, y, Z2,label="Uplifts with prices from Row-Generation", title="Uplifts according to price of hour 1 and hour 2" ,fill=true, levels = 0:1:700,xlabel="Price hour 1 [€/MWh]",ylabel="Price hour 2 [€/MWh]",xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,zguidefontsize=20,legendfontsize=15,titlefontsize=20); # , fill = true , levels = 100:1:101
	Plots.scatter!(p10, [price_EF[1]], [price_EF[2]], markershape=:circle,markersize=10,label="EF prices", color="blue");
	Plots.scatter!(p10, [price_LP[1]], [price_LP[2]], markershape=:rect,markersize=10,label="LP prices"); # , color="blue"
	Plots.scatter!(p10, price_RG_hour1, price_RG_hour2, markershape=:x,markersize=10,label="Row Generation iterations", color="white");
	Plots.scatter!(p10, [prices_CG_hour1[3]], [prices_CG_hour2[3]], markershape=:star5,markersize=18,label="Column Generation iterations", color="orange");
	Plots.scatter!(p10, prices_CRG_hour1, prices_CRG_hour2, markershape=:+,markersize=10,label="Column-and-Row Generation iterations", color="red");
	Plots.plot!(p10, price_RG_hour1, price_RG_hour2,color="white",label=false,linestyle=:dot);
	Plots.plot!(p10, prices_CRG_hour1, prices_CRG_hour2,color="red",label=false,linestyle=:dot);
	display(p10);
end

# Small_Example_Manuscript()


function Smoke_Stack_Example_1_Difference_RG_CRG()
	### Data of the example
	# Generators
	nb_gen = 1;
	MinRunCapacity 		= [6];
	MaxRunCapacity 		= [16];
	RampUp				= [5];
	RampDown 			= [5];
	UT 					= [1];
	DT		 			= [1];
	NoLoadConsumption 	= [10];
	F 					= [53];
	C 					= [30];
	SU = MinRunCapacity;
	SD = MinRunCapacity;
	# Demand/Load
	VOLL = 3000;
	# data_demand = [6 11 16 11 13 15]; # 6 10 16 9 13 15 -- 6 10 16 11 13 15 -- 6 11 16 11 13 15
	data_demand = [6 11 16 11];
	T_max = length(data_demand);
	
	### Optimization programs
	u_prior 		= zeros(nb_gen);
	v_prior 		= zeros(nb_gen);
	w_prior 		= zeros(nb_gen);
	p_prior 		= zeros(nb_gen);
	pbar_prior 		= zeros(nb_gen);
	u_posterior 	= zeros(nb_gen);
	v_posterior 	= zeros(nb_gen);
	w_posterior 	= zeros(nb_gen);
	p_posterior 	= zeros(nb_gen);
	pbar_posterior 	= zeros(nb_gen);

	Printer = false;
	G_c = Int64[g for g=1:nb_gen]
	# gather_gen, indices_gen = Gather_Generators(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, G_c);
	gather_gen = Int64[i for i=1:nb_gen]; indices_gen = Int64[i for i=1:nb_gen];
	delta_criterion = 0;

	(iter_max, price_BD, obj_vec_BD, p_BD, pbar_BD, u_BD, v_BD, w_BD, l_BD, obj_real_BD, Solving_time, nb_cuts_tot, nb_salves_solve, prices_vect, p_vect, u_vect, v_vect, w_vect) = Bender_Decomposition(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, gather_gen, indices_gen, G_c, delta_criterion, Printer);
	(p_time_restricted, pbar_time_restricted, u_restricted, v_restricted, w_restricted, price_restricted, l_restricted, obj_restricted, obj_pricing, Solving_time_CG, obj_restricted_vector, obj_pricing_vector, prices_vect_cr) = Column_And_Row_Generation(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer);
	(p_matching, pbar_matching, u_matching, v_matching, w_matching, l_matching, obj_matching, Solving_time_matching) = Matching(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior)
	println();println();
	println("price_BD : $(price_BD)");
	println("price_restricted : $(price_restricted)");
	println("")
	println("p_vect : $(p_vect)");
	println("u_vect : $(u_vect)");
	println("v_vect : $(v_vect)");
	println("w_vect : $(w_vect)");
	println("")
	println();
	println("obj_real_BD    : $(obj_real_BD)");
	println("obj_restricted : $(obj_restricted)");
	println("obj_pricing    : $(obj_pricing)");

	#=
	pyplot();
	#figure();
	p1 = Plots.scatter(1:length(p_vect[1][2:end-1]),p_vect[1][2:end-1],label="Iter 1"); # ,linestyle=:dot
	for (i,p) in enumerate(p_vect[2:end])
		p = p[2:end-1]
		Plots.scatter!(p1,1:length(p),p,label="Iter $(i+1)");
	end
	display(p1);

	figure();
	p2 = Plots.scatter(1:length(u_vect[1][2:end-1]),u_vect[1][2:end-1],label="Iter 1"); # ,linestyle=:dot
	for (i,u) in enumerate(u_vect[2:end])
		u = u[2:end-1]
		Plots.scatter!(p2,1:length(u),u,label="Iter $(i+1)");
	end
	display(p2);

	figure();
	p3 = Plots.scatter(1:length(v_vect[1][2:end-1]),v_vect[1][2:end-1],label="Iter 1"); # ,linestyle=:dot
	for (i,v) in enumerate(v_vect[2:end])
		v = v[2:end-1]
		Plots.scatter!(p3,1:length(v),v,label="Iter $(i+1)");
	end
	display(p3);

	figure();
	p4 = Plots.scatter(1:length(w_vect[1][2:end-1]),w_vect[1][2:end-1],label="Iter 1"); # ,linestyle=:dot
	for (i,w) in enumerate(w_vect[2:end])
		w = w[2:end-1]
		Plots.scatter!(p4,1:length(w),w,label="Iter $(i+1)");
	end
	display(p4);
	=#
	println("obj_vec_BD : $(obj_vec_BD)");
	append!(obj_vec_BD, obj_matching);
	println("obj_vec_BD : $(obj_vec_BD)");
	println("obj_matching : $(obj_matching)")
	figure();

	
	p5 = Plots.bar([1],[abs(obj_vec_BD[2]-obj_vec_BD[1])],markersize=7,label="Iter 1");
	obj_vec_BD_look = obj_vec_BD[2:end];
	for (i,obj) in enumerate(obj_vec_BD_look)
		if i+1<=length(obj_vec_BD_look)
			Plots.bar!(p5,[i+1],[abs(obj_vec_BD_look[i+1]-obj_vec_BD_look[i])],markersize=7,label="Iter $(i+1)");
		else
			break
		end
	end
	display(p5);

	append!(obj_restricted_vector,obj_matching)
	println("obj_restricted_vector : $(obj_restricted_vector)");
	figure();
	p6 = Plots.bar([1],[abs(obj_restricted_vector[1] - obj_restricted_vector[2])],markersize=7,label="Iter 1");
	obj_restricted_vector_look = obj_restricted_vector[2:end];
	for (i,obj) in enumerate(obj_restricted_vector_look)
		if i+1<=length(obj_restricted_vector_look)
			Plots.bar!(p6,[i+1],[abs(obj_restricted_vector_look[i] - obj_restricted_vector_look[i+1])],markersize=7,label="Iter $(i+1)");
		else
			break
		end
	end
	display(p6);

end

# Smoke_Stack_Example_1_Difference_RG_CRG()

function Real_World_Example_Difference_RG_CRG(save_fig=0)
	
	### Data of the example
	df1 = CSV.read("D:/Documents/Mémoire Papavasiliou/Unit commitment Toy/UCDataBE/Demand.csv", DataFrame);
	Demand = convert(Matrix, df1);
	df2 = CSV.read("D:/Documents/Mémoire Papavasiliou/Unit commitment Toy/UCDataBE/Generators.csv", DataFrame);
	Generators = convert(Matrix, df2);
	######## Demand ########
	Data_types = Demand[:,1];
	(a,b) = size(Demand)
	Demand = convert(Matrix{Float64}, Demand[:,2:b]);
	AutumnWE = Demand[5,:];
	WinterWE = Demand[8,:];
	SpringWE = Demand[6,:];
	SummerWE = Demand[7,:];
	AutumnWD = Demand[1,:];
	WinterWD = Demand[4,:];
	SpringWD = Demand[2,:];
	SummerWD = Demand[3,:];
	######## Generators ########
	# nb_gen = 68; # There are normally 68 generators
	nb_gen = 68;

	(a,b) = size(Generators);
	dt = 1; # [hour]
	Generators = convert(Matrix{Float64}, Generators[1:nb_gen,2:b]);
	MinRunCapacity 		= Generators[:,1];								# [MW]
	MaxRunCapacity 		= Generators[:,2];
	RampUp				= (60*dt)*Generators[:,3];						# [MW/min] must multiply by 15 to have [MW/15-min]
	RampDown 			= (60*dt)*Generators[:,4];						# [MW/min] must multiply by 15
	UT_1 				= convert(Array{Int64}, (1/dt)*Generators[:,5]);# [15 min]
	DT_1		 		= convert(Array{Int64}, (1/dt)*Generators[:,6]);# [15 min]
	NoLoadConsumption 	= Generators[:,7];								# [MW], NoLoadCost = NoLoadConsumption * MargCost
	F 					= Generators[:,8];								# [€]
	C 					= dt*Generators[:,9];							# [€/MWh]
	SU = copy(MinRunCapacity);
	SD = copy(MinRunCapacity);

	VOLL = 3000;
	# Demand in hours
	data_demand_vec = WinterWE;
	data_demand = zeros(24);
	c = 1;
	for i=1:24
		data_demand[i] = data_demand_vec[c];
		c+=4;
	end
	T_max = 24;
	UT = check_UT_DT(UT_1, T_max);
	DT = check_UT_DT(DT_1, T_max);
	
	### Optimization programs
	u_prior 		= zeros(nb_gen);
	v_prior 		= zeros(nb_gen);
	w_prior 		= zeros(nb_gen);
	p_prior 		= zeros(nb_gen);
	pbar_prior 		= zeros(nb_gen);
	u_posterior 	= zeros(nb_gen);
	v_posterior 	= zeros(nb_gen);
	w_posterior 	= zeros(nb_gen);
	p_posterior 	= zeros(nb_gen);
	pbar_posterior 	= zeros(nb_gen);

	Printer = false;
	G_c = Int64[g for g=1:nb_gen];gather_gen = Int64[i for i=1:nb_gen]; indices_gen = Int64[i for i=1:nb_gen];
	println("G_c : $(G_c)");
	println("gather_gen : $(gather_gen)");
	println("indices_gen : $(indices_gen)");

	delta_criterion = 0;
	
	#=(iter_max, price_BD, obj_vec_BD, p_BD, pbar_BD, u_BD, v_BD, w_BD, l_BD, obj_real_BD, Solving_time, nb_cuts_tot, nb_salves_solve, prices_vect, p_vect, u_vect, v_vect, w_vect) = Bender_Decomposition(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, gather_gen, indices_gen, G_c, delta_criterion, Printer);
	(obj_master, price_CG, iter_stop, Solving_time, solving_time_master, solving_time_slaves, obj_vect, prices_iterates) = Column_Generation_DW(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer);
	(p_time_restricted, pbar_time_restricted, u_restricted, v_restricted, w_restricted, price_restricted, l_restricted, obj_restricted, obj_pricing, Solving_time_CG, obj_restricted_vector, obj_pricing_vector, prices_vect_cr) = Column_And_Row_Generation(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer);
	(p_matching, pbar_matching, u_matching, v_matching, w_matching, l_matching, obj_matching, Solving_time_matching) = Matching(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior)
	(p_EF, pbar_EF,gamma_EF, p_time_EF, pbar_time_EF, u_EF, v_EF, w_EF, price_EF, l_EF, obj_EF, Solving_time_EF) = Extended_Formulation(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
	println();println();
	println("price_BD : $(price_BD)");
	println("price_restricted : $(price_restricted)");
	println("price_CG : $(price_CG)");
	println("obj_real_BD    : $(obj_real_BD)");
	println("obj_restricted : $(obj_restricted)");
	println("obj_pricing    : $(obj_pricing)");
	println("obj_vec_BD     : $(obj_vec_BD)");
	println("obj_EF : $(obj_EF)");
	println("obj_vect : $(obj_vect)");=#


	obj_vec_BD = [-6.07436e8, -6.07434e8, -6.07431e8, -6.07415e8, -6.07386e8, -6.07377e8, -6.07358e8, -6.07339e8, -6.0733e8, -6.07323e8, -6.07315e8, -6.07305e8, -6.07302e8, -6.07296e8, -6.07294e8, -6.07291e8, -6.07289e8, -6.07286e8, -6.07281e8, -6.07275e8, -6.07272e8, -6.0727e8, -6.07268e8, -6.07263e8, -6.07259e8, -6.07256e8, -6.07254e8, -6.07251e8, -6.0725e8, -6.07246e8, -6.07245e8];
	obj_restricted_vector = [-6.07224e8, -6.07245e8, -6.07245e8];
	price_BD = [3000.0, 3000.0, 68.5876, 68.5876, 64.4349, 64.4349, 68.5876, 68.5876, 68.5876, 73.4634, 168.091, 3000.0, 3000.0, 3000.0, 244.542, 87.585, 69.4517, 69.4517, 3000.0, 3000.0, 3000.0, 3000.0, 3000.0, 3000.0]
	price_restricted = [3000.0, 3000.0, 68.5876, 68.5876, 64.4349, 64.4349, 68.5876, 68.5876, 68.5876, 73.4634, 168.091, 3000.0, 3000.0, 3000.0, 244.542, 87.585, 69.4517, 69.4517, 3000.0, 3000.0, 3000.0, 3000.0, 3000.0, 3000.0]
	price_CG = [3000.0, 3000.0, 68.5876, 68.5876, 64.4349, 64.4349, 68.5876, 68.5876, 68.5876, 73.4634, 168.091, 3000.0, 3000.0, 3000.0, 244.542, 87.585, 69.4517, 69.4517, 3000.0, 3000.0, 3000.0, 3000.0, 3000.0, 3000.0]
	obj_real_BD    = -6.072454773401365e8;
	obj_restricted = -6.072454773401365e8;
	obj_pricing    = -6.072454773401369e8;
	obj_EF = -6.072454773401365e8;
	obj_vect = [-6.07224e8, -6.07224e8, -6.07224e8, -6.07224e8, -6.07224e8, -6.07224e8, -6.07224e8, -6.07224e8, -6.07224e8, -6.07224e8, -6.07226e8, -6.07235e8, -6.07243e8, -6.07243e8, -6.07244e8, -6.07245e8, -6.07245e8, -6.07245e8, -6.07245e8, -6.07245e8, -6.07245e8, -6.07245e8, -6.07245e8, -6.07245e8, -6.07245e8, -6.07245e8, -6.07245e8, -6.07245e8, -6.07245e8];
	
	path = "D:/Documents/Mémoire Papavasiliou/Master Thesis - Manuscript/Temporary Images/Difference_RG_CRG"

	pyplot();
	pygui(true);
	p7 = Plots.plot(1:length(obj_vec_BD),obj_vec_BD,line=(true, 2), marker=:rect,title="Evolution Row Generation Method",ylabel="Objective value[€]",xlabel="Iteration",xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=20,titlefontsize=20,label=L"$RG$",legend=:bottomright,color=:blue);
	Plots.plot!(p7, 1:length(obj_vec_BD), obj_EF.*ones(length(obj_vec_BD)),linestyle=:dash,label=L"$EF$",color="green",linewidth=2);
	#Plots.plot!(p7, 1:length(obj_vec_BD), obj_matching.*ones(length(obj_vec_BD)),linestyle=:dot,label=false)
	#append!(obj_vec_BD, obj_matching);
	p5 = Plots.bar([1],[abs(obj_vec_BD[2]-obj_vec_BD[1])],markersize=7,label=false,legend=false,xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=20,xlabel="Iteration",ylabel="Enhancement objective value",color=:blue,bar_width=0.7);
	obj_vec_BD_look = obj_vec_BD[2:end];
	for (i,obj) in enumerate(obj_vec_BD_look)
		if i+1<=length(obj_vec_BD_look)
			Plots.bar!(p5,[i+1],[abs(obj_vec_BD_look[i+1]-obj_vec_BD_look[i])],markersize=7,label=false,color=:blue,bar_width=0.7);
		else
			break
		end
	end
	p8 = Plots.plot(p7, p5, layout=(2,1));
	if save_fig==1
		Plots.savefig(p8, "$(path)/row_generation.png");
	end
	figure();
	display(p8);
	
	p9 = Plots.plot(1:length(obj_restricted_vector),obj_restricted_vector,line=(true, 5), marker=:circle,title="Evolution Column-and-Row Generation Method",ylabel="Objective value[€]",xlabel="Iteration",xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=20,titlefontsize=20,legend=false);
	#Plots.plot!(p9, 1:length(obj_restricted_vector), obj_matching.*ones(length(obj_restricted)),linestyle=:dot,label=false);

	#append!(obj_restricted_vector,obj_matching)
	println("obj_restricted_vector : $(obj_restricted_vector)");
	p6 = Plots.bar([1],[abs(obj_restricted_vector[1] - obj_restricted_vector[2])],markersize=7,label=false,legend=false,xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=20,xlabel="Iteration",ylabel="Enhancement objective value");
	obj_restricted_vector_look = obj_restricted_vector[2:end];
	for (i,obj) in enumerate(obj_restricted_vector_look)
		if i+1<=length(obj_restricted_vector_look)
			Plots.bar!(p6,[i+1],[abs(obj_restricted_vector_look[i] - obj_restricted_vector_look[i+1])],markersize=7,label=false);
		else
			break
		end
	end
	p10 = Plots.plot(p9, p6, layout=(2,1));
	if save_fig==1
		Plots.savefig(p10, "$(path)/column_and_row_generation.png");
	end
	figure();
	display(p10);

	p11 = Plots.plot(1:length(obj_vec_BD), obj_vec_BD[end].*ones(length(obj_vec_BD)), linestyle=:dash,label="Max Dual Lagrangian",linewidth=2,xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=20,titlefontsize=20,color="green");
	Plots.plot!(p11,1:length(obj_vec_BD),obj_vec_BD,line=(true, 2), marker=:rect,title="Evolution Objective Value",ylabel="Objective value[€]",xlabel="Iteration",label="Row Generation",legend=:best,color=:red); # color=:red
	Plots.plot!(p11, 1:length(obj_restricted_vector),obj_restricted_vector,line=(true, 2), marker=:circle,ylabel="Objective value[€]",xlabel="Iteration",label="Column-and-Row Generation",color=:blue);
	if save_fig==1
		Plots.savefig(p11, "$(path)/both_generations.png");
	end

	obj_restricted_vector_2 = zeros(length(obj_vec_BD));
	obj_restricted_vector_2[1:length(obj_restricted_vector)] = obj_restricted_vector;
	
	#obj_restricted_vector_look_2 = zeros(length(obj_restricted_vector_2)-1);
	obj_look_2 = zeros(length(obj_restricted_vector_2)-1,2);
	println("obj_look_2 : $(size(obj_look_2))");
	ctg = repeat(["Row Generation", "Column-and-Row Generation"], inner=length(obj_restricted_vector_2)-1);
	for i=1:length(obj_vec_BD)-1
		obj_look_2[i,1] = abs(obj_vec_BD[i+1]-obj_vec_BD[i]);
		if i<=length(obj_restricted_vector)-1
			obj_look_2[i,2] = abs(obj_restricted_vector_2[i] - obj_restricted_vector_2[i+1]);
		end
	end
	x = repeat(1:length(obj_restricted_vector_2)-1,outer=2);
	println("x : $(size(x))");
	println("obj_look_2 : $(size(obj_look_2))");
	#p12 = StatsPlots.groupedbar(x,obj_look_2); # ,xlabel="Iterations" ,group=ctg, bar_width = 0.67,lw = 0, framestyle=:box
	p12 = StatsPlots.groupedbar(obj_look_2,bar_position=:dodge,bar_width=0.7,xlabel="Iteration",group=ctg,color=[:blue :red],xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=20,ylabel="Enhancement objective value");
	
	p13 = Plots.plot(p11, p12, layout=(2,1));
	figure();
	display(p13);

	p14 = Plots.plot(1:length(obj_vec_BD), obj_vec_BD[end].*ones(length(obj_vec_BD)), linestyle=:dash,label=L"$EF$",linewidth=2,xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=20,titlefontsize=20,color="green");
	Plots.plot!(p14,1:length(obj_vec_BD),obj_vec_BD,line=(true, 2),marker=:rect,title="Evolution Objective Value",ylabel="Objective value[€]",xlabel="Iteration",label="RG",legend=:best,color=:blue); # color=:red
	Plots.plot!(p14, 1:length(obj_vect),obj_vect,line=(true, 2),marker=:circle,ylabel="Objective value[€]",xlabel="Iteration",label=L"CG",color=:red);

	obj_vect_CG2 = zeros(length(obj_vec_BD));
	obj_vect_CG2[1:length(obj_vect)] = obj_vect;
	obj_look_2 = zeros(length(obj_vect_CG2)-1,2);
	ctg = repeat(["Row Generation", "Column Generation"], inner=length(obj_vect_CG2)-1);
	for i=1:length(obj_vec_BD)-1
		obj_look_2[i,1] = abs(obj_vec_BD[i+1]-obj_vec_BD[i]);
		if i<=length(obj_vect)-1
			obj_look_2[i,2] = abs(obj_vect_CG2[i] - obj_vect_CG2[i+1]);
		end
	end
	println("obj_look_2 : $(obj_look_2)");
	p15 = StatsPlots.groupedbar(obj_look_2,bar_position=:dodge,bar_width=0.7,xlabel="Iteration",group=ctg,color=[:red :blue],xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=20,ylabel="Enhancement objective value");
	
	p16 = Plots.plot(p14, p15, layout=(2,1));
	figure();
	display(p16);


	p14 = Plots.plot(1:length(obj_vec_BD), obj_vec_BD[end].*ones(length(obj_vec_BD)), linestyle=:dash,label=L"$EF$",linewidth=2,xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=20,titlefontsize=20,color="green");
	Plots.plot!(p14,1:length(obj_vec_BD),obj_vec_BD,line=(true, 2),marker=:rect,title="Evolution Objective Value",ylabel="Objective value[€]",xlabel="Iteration",label=L"RG",legend=:best,color=:blue); # color=:red
	Plots.plot!(p14, 1:length(obj_vect),obj_vect,line=(true, 2),marker=:circle,ylabel="Objective value[€]",xlabel="Iteration",label=L"CG",color=:red);

	obj_vect_CG2 = zeros(length(obj_vec_BD));
	obj_vect_CG2[1:length(obj_vect)] = obj_vect;
	obj_look_2 = zeros(length(obj_vect_CG2)-1,2);
	ctg = repeat(["Row Generation", "Column Generation"], inner=length(obj_vect_CG2)-1);
	for i=1:length(obj_vec_BD)-1
		obj_look_2[i,1] = abs(obj_vec_BD[i+1]-obj_vec_BD[i]);
		if i<=length(obj_vect)-1
			obj_look_2[i,2] = abs(obj_vect_CG2[i] - obj_vect_CG2[i+1]);
		end
	end
	println("obj_look_2 : $(obj_look_2)");
	p15 = StatsPlots.groupedbar(obj_look_2,bar_position=:dodge,bar_width=0.7,xlabel="Iteration",group=ctg,color=[:red :blue],xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=20,ylabel="Enhancement objective value");
	
	p16 = Plots.plot(p14, p15, layout=(2,1));
	figure();
	display(p16);

	##### RG, CG, CRG
	p17 = Plots.plot(1:length(obj_vec_BD), obj_vec_BD[end].*ones(length(obj_vec_BD)), linestyle=:dash,label=L"$EF$",linewidth=2,xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=20,titlefontsize=20,color="green");
	Plots.plot!(p17,1:length(obj_vec_BD),obj_vec_BD,line=(true, 2),marker=:rect,title="Evolution Objective Value",ylabel="Objective value[€]",xlabel="Iteration",label=L"RG",legend=:best,color=:blue); # color=:red
	Plots.plot!(p17, 1:length(obj_vect),obj_vect,line=(true, 2),marker=:circle,ylabel="Objective value[€]",xlabel="Iteration",label=L"CG",color=:red);
	Plots.plot!(p17, 1:length(obj_restricted_vector),obj_restricted_vector,line=(true, 2), marker=:utriangle,ylabel="Objective value[€]",xlabel="Iteration",label=L"CRG",color=:orange);

	obj_vect_CG2 = zeros(length(obj_vec_BD));
	obj_vect_CG2[1:length(obj_vect)] = obj_vect;
	obj_restricted_vector_2 = zeros(length(obj_vec_BD));
	obj_restricted_vector_2[1:length(obj_restricted_vector)] = obj_restricted_vector;
	obj_look_3 = zeros(length(obj_vect_CG2)-1,3);
	ctg3 = repeat(["Row Generation", "Column Generation", "Column-and-Row Generation"], inner=length(obj_vect_CG2)-1);
	for i=1:length(obj_vec_BD)-1
		obj_look_3[i,1] = abs(obj_vec_BD[i+1]-obj_vec_BD[i]);
		if i<=length(obj_vect)-1
			obj_look_3[i,2] = abs(obj_vect_CG2[i] - obj_vect_CG2[i+1]);
		end
		if i<=length(obj_restricted_vector)-1
			obj_look_3[i,3] = abs(obj_restricted_vector_2[i] - obj_restricted_vector_2[i+1]);
		end
	end
	p18 = StatsPlots.groupedbar(obj_look_3,bar_position=:dodge,bar_width=0.7,xlabel="Iteration",group=ctg3,color=[:red :orange :blue],xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=20,ylabel="Enhancement objective value");
	p19 = Plots.plot(p17, p18, layout=(2,1));
	figure();
	display(p19);
end

# Real_World_Example_Difference_RG_CRG()


function dual_lagrangian_example_price(price, Pmin, Pmax, L, VOLL, C, F)
	m = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
	JuMP.set_optimizer_attribute(m,"OutputFlag",0);
	@variable(m, p[1:2], lower_bound = 0);
	@variable(m, u[1:2], Bin);
	@variable(m, l, lower_bound = 0);
	@constraint(m, [g=1:2], u[g]Pmin[g]<=p[g]);
	@constraint(m, [g=1:2], p[g]<=u[g]Pmax[g]);
	@constraint(m, l<=L[1]);
	@objective(m, Min, sum( C[g]*p[g] + u[g]*F[g] for g=1:2) - VOLL*l - price*(sum(p[g] for g=1:2) - l) );
	optimize!(m);
	return objective_value(m);
end

function unit_commitment_example_price(Pmin, Pmax, L, VOLL, C, F)
	m = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
	JuMP.set_optimizer_attribute(m,"OutputFlag",0);
	@variable(m, p[1:2], lower_bound = 0);
	@variable(m, u[1:2], Bin);
	@variable(m, l, lower_bound = 0);
	@constraint(m, [g=1:2], u[g]Pmin[g]<=p[g]);
	@constraint(m, [g=1:2], p[g]<=u[g]Pmax[g]);
	@constraint(m, l<=L[1]);
	@constraint(m, sum( p[g] for g=1:2) - l==0);
	@objective(m, Min, sum( C[g]*p[g] + u[g]*F[g] for g=1:2) - VOLL*l);
	optimize!(m);
	return objective_value(m);
end

function examples_price(save_fig=1)
	### EXAMPLE 1
	Pmin = [0 6];Pmax = [15 10];L = [20]; VOLL = 60; C = [10 40];F = [0 0];
	price_max = VOLL+10;
	# price_max = 5;
	price_vec = 0:0.1:price_max
	uplift = zeros(length(price_vec));
	dual_lagrangian = zeros(length(price_vec));
	obj_uc = unit_commitment_example_price(Pmin, Pmax, L, VOLL, C, F);
	for (i,price) in enumerate(price_vec)
		obj_dl = dual_lagrangian_example_price(price, Pmin, Pmax, L, VOLL, C, F);
		uplift[i] = obj_uc - obj_dl;
		dual_lagrangian[i] = obj_dl;
	end

	price_ex1 = 40;
	dual_lagrangian_marginal_cost_pricing = dual_lagrangian_example_price(price_ex1, Pmin, Pmax, L, VOLL, C, F);
	uplift_marginal_cost_pricing = obj_uc - dual_lagrangian_marginal_cost_pricing;

	path = "D:/Documents/Mémoire Papavasiliou/Master Thesis - Manuscript/Images/Minimum Uplifts"
	pyplot();
	#figure();
	p1_1 = Plots.plot(price_vec, uplift, label="Uplifts", xlabel="price [€/MW]", ylabel="Sum Uplifts [€]", xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=15,titlefontsize=20,legend=:outerbottom); # legend=:outerbottom :outerright
	Plots.plot!(p1_1, price_vec, dual_lagrangian, label="Dual Lagrangian");
	Plots.scatter!(p1_1, [price_ex1], [uplift_marginal_cost_pricing], markersize=6,label=L"Uplift($\lambda=40$)");
	Plots.scatter!(p1_1, [price_ex1], [dual_lagrangian_marginal_cost_pricing],markersize=6,label=L"Dual Lagrangian($\lambda=40$)");
	Plots.plot!(p1_1, price_vec, obj_uc.*ones(length(price_vec)),linestyle=:dash,label="Unit Commitment");


	p3_1 = Plots.plot(price_vec, uplift, label="Uplifts", xlabel="price [€/MW]", ylabel="Sum Uplifts [€]", xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=15,titlefontsize=20,legend=:outerright); # legend=:outerbottom :outerright
	Plots.plot!(p3_1, price_vec, dual_lagrangian, label="Dual Lagrangian");
	Plots.scatter!(p3_1, [price_ex1], [uplift_marginal_cost_pricing], markersize=6,label=L"Uplift($\lambda=40$)");
	Plots.scatter!(p3_1, [price_ex1], [dual_lagrangian_marginal_cost_pricing],markersize=6,label=L"Dual Lagrangian($\lambda=40$)");
	Plots.plot!(p3_1, price_vec, obj_uc.*ones(length(price_vec)),linestyle=:dash,label="Unit Commitment");

	p2 = Plots.plot(price_vec, dual_lagrangian, label=false, xlabel="price [€/MW]", ylabel="Sum Uplifts [€]", xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=15,titlefontsize=20);
	if save_fig==1
		Plots.savefig(p1_1, "$(path)/uplifts_example1.png");
		Plots.savefig(p2, "$(path)/dual_lagrangian_example1.png");
		Plots.savefig(p3_1, "$(path)/uplifts_example1_2.png");
	end
	
	### EXAMPLE 2
	Pmin = [0 10];Pmax = [10 20];L = [25]; VOLL = 500; C = [10 30];F = [400 0];
	price_max = VOLL+10;
	obj_uc = unit_commitment_example_price(Pmin, Pmax, L, VOLL, C, F);
	# price_max = 5;
	price_vec = 0:0.1:price_max
	#price_vec = 0:10:price_max
	uplift = zeros(length(price_vec));
	dual_lagrangian = zeros(length(price_vec));
	for (i,price) in enumerate(price_vec)
		obj_dl = dual_lagrangian_example_price(price, Pmin, Pmax, L, VOLL, C, F);
		uplift[i] = obj_uc - obj_dl;
		dual_lagrangian[i] = obj_dl;
	end
	ind = argmin(uplift);
	println("price leading to smallest uplifts is ",price_vec[ind]," for uplifts ",uplift[ind]);
	price_ex2 = 30;
	dual_lagrangian_marginal_cost_pricing = dual_lagrangian_example_price(price_ex2, Pmin, Pmax, L, VOLL, C, F);
	uplift_marginal_cost_pricing = obj_uc - dual_lagrangian_marginal_cost_pricing;
	path = "D:/Documents/Mémoire Papavasiliou/Master Thesis - Manuscript/Images/Minimum Uplifts"
	
	#=println("price_vec = $(price_vec);");
	println("uplift = $(uplift);");
	println("dual_lagrangian = $(dual_lagrangian);");
	println("price_ex2 = $(price_ex2)");
	println("uplift_marginal_cost_pricing = $(uplift_marginal_cost_pricing);");
	println("dual_lagrangian_marginal_cost_pricing = $(dual_lagrangian_marginal_cost_pricing);")
	println("obj_uc = $(obj_uc);");=#
	println("Dual Lagrangian with price=50 : $(dual_lagrangian_example_price(50, Pmin, Pmax, L, VOLL, C, F))");

	pyplot();
	figure();
	p1 = Plots.plot(price_vec, uplift, label="Uplifts", xlabel="price [€/MW]", ylabel="Sum Uplifts [€]", xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=25,titlefontsize=20,legend=:best, linewidth=3); # legend=:outerbottom :outerright
	Plots.plot!(p1, price_vec, dual_lagrangian, label="Dual Lagrangian",linewidth=3);
	Plots.scatter!(p1, [price_ex2], [uplift_marginal_cost_pricing], markersize=8,label=L"Uplift($\lambda=30$)");
	Plots.scatter!(p1, [price_ex2], [dual_lagrangian_marginal_cost_pricing], markersize=8,label=L"Dual Lagrangian($\lambda=30$)");
	Plots.plot!(p1, price_vec, obj_uc.*ones(length(price_vec)),linestyle=:dash,label="Unit Commitment",linewidth=3);

	p2 = Plots.plot(price_vec, dual_lagrangian, label=false, xlabel="Price [€/MW]", ylabel="Sum Uplifts [€]", xtickfontsize=20,ytickfontsize=20,xguidefontsize=20,yguidefontsize=20,legendfontsize=25,titlefontsize=20);
	figure();
	display(p1);
	if save_fig==1
		Plots.savefig(p2, "$(path)/dual_lagrangian_example2.png");
		Plots.savefig(p1, "$(path)/uplifts_example2.png");
	end

end

examples_price()
