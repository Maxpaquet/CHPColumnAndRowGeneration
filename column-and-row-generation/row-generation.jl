##################################################################################################
##################################################################################################
########################################### ROW METHOD ###########################################
##################################################################################################
##################################################################################################

"""
Bender_Decomposition_1 does NOT consider "hot-start" models.
This function must be used to understand how row generation is working, not to measure any performances.

This function uses the auxilirary functions Cut_Generating_Linear_Program and LP_Relaxation_3bin_C describe below.
"""
function Bender_Decomposition_1(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, delat_criterion=0, Printer=false)
	iter_max = 500;
	obj_vec = zeros(iter_max);
	# Initialisation of Delta, Epsilon, Mu, Xi, Alpha, Sigma
	Delta = zeros(1+T_max,1);
	Epsilon = zeros(1+T_max,1);
	Mu = zeros(1+T_max,1);
	Xi = zeros(1+T_max,1);
	Alpha = zeros(1+T_max,1);
	Sigma = zeros(1+T_max,1);
	Solving_time = 0;

	nb_cuts_tot = 0;
	price_opt = 0;
	total_welfare_opt = 0;
	no_cuts_added = 0;
	obj_prev = 0;

	p_time_opt = 0; pbar_time_opt = 0; u_opt = 0; v_opt = 0; w_opt = 0;
	l_BD = zeros(T_max);
	obj_real_BD = 0;

	#start_time = time();
	#slave_mean_con = 0;
	#slave_mean_var = 0;

	for iter=1:iter_max
        if Printer
            println("\n");
            println("[Bender_Decomposition_1] iter : ", iter);
        end
		## Find solution with relaxation of 3-bin formulation
		(p_time_opt, pbar_time_opt, u_opt, v_opt, w_opt, price, l, obj_real, time_3bin_C) = LP_Relaxation_3bin_C(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, Delta, Epsilon, Mu, Xi, Alpha, Sigma, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer);
		Solving_time+=time_3bin_C;
		price_opt = price;
		obj_vec[iter] = obj_real;
		l_BD = l;
		obj_real_BD = obj_real;
		if Printer
            println("[Bender_Decomposition_1] Objective value of the Linear Relaxation (Master program) : $(obj_real)");
			println("[Bender_Decomposition_1] price = $(price)");
        end
		### Check for generating cut
		cut_added = false;
		for g=1:nb_gen
			if length(u_prior) > 0 && length(u_posterior)>0
				#start_time = time();
				(obj, z_opt, delta_opt, epsilon_opt, mu_opt, xi_opt, alpha_opt, sigma_opt, time_cut, nb_var, nb_con) = Cut_Generating_Linear_Program(MinRunCapacity[g], MaxRunCapacity[g], RU[g], RD[g], UT[g], DT[g], SU[g], SD[g], T_max, nb_gen, p_time_opt[g,:], pbar_time_opt[g,:], u_opt[g,:], v_opt[g,:], w_opt[g,:], u_prior[g], v_prior[g], w_prior[g], p_prior[g], pbar_prior[g], u_posterior[g], v_posterior[g], w_posterior[g], p_posterior[g], pbar_posterior[g]);
				#=slave_mean_con+=nb_con;
				slave_mean_var+=nb_var;
				elapsed = time()-start_time;
				println("################################# elapsed slave = $(elapsed)");=#
			elseif length(u_prior) > 0 && length(u_posterior)==0
				(obj, z_opt, delta_opt, epsilon_opt, mu_opt, xi_opt, alpha_opt, sigma_opt, time_cut) = Cut_Generating_Linear_Program(MinRunCapacity[g], MaxRunCapacity[g], RU[g], RD[g], UT[g], DT[g], SU[g], SD[g], T_max, nb_gen, p_time_opt[g,:], pbar_time_opt[g,:], u_opt[g,:], v_opt[g,:], w_opt[g,:], u_prior[g], v_prior[g], w_prior[g], p_prior[g], pbar_prior[g], [], [], [], [], []);
			elseif length(u_prior)==0 && length(u_posterior)>0
				(obj, z_opt, delta_opt, epsilon_opt, mu_opt, xi_opt, alpha_opt, sigma_opt, time_cut) = Cut_Generating_Linear_Program(MinRunCapacity[g], MaxRunCapacity[g], RU[g], RD[g], UT[g], DT[g], SU[g], SD[g], T_max, nb_gen, p_time_opt[g,:], pbar_time_opt[g,:], u_opt[g,:], v_opt[g,:], w_opt[g,:], [], [], [], [], [], u_posterior[g], v_posterior[g], w_posterior[g], p_posterior[g], pbar_posterior[g]);
			elseif length(u_prior)==0 && length(u_posterior)==0
				(obj, z_opt, delta_opt, epsilon_opt, mu_opt, xi_opt, alpha_opt, sigma_opt, time_cut) = Cut_Generating_Linear_Program(MinRunCapacity[g], MaxRunCapacity[g], RU[g], RD[g], UT[g], DT[g], SU[g], SD[g], T_max, nb_gen, p_time_opt[g,:], pbar_time_opt[g,:], u_opt[g,:], v_opt[g,:], w_opt[g,:], [], [], [], [], [], [], [], [], [], []);
			end
			Solving_time+=time_cut;
			if z_opt>delat_criterion
				cut_added = true;
                if Printer
					println("[Bender_Decomposition_1] g = $(g) -> z_opt = $(z_opt)");
                    println("[Bender_Decomposition_1] Cut added to master for generator $(g)");
                end
				### Add a cut to the master problem
				# Add vector delta_opt, epsilon_opt, mu_opt, xi_opt, alpha_opt, sigma_opt
				if(nb_cuts_tot == 0)
					Delta = Add_vector_first(Delta, delta_opt, g);
					Epsilon = Add_vector_first(Epsilon, epsilon_opt, g);
					Mu = Add_vector_first(Mu, mu_opt, g);
					Xi = Add_vector_first(Xi, xi_opt, g);
					Alpha = Add_vector_first(Alpha, alpha_opt,g);
					Sigma = Add_vector_first(Sigma, sigma_opt, g);
					nb_cuts_tot+=1;
				else
					Delta = Add_vector(Delta, delta_opt, g);
					Epsilon = Add_vector(Epsilon, epsilon_opt, g);
					Mu = Add_vector(Mu, mu_opt, g);
					Xi = Add_vector(Xi, xi_opt, g);
					Alpha = Add_vector(Alpha, alpha_opt,g);
					Sigma = Add_vector(Sigma, sigma_opt, g);
					nb_cuts_tot+=1;
				end
			end
			#elapsed_tmp = time()-start_time;
			#println("Time need for one generator on average : $(elapsed_tmp/g)");
			#println("slave_mean_var on average : $(slave_mean_var/g)");
			#println("slave_mean_con on average : $(slave_mean_con/g)");
		end

		if (cut_added==false) # Added no cut to none of the generator
			#=println("Stop at iteration ", iter);
			time_1iter = time()-start_time;
			println("==============================================");
			println("time_1iter : $(time_1iter)");
			println("slave_mean_var = $(slave_mean_var/nb_gen)");
			println("slave_mean_con = $(slave_mean_con/nb_gen)");
			println("==============================================");=#
			return (iter, price_opt, obj_vec, p_time_opt, pbar_time_opt, u_opt, v_opt, w_opt, l_BD, obj_real_BD, Solving_time);
		end
	end
	#=time_1iter = time()-start_time;
	println("==============================================");
	println("time_1iter : $(time_1iter)");
	println("slave_mean_var = $(slave_mean_var/nb_gen)");
	println("slave_mean_con = $(slave_mean_con/nb_gen)");
	println("==============================================");=#


	println("Stop at maximum iterations ", iter_max);
	return (iter_max, price_opt, obj_vec, p_time_opt, pbar_time_opt, u_opt, v_opt, w_opt, l_BD, obj_real_BD, Solving_time);
end

function Cut_Generating_Linear_Program(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, p_time_opt, pbar_time_opt, u_opt, v_opt, w_opt, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior)
	(A,B,nb_intervals_gen) = set_T(UT, T_max, 1);
	A = A[1,:];
	B = B[1,:];
	println("A = $(A)");
	println("B = $(B)");
	println("nb_intervals_gen = $(nb_intervals_gen)");
	nb_intervals_gen = nb_intervals_gen[1];
	if (!(nb_intervals_gen > 0))
		#println("	Generator not taken into account.")
		return (0, 0, zeros(T_max), zeros(T_max), zeros(T_max), zeros(T_max), zeros(T_max), zeros(T_max));
	end
	
	# m_cut = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
	m_cut = JuMP.direct_model(Gurobi.Optimizer(GUROBI_ENV)); # OutputFlag=0
	JuMP.set_optimizer_attribute(m_cut,"OutputFlag",0);
	# JuMP.set_optimizer_attribute(m_cut,"FeasibilityTol",1e-3);

	# Variables
	@variable(m_cut, 0 <= z);
	@variable(m_cut, p[i=1:nb_intervals_gen, t=0:T_max+1], lower_bound = 0);		# Power output vector of the generators		p[generator g, interval(a,b) i, time t]
	@variable(m_cut, pbar[i=1:nb_intervals_gen, t=0:T_max+1], lower_bound = 0);		# Maximum power available from generators	pbar[generator g, interval(a,b) i, time t]
	@variable(m_cut, gamma[i=1:nb_intervals_gen], lower_bound = 0, upper_bound = 1);# Variable that indicates whether a vertex (interval) is in the packing or not   gamma[generator g, interval(a,b) i]
	@variable(m_cut, p_time[t=0:T_max+1], lower_bound = 0);							# (to simplify computations) production of each generator on each time period
	@variable(m_cut, pbar_time[t=0:T_max+1], lower_bound = 0);						# (to simplify computations) maximum power available of each generator on each time period
	@variable(m_cut, u[t=0:T_max+1], lower_bound = 0, upper_bound = 1); 			# Commitment status variables
	@variable(m_cut, v[t=0:T_max+1], lower_bound = 0, upper_bound = 1); 			# Startup status variables
	@variable(m_cut, w[t=0:T_max+1], lower_bound = 0, upper_bound = 1); 			# Shutdown status variables
	### Feasible Dispatch Polytope (equations 13 from source)
	# For interval i of generator g, a = A[g,i] and b=A[g,i]
	# 13a and 13b
	for i=1:nb_intervals_gen
		for t=0:T_max+1
			if t<A[i] || t>B[i]
				@constraint(m_cut, p[i,t] <= 0 + z);
				@constraint(m_cut, pbar[i,t] <= 0 + z);
			end
		end
	end
	# 13c
	@constraint(m_cut, min_output[i=1:nb_intervals_gen, t=A[i]:B[i] ], -p[i,t] <= -gamma[i] * MinRunCapacity + z);
	# 13d
	@constraint(m_cut, max_output[i=1:nb_intervals_gen, t=A[i]:B[i] ], p[i,t] - pbar[i,t] <= 0 + z);
	# 13e
	@constraint(m_cut, upper_bound_on_power_output_1[i=1:nb_intervals_gen, t=A[i]:B[i] ] , pbar[i,t] <= gamma[i] * MaxRunCapacity + z);
	@constraint(m_cut, upper_bound_on_power_output_2[i=1:nb_intervals_gen, t=A[i]:B[i] ] , pbar[i,t] <= gamma[i] * (SU + (t - A[i])*RU) + z);
	@constraint(m_cut, upper_bound_on_power_output_3[i=1:nb_intervals_gen, t=A[i]:B[i] ] , pbar[i,t] <= gamma[i] * (SD + (B[i] - t)*RD) + z);
	# 13f For this constraint assume that a generator can already produces before the set [T_max] of time periods.
	@constraint(m_cut, limit_power_jumps_up_1[i=1:nb_intervals_gen, t=A[i]+1:B[i] ], pbar[i,t] - p[i,t-1] <= gamma[i] * RU + z);
	@constraint(m_cut, limit_power_jumps_up_2[i=1:nb_intervals_gen, t=A[i]+1:B[i] ], pbar[i,t] - p[i,t-1] <= gamma[i] * (SD + (B[i] - t)*RD - MinRunCapacity) + z);
	# 13g
	@constraint(m_cut, limit_power_jumps_down_1[i=1:nb_intervals_gen, t=A[i]+1:B[i] ], pbar[i,t-1] - p[i,t] <= gamma[i] * RD + z);
	@constraint(m_cut, limit_power_jumps_down_2[i=1:nb_intervals_gen, t=A[i]+1:B[i] ], pbar[i,t-1] - p[i,t] <= gamma[i] * (SU + (t - A[i])*RU - MinRunCapacity) + z);
	### Packing Dispatch Polytopes
	# 15e : to eliminate the impossible combinations
	@constraint(m_cut, minimum_down_time_constraint[t=1:T_max], sum(gamma[i] for i=1:nb_intervals_gen if (A[i] <= t && t <= B[i]+DT)) <= 1.0 + z);
	# 15b
	@constraint(m_cut, power_output_global[t=1:T_max], sum(p[i,t]  for i=1:nb_intervals_gen) == p_time_opt[t+1]); # Master program does not go until T_max+1
	#@constraint(m_cut, power_output_global_2[t=0,T_max+1], sum(p[i,t]  for i=1:nb_intervals_gen) == p_time_opt[t+1]);
	# 15c
	println("pbar_time_opt = $(pbar_time_opt)");
	@constraint(m_cut, maximum_power_available[t=1:T_max], sum(pbar[i,t]  for i=1:nb_intervals_gen) == pbar_time_opt[t+1]); # Master program does not go until T_max+1
	@constraint(m_cut, maximum_power_available_2[t=0,T_max+1], sum(pbar[i,t]  for i=1:nb_intervals_gen) == pbar_time_opt[t+1]);
	### Binary variable of UC (commitment status(on/off), startup status, shutdown status)
	@constraint(m_cut, commitment_status[t=1:T_max], sum(gamma[i] for i=1:nb_intervals_gen if (A[i] <= t && t <= B[i])) == u_opt[t+1]); # Master program does not go until T_max+
	@constraint(m_cut, commitment_status_2[t=0,T_max+1], sum(gamma[i] for i=1:nb_intervals_gen if (A[i] <= t && t <= B[i])) == u_opt[t+1]);

	@constraint(m_cut, startup_status[t=1:T_max], sum(gamma[i]  for i=1:nb_intervals_gen if (t == A[i])) == v_opt[t+1]); # Master program does not go until T_max+1
	@constraint(m_cut, startup_status_2[t=0,T_max+1], sum(gamma[i]  for i=1:nb_intervals_gen if (t == A[i])) == v_opt[t+1]);

	@constraint(m_cut, shutdown_status[t=1:T_max], sum(gamma[i] for i=1:nb_intervals_gen if (t == B[i]+1)) == w_opt[t+1]); # Master program does not go until T_max+1
	@constraint(m_cut, shutdown_status_2[t=0,T_max+1], sum(gamma[i] for i=1:nb_intervals_gen if (t == B[i]+1)) == w_opt[t+1]);
	### Prior data
	if length(u_prior) > 0
		JuMP.fix(u[0],u_prior; force = true);
		JuMP.fix(v[0],v_prior; force = true);
		JuMP.fix(w[0],w_prior; force = true);
		JuMP.fix(p_time[0],p_prior; force = true);
		JuMP.fix(pbar_time[0],pbar_prior; force = true);
	end
	### Posterior data
	if length(u_posterior) > 0
		JuMP.fix(u[T_max+1],u_posterior; force = true);
		JuMP.fix(v[T_max+1],v_posterior; force = true);
		JuMP.fix(w[T_max+1],w_posterior; force = true);
		JuMP.fix(p_time[T_max+1],p_posterior; force = true);
		JuMP.fix(pbar_time[T_max+1],pbar_posterior; force = true);
	end
	### Objective function
	@objective(m_cut, Min, z);

	optimize!(m_cut);
	delta_opt = dual.(minimum_down_time_constraint);
	epsilon_opt = dual.(power_output_global);
	mu_opt = dual.(maximum_power_available);
	xi_opt = dual.(commitment_status);
	alpha_opt = dual.(startup_status);
	sigma_opt = dual.(shutdown_status);

	#println("################################# num variables slave = $(num_variables(m_cut))");
	#println("################################# num constraints $(sum([num_constraints(m_cut, i, j) for (i,j) in list_of_constraint_types(m_cut)]))");
	nb_var = num_variables(m_cut);
	nb_con = sum([num_constraints(m_cut, i, j) for (i,j) in list_of_constraint_types(m_cut)]);
	# return (objective_value(m_cut), value(z), delta_opt, epsilon_opt, mu_opt, xi_opt, alpha_opt, sigma_opt, MOI.get(m_cut, MOI.SolveTime()), nb_var, nb_con);
	return (objective_value(m_cut), value(z), delta_opt, epsilon_opt, mu_opt, xi_opt, alpha_opt, sigma_opt, 0, nb_var, nb_con);
end

function LP_Relaxation_3bin_C(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD,T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, Delta, Epsilon, Mu, Xi, Alpha, Sigma, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer=false)
	# m_relax = JuMP.direct_model(Gurobi.Optimizer());
	m_relax = JuMP.direct_model(Gurobi.Optimizer(GUROBI_ENV));
	JuMP.set_optimizer_attribute(m_relax,"OutputFlag",0);

	@variable(m_relax, 0 <= p[g=1:nb_gen, t=0:T_max+1]);		# power output of each generator g at time period t
	@variable(m_relax, 0 <= pbar[g=1:nb_gen, t=0:T_max+1]);		# maximum power available of each generator g at time period t
	@variable(m_relax, 0 <= u[g=1:nb_gen, t=0:T_max+1]);		# commitment status of geenrator g at time period t
	@variable(m_relax, 0 <= v[g=1:nb_gen, t=0:T_max+1]); 		# startup status of generator g at time period t
	@variable(m_relax, 0 <= w[g=1:nb_gen, t=0:T_max+1]);		# shutdown status of generator g at time period t
	
	@variable(m_relax, 0 <= l[t=1:T_max]);
	@constraint(m_relax, [t=1:T_max], l[t] <= data_demand[t]);			# upper bound for elasticity of the demand

	@constraint(m_relax, LP_relax_u[g=1:nb_gen, t=0:T_max+1], u[g,t] <= 1);	# LP relaxation of the binary variable u
	@constraint(m_relax, LP_relax_v[g=1:nb_gen, t=0:T_max+1], v[g,t] <= 1);	# LP relaxation of the binary variable v
	@constraint(m_relax, LP_relax_w[g=1:nb_gen, t=0:T_max+1], w[g,t] <= 1);	# LP relaxation of the binary variable w

	### The following constriants described the feasible set for technical constraints, the relaxation of 3-bin space
	@constraint(m_relax, logical_constraint[g=1:nb_gen, t=1:T_max+1], u[g,t] - u[g,t-1] == v[g,t] - w[g,t]);				# Constraint 2, logical constraint on the different status (commitment status, startup status, shutdown status)
	@constraint(m_relax, minimum_up_time[g=1:nb_gen, t=UT[g]:T_max+1], sum(v[g,i]  for i=t-UT[g]+1:t) <= u[g,t]);			# Constraint 3, minimum up time constraints
	@constraint(m_relax, minimum_down_time[g=1:nb_gen, t=DT[g]:T_max+1], sum(w[g,i]  for i=t-DT[g]+1:t) <= (1 - u[g,t])); 	# Constraint 4, minimum down time constraints
	
	@constraint(m_relax, generation_limits_1[g=1:nb_gen, t=0:T_max+1], MinRunCapacity[g]*u[g,t] <= p[g,t]);					# Constraint 5_1, generation limits constraints
	@constraint(m_relax, generation_limits_2[g=1:nb_gen, t=0:T_max+1], p[g,t] <= pbar[g,t]);								# Constraint 5_2, generation limits constraints
	@constraint(m_relax, generation_limits_3[g=1:nb_gen, t=0:T_max+1], pbar[g,t] <= MaxRunCapacity[g]*u[g,t]);				# Constraint 5_3, generation limits constraints
	
	@constraint(m_relax, ramp_up_constraint[g=1:nb_gen, t=1:T_max+1], pbar[g,t] - p[g,t-1] <= RU[g]*u[g,t-1] + SU[g]*v[g,t]);	# Constraint 6, ramp-up constraints and start-up mode
	@constraint(m_relax, ramp_down_consraint[g=1:nb_gen, t=1:T_max+1], pbar[g,t-1] - p[g,t] <= RD[g]*u[g,t] + SD[g]*w[g,t]);	# Constraint 7, ramp-down constraints and shutdown mode
	
	### Production >= Demand
	@constraint(m_relax, loads[t=1:T_max], sum( p[g,t] for g=1:nb_gen) == l[t]);
	### Prior data
	if length(u_prior) > 0
		for g=1:nb_gen
			JuMP.fix(u[g,0], u_prior[g]; force = true);
			JuMP.fix(v[g,0], v_prior[g]; force = true);
			JuMP.fix(w[g,0], w_prior[g]; force = true);
			JuMP.fix(p[g,0], p_prior[g]; force = true);
			JuMP.fix(pbar[g,0], pbar_prior[g]; force = true);
		end
	end
	### Posterior data
	if length(u_posterior) > 0
		for g=1:nb_gen
			JuMP.fix(u[g,T_max+1],u_posterior[g]; force = true);
			JuMP.fix(v[g,T_max+1],v_posterior[g]; force = true);
			JuMP.fix(w[g,T_max+1],w_posterior[g]; force = true);
			JuMP.fix(p[g,T_max+1],p_posterior[g]; force = true);
			JuMP.fix(pbar[g,T_max+1],pbar_posterior[g]; force = true);
		end
	end

	### Generating cut : delta^T * ones + epsilon^T * p + mu^T * pbar + xi^T * u + alpha^T * v + sigma^T*w <= 0
	if Printer
		println("Delta : ",Delta);
		println("Epsilon : ",Epsilon);
		println("Mu : ",Mu);
		println("Xi : ",Xi);
		println("Alpha : ",Alpha);
		println("Sigma : ",Sigma);
	end
	(nb_rows, nb_cuts) = size(Delta); # The number of cuts is the number of columns
	for cut=1:nb_cuts
		g = convert(Int64, Delta[1,cut]);
		if (g!=0)
			@constraint(m_relax, sum(Delta[t,cut] for t=2:nb_rows) + sum(Epsilon[t,cut]*p[g,t-1] for t=2:nb_rows) + sum( Mu[t,cut]*pbar[g,t-1] for t=2:nb_rows)  + sum( Xi[t,cut]*u[g,t-1] for t=2:nb_rows) + sum( Alpha[t,cut]*v[g,t-1] for t=2:nb_rows) + sum(Sigma[t,cut] *w[g,t-1] for t=2:nb_rows) <= 0);
			#@constraint(m_relax, constraint_defining_C[cut=1:nb_cuts], sum(Delta[t,cut] for t=2:nb_rows) + sum(Epsilon[t,cut]*p[Delta[1,cut],t-1] for t=2:nb_rows) + sum( Mu[t,cut]*pbar[Delta[1,cut],t-1] for t=2:nb_rows)  + sum( Xi[t,cut]*u[Delta[1,cut],t-1] for t=2:nb_rows) + sum( Alpha[t,cut]*v[Delta[1,cut],t-1] for t=2:nb_rows) + sum(Sigma[t,cut] *w[Delta[1,cut],t-1] for t=2:nb_rows) <= 0);
		end
	end

	### Objective function
	@objective(m_relax, Min, sum( sum( NoLoadConsumption[g]*C[g]*u[g,t] + F[g]*v[g,t] + C[g]*p[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum( l[t] for t=1:T_max));
	optimize!(m_relax);
	price = dual.(loads);
	u_val = value.(u); v_val = value.(v); p_val = value.(p);
	if Printer
		println("[LP_Relaxation_3bin_C] objective of master program is $(objective_value(m_relax))");
		println("[LP_Relaxation_3bin_C] price is $(price)")
	end
	# return (value.(p).data, value.(pbar).data, value.(u).data, value.(v).data, value.(w).data, price, value.(l), objective_value(m_relax), MOI.get(m_relax, MOI.SolveTime()));
	return (value.(p).data, value.(pbar).data, value.(u).data, value.(v).data, value.(w).data, price, value.(l), objective_value(m_relax), 0);
end


"""
Need to specify which generator we are using.
G_c = Int64[g for g=1:nb_gen]
# gather_gen, indices_gen = Gather_Generators(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, G_c);
gather_gen = Int64[i for i=1:nb_gen]; indices_gen = Int64[i for i=1:nb_gen];
"""
function Bender_Decomposition(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, gather_gen, indices_gen, G_c, delta_criterion, Printer=false)
	iter_max = 500;

	nb_cuts_tot 	= 0;
	nb_salves_solve = 0;
	price_opt 		= 0;
	p_time_opt 		= 0;
	pbar_time_opt 	= 0;
	u_opt 			= 0;
	v_opt 			= 0;
	w_opt 			= 0;
	l_BD 			= zeros(T_max);
	obj_real_BD 	= 0;

	(A_vec,B_vec,nb_intervals_gen_vec) = set_T(UT, T_max, nb_gen);
	Solving_time = 0;
	prices_vect = [];
	p_vect = [];
	u_vect = [];
	v_vect = [];
	w_vect = [];
	obj_vec = [];
	tab_CGLP = Array{ModelCGLP}(undef, length(G_c));
	count = 1;
	### Definition of Cut-Generating Linear Program (slave programs) for all generators in G_c
	for g in G_c
		if Printer
			println("Building slave for generator $(g)");
		end
		### Data used for the problem
		nb_intervals_gen = nb_intervals_gen_vec[g];
		A = A_vec[g,:];
		B = B_vec[g,:];
		### Model variable
		m_cut = JuMP.direct_model(Gurobi.Optimizer());
		JuMP.set_optimizer_attribute(m_cut,"OutputFlag",0);
		# JuMP.set_optimizer_attribute(m_cut,"FeasibilityTol",1e-3); If not feasible, de-comment this line
		### Variables
		@variable(m_cut, z, lower_bound = 0);
		@variable(m_cut, p_cut[i=1:nb_intervals_gen, t=0:T_max+1], lower_bound = 0);	# Power output vector of the generators		p[generator g, interval(a,b) i, time t]
		@variable(m_cut, pbar_cut[i=1:nb_intervals_gen, t=0:T_max+1], lower_bound = 0);	# Maximum power available from generators	pbar[generator g, interval(a,b) i, time t]
		# @variable(m_cut, gamma_cut[i=1:nb_intervals_gen]);				# Variable that indicates whether a vertex (interval) is in the packing or not   gamma[generator g, interval(a,b) i]
		@variable(m_cut, gamma_cut[i=1:nb_intervals_gen], lower_bound = 0, upper_bound = 1);
		@variable(m_cut, p_time_cut[t=0:T_max+1], lower_bound = 0);						# (to simplify computations) production of each generator on each time period
		@variable(m_cut, pbar_time_cut[t=0:T_max+1], lower_bound = 0);					# (to simplify computations) maximum power available of each generator on each time period
		@variable(m_cut, u_cut[t=0:T_max+1], lower_bound = 0, upper_bound = 1); 							# Commitment status variables
		@variable(m_cut, v_cut[t=0:T_max+1], lower_bound = 0, upper_bound = 1); 							# Startup status variables
		@variable(m_cut, w_cut[t=0:T_max+1], lower_bound = 0, upper_bound = 1); 							# Shutdown status variables		
		@objective(m_cut, Min, z);
		### Feasible Dispatch Polytope (equations 13 from source)
		for i=1:nb_intervals_gen
			for t=0:T_max+1
				if t<A[i] || t>B[i]
					@constraint(m_cut, p_cut[i,t] <= 0 + z);
					@constraint(m_cut, pbar_cut[i,t] <= 0 + z);
				end
			end
		end
		min_output 						= @constraint(m_cut, [i=1:nb_intervals_gen, t=A[i]:B[i] ], -p_cut[i,t] <= -gamma_cut[i] * MinRunCapacity[g] + z);
		max_output 						= @constraint(m_cut, [i=1:nb_intervals_gen, t=A[i]:B[i] ], p_cut[i,t] - pbar_cut[i,t] <= 0 + z);
		upper_bound_on_power_output_1 	= @constraint(m_cut, [i=1:nb_intervals_gen, t=A[i]:B[i] ] , pbar_cut[i,t] <= gamma_cut[i] * MaxRunCapacity[g] + z);
		upper_bound_on_power_output_2 	= @constraint(m_cut, [i=1:nb_intervals_gen, t=A[i]:B[i] ] , pbar_cut[i,t] <= gamma_cut[i] * (SU[g] + (t - A[i])*RU[g]) + z);
		upper_bound_on_power_output_3 	= @constraint(m_cut, [i=1:nb_intervals_gen, t=A[i]:B[i] ] , pbar_cut[i,t] <= gamma_cut[i] * (SD[g] + (B[i] - t)*RD[g]) + z);
		limit_power_jumps_up_1 			= @constraint(m_cut, [i=1:nb_intervals_gen, t=A[i]+1:B[i] ], pbar_cut[i,t] - p_cut[i,t-1] <= gamma_cut[i] * RU[g] + z);
		limit_power_jumps_up_2 			= @constraint(m_cut, [i=1:nb_intervals_gen, t=A[i]+1:B[i] ], pbar_cut[i,t] - p_cut[i,t-1] <= gamma_cut[i] * (SD[g] + (B[i] - t)*RD[g] - MinRunCapacity[g]) + z);
		limit_power_jumps_down_1 		= @constraint(m_cut, [i=1:nb_intervals_gen, t=A[i]+1:B[i] ], pbar_cut[i,t-1] - p_cut[i,t] <= gamma_cut[i] * RD[g] + z);
		limit_power_jumps_down_2 		= @constraint(m_cut, [i=1:nb_intervals_gen, t=A[i]+1:B[i] ], pbar_cut[i,t-1] - p_cut[i,t] <= gamma_cut[i] * (SU[g] + (t - A[i])*RU[g] - MinRunCapacity[g]) + z);
		minimum_down_time_constraint 	= @constraint(m_cut, [t=1:T_max], sum(gamma_cut[i] for i=1:nb_intervals_gen if (A[i] <= t && t <= B[i]+DT[g])) <= 1.0 + z);
		### Will modify rhs
		power_output_global 	= @constraint(m_cut, [t=1:T_max], sum(p_cut[i,t]  for i=1:nb_intervals_gen) == 0); # Master program does not go until T_max+1
		maximum_power_available = @constraint(m_cut, [t=1:T_max], sum(pbar_cut[i,t]  for i=1:nb_intervals_gen) == 0); # Master program does not go until T_max+1
		commitment_status 		= @constraint(m_cut, [t=1:T_max], sum(gamma_cut[i] for i=1:nb_intervals_gen if (A[i] <= t && t <= B[i])) == 0); # Master program does not go until T_max+
		startup_status 			= @constraint(m_cut, [t=1:T_max], sum(gamma_cut[i]  for i=1:nb_intervals_gen if (t == A[i])) == 0); # Master program does not go until T_max+1
		shutdown_status 		= @constraint(m_cut, [t=1:T_max], sum(gamma_cut[i] for i=1:nb_intervals_gen if (t == B[i]+1)) == 0); # Master program does not go until T_max+1
		k_vec = [0,T_max+1];
		power_output_global_2		= @constraint(m_cut, [k=1:2], sum(p_cut[i,k_vec[k]]  for i=1:nb_intervals_gen) == 0);
		maximum_power_available_2 	= @constraint(m_cut, [k=1:2], sum(pbar_cut[i,k_vec[k]]  for i=1:nb_intervals_gen) == 0);
		commitment_status_2 		= @constraint(m_cut, [k=1:2], sum(gamma_cut[i] for i=1:nb_intervals_gen if (A[i] <= k_vec[k] && k_vec[k] <= B[i])) == 0);
		startup_status_2 			= @constraint(m_cut, [k=1:2], sum(gamma_cut[i]  for i=1:nb_intervals_gen if (k_vec[k] == A[i])) == 0);
		shutdown_status_2 			= @constraint(m_cut, [k=1:2], sum(gamma_cut[i] for i=1:nb_intervals_gen if (k_vec[k] == B[i]+1)) == 0);

		#obj = ModelCGLP(m_cut, z, p_cut, pbar_cut, p_time_cut, pbar_time_cut, u_cut, v_cut, w_cut, gamma_cut, min_output, max_output, upper_bound_on_power_output_1, upper_bound_on_power_output_2, upper_bound_on_power_output_3, limit_power_jumps_up_1, limit_power_jumps_up_2, limit_power_jumps_down_1, limit_power_jumps_down_2, minimum_down_time_constraint, power_output_global, maximum_power_available, maximum_power_available_2, commitment_status, commitment_status_2, startup_status, startup_status_2, shutdown_status, shutdown_status_2, con_added_u, con_added_v, con_added_w);
		obj = ModelCGLP(m_cut, z, p_cut, pbar_cut, p_time_cut, pbar_time_cut, u_cut, v_cut, w_cut, gamma_cut, min_output, max_output, upper_bound_on_power_output_1, upper_bound_on_power_output_2, upper_bound_on_power_output_3, limit_power_jumps_up_1, limit_power_jumps_up_2, limit_power_jumps_down_1, limit_power_jumps_down_2, minimum_down_time_constraint, power_output_global, power_output_global_2, maximum_power_available, maximum_power_available_2, commitment_status, commitment_status_2, startup_status, startup_status_2, shutdown_status, shutdown_status_2);
		tab_CGLP[count] = obj;
		count+=1;
		
		if length(u_prior) > 0
			JuMP.fix(obj.u[0], u_prior[g]; force = true);
			JuMP.fix(obj.v[0], v_prior[g]; force = true);
			JuMP.fix(obj.w[0], w_prior[g]; force = true);
			JuMP.fix(obj.p_time[0], p_prior[g]; force = true);
			JuMP.fix(obj.pbar_time[0], pbar_prior[g]; force = true);
		end
		### Posterior data
		if length(u_posterior) > 0
			JuMP.fix(obj.u[T_max+1], u_posterior[g]; force = true);
			JuMP.fix(obj.v[T_max+1], v_posterior[g]; force = true);
			JuMP.fix(obj.w[T_max+1], w_posterior[g]; force = true);
			JuMP.fix(obj.p_time[T_max+1], p_posterior[g]; force = true);
			JuMP.fix(obj.pbar_time[T_max+1], pbar_posterior[g]; force = true);
		end
	end
	### Definition of master program
	m_relax = JuMP.direct_model(Gurobi.Optimizer());
	JuMP.set_optimizer_attribute(m_relax,"OutputFlag",0);

	@variable(m_relax, p[g=1:nb_gen, t=0:T_max+1], lower_bound=0);					# power output of each generator g at time period t
	@variable(m_relax, pbar[g=1:nb_gen, t=0:T_max+1], lower_bound=0);				# maximum power available of each generator g at time period t
	@variable(m_relax, u[g=1:nb_gen, t=0:T_max+1], lower_bound=0, upper_bound=1);	# commitment status of geenrator g at time period t
	@variable(m_relax, v[g=1:nb_gen, t=0:T_max+1], lower_bound=0, upper_bound=1); 	# startup status of generator g at time period t
	@variable(m_relax, w[g=1:nb_gen, t=0:T_max+1], lower_bound=0, upper_bound=1);	# shutdown status of generator g at time period t

	@variable(m_relax, 0 <= l[t=1:T_max]);
	@constraint(m_relax, [t=1:T_max], l[t] <= data_demand[t]);	# upper bound for elasticity of the demand
	### The following constriants described the feasible set for technical constraints, the relaxation of 3-bin space
	@constraint(m_relax, logical_constraint[g=1:nb_gen, t=1:T_max+1], u[g,t] - u[g,t-1] == v[g,t] - w[g,t]);				# Constraint 2, logical constraint on the different status (commitment status, startup status, shutdown status)
	@constraint(m_relax, minimum_up_time[g=1:nb_gen, t=UT[g]:T_max+1], sum(v[g,i]  for i=t-UT[g]+1:t) <= u[g,t]);			# Constraint 3, minimum up time constraints
	@constraint(m_relax, minimum_down_time[g=1:nb_gen, t=DT[g]:T_max+1], sum(w[g,i]  for i=t-DT[g]+1:t) <= (1 - u[g,t])); 	# Constraint 4, minimum down time constraints
	@constraint(m_relax, generation_limits_1[g=1:nb_gen, t=0:T_max+1], MinRunCapacity[g]*u[g,t] <= p[g,t]);					# Constraint 5_1, generation limits constraints
	@constraint(m_relax, generation_limits_2[g=1:nb_gen, t=0:T_max+1], p[g,t] <= pbar[g,t]);								# Constraint 5_2, generation limits constraints
	@constraint(m_relax, generation_limits_3[g=1:nb_gen, t=0:T_max+1], pbar[g,t] <= MaxRunCapacity[g]*u[g,t]);				# Constraint 5_3, generation limits constraints
	@constraint(m_relax, ramp_up_constraint[g=1:nb_gen, t=1:T_max+1], pbar[g,t] - p[g,t-1] <= RU[g]*u[g,t-1] + SU[g]*v[g,t]);	# Constraint 6, ramp-up constraints and start-up mode
	@constraint(m_relax, ramp_down_consraint[g=1:nb_gen, t=1:T_max+1], pbar[g,t-1] - p[g,t] <= RD[g]*u[g,t] + SD[g]*w[g,t]);	# Constraint 7, ramp-down constraints and shutdown mode
	### Production >= Demand
	@constraint(m_relax, loads[t=1:T_max], sum( p[g,t] for g=1:nb_gen) == l[t]);
	### Prior data
	if length(u_prior) > 0
		for g=1:nb_gen
			JuMP.fix(u[g,0], u_prior[g]; force = true);
			JuMP.fix(v[g,0], v_prior[g]; force = true);
			JuMP.fix(w[g,0], w_prior[g]; force = true);
			JuMP.fix(p[g,0], p_prior[g]; force = true);
			JuMP.fix(pbar[g,0], pbar_prior[g]; force = true);
		end
	end
	### Posterior data
	if length(u_posterior) > 0
		for g=1:nb_gen
			JuMP.fix(u[g,T_max+1], u_posterior[g]; force = true);
			JuMP.fix(v[g,T_max+1], v_posterior[g]; force = true);
			JuMP.fix(w[g,T_max+1], w_posterior[g]; force = true);
			JuMP.fix(p[g,T_max+1], p_posterior[g]; force = true);
			JuMP.fix(pbar[g,T_max+1], pbar_posterior[g]; force = true);
		end
	end
	### Objective function of master program
	@objective(m_relax, Min, sum( sum( NoLoadConsumption[g]*C[g]*u[g,t] + F[g]*v[g,t] + C[g]*p[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum( l[t] for t=1:T_max));
	slave_ST_mean = 0;
	master_ST_mean = 0;
	### Start iterations of Bender decomposition
	for iter=1:iter_max
		### Generating cut : delta^T * ones + epsilon^T * p + mu^T * pbar + xi^T * u + alpha^T * v + sigma^T*w <= 0
		optimize!(m_relax);
		# Solving_time+=MOI.get(m_relax, MOI.SolveTime());
		Solving_time+=0;
		# master_ST_mean+=MOI.get(m_relax, MOI.SolveTime());
		master_ST_mean+=0;
		price 			= dual.(loads);
		push!(prices_vect,price);
		p_time_opt 		= value.(p).data;
		push!(p_vect,p_time_opt);
		pbar_time_opt 	= value.(pbar).data;
		u_opt 			= value.(u).data;
		push!(u_vect, u_opt);
		v_opt 			= value.(v).data;
		push!(v_vect, v_opt);
		w_opt 			= value.(w).data;
		push!(w_vect, w_opt);
		obj_real 		= objective_value(m_relax);
		push!(obj_vec,obj_real);
		if Printer
			println();
			println("iter ",iter);
			println("obj : ",objective_value(m_relax));
			println("price : $(price)");
			println("obj_vec : $(obj_vec)");
			println("p_time_opt : $(p_time_opt)");
			println("pbar_time_opt : $(pbar_time_opt)");
			println("u_opt : $(u_opt)");
			println("v_opt : $(v_opt)");
			println("w_opt : $(w_opt)");
			println();
		end
		price_opt 		= price;
		l_BD 			= value.(l);
		obj_real_BD 	= obj_real;
		### Check for generating cut
		cut_added = false;
		for g_index=1:length(G_c)
			g_considered 	= G_c[g_index];
			elem_CGLP 		= tab_CGLP[g_index];
			if Printer
				println();println("g_considered : ",g_considered);
			end
			for t=1:T_max
				set_normalized_rhs(elem_CGLP.power_output_global[t], p_time_opt[g_considered,t+1]);
				set_normalized_rhs(elem_CGLP.maximum_power_available[t], pbar_time_opt[g_considered,t+1]);
				set_normalized_rhs(elem_CGLP.commitment_status[t], u_opt[g_considered,t+1]);
				set_normalized_rhs(elem_CGLP.startup_status[t], v_opt[g_considered,t+1]);
				set_normalized_rhs(elem_CGLP.shutdown_status[t], w_opt[g_considered,t+1]);
			end
			k_vec = [0,T_max+1]
			for k=1:length(k_vec)
				set_normalized_rhs(elem_CGLP.power_output_global_2[k], p_time_opt[g_considered,k_vec[k]+1]);
				set_normalized_rhs(elem_CGLP.maximum_power_available_2[k], pbar_time_opt[g_considered,k_vec[k]+1]);
				set_normalized_rhs(elem_CGLP.commitment_status_2[k], u_opt[g_considered,k_vec[k]+1]);
				set_normalized_rhs(elem_CGLP.startup_status_2[k], v_opt[g_considered,k_vec[k]+1]);
				set_normalized_rhs(elem_CGLP.shutdown_status_2[k], w_opt[g_considered,k_vec[k]+1]);
			end

			optimize!(elem_CGLP.model);
			# Solving_time+=MOI.get(elem_CGLP.model, MOI.SolveTime());
			Solving_time+=0;
			# slave_ST_mean+=MOI.get(elem_CGLP.model, MOI.SolveTime());
			slave_ST_mean+=0;
			z_opt = value.(elem_CGLP.z);
			delta_opt 	= dual.(elem_CGLP.minimum_down_time_constraint);
			epsilon_opt = dual.(elem_CGLP.power_output_global);
			mu_opt 		= dual.(elem_CGLP.maximum_power_available);
			xi_opt 		= dual.(elem_CGLP.commitment_status);
			alpha_opt 	= dual.(elem_CGLP.startup_status);
			sigma_opt 	= dual.(elem_CGLP.shutdown_status);

			if Printer
				println("z_opt : $(z_opt)");
				println("delta_opt (ones) : $(delta_opt)");
				println("epsilon_opt (p) : $(epsilon_opt)");
				println("mu_opt (pbar) : $(mu_opt)");
				println("xi_opt (u) : $(xi_opt)");
				println("alpha_opt (v) : $(alpha_opt)");
				println("sigma_opt (w) : $(sigma_opt)");
			end

			if z_opt>delta_criterion
				if Printer
					println("Cut added due to generator $(g_considered)");
				end
				cut_added = true;
				nb_salves_solve+=1;
				for i in gather_gen[indices_gen[g_considered]]
					nb_cuts_tot+=1;
					# println("cut : $(sum(delta_opt[t] + epsilon_opt[t]*p[i,t] + mu_opt[t]*pbar[i,t] + xi_opt[t]*u[i,t] + alpha_opt[t]*v[i,t] + sigma_opt[t]*w[i,t] for t=1:T_max)) <= 0");
					@constraint(m_relax, sum(delta_opt[t] + epsilon_opt[t]*p[i,t] + mu_opt[t]*pbar[i,t] + xi_opt[t]*u[i,t] + alpha_opt[t]*v[i,t] + sigma_opt[t]*w[i,t] for t=1:T_max) <= 0);
				end
				# @constraint(m_relax, sum(delta_opt[t] + epsilon_opt[t]*p[g_considered,t] + mu_opt[t]*pbar[g_considered,t] + xi_opt[t]*u[g_considered,t] + alpha_opt[t]*v[g_considered,t] + sigma_opt[t]*w[g_considered,t] for t=1:T_max) <= 0);
			end
		end
		if (cut_added==false) # Added no cut to none of the generator
			println("Stop at iteration ", iter);
			return (iter, price_opt, obj_vec, p_time_opt, pbar_time_opt, u_opt, v_opt, w_opt, l_BD, obj_real_BD, Solving_time, nb_cuts_tot, nb_salves_solve, prices_vect, p_vect, u_vect, v_vect, w_vect, slave_ST_mean, master_ST_mean);
		end
	end
	println("Stop at maximum iterations ", iter_max);
	return (iter_max, price_opt, obj_vec, p_time_opt, pbar_time_opt, u_opt, v_opt, w_opt, l_BD, obj_real_BD, Solving_time, nb_cuts_tot, nb_salves_solve, prices_vect, p_vect, u_vect, v_vect, w_vect, slave_ST_mean, master_ST_mean);
end




function Bender_Decomposition_NO_Posterior(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, gather_gen, indices_gen, G_c, delta_criterion, Printer=false)
	iter_max = 500;
	obj_vec = zeros(iter_max);

	nb_cuts_tot 	= 0;
	nb_salves_solve = 0;
	price_opt 		= 0;
	p_time_opt 		= 0;
	pbar_time_opt 	= 0;
	u_opt 			= 0;
	v_opt 			= 0;
	w_opt 			= 0;
	l_BD 			= zeros(T_max);
	obj_real_BD 	= 0;

	(A_vec,B_vec,nb_intervals_gen_vec) = set_T(UT, T_max, nb_gen);
	Solving_time = 0;
	prices_vect = [];
	p_vect = [];
	tab_CGLP = Array{ModelCGLP}(undef, length(G_c));
	count = 1;
	### Definition of Cut-Generating Linear Program (slave programs) for all generators in G_c
	for g in G_c
		### Data used for the problem
		nb_intervals_gen = nb_intervals_gen_vec[g];
		A = A_vec[g,:];
		B = B_vec[g,:];
		### Model variable
		m_cut = JuMP.direct_model(Gurobi.Optimizer());
		JuMP.set_optimizer_attribute(m_cut,"OutputFlag",0);
		# JuMP.set_optimizer_attribute(m_cut,"FeasibilityTol",1e-3); # If not feasible, de-comment this line
		### Variables
		@variable(m_cut, z, lower_bound = 0);
		@variable(m_cut, p_cut[i=1:nb_intervals_gen, t=0:T_max], lower_bound = 0);	# Power output vector of the generators		p[generator g, interval(a,b) i, time t]
		@variable(m_cut, pbar_cut[i=1:nb_intervals_gen, t=0:T_max], lower_bound = 0);	# Maximum power available from generators	pbar[generator g, interval(a,b) i, time t]
		@variable(m_cut, gamma_cut[i=1:nb_intervals_gen], lower_bound = 0); # Variable that indicates whether a vertex (interval) is in the packing or not   gamma[generator g, interval(a,b) i]   , upper_bound = 1
		@variable(m_cut, p_time_cut[t=0:T_max], lower_bound = 0);						# (to simplify computations) production of each generator on each time period
		@variable(m_cut, pbar_time_cut[t=0:T_max], lower_bound = 0);					# (to simplify computations) maximum power available of each generator on each time period
		@variable(m_cut, u_cut[t=0:T_max], lower_bound = 0, upper_bound = 1); 							# Commitment status variables
		@variable(m_cut, v_cut[t=0:T_max], lower_bound = 0, upper_bound = 1); 							# Startup status variables
		@variable(m_cut, w_cut[t=0:T_max], lower_bound = 0, upper_bound = 1); 							# Shutdown status variables		
		@objective(m_cut, Min, z);
		### Feasible Dispatch Polytope (equations 13 from source)
		for i=1:nb_intervals_gen
			for t=0:T_max
				if t<A[i] || t>B[i]
					@constraint(m_cut, p_cut[i,t] <= 0 + z);
					@constraint(m_cut, pbar_cut[i,t] <= 0 + z);
				end
			end
		end
		min_output 						= @constraint(m_cut, [i=1:nb_intervals_gen, t=A[i]:B[i] ], -p_cut[i,t] <= -gamma_cut[i] * MinRunCapacity[g] + z);
		max_output 						= @constraint(m_cut, [i=1:nb_intervals_gen, t=A[i]:B[i] ], p_cut[i,t] - pbar_cut[i,t] <= 0 + z);
		upper_bound_on_power_output_1 	= @constraint(m_cut, [i=1:nb_intervals_gen, t=A[i]:B[i] ] , pbar_cut[i,t] <= gamma_cut[i] * MaxRunCapacity[g] + z);
		upper_bound_on_power_output_2 	= @constraint(m_cut, [i=1:nb_intervals_gen, t=A[i]:B[i] ] , pbar_cut[i,t] <= gamma_cut[i] * (SU[g] + (t - A[i])*RU[g]) + z);
		upper_bound_on_power_output_3 	= @constraint(m_cut, [i=1:nb_intervals_gen, t=A[i]:B[i] ] , pbar_cut[i,t] <= gamma_cut[i] * (SD[g] + (B[i] - t)*RD[g]) + z);
		limit_power_jumps_up_1 			= @constraint(m_cut, [i=1:nb_intervals_gen, t=A[i]+1:B[i] ], pbar_cut[i,t] - p_cut[i,t-1] <= gamma_cut[i] * RU[g] + z);
		limit_power_jumps_up_2 			= @constraint(m_cut, [i=1:nb_intervals_gen, t=A[i]+1:B[i] ], pbar_cut[i,t] - p_cut[i,t-1] <= gamma_cut[i] * (SD[g] + (B[i] - t)*RD[g] - MinRunCapacity[g]) + z);
		limit_power_jumps_down_1 		= @constraint(m_cut, [i=1:nb_intervals_gen, t=A[i]+1:B[i] ], pbar_cut[i,t-1] - p_cut[i,t] <= gamma_cut[i] * RD[g] + z);
		limit_power_jumps_down_2 		= @constraint(m_cut, [i=1:nb_intervals_gen, t=A[i]+1:B[i] ], pbar_cut[i,t-1] - p_cut[i,t] <= gamma_cut[i] * (SU[g] + (t - A[i])*RU[g] - MinRunCapacity[g]) + z);
		minimum_down_time_constraint 	= @constraint(m_cut, [t=1:T_max], sum(gamma_cut[i] for i=1:nb_intervals_gen if (A[i] <= t && t <= B[i]+DT[g])) <= 1.0 + z);
		### Will modify rhs
		power_output_global 	= @constraint(m_cut, [t=1:T_max], sum(p_cut[i,t]  for i=1:nb_intervals_gen) == 0); # Master program does not go until T_max+1
		maximum_power_available = @constraint(m_cut, [t=1:T_max], sum(pbar_cut[i,t]  for i=1:nb_intervals_gen) == 0); # Master program does not go until T_max+1
		commitment_status 		= @constraint(m_cut, [t=1:T_max], sum(gamma_cut[i] for i=1:nb_intervals_gen if (A[i] <= t && t <= B[i])) == 0); # Master program does not go until T_max+
		startup_status 			= @constraint(m_cut, [t=1:T_max], sum(gamma_cut[i]  for i=1:nb_intervals_gen if (t == A[i])) == 0); # Master program does not go until T_max+1
		shutdown_status 		= @constraint(m_cut, [t=1:T_max], sum(gamma_cut[i] for i=1:nb_intervals_gen if (t == B[i]+1)) == 0); # Master program does not go until T_max+1
		k_vec = [0];
		power_output_global_2		= @constraint(m_cut, [k=1], sum(p_cut[i,k_vec[k]]  for i=1:nb_intervals_gen) == p_prior[g]);
		maximum_power_available_2 	= @constraint(m_cut, [k=1], sum(pbar_cut[i,k_vec[k]]  for i=1:nb_intervals_gen) == pbar_prior[g]);
		commitment_status_2 		= @constraint(m_cut, [k=1], sum(gamma_cut[i] for i=1:nb_intervals_gen if (A[i] <= k_vec[k] && k_vec[k] <= B[i])) == u_prior[g]);
		startup_status_2 			= @constraint(m_cut, [k=1], sum(gamma_cut[i] for i=1:nb_intervals_gen if (k_vec[k] == A[i])) == v_prior[g]);
		shutdown_status_2 			= @constraint(m_cut, [k=1], sum(gamma_cut[i] for i=1:nb_intervals_gen if (k_vec[k] == B[i]+1)) == w_prior[g]);
		
		if length(u_prior) > 0
			JuMP.fix(u_cut[0], u_prior[g]; force = true);
			JuMP.fix(v_cut[0], v_prior[g]; force = true);
			JuMP.fix(w_cut[0], w_prior[g]; force = true);
			JuMP.fix(p_time_cut[0], p_prior[g]; force = true);
			JuMP.fix(pbar_time_cut[0], pbar_prior[g]; force = true);
		end

		obj = ModelCGLP(m_cut, z, p_cut, pbar_cut, p_time_cut, pbar_time_cut, u_cut, v_cut, w_cut, gamma_cut, min_output, max_output, upper_bound_on_power_output_1, upper_bound_on_power_output_2, upper_bound_on_power_output_3, limit_power_jumps_up_1, limit_power_jumps_up_2, limit_power_jumps_down_1, limit_power_jumps_down_2, minimum_down_time_constraint, power_output_global, power_output_global_2, maximum_power_available, maximum_power_available_2, commitment_status, commitment_status_2, startup_status, startup_status_2, shutdown_status, shutdown_status_2);
		tab_CGLP[count] = obj;
		count+=1;
	end
	### Definition of master program
	m_relax = JuMP.direct_model(Gurobi.Optimizer());
	# JuMP.set_optimizer_attribute(m_relax,"OutputFlag",0);
	@variable(m_relax, p[g=1:nb_gen, t=0:T_max], lower_bound=0);					# power output of each generator g at time period t
	@variable(m_relax, pbar[g=1:nb_gen, t=0:T_max], lower_bound=0);				# maximum power available of each generator g at time period t
	@variable(m_relax, u[g=1:nb_gen, t=0:T_max], lower_bound=0, upper_bound=1);	# commitment status of geenrator g at time period t
	@variable(m_relax, v[g=1:nb_gen, t=0:T_max], lower_bound=0, upper_bound=1); 	# startup status of generator g at time period t
	@variable(m_relax, w[g=1:nb_gen, t=0:T_max], lower_bound=0, upper_bound=1);	# shutdown status of generator g at time period t

	@variable(m_relax, 0 <= l[t=1:T_max]);
	@constraint(m_relax, [t=1:T_max], l[t] <= data_demand[t]);	# upper bound for elasticity of the demand
	### The following constriants described the feasible set for technical constraints, the relaxation of 3-bin space
	@constraint(m_relax, logical_constraint[g=1:nb_gen, t=1:T_max], u[g,t] - u[g,t-1] == v[g,t] - w[g,t]);				# Constraint 2, logical constraint on the different status (commitment status, startup status, shutdown status)
	@constraint(m_relax, minimum_up_time[g=1:nb_gen, t=UT[g]:T_max], sum(v[g,i]  for i=t-UT[g]+1:t) <= u[g,t]);			# Constraint 3, minimum up time constraints
	@constraint(m_relax, minimum_down_time[g=1:nb_gen, t=DT[g]:T_max], sum(w[g,i]  for i=t-DT[g]+1:t) <= (1 - u[g,t])); 	# Constraint 4, minimum down time constraints
	@constraint(m_relax, generation_limits_1[g=1:nb_gen, t=0:T_max], MinRunCapacity[g]*u[g,t] <= p[g,t]);					# Constraint 5_1, generation limits constraints
	@constraint(m_relax, generation_limits_2[g=1:nb_gen, t=0:T_max], p[g,t] <= pbar[g,t]);								# Constraint 5_2, generation limits constraints
	@constraint(m_relax, generation_limits_3[g=1:nb_gen, t=0:T_max], pbar[g,t] <= MaxRunCapacity[g]*u[g,t]);				# Constraint 5_3, generation limits constraints
	@constraint(m_relax, ramp_up_constraint[g=1:nb_gen, t=1:T_max], pbar[g,t] - p[g,t-1] <= RU[g]*u[g,t-1] + SU[g]*v[g,t]);	# Constraint 6, ramp-up constraints and start-up mode
	@constraint(m_relax, ramp_down_consraint[g=1:nb_gen, t=1:T_max], pbar[g,t-1] - p[g,t] <= RD[g]*u[g,t] + SD[g]*w[g,t]);	# Constraint 7, ramp-down constraints and shutdown mode
	### Production >= Demand
	@constraint(m_relax, loads[t=1:T_max], sum( p[g,t] for g=1:nb_gen) == l[t]);
	### Prior data
	if length(u_prior) > 0
		for g=1:nb_gen
			JuMP.fix(u[g,0], u_prior[g]; force = true);
			JuMP.fix(v[g,0], v_prior[g]; force = true);
			JuMP.fix(w[g,0], w_prior[g]; force = true);
			JuMP.fix(p[g,0], p_prior[g]; force = true);
			JuMP.fix(pbar[g,0], pbar_prior[g]; force = true);
		end
	end
	### Objective function of master program
	@objective(m_relax, Min, sum( sum( NoLoadConsumption[g]*C[g]*u[g,t] + F[g]*v[g,t] + C[g]*p[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum( l[t] for t=1:T_max));
	### Start iterations of Bender decomposition
	for iter=1:iter_max
		### Generating cut : delta^T * ones + epsilon^T * p + mu^T * pbar + xi^T * u + alpha^T * v + sigma^T*w <= 0
		optimize!(m_relax);
		# Solving_time+=MOI.get(m_relax, MOI.SolveTime());
		Solving_time+=0;
		price 			= dual.(loads);
		push!(prices_vect,price);
		p_time_opt 		= value.(p).data;
		push!(p_vect,p_time_opt);
		pbar_time_opt 	= value.(pbar).data;
		u_opt 			= value.(u).data;
		v_opt 			= value.(v).data;
		w_opt 			= value.(w).data;
		obj_real 		= objective_value(m_relax);
		if Printer
			println();
			println("iter ",iter);
			println("obj : ",objective_value(m_relax));
			println("price : $(price)");
			println();
		end
		price_opt 		= price;
		obj_vec[iter] 	= obj_real;
		l_BD 			= value.(l);
		obj_real_BD 	= obj_real;

		### Check for generating cut
		cut_added = false;
		for g_index=1:length(G_c)
			g_considered 	= G_c[g_index];
			elem_CGLP 		= tab_CGLP[g_index];
			println();println("g_considered : ",g_considered);
			println("p_time_opt : $(p_time_opt[g_considered,:])");
			println("u_opt      : $(u_opt[g_considered,:])");
			println("v_opt      : $(v_opt[g_considered,:])");
			println("w_opt      : $(w_opt[g_considered,:])");
			for t=2:T_max+1
				println("$(u_opt[g_considered,t] - u_opt[g_considered,t-1]) = $(v_opt[g_considered,t] - w_opt[g_considered,t])");
			end
			
			for t=1:T_max
				set_normalized_rhs(elem_CGLP.power_output_global[t], p_time_opt[g_considered,t+1]);
				set_normalized_rhs(elem_CGLP.maximum_power_available[t], pbar_time_opt[g_considered,t+1]);
				set_normalized_rhs(elem_CGLP.commitment_status[t], u_opt[g_considered,t+1]);
				set_normalized_rhs(elem_CGLP.startup_status[t], v_opt[g_considered,t+1]);
				set_normalized_rhs(elem_CGLP.shutdown_status[t], w_opt[g_considered,t+1]);
			end
			### Normally no need to timpose those right-hand-side
			k_vec = [0]
			for k=1:length(k_vec)
				set_normalized_rhs(elem_CGLP.power_output_global_2[k], p_time_opt[g_considered,k_vec[k]+1]);
				set_normalized_rhs(elem_CGLP.maximum_power_available_2[k], pbar_time_opt[g_considered,k_vec[k]+1]);
				set_normalized_rhs(elem_CGLP.commitment_status_2[k], u_opt[g_considered,k_vec[k]+1]);
				set_normalized_rhs(elem_CGLP.startup_status_2[k], v_opt[g_considered,k_vec[k]+1]);
				set_normalized_rhs(elem_CGLP.shutdown_status_2[k], w_opt[g_considered,k_vec[k]+1]);
			end

			optimize!(elem_CGLP.model);
			# Solving_time+=MOI.get(elem_CGLP.model, MOI.SolveTime());
			Solving_time+=0;
			z_opt = value.(elem_CGLP.z);
			delta_opt 	= dual.(elem_CGLP.minimum_down_time_constraint);
			epsilon_opt = dual.(elem_CGLP.power_output_global);
			mu_opt 		= dual.(elem_CGLP.maximum_power_available);
			xi_opt 		= dual.(elem_CGLP.commitment_status);
			alpha_opt 	= dual.(elem_CGLP.startup_status);
			sigma_opt 	= dual.(elem_CGLP.shutdown_status);
			if z_opt>delta_criterion
				if Printer==true
					println("Add a cut. Master infeasible for $(g_considered)");
				end
				cut_added = true;
				nb_salves_solve+=1;
				for i in gather_gen[indices_gen[g_considered]]
					nb_cuts_tot+=1;
					@constraint(m_relax, sum(delta_opt[t] + epsilon_opt[t]*p[i,t] + mu_opt[t]*pbar[i,t] + xi_opt[t]*u[i,t] + alpha_opt[t]*v[i,t] + sigma_opt[t]*w[i,t] for t=1:T_max) <= 0);
				end
			end
		end
		if (cut_added==false) # Added no cut to none of the generator
			println("Stop at iteration ", iter);
			return (iter, price_opt, obj_vec, p_time_opt, pbar_time_opt, u_opt, v_opt, w_opt, l_BD, obj_real_BD, Solving_time, nb_cuts_tot, nb_salves_solve, prices_vect, p_vect);
		end
	end
	println("Stop at maximum iterations ", iter_max);
	return (iter_max, price_opt, obj_vec, p_time_opt, pbar_time_opt, u_opt, v_opt, w_opt, l_BD, obj_real_BD, Solving_time, nb_cuts_tot, nb_salves_solve, prices_vect, p_vect);
end