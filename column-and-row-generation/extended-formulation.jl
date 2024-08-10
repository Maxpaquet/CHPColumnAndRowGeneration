############################################################################################################
############################################################################################################
########################################### Extended Formulation ###########################################
############################################################################################################
############################################################################################################

"""
@in : 	- MinRunCapacity (array of size 1 x nb_gen) : minimum run capacity of each generator
		- MaxRunCapacity (array of size 1 x nb_gen) : maximum run capacity of each generator
		- RU (array of size 1 x nb_gen) : maximum ramp-up constraint
		- RD (array of size 1 x nb_gen) : maximum ramp-down constraint
		- UT (array of size 1 x nb_gen) : minimum up time constraint
		- DT (array of size 1 x nb_gen) : minimum down time cosntraint
		- SU (array of size 1 x nb_gen) : start-up levels
		- SD (array of size 1 x nb_gen) : shut-down levels
		- T_max (scalar) : number of considered periods
		- nb_gen (scalar) : number of considered generators
		- F (array of size 1 x nb_gen) : fixed start up costs
		- C (array of size 1 x nb_gen) : marginal cost
		- NoLoadConsumption (array of size 1 x nb_gen) : fixed cost to pay when generator g produces something
		- data_demand (array of size 1 x T_max) : demand that production has to meet
		- VOLL (scalar) : Value Of Lost Load (normally it should be 3000[€/MWh])
		- u_prior (array of size 1 x nb_gen) : Commitment status prior to planning period
		- v_prior (array of size 1 x nb_gen) : Startup status prior to planning period
		- w_prior (array of size 1 x nb_gen) : Shutdown status prior to planning period
		- p_prior (array of size 1 x nb_gen) : Power output vector of the generators prior to planning period
		- pbar_prior (array of size 1 x nb_gen) : Maximum power available from generators prior to planning period

@out :	- p (array of size nb_gen x nb_intervals_gen[g] x T_max) : power output vector of the generators g on interval(a,b) i at time t
		- pbar (array of size nb_gen x nb_intervals_gen[g] x T_max) : maximum power available from generators g on interval(a,b) i at time t
		- gamma (array of size nb_gen x nb_intervals_gen[g]) : variable that indicates whether a vertex (interval(a,b) i) is in the packing or not
		- p_global (array of size 1 x nb_gen) : production global of each generator on all the time periods
		- pbar_global (array of size 1 x nb_gen) : maximum power available of each generator on all the time periods
		- p_time (array of size nb_gen x nb_intervals_gen[g]) : production of each generator g on each time period t
		- pbar_time (array of size nb_gen x nb_intervals_gen[g]) : maximum power available of each generator g on each time period t
		- u (array of size nb_gen x T_max) : commitment status variable
		- v (array of size nb_gen x T_max) : startup status variable
		- w (array of size nb_gen x T_mas) : shutdown status variable
"""
function Extended_Formulation(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer=false)
	(A,B,nb_intervals_gen) = set_T(UT, T_max, nb_gen);
	# m = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
	# m = JuMP.direct_model(Gurobi.Optimizer(GUROBI_ENV)); # OutputFlag=0
	# JuMP.set_optimizer_attribute(m,"OutputFlag",0);
	model = Model(HiGHS.Optimizer; add_bridges = false);
	set_attribute(model, "output_flag", false)

	#JuMP.set_optimizer_attribute(m,"FeasibilityTol",1e-9);
	#JuMP.set_optimizer_attribute(m,"OptimalityTol",1e-9);
	#JuMP.set_optimizer_attribute(m,"BarConvTol",1e-8);

	# Variables
	# Power output vector of the generators		p[generator g, interval(a,b) i, time t]
	@variable(model, p[g=1:nb_gen, i=1:nb_intervals_gen[g], t=0:T_max+1], lower_bound = 0);
	@variable(model, pbar[g=1:nb_gen, i=1:nb_intervals_gen[g], t=0:T_max+1], lower_bound = 0);# Maximum power available from generators	pbar[generator g, interval(a,b) i, time t]
	@variable(model, gamma[g=1:nb_gen, i=1:nb_intervals_gen[g]], lower_bound = 0, upper_bound = 1);	# Variable that indicates whether a vertex (interval) is in the packing or not   gamma[generator g, interval(a,b) i]
	@variable(model, p_time[g=1:nb_gen, t=0:T_max+1], lower_bound = 0);	# (to simplify computations) production of each generator on each time period
	@variable(model, pbar_time[g=1:nb_gen, t=0:T_max+1], lower_bound = 0); # (to simplify computations) maximum power available of each generator on each time period
	@variable(model, u[g=1:nb_gen, t=0:T_max+1], lower_bound = 0, upper_bound = 1); # Commitment status variables
	@variable(model, v[g=1:nb_gen, t=0:T_max+1], lower_bound = 0, upper_bound = 1); # Startup status variables
	@variable(model, w[g=1:nb_gen, t=0:T_max+1], lower_bound = 0, upper_bound = 1); 							# Shutdown status variables
	@variable(model, l[t=1:T_max], lower_bound = 0); # Modelling an inelastic demand
	@constraint(model, [t=1:T_max], l[t] <= data_demand[t]);
	### Feasible Dispatch Polytope (equations 13 from source)
	# For interval i of generator g, a = A[g,i] and b=A[g,i]
	# 13a and 13b
	for g=1:nb_gen
		for i=1:nb_intervals_gen[g]
			for t=0:T_max+1
				if t<A[g,i] || t>B[g,i]
					@constraint(model, p[g,i,t] <= 0);
					@constraint(model, pbar[g,i,t] <= 0);
				end
			end
		end
	end
	# 13c
	@constraint(model, min_output[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]:B[g,i] ], -p[g,i,t] <= -gamma[g,i] * MinRunCapacity[g]);
	# 13d
	@constraint(model, max_output[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]:B[g,i] ], p[g,i,t] - pbar[g,i,t] <= 0);
	# 13e
	@constraint(model, upper_bound_on_power_output_1[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]:B[g,i] ] , pbar[g,i,t] <= gamma[g,i] * MaxRunCapacity[g]);
	@constraint(model, upper_bound_on_power_output_2[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]:B[g,i] ] , pbar[g,i,t] <= gamma[g,i] * (SU[g] + (t - A[g,i])*RU[g] ));
	@constraint(model, upper_bound_on_power_output_3[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]:B[g,i] ] , pbar[g,i,t] <= gamma[g,i] * (SD[g] + (B[g,i] - t)*RD[g] ));
	### The followin loop (can be) is useful when considering generators already on for prior and posterior intervals.
	#=for g=1:nb_gen
		for i=1:nb_intervals_gen[g]
			if A[g,i]!=0
				for t=A[g,i]:B[g,i]
					@constraint(m, pbar[g,i,t] <= gamma[g,i] * (SU[g] + (t - A[g,i])*RU[g])); # 13e.2
				end
			end
			if B[g,i]!=T_max+1
				for t=A[g,i]:B[g,i]
					@constraint(m, pbar[g,i,t] <= gamma[g,i] * (SD[g] + (B[g,i] - t)*RD[g])); # 13e.3
				end
			end
		end
	end=#
	# 13f For this constraint assume that a generator can already produces before the set [T_max] of time periods.
	@constraint(model, limit_power_jumps_up_1[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], pbar[g,i,t] - p[g,i,t-1] <= gamma[g,i] * RU[g]);
	@constraint(model, limit_power_jumps_up_2[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], pbar[g,i,t] - p[g,i,t-1] <= gamma[g,i] * (SD[g] + (B[g,i] - t)*RD[g] - MinRunCapacity[g]));
	### The followin loop (can be) is useful when considering generators already on for prior and posterior intervals.
	#=for g=1:nb_gen
		for i=1:nb_intervals_gen[g]
			if B[g,i]!=T_max+1
				for t=A[g,i]+1:B[g,i]
					@constraint(m, pbar[g,i,t-1] - p[g,i,t] <= gamma[g,i] * (SU[g] + (t - A[g,i])*RU[g] - MinRunCapacity[g])); # 13f.2
				end
			end
		end
	end=#
	# 13g
	@constraint(model, limit_power_jumps_down_1[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], pbar[g,i,t-1] - p[g,i,t] <= gamma[g,i] * RD[g]);
	@constraint(model, limit_power_jumps_down_2[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], pbar[g,i,t-1] - p[g,i,t] <= gamma[g,i] * (SU[g] + (t - A[g,i])*RU[g] - MinRunCapacity[g]));
	### The followin loop (can be) is useful when considering generators already on for prior and posterior intervals.
	#=for g=1:nb_gen
		for i=1:nb_intervals_gen[g]
			if A[g,i]!=0
				for t=A[g,i]+1:B[g,i]
					@constraint(m, pbar[g,i,t-1] - p[g,i,t] <= gamma[g,i] * (SU[g] + (t - A[g,i])*RU[g] - MinRunCapacity[g])); # 13g.2
				end
			end
		end
	end=#
	### Packing Dispatch Polytopes
	# 15e : to eliminate the impossible combinations (to have the convex hull)
	@constraint(model, minimum_down_time_constraint[g=1:nb_gen, t=1:T_max], sum(gamma[g,i] for i=1:nb_intervals_gen[g] if (A[g,i] <= t && t <= B[g,i]+DT[g])) <= 1);
	# 15b : to be able to work with production on time instead of intervals in the objective function
	@constraint(model, production_time[g=1:nb_gen, t=1:T_max], p_time[g,t] == sum(p[g,i,t]  for i=1:nb_intervals_gen[g]));
	# 15c : to be able to work with production on time instead of intervals in the objective function
	@constraint(model, productionbar_time[g=1:nb_gen, t=1:T_max], pbar_time[g,t] == sum(pbar[g,i,t]  for i=1:nb_intervals_gen[g]));
	### Binary variable of UC (commitment status(on/off), startup status, shutdown status)
	@constraint(model, commitment_status[g=1:nb_gen,t=0:T_max+1], sum(gamma[g,i] for i=1:nb_intervals_gen[g] if (A[g,i] <= t && t <= B[g,i])) == u[g,t]);
	@constraint(model, startup_status[g=1:nb_gen, t=0:T_max+1], sum(gamma[g,i]  for i=1:nb_intervals_gen[g] if (t == A[g,i])) == v[g,t]);
	@constraint(model, shutdown_status[g=1:nb_gen, t=0:T_max+1], sum(gamma[g,i] for i=1:nb_intervals_gen[g] if (t == B[g,i]+1)) == w[g,t]);

	### Production >= Demand
	@constraint(model, loads[t=1:T_max], sum( p_time[g,t] for g=1:nb_gen) == l[t]);
	### Prior data
	if length(u_prior) > 0
		for g=1:nb_gen
			JuMP.fix(u[g,0], u_prior[g]; force = true);
			JuMP.fix(v[g,0], v_prior[g]; force = true);
			JuMP.fix(w[g,0], w_prior[g]; force = true);
			JuMP.fix(p_time[g,0], p_prior[g]; force = true);
			JuMP.fix(pbar_time[g,0], pbar_prior[g]; force = true);
		end
	end
	### Posterior data
	if length(u_posterior) > 0
		for g=1:nb_gen
			JuMP.fix(u[g,T_max+1],u_posterior[g]; force = true);
			JuMP.fix(v[g,T_max+1],v_posterior[g]; force = true);
			JuMP.fix(w[g,T_max+1],w_posterior[g]; force = true);
			JuMP.fix(p_time[g,T_max+1],p_posterior[g]; force = true);
			JuMP.fix(pbar_time[g,T_max+1],pbar_posterior[g]; force = true);
		end
	end
	### Objective function
	@objective(model, Min, sum( sum( NoLoadConsumption[g]*C[g]*u[g,t] + F[g]*v[g,t] + C[g]*p_time[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum(l[t] for t=1:T_max));
	optimize!(model);
	# Solving_time = MOI.get(m, MOI.SolveTime());
	Solving_time = 0;
	price = dual.(loads);
	u_val = value.(u);
	v_val = value.(v);
	p_time_val = value.(p_time);
	l_val = value.(l);

	if Printer
		nb_con = sum([num_constraints(model, i, j) for (i,j) in list_of_constraint_types(model)]);
		nb_var = num_variables(model);
		println("[EF] nb_con : $(nb_con) and nb_var : $(nb_var)");
	end
	nb_con = sum([num_constraints(model, i, j) for (i,j) in list_of_constraint_types(model)]);
	nb_var = num_variables(model);

	return (value.(p).data, value.(pbar).data, value.(gamma), value.(p_time).data, value.(pbar_time).data, value.(u).data, value.(v).data, value.(w).data, price, l_val, objective_value(model), Solving_time, nb_var, nb_con);
end


function Extended_Formulation_NO_Posterior(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior)
	(A,B,nb_intervals_gen) = set_T(UT, T_max, nb_gen);
	m = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
	JuMP.set_optimizer_attribute(m,"OutputFlag",0);
	#JuMP.set_optimizer_attribute(m,"FeasibilityTol",1e-9);
	#JuMP.set_optimizer_attribute(m,"OptimalityTol",1e-9);
	#JuMP.set_optimizer_attribute(m,"BarConvTol",1e-8);
	# Variables
	@variable(m, p[g=1:nb_gen, i=1:nb_intervals_gen[g], t=0:T_max], lower_bound = 0);	# Power output vector of the generators		p[generator g, interval(a,b) i, time t]
	@variable(m, pbar[g=1:nb_gen, i=1:nb_intervals_gen[g], t=0:T_max], lower_bound = 0);# Maximum power available from generators	pbar[generator g, interval(a,b) i, time t]
	@variable(m, gamma[g=1:nb_gen, i=1:nb_intervals_gen[g]], lower_bound = 0, upper_bound = 1);			# Variable that indicates whether a vertex (interval) is in the packing or not   gamma[generator g, interval(a,b) i]
	@variable(m, p_time[g=1:nb_gen, t=0:T_max], lower_bound = 0);						# (to simplify computations) production of each generator on each time period
	@variable(m, pbar_time[g=1:nb_gen, t=0:T_max], lower_bound = 0);						# (to simplify computations) maximum power available of each generator on each time period
	@variable(m, u[g=1:nb_gen, t=0:T_max], lower_bound = 0, upper_bound = 1); 							# Commitment status variables
	@variable(m, v[g=1:nb_gen, t=0:T_max], lower_bound = 0, upper_bound = 1); 							# Startup status variables
	@variable(m, w[g=1:nb_gen, t=0:T_max], lower_bound = 0, upper_bound = 1); 							# Shutdown status variables
	@variable(m, l[t=1:T_max], lower_bound = 0); # Modelling an inelastic demand
	@constraint(m, [t=1:T_max], l[t] <= data_demand[t]);
	### Feasible Dispatch Polytope (equations 13 from source)
	# For interval i of generator g, a = A[g,i] and b=A[g,i]
	# 13a and 13b
	for g=1:nb_gen
		for i=1:nb_intervals_gen[g]
			for t=0:T_max
				if t<A[g,i] || t>B[g,i]
					@constraint(m, p[g,i,t] <= 0);
					@constraint(m, pbar[g,i,t] <= 0);
				end
			end
		end
	end
	# 13c
	@constraint(m, min_output[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]:B[g,i] ], -p[g,i,t] <= -gamma[g,i] * MinRunCapacity[g]);
	# 13d
	@constraint(m, max_output[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]:B[g,i] ], p[g,i,t] - pbar[g,i,t] <= 0);
	# 13e
	@constraint(m, upper_bound_on_power_output_1[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]:B[g,i] ] , pbar[g,i,t] <= gamma[g,i] * MaxRunCapacity[g]);
	@constraint(m, upper_bound_on_power_output_2[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]:B[g,i] ] , pbar[g,i,t] <= gamma[g,i] * (SU[g] + (t - A[g,i])*RU[g] ));
	@constraint(m, upper_bound_on_power_output_3[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]:B[g,i] ] , pbar[g,i,t] <= gamma[g,i] * (SD[g] + (B[g,i] - t)*RD[g] ));
	# 13f For this constraint assume that a generator can already produces before the set [T_max] of time periods.
	@constraint(m, limit_power_jumps_up_1[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], pbar[g,i,t] - p[g,i,t-1] <= gamma[g,i] * RU[g]);
	@constraint(m, limit_power_jumps_up_2[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], pbar[g,i,t] - p[g,i,t-1] <= gamma[g,i] * (SD[g] + (B[g,i] - t)*RD[g] - MinRunCapacity[g]));
	# 13g
	@constraint(m, limit_power_jumps_down_1[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], pbar[g,i,t-1] - p[g,i,t] <= gamma[g,i] * RD[g]);
	@constraint(m, limit_power_jumps_down_2[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], pbar[g,i,t-1] - p[g,i,t] <= gamma[g,i] * (SU[g] + (t - A[g,i])*RU[g] - MinRunCapacity[g]));
	### Packing Dispatch Polytopes
	# 15e : to eliminate the impossible combinations (to have the convex hull)
	@constraint(m, minimum_down_time_constraint[g=1:nb_gen, t=1:T_max], sum(gamma[g,i] for i=1:nb_intervals_gen[g] if (A[g,i] <= t && t <= B[g,i]+DT[g])) <= 1);
	# 15b : to be able to work with production on time instead of intervals in the objective function
	@constraint(m, production_time[g=1:nb_gen, t=1:T_max], p_time[g,t] == sum(p[g,i,t]  for i=1:nb_intervals_gen[g]));
	# 15c : to be able to work with production on time instead of intervals in the objective function
	@constraint(m, productionbar_time[g=1:nb_gen, t=1:T_max], pbar_time[g,t] == sum(pbar[g,i,t]  for i=1:nb_intervals_gen[g]));
	
	### Binary variable of UC (commitment status(on/off), startup status, shutdown status)
	@constraint(m, commitment_status[g=1:nb_gen,t=0:T_max], sum(gamma[g,i] for i=1:nb_intervals_gen[g] if (A[g,i] <= t && t <= B[g,i])) == u[g,t]);
	@constraint(m, startup_status[g=1:nb_gen, t=0:T_max], sum(gamma[g,i]  for i=1:nb_intervals_gen[g] if (t == A[g,i])) == v[g,t]);
	@constraint(m, shutdown_status[g=1:nb_gen, t=0:T_max], sum(gamma[g,i] for i=1:nb_intervals_gen[g] if (t == B[g,i]+1)) == w[g,t]);

	### Production >= Demand
	@constraint(m, loads[t=1:T_max], sum( p_time[g,t] for g=1:nb_gen) == l[t]);
	### Prior data
	if length(u_prior) > 0
		for g=1:nb_gen
			JuMP.fix(u[g,0], u_prior[g]; force = true);
			JuMP.fix(v[g,0], v_prior[g]; force = true);
			JuMP.fix(w[g,0], w_prior[g]; force = true);
			JuMP.fix(p_time[g,0], p_prior[g]; force = true);
			JuMP.fix(pbar_time[g,0], pbar_prior[g]; force = true);
		end
	end
	### Objective function
	@objective(m, Min, sum( sum( NoLoadConsumption[g]*C[g]*u[g,t] + F[g]*v[g,t] + C[g]*p_time[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum(l[t] for t=1:T_max));
	optimize!(m);
	price = dual.(loads);
	u_val = value.(u); v_val = value.(v); p_time_val = value.(p_time); l_val = value.(l);
	# Solving_time = MOI.get(m, MOI.SolveTime());
	Solving_time = 0
	return (value.(p).data, value.(pbar).data, value.(gamma), value.(p_time).data, value.(pbar_time).data, value.(u).data, value.(v).data, value.(w).data, price, l_val, objective_value(m), Solving_time);
end