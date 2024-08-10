

"""
Test LP_Relaxation when dualize balance constraint.
"""
function Sanity_LP_Relaxation(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price)
	m_relax = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
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
	@constraint(m_relax, generation_limits_2[g=1:nb_gen, t=0:T_max+1], p[g,t] - pbar[g,t] <= 0);							# Constraint 5_2, generation limits constraints
	@constraint(m_relax, generation_limits_3[g=1:nb_gen, t=0:T_max+1], pbar[g,t] <= MaxRunCapacity[g]*u[g,t]);				# Constraint 5_3, generation limits constraints
	
	@constraint(m_relax, ramp_up_constraint[g=1:nb_gen, t=1:T_max+1], pbar[g,t] - p[g,t-1] <= RU[g]*u[g,t-1] + SU[g]*v[g,t]);	# Constraint 6, ramp-up constraints and start-up mode
	@constraint(m_relax, ramp_down_consraint[g=1:nb_gen, t=1:T_max+1], pbar[g,t-1] - p[g,t] <= RD[g]*u[g,t] + SD[g]*w[g,t]);	# Constraint 7, ramp-down constraints and shutdown mode
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
	### Objective function
	#@objective(m_relax, Min, sum( sum( NoLoadConsumption[g]*C[g]*u[g,t] + F[g]*v[g,t] + C[g]*p[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum( x[t]*data_demand[t] for t=1:T_max) - sum( price[t]*( sum(p[g,t] for g=1:nb_gen) - data_demand[t]*x[t] ) for t=1:T_max));
	@objective(m_relax, Min, sum( sum( NoLoadConsumption[g]*C[g]*u[g,t] + F[g]*v[g,t] + C[g]*p[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum( l[t] for t=1:T_max) - sum( price[t]*( sum(p[g,t] for g=1:nb_gen) - l[t] ) for t=1:T_max));
	optimize!(m_relax);
	return (objective_value(m_relax));
end

"""
Test Extended Formulation when dualize balance constraint.
"""
function Sanity_Extended_Formulation(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price)
	(A,B,nb_intervals_gen) = set_T(UT, T_max, nb_gen);
	m = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
	JuMP.set_optimizer_attribute(m,"OutputFlag",0);
	# Variables
	@variable(m, 0 <= p[g=1:nb_gen, i=1:nb_intervals_gen[g], t=0:T_max+1]);	# Power output vector of the generators		p[generator g, interval(a,b) i, time t]
	@variable(m, 0 <= pbar[g=1:nb_gen, i=1:nb_intervals_gen[g], t=0:T_max+1]);# Maximum power available from generators	pbar[generator g, interval(a,b) i, time t]
	@variable(m, 0 <= gamma[g=1:nb_gen, i=1:nb_intervals_gen[g]]);			# Variable that indicates whether a vertex (interval) is in the packing or not   gamma[generator g, interval(a,b) i]
	@variable(m, 0 <= p_time[g=1:nb_gen, t=0:T_max+1]);						# (to simplify computations) production of each generator on each time period
	@variable(m, 0 <= pbar_time[g=1:nb_gen, t=0:T_max+1]);						# (to simplify computations) maximum power available of each generator on each time period
	@variable(m, 0 <= u[g=1:nb_gen, t=0:T_max+1]); 							# Commitment status variables
	@variable(m, 0 <= v[g=1:nb_gen, t=0:T_max+1]); 							# Startup status variables
	@variable(m, 0 <= w[g=1:nb_gen, t=0:T_max+1]); 							# Shutdown status variables
	
	@variable(m, 0 <= l[t=1:T_max]); # Modelling an elastic demand
	@constraint(m, [t=1:T_max], l[t] <= data_demand[t]);
	@constraint(m, [g=1:nb_gen, i=1:nb_intervals_gen[g]], gamma[g,i] <=1); # By security
	### Feasible Dispatch Polytope (equations 13 from source)
	# For interval i of generator g, a = A[g,i] and b=A[g,i]
	# 13a and 13b
	for g=1:nb_gen
		for i=1:nb_intervals_gen[g]
			for t=0:T_max+1
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
	### The followin loop (can be) is useful when considering generators already on for prior and posterior intervals.
	# 13f For this constraint assume that a generator can already produces before the set [T_max] of time periods.
	@constraint(m, limit_power_jumps_up_1[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], pbar[g,i,t] - p[g,i,t-1] <= gamma[g,i] * RU[g]);
	@constraint(m, limit_power_jumps_up_2[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], pbar[g,i,t] - p[g,i,t-1] <= gamma[g,i] * (SD[g] + (B[g,i] - t)*RD[g] - MinRunCapacity[g]));
	### The followin loop (can be) is useful when considering generators already on for prior and posterior intervals.
	# 13g
	@constraint(m, limit_power_jumps_down_1[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], pbar[g,i,t-1] - p[g,i,t] <= gamma[g,i] * RD[g]);
	@constraint(m, limit_power_jumps_down_2[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], pbar[g,i,t-1] - p[g,i,t] <= gamma[g,i] * (SU[g] + (t - A[g,i])*RU[g] - MinRunCapacity[g]));
	### The followin loop (can be) is useful when considering generators already on for prior and posterior intervals.
	### Packing Dispatch Polytopes
	# 15e : to eliminate the impossible combinations (to have the convex hull)
	@constraint(m, minimum_down_time_constraint[g=1:nb_gen, t=1:T_max], sum(gamma[g,i] for i=1:nb_intervals_gen[g] if (A[g,i] <= t && t <= B[g,i]+DT[g])) <= 1);
	# 15b : to be able to work with production on time instead of intervals in the objective function
	@constraint(m, production_time[g=1:nb_gen, t=1:T_max], p_time[g,t] == sum(p[g,i,t]  for i=1:nb_intervals_gen[g]));
	# 15c : to be able to work with production on time instead of intervals in the objective function
	@constraint(m, productionbar_time[g=1:nb_gen, t=1:T_max], pbar_time[g,t] == sum(pbar[g,i,t]  for i=1:nb_intervals_gen[g]));
	### Binary variable of UC (commitment status(on/off), startup status, shutdown status)
	@constraint(m, commitment_status[g=1:nb_gen,t=0:T_max+1], sum(gamma[g,i] for i=1:nb_intervals_gen[g] if (A[g,i] <= t && t <= B[g,i])) == u[g,t]);
	@constraint(m, startup_status[g=1:nb_gen, t=0:T_max+1], sum(gamma[g,i]  for i=1:nb_intervals_gen[g] if (t == A[g,i])) == v[g,t]);
	@constraint(m, shutdown_status[g=1:nb_gen, t=0:T_max+1], sum(gamma[g,i] for i=1:nb_intervals_gen[g] if (t == B[g,i]+1)) == w[g,t]);
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
	# @objective(m, Min, sum( sum( NoLoadConsumption[g]*C[g]*u[g,t] + F[g]*v[g,t] + C[g]*p_time[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum(x[t]*data_demand[t] for t=1:T_max) - sum( price[t]*(sum( p_time[g,t] for g=1:nb_gen) - data_demand[t]*x[t]) for t=1:T_max));
	@objective(m, Min, sum( sum( NoLoadConsumption[g]*C[g]*u[g,t] + F[g]*v[g,t] + C[g]*p_time[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum(l[t] for t=1:T_max) - sum( price[t]*(sum( p_time[g,t] for g=1:nb_gen) - l[t]) for t=1:T_max));
	optimize!(m);
	return (objective_value(m));
end


function sanity_check_MaxProfit(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price)
	m = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
	JuMP.set_optimizer_attribute(m,"OutputFlag",0);

	@variable(m, 0 <= p[g=1:nb_gen, t=0:T_max+1]);		# power output of each generator g at time period t
	@variable(m, 0 <= pbar[g=1:nb_gen, t=0:T_max+1]);	# maximum power available of each generator g at time period t
	@variable(m, u[g=1:nb_gen, t=0:T_max+1], Bin);		# commitment status of geenrator g at time period t
	@variable(m, v[g=1:nb_gen, t=0:T_max+1], Bin); 		# startup status of generator g at time period t
	@variable(m, w[g=1:nb_gen, t=0:T_max+1], Bin);		# shutdown status of generator g at time period t
	@variable(m, 0 <= l[t=1:T_max]);
	@constraint(m, [t=1:T_max], l[t] <= data_demand[t]);				# upper bound for elasticity of the demand
	### The following constriants described the feasible set for technical constraints, the relaxation of 3-bin space
	@constraint(m, logical_constraint[g=1:nb_gen, t=1:T_max+1], u[g,t] - u[g,t-1] == v[g,t] - w[g,t]);					# Constraint 2, logical constraint on the different status (commitment status, startup status, shutdown status)
	@constraint(m, minimum_up_time[g=1:nb_gen, t=UT[g]:T_max+1], sum(v[g,i]  for i=t-UT[g]+1:t) <= u[g,t]);			# Constraint 3, minimum up time constraints
	@constraint(m, minimum_down_time[g=1:nb_gen, t=DT[g]:T_max+1], sum(w[g,i]  for i=t-DT[g]+1:t) <= (1 - u[g,t])); 	# Constraint 4, minimum down time constraints
	
	@constraint(m, generation_limits_1[g=1:nb_gen, t=0:T_max+1], MinRunCapacity[g]*u[g,t] <= p[g,t]);					# Constraint 5_1, generation limits constraints
	@constraint(m, generation_limits_2[g=1:nb_gen, t=0:T_max+1], p[g,t] <= pbar[g,t]);									# Constraint 5_2, generation limits constraints
	@constraint(m, generation_limits_3[g=1:nb_gen, t=0:T_max+1], pbar[g,t] <= MaxRunCapacity[g]*u[g,t]);				# Constraint 5_3, generation limits constraints
	@constraint(m, ramp_up_constraint[g=1:nb_gen, t=1:T_max+1], pbar[g,t] - p[g,t-1] <= RU[g]*u[g,t-1] + SU[g]*v[g,t]);# Constraint 6, ramp-up constraints and start-up mode
	@constraint(m, ramp_down_consraint[g=1:nb_gen, t=1:T_max+1], pbar[g,t-1] - p[g,t] <= RD[g]*u[g,t] + SD[g]*w[g,t]);	# Constraint 7, ramp-down constraints and shutdown mode
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
	#@objective(m, Max, sum( sum( price[t]*p[g,t] - (NoLoadConsumption[g]*C[g]*u[g,t] + F[g]*v[g,t] + C[g]*p[g,t]) for t=1:T_max) for g=1:nb_gen) + sum( (VOLL-price[t])*x[t]*data_demand[t] for t=1:T_max) );
	@objective(m, Max, sum( sum( price[t]*p[g,t] - (NoLoadConsumption[g]*C[g]*u[g,t] + F[g]*v[g,t] + C[g]*p[g,t]) for t=1:T_max) for g=1:nb_gen) + sum( (VOLL-price[t])*l[t] for t=1:T_max) );
	optimize!(m);
	return objective_value(m);
end


