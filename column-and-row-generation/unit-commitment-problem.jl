"""
Compute set of generators G_c as described in Kneuven's paper, i.e. MaxRunCapacity[g]-MinRunCapacity[g]>minimum([RD[g] RU[g]]) && UT[g]>=2
G_c = Compute_G_c(MinRunCapacity, MaxRunCapacity, RU, RD, UT)
"""
function Compute_G_c(MinRunCapacity, MaxRunCapacity, RU, RD, UT)
	G_c = Int64[];
	for g=1:length(MinRunCapacity)
		if (MaxRunCapacity[g]-MinRunCapacity[g])>minimum([RD[g] RU[g]]) && UT[g]>=2
			append!(G_c,g)
		end
	end
	return G_c
end

"""
Matching program from Unit Commitment problem. We consider all generators.

This function is used to compute the uplifts of the generators.
"""
function Matching(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer=false)
	# m_matching = JuMP.direct_model(Gurobi.Optimizer(GUROBI_ENV)); # OutputFlag=0
	# m_matching = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
	# JuMP.set_optimizer_attribute(m_matching,"OutputFlag",0);
	m_matching = Model(HiGHS.Optimizer; add_bridges = false);
	set_attribute(m_matching, "output_flag", false)

	@variable(m_matching, 0 <= p[g=1:nb_gen, t=0:T_max+1]);		# power output of each generator g at time period t
	@variable(m_matching, 0 <= pbar[g=1:nb_gen, t=0:T_max+1]);	# maximum power available of each generator g at time period t
	@variable(m_matching, u[g=1:nb_gen, t=0:T_max+1], Bin);		# commitment status of geenrator g at time period t
	@variable(m_matching, v[g=1:nb_gen, t=0:T_max+1], Bin); 	# startup status of generator g at time period t
	@variable(m_matching, w[g=1:nb_gen, t=0:T_max+1], Bin);		# shutdown status of generator g at time period t

	@variable(m_matching, 0 <= l[t=1:T_max]);
	@constraint(m_matching, [t=1:T_max], l[t] <= data_demand[t]);			# upper bound for elasticity of the demand

	### The following constriants described the feasible set for technical constraints, the relaxation of 3-bin space
	@constraint(m_matching, logical_constraint[g=1:nb_gen, t=1:T_max+1], u[g,t] - u[g,t-1] == v[g,t] - w[g,t]);					# Constraint 2, logical constraint on the different status (commitment status, startup status, shutdown status)
	@constraint(m_matching, minimum_up_time[g=1:nb_gen, t=UT[g]:T_max+1], sum(v[g,i]  for i=t-UT[g]+1:t) <= u[g,t]);			# Constraint 3, minimum up time constraints
	@constraint(m_matching, minimum_down_time[g=1:nb_gen, t=DT[g]:T_max+1], sum(w[g,i]  for i=t-DT[g]+1:t) <= (1 - u[g,t])); 	# Constraint 4, minimum down time constraints
	
	@constraint(m_matching, generation_limits_1[g=1:nb_gen, t=0:T_max+1], MinRunCapacity[g]*u[g,t] <= p[g,t]);					# Constraint 5_1, generation limits constraints
	@constraint(m_matching, generation_limits_2[g=1:nb_gen, t=0:T_max+1], p[g,t] <= pbar[g,t]);									# Constraint 5_2, generation limits constraints
	@constraint(m_matching, generation_limits_3[g=1:nb_gen, t=0:T_max+1], pbar[g,t] <= MaxRunCapacity[g]*u[g,t]);				# Constraint 5_3, generation limits constraints
	@constraint(m_matching, ramp_up_constraint[g=1:nb_gen, t=1:T_max+1], pbar[g,t] - p[g,t-1] <= RU[g]*u[g,t-1] + SU[g]*v[g,t]);# Constraint 6, ramp-up constraints and start-up mode
	@constraint(m_matching, ramp_down_consraint[g=1:nb_gen, t=1:T_max+1], pbar[g,t-1] - p[g,t] <= RD[g]*u[g,t] + SD[g]*w[g,t]);	# Constraint 7, ramp-down constraints and shutdown mode
	
	### Production >= Demand
	@constraint(m_matching, loads[t=1:T_max], sum( p[g,t] for g=1:nb_gen) == l[t]);
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
	@objective(m_matching, Min, sum( sum( NoLoadConsumption[g]*C[g]*u[g,t] + F[g]*v[g,t] + C[g]*p[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum( l[t] for t=1:T_max));
	optimize!(m_matching);
	u_val = value.(u); v_val = value.(v); w_val = value.(w); p_val = value.(p); pbar_val = value.(pbar);
	l_val = value.(l);
	# Solving_time = MOI.get(m_matching, MOI.SolveTime());
	Solving_time = 0

	if Printer
		nb_con = sum([num_constraints(m_matching, i, j) for (i,j) in list_of_constraint_types(m_matching)]);
		nb_var = num_variables(m_matching);
		println("[Matching] nb_con : $(nb_con) and nb_var : $(nb_var)");
	end

	return (p_val.data, pbar_val.data, u_val.data, v_val.data, w_val.data, l_val, objective_value(m_matching), Solving_time);
end


function Dual_Lagrangian(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price)
	# m_matching = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
	m_matching = JuMP.direct_model(Gurobi.Optimizer(GUROBI_ENV)); # OutputFlag=0
	JuMP.set_optimizer_attribute(m_matching,"OutputFlag",0);
	# m_matching = Model(HiGHS.Optimizer);

	@variable(m_matching, 0 <= p[g=1:nb_gen, t=0:T_max+1]);		# power output of each generator g at time period t
	@variable(m_matching, 0 <= pbar[g=1:nb_gen, t=0:T_max+1]);	# maximum power available of each generator g at time period t
	@variable(m_matching, u[g=1:nb_gen, t=0:T_max+1], Bin);		# commitment status of geenrator g at time period t
	@variable(m_matching, v[g=1:nb_gen, t=0:T_max+1], Bin); 	# startup status of generator g at time period t
	@variable(m_matching, w[g=1:nb_gen, t=0:T_max+1], Bin);		# shutdown status of generator g at time period t

	@variable(m_matching, 0 <= l[t=1:T_max]);
	@constraint(m_matching, [t=1:T_max], l[t] <= data_demand[t]);			# upper bound for elasticity of the demand

	### The following constriants described the feasible set for technical constraints, the relaxation of 3-bin space
	@constraint(m_matching, logical_constraint[g=1:nb_gen, t=1:T_max+1], u[g,t] - u[g,t-1] == v[g,t] - w[g,t]);					# Constraint 2, logical constraint on the different status (commitment status, startup status, shutdown status)
	@constraint(m_matching, minimum_up_time[g=1:nb_gen, t=UT[g]:T_max+1], sum(v[g,i]  for i=t-UT[g]+1:t) <= u[g,t]);			# Constraint 3, minimum up time constraints
	@constraint(m_matching, minimum_down_time[g=1:nb_gen, t=DT[g]:T_max+1], sum(w[g,i]  for i=t-DT[g]+1:t) <= (1 - u[g,t])); 	# Constraint 4, minimum down time constraints
	
	@constraint(m_matching, generation_limits_1[g=1:nb_gen, t=0:T_max+1], MinRunCapacity[g]*u[g,t] <= p[g,t]);					# Constraint 5_1, generation limits constraints
	@constraint(m_matching, generation_limits_2[g=1:nb_gen, t=0:T_max+1], p[g,t] <= pbar[g,t]);									# Constraint 5_2, generation limits constraints
	@constraint(m_matching, generation_limits_3[g=1:nb_gen, t=0:T_max+1], pbar[g,t] <= MaxRunCapacity[g]*u[g,t]);				# Constraint 5_3, generation limits constraints
	@constraint(m_matching, ramp_up_constraint[g=1:nb_gen, t=1:T_max+1], pbar[g,t] - p[g,t-1] <= RU[g]*u[g,t-1] + SU[g]*v[g,t]);# Constraint 6, ramp-up constraints and start-up mode
	@constraint(m_matching, ramp_down_consraint[g=1:nb_gen, t=1:T_max+1], pbar[g,t-1] - p[g,t] <= RD[g]*u[g,t] + SD[g]*w[g,t]);	# Constraint 7, ramp-down constraints and shutdown mode
	
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
	@objective(m_matching, Min, sum( sum( NoLoadConsumption[g]*C[g]*u[g,t] + F[g]*v[g,t] + C[g]*p[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum( l[t] for t=1:T_max) - sum( price[t]*( sum( p[g,t] for g=1:nb_gen) - l[t] ) for t=1:T_max));
	optimize!(m_matching);
	u_val = value.(u); v_val = value.(v); w_val = value.(w); p_val = value.(p); pbar_val = value.(pbar);
	l_val = value.(l);
	# Solving_time = MOI.get(m_matching, MOI.SolveTime());
	Solving_time = 0
	return (p_val.data, pbar_val.data, u_val.data, v_val.data, w_val.data, l_val, objective_value(m_matching), Solving_time);
end



function Matching_NO_Posterior(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior)
	m_matching = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
	JuMP.set_optimizer_attribute(m_matching,"OutputFlag",0);

	@variable(m_matching, p[g=1:nb_gen, t=0:T_max], lower_bound = 0);		# power output of each generator g at time period t
	@variable(m_matching, pbar[g=1:nb_gen, t=0:T_max], lower_bound = 0);	# maximum power available of each generator g at time period t
	@variable(m_matching, u[g=1:nb_gen, t=0:T_max], Bin);		# commitment status of geenrator g at time period t
	@variable(m_matching, v[g=1:nb_gen, t=0:T_max], Bin); 	# startup status of generator g at time period t
	@variable(m_matching, w[g=1:nb_gen, t=0:T_max], Bin);		# shutdown status of generator g at time period t

	@variable(m_matching, 0 <= l[t=1:T_max]);
	@constraint(m_matching, [t=1:T_max], l[t] <= data_demand[t]);			# upper bound for elasticity of the demand

	### The following constriants described the feasible set for technical constraints, the relaxation of 3-bin space
	@constraint(m_matching, logical_constraint[g=1:nb_gen, t=1:T_max], u[g,t] - u[g,t-1] == v[g,t] - w[g,t]);					# Constraint 2, logical constraint on the different status (commitment status, startup status, shutdown status)
	@constraint(m_matching, minimum_up_time[g=1:nb_gen, t=UT[g]:T_max], sum(v[g,i]  for i=t-UT[g]+1:t) <= u[g,t]);			# Constraint 3, minimum up time constraints
	@constraint(m_matching, minimum_down_time[g=1:nb_gen, t=DT[g]:T_max], sum(w[g,i]  for i=t-DT[g]+1:t) <= (1 - u[g,t])); 	# Constraint 4, minimum down time constraints
	
	@constraint(m_matching, generation_limits_1[g=1:nb_gen, t=0:T_max], MinRunCapacity[g]*u[g,t] <= p[g,t]);					# Constraint 5_1, generation limits constraints
	@constraint(m_matching, generation_limits_2[g=1:nb_gen, t=0:T_max], p[g,t] <= pbar[g,t]);									# Constraint 5_2, generation limits constraints
	@constraint(m_matching, generation_limits_3[g=1:nb_gen, t=0:T_max], pbar[g,t] <= MaxRunCapacity[g]*u[g,t]);				# Constraint 5_3, generation limits constraints
	@constraint(m_matching, ramp_up_constraint[g=1:nb_gen, t=1:T_max], pbar[g,t] - p[g,t-1] <= RU[g]*u[g,t-1] + SU[g]*v[g,t]);# Constraint 6, ramp-up constraints and start-up mode
	@constraint(m_matching, ramp_down_consraint[g=1:nb_gen, t=1:T_max], pbar[g,t-1] - p[g,t] <= RD[g]*u[g,t] + SD[g]*w[g,t]);	# Constraint 7, ramp-down constraints and shutdown mode
	
	### Production >= Demand
	@constraint(m_matching, loads[t=1:T_max], sum( p[g,t] for g=1:nb_gen) == l[t]);
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
	### Objective function
	@objective(m_matching, Min, sum( sum( NoLoadConsumption[g]*C[g]*u[g,t] + F[g]*v[g,t] + C[g]*p[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum( l[t] for t=1:T_max));
	optimize!(m_matching);
	u_val = value.(u); v_val = value.(v); w_val = value.(w); p_val = value.(p); pbar_val = value.(pbar);
	l_val = value.(l);
	# Solving_time = MOI.get(m_matching, MOI.SolveTime());
	Solving_time = 0
	return (p_val.data, pbar_val.data, u_val.data, v_val.data, w_val.data, l_val, objective_value(m_matching), Solving_time);
end


"""
Maximum Profit for each generator g. We consider one generator at a time.

This function is used to compute the uplifts of the generator g.
"""
function MaxProfit_Producer(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, F, C, NoLoadConsumption, price, u_pri, v_pri, w_pri, p_pri, pbar_pri, u_post, v_post, w_post, p_post, pbar_post, g)
	if length(u_pri) > 0
		u_prior = u_pri[g];
		v_prior = v_pri[g];
		w_prior = w_pri[g];
		p_prior = p_pri[g];
		pbar_prior = pbar_pri[g];
	end
	if length(u_post) > 0
		u_posterior = u_post[g];
		v_posterior = v_post[g];
		w_posterior = w_post[g];
		p_posterior = p_post[g];
		pbar_posterior = pbar_post[g];
	end
	# m_max_profit = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
	m_max_profit = JuMP.direct_model(Gurobi.Optimizer(GUROBI_ENV)); # OutputFlag=0
	JuMP.set_optimizer_attribute(m_max_profit,"OutputFlag",0);

	@variable(m_max_profit, p[t=0:T_max+1], lower_bound = 0);		# power output of each generator g at time period t
	@variable(m_max_profit, pbar[t=0:T_max+1], lower_bound = 0);	# maximum power available of each generator g at time period t
	
	@variable(m_max_profit, u[t=0:T_max+1], Bin);		# commitment status of geenrator g at time period t
	@variable(m_max_profit, v[t=0:T_max+1], Bin); 		# startup status of generator g at time period t
	@variable(m_max_profit, w[t=0:T_max+1], Bin);		# shutdown status of generator g at time period t
	
	### The following constriants described the feasible set for technical constraints, the relaxation of 3-bin space
	@constraint(m_max_profit, logical_constraint[t=1:T_max+1], u[t] - u[t-1] == v[t] - w[t]);				# Constraint 2, logical constraint on the different status (commitment status, startup status, shutdown status)
	@constraint(m_max_profit, minimum_up_time[t=UT:T_max+1], sum(v[i] for i=t-UT+1:t) <= u[t]);				# Constraint 3, minimum up time constraints
	@constraint(m_max_profit, minimum_down_time[t=DT:T_max+1], sum(w[i] for i=t-DT+1:t) <= (1 - u[t])); 	# Constraint 4, minimum down time constraints
	
	@constraint(m_max_profit, generation_limits_1[t=0:T_max+1], MinRunCapacity*u[t] <= p[t]);				# Constraint 5_1, generation limits constraints
	@constraint(m_max_profit, generation_limits_2[t=0:T_max+1], p[t] - pbar[t] <= 0);						# Constraint 5_2, generation limits constraints
	@constraint(m_max_profit, generation_limits_3[t=0:T_max+1], pbar[t] <= MaxRunCapacity*u[t]);			# Constraint 5_3, generation limits constraints
	
	@constraint(m_max_profit, ramp_up_constraint[t=1:T_max+1], pbar[t] - p[t-1] <= RU*u[t-1] + SU*v[t]);	# Constraint 6, ramp-up constraints and start-up mode
	@constraint(m_max_profit, ramp_down_consraint[t=1:T_max+1], pbar[t-1] - p[t] <= RD*u[t] + SD*w[t]);		# Constraint 7, ramp-down constraints and shutdown mode

	### Prior data
	if length(u_pri) > 0
		JuMP.fix(u[0], u_prior; force = true);
		JuMP.fix(v[0], v_prior; force = true);
		JuMP.fix(w[0], w_prior; force = true);
		JuMP.fix(p[0], p_prior; force = true);
		JuMP.fix(pbar[0], pbar_prior; force = true);
	end
	### Posterior data
	if length(u_post) > 0
		JuMP.fix(u[T_max+1], u_posterior; force = true);
		JuMP.fix(v[T_max+1], v_posterior; force = true);
		JuMP.fix(w[T_max+1], w_posterior; force = true);
		JuMP.fix(p[T_max+1], p_posterior; force = true);
		JuMP.fix(pbar[T_max+1], pbar_posterior; force = true);
	end

	@objective(m_max_profit, Max, sum(price[t]*p[t] -  (NoLoadConsumption*C*u[t] + F*v[t] + C*p[t]) for t=1:T_max));
	optimize!(m_max_profit);
	p_val =  value.(p).data; pbar_val = value.(pbar).data; u_val = value.(u).data; v_val = value.(v).data; w_val = value.(w).data;
	return (objective_value(m_max_profit), p_val, pbar_val, u_val, v_val, w_val);
end

"""
Maximum Profit for each generator g. We consider one generator at a time.

This function is used to compute the uplifts of the generators.
"""
function MaxProfit_Producer_LP(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, F, C, NoLoadConsumption, price, u_pri, v_pri, w_pri, p_pri, pbar_pri, u_post, v_post, w_post, p_post, pbar_post, g)
	if length(u_pri) > 0
		u_prior = u_pri[g];
		v_prior = v_pri[g];
		w_prior = w_pri[g];
		p_prior = p_pri[g];
		pbar_prior = pbar_pri[g];
	end
	if length(u_post) > 0
		u_posterior = u_post[g];
		v_posterior = v_post[g];
		w_posterior = w_post[g];
		p_posterior = p_post[g];
		pbar_posterior = pbar_post[g];
	end
	m_max_profit = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
	JuMP.set_optimizer_attribute(m_max_profit,"OutputFlag",0);

	@variable(m_max_profit, p[t=0:T_max+1], lower_bound = 0);		# power output of each generator g at time period t
	@variable(m_max_profit, pbar[t=0:T_max+1], lower_bound = 0);	# maximum power available of each generator g at time period t
	@variable(m_max_profit, u[t=0:T_max+1], lower_bound = 0, upper_bound = 1);		# commitment status of geenrator g at time period t
	@variable(m_max_profit, v[t=0:T_max+1], lower_bound = 0, upper_bound = 1); 		# startup status of generator g at time period t
	@variable(m_max_profit, w[t=0:T_max+1], lower_bound = 0, upper_bound = 1);		# shutdown status of generator g at time period t
	### The following constriants described the feasible set for technical constraints, the relaxation of 3-bin space
	@constraint(m_max_profit, logical_constraint[t=1:T_max+1], u[t] - u[t-1] == v[t] - w[t]);				# Constraint 2, logical constraint on the different status (commitment status, startup status, shutdown status)
	@constraint(m_max_profit, minimum_up_time[t=UT:T_max+1], sum(v[i] for i=t-UT+1:t) <= u[t]);				# Constraint 3, minimum up time constraints
	@constraint(m_max_profit, minimum_down_time[t=DT:T_max+1], sum(w[i] for i=t-DT+1:t) <= (1 - u[t])); 	# Constraint 4, minimum down time constraints
	@constraint(m_max_profit, generation_limits_1[t=0:T_max+1], MinRunCapacity*u[t] <= p[t]);				# Constraint 5_1, generation limits constraints
	@constraint(m_max_profit, generation_limits_2[t=0:T_max+1], p[t] - pbar[t] <= 0);						# Constraint 5_2, generation limits constraints
	@constraint(m_max_profit, generation_limits_3[t=0:T_max+1], pbar[t] <= MaxRunCapacity*u[t]);			# Constraint 5_3, generation limits constraints
	@constraint(m_max_profit, ramp_up_constraint[t=1:T_max+1], pbar[t] - p[t-1] <= RU*u[t-1] + SU*v[t]);	# Constraint 6, ramp-up constraints and start-up mode
	@constraint(m_max_profit, ramp_down_consraint[t=1:T_max+1], pbar[t-1] - p[t] <= RD*u[t] + SD*w[t]);		# Constraint 7, ramp-down constraints and shutdown mode
	### Prior data
	if length(u_pri) > 0
		#@constraint(m_max_profit, u_prior_constraint, u[0]==u_prior);@constraint(m_max_profit, v_prior_constraint, v[0]==v_prior);@constraint(m_max_profit, w_prior_constraint, w[0]==w_prior);@constraint(m_max_profit, p_prior_constraint, p[0]==p_prior);@constraint(m_max_profit, pbar_prior_constraint, pbar[0]==pbar_prior);
		JuMP.fix(u[0], u_prior; force = true);
		JuMP.fix(v[0], v_prior; force = true);
		JuMP.fix(w[0], w_prior; force = true);
		JuMP.fix(p[0], p_prior; force = true);
		JuMP.fix(pbar[0], pbar_prior; force = true);
	end
	### Posterior data
	if length(u_post) > 0
		#@constraint(m_max_profit, u_posterior_constraint, u[T_max+1]==u_posterior);@constraint(m_max_profit, v_posterior_constraint, v[T_max+1]==v_posterior);@constraint(m_max_profit, w_posterior_constraint, w[T_max+1]==w_posterior);@constraint(m_max_profit, p_posterior_constraint, p[T_max+1]==p_posterior);@constraint(m_max_profit, pbar_posterior_constraint, pbar[T_max+1]==pbar_posterior);
		JuMP.fix(u[T_max+1], u_posterior; force = true);
		JuMP.fix(v[T_max+1], v_posterior; force = true);
		JuMP.fix(w[T_max+1], w_posterior; force = true);
		JuMP.fix(p[T_max+1], p_posterior; force = true);
		JuMP.fix(pbar[T_max+1], pbar_posterior; force = true);
	end
	@objective(m_max_profit, Max, sum(price[t]*p[t] -  (NoLoadConsumption*C*u[t] + F*v[t] + C*p[t]) for t=1:T_max));
	optimize!(m_max_profit);
	p_val =  value.(p).data; pbar_val = value.(pbar).data; u_val = value.(u).data; v_val = value.(v).data; w_val = value.(w).data;
	return (objective_value(m_max_profit), p_val, pbar_val, u_val, v_val, w_val);
end


function MaxProfit_Consumer(data_demand, T_max, price, VOLL)
	# m_max_profit_consumer = JuMP.direct_model(Gurobi.Optimizer());
	m_max_profit_consumer = JuMP.direct_model(Gurobi.Optimizer(GUROBI_ENV));
	JuMP.set_optimizer_attribute(m_max_profit_consumer,"OutputFlag",0);
	@variable(m_max_profit_consumer, l[t=1:T_max], lower_bound = 0);
	@constraint(m_max_profit_consumer, upper_bound[t=1:T_max], l[t]<=data_demand[t]);
	@objective(m_max_profit_consumer, Max, sum( (VOLL - price[t])*l[t] for t=1:T_max));
	optimize!(m_max_profit_consumer);
	return (objective_value(m_max_profit_consumer),value.(l));
end
