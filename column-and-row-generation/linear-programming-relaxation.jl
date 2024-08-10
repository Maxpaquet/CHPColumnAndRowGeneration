#####################################################################################################
#####################################################################################################
########################################### LP Relaxation ###########################################
#####################################################################################################
#####################################################################################################

function LP_Relaxation(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer=false)
	m_relax = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
	JuMP.set_optimizer_attribute(m_relax,"OutputFlag",0);
	@variable(m_relax, p[g=1:nb_gen, t=0:T_max+1], lower_bound = 0);		# power output of each generator g at time period t
	@variable(m_relax, pbar[g=1:nb_gen, t=0:T_max+1], lower_bound = 0);		# maximum power available of each generator g at time period t
	@variable(m_relax, u[g=1:nb_gen, t=0:T_max+1], lower_bound = 0, upper_bound = 1);		# commitment status of geenrator g at time period t
	@variable(m_relax, v[g=1:nb_gen, t=0:T_max+1], lower_bound = 0, upper_bound = 1); 		# startup status of generator g at time period t
	@variable(m_relax, w[g=1:nb_gen, t=0:T_max+1], lower_bound = 0, upper_bound = 1);		# shutdown status of generator g at time period t

	@variable(m_relax, 0 <= l[t=1:T_max]);
	@constraint(m_relax, [t=1:T_max], l[t] <= data_demand[t]);			# upper bound for elasticity of the Demand

	### The following constriants described the feasible set for technical constraints, the relaxation of 3-bin space
	@constraint(m_relax, logical_constraint[g=1:nb_gen, t=1:T_max+1], u[g,t] - u[g,t-1] == v[g,t] - w[g,t]);				# Constraint 2, logical constraint on the different status (commitment status, startup status, shutdown status)
	@constraint(m_relax, minimum_up_time[g=1:nb_gen, t=UT[g]:T_max+1], sum(v[g,i]  for i=t-UT[g]+1:t) <= u[g,t]);			# Constraint 3, minimum up time constraints
	@constraint(m_relax, minimum_down_time[g=1:nb_gen, t=DT[g]:T_max+1], sum(w[g,i]  for i=t-DT[g]+1:t) <= (1 - u[g,t])); 	# Constraint 4, minimum down time constraints

	@constraint(m_relax, generation_limits_1[g=1:nb_gen, t=0:T_max+1], MinRunCapacity[g]*u[g,t] <= p[g,t]);					# Constraint 5_1, generation limits constraints
	@constraint(m_relax, generation_limits_2[g=1:nb_gen, t=0:T_max+1], p[g,t] - pbar[g,t] <= 0);							# Constraint 5_2, generation limits constraints
	@constraint(m_relax, generation_limits_3[g=1:nb_gen, t=0:T_max+1], pbar[g,t] <= MaxRunCapacity[g]*u[g,t]);				# Constraint 5_3, generation limits constraints
	
	@constraint(m_relax, ramp_up_constraint[g=1:nb_gen, t=1:T_max+1], pbar[g,t] - p[g,t-1] <= RU[g]*u[g,t-1] + SU[g]*v[g,t]);	# Constraint 6, ramp-up constraints and start-up mode
	@constraint(m_relax, ramp_down_consraint[g=1:nb_gen, t=1:T_max+1], pbar[g,t-1] - p[g,t] <= RD[g]*u[g,t] + SD[g]*w[g,t]);	# Constraint 7, ramp-down constraints and shutdown mode
	### Production >= Demand
	@constraint(m_relax, loads[t=1:T_max], sum(p[g,t] for g=1:nb_gen) == l[t]);
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
	@objective(m_relax, Min, sum( sum( NoLoadConsumption[g]*C[g]*u[g,t] + F[g]*v[g,t] + C[g]*p[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum( l[t] for t=1:T_max));
	optimize!(m_relax);
	# Solving_time = MOI.get(m_relax, MOI.SolveTime());
	Solving_time = 0;
	price_UB = dual.(loads);
	u_val = value.(u); v_val = value.(v); p_val = value.(p);

	if Printer
		nb_con = sum([num_constraints(m_relax, i, j) for (i,j) in list_of_constraint_types(m_relax)]);
		nb_var = num_variables(m_relax);
		println("[LP] nb_con : $(nb_con) and nb_var : $(nb_var)");
	end

	return (value.(p).data, value.(pbar).data, value.(u).data, value.(v).data, value.(w).data, price_UB, value.(l), objective_value(m_relax), Solving_time);
end

function LP_Relaxation_NO_Posterior(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior)
	m_relax = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
	JuMP.set_optimizer_attribute(m_relax,"OutputFlag",0);
	@variable(m_relax, p[g=1:nb_gen, t=0:T_max], lower_bound = 0);					# power output of each generator g at time period t
	@variable(m_relax, pbar[g=1:nb_gen, t=0:T_max], lower_bound = 0);				# maximum power available of each generator g at time period t
	@variable(m_relax, u[g=1:nb_gen, t=0:T_max], lower_bound = 0, upper_bound = 1);# commitment status of geenrator g at time period t
	@variable(m_relax, v[g=1:nb_gen, t=0:T_max], lower_bound = 0, upper_bound = 1);# startup status of generator g at time period t
	@variable(m_relax, w[g=1:nb_gen, t=0:T_max], lower_bound = 0, upper_bound = 1);# shutdown status of generator g at time period t

	@variable(m_relax, l[t=1:T_max], lower_bound = 0);
	@constraint(m_relax, [t=1:T_max], l[t] <= data_demand[t]);			# upper bound for elasticity of the demand

	### The following constriants described the feasible set for technical constraints, the relaxation of 3-bin space
	@constraint(m_relax, logical_constraint[g=1:nb_gen, t=1:T_max], u[g,t] - u[g,t-1] == v[g,t] - w[g,t]);					# Constraint 2, logical constraint on the different status (commitment status, startup status, shutdown status)
	@constraint(m_relax, minimum_up_time[g=1:nb_gen, t=UT[g]:T_max], sum(v[g,i]  for i=t-UT[g]+1:t) <= u[g,t]);				# Constraint 3, minimum up time constraints
	@constraint(m_relax, minimum_down_time[g=1:nb_gen, t=DT[g]:T_max], sum(w[g,i]  for i=t-DT[g]+1:t) <= (1 - u[g,t])); 	# Constraint 4, minimum down time constraints

	@constraint(m_relax, generation_limits_1[g=1:nb_gen, t=0:T_max], MinRunCapacity[g]*u[g,t] <= p[g,t]);					# Constraint 5_1, generation limits constraints
	@constraint(m_relax, generation_limits_2[g=1:nb_gen, t=0:T_max], p[g,t] - pbar[g,t] <= 0);								# Constraint 5_2, generation limits constraints
	@constraint(m_relax, generation_limits_3[g=1:nb_gen, t=0:T_max], pbar[g,t] <= MaxRunCapacity[g]*u[g,t]);				# Constraint 5_3, generation limits constraints
	
	@constraint(m_relax, ramp_up_constraint[g=1:nb_gen, t=1:T_max], pbar[g,t] - p[g,t-1] <= RU[g]*u[g,t-1] + SU[g]*v[g,t]);	# Constraint 6, ramp-up constraints and start-up mode
	@constraint(m_relax, ramp_down_consraint[g=1:nb_gen, t=1:T_max], pbar[g,t-1] - p[g,t] <= RD[g]*u[g,t] + SD[g]*w[g,t]);	# Constraint 7, ramp-down constraints and shutdown mode
	### Production >= Demand
	@constraint(m_relax, loads[t=1:T_max], sum(p[g,t] for g=1:nb_gen) == l[t]);
	### Prior data
	if length(u_prior) > 0
		# @constraint(m_relax, u_prior_constraint[g=1:nb_gen], u[g,0]==u_prior[g]);@constraint(m_relax, v_prior_constraint[g=1:nb_gen], v[g,0]==v_prior[g]);@constraint(m_relax, w_prior_constraint[g=1:nb_gen], w[g,0]==w_prior[g]);@constraint(m_relax, p_prior_constraint[g=1:nb_gen], p[g,0]==p_prior[g]);@constraint(m_relax, pbar_prior_constraint[g=1:nb_gen], pbar[g,0]==pbar_prior[g]);
		for g=1:nb_gen
			JuMP.fix(u[g,0], u_prior[g]; force = true);
			JuMP.fix(v[g,0], v_prior[g]; force = true);
			JuMP.fix(w[g,0], w_prior[g]; force = true);
			JuMP.fix(p[g,0], p_prior[g]; force = true);
			JuMP.fix(pbar[g,0], pbar_prior[g]; force = true);
		end
	end
	### Objective function
	@objective(m_relax, Min, sum( sum( NoLoadConsumption[g]*C[g]*u[g,t] + F[g]*v[g,t] + C[g]*p[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum( l[t] for t=1:T_max));
	optimize!(m_relax);
	# Solving_time = MOI.get(m_relax, MOI.SolveTime());
	Solving_time = 0;
	price_UB = dual.(loads);
	u_val = value.(u); v_val = value.(v); p_val = value.(p);
	return (value.(p).data, value.(pbar).data, value.(u).data, value.(v).data, value.(w).data, price_UB, value.(l), objective_value(m_relax), Solving_time);
end