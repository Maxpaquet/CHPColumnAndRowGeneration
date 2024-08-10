#######################################################################################################
#######################################################################################################
############################################ DW METHODS ###############################################
#######################################################################################################
#######################################################################################################

function isInFunction(p,p_dict)
	T_max = length(p);
	for p2 in values(p_dict)
		println("p2 : $(p2)")
		for t=1:T_max
			if p[t]!=p2[t]
				return false
			end
		end
	end
	return true
end


function Restricted_Master_Program_Column_Generation(data_demand, counter_schedules, p_schedules, cost_schedules, nb_gen, T_max, VOLL)
	### Build Master Restricted Problem
	# m_restricted = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
	m_restricted = JuMP.direct_model(Gurobi.Optimizer(GUROBI_ENV)); # OutputFlag=0
	JuMP.set_optimizer_attribute(m_restricted,"OutputFlag",0);
	JuMP.set_optimizer_attribute(m_restricted,"FeasibilityTol",1e-3);

	@variable(m_restricted, z[g=1:nb_gen, i=1:counter_schedules[g] ], lower_bound = 0, upper_bound=1);
	@variable(m_restricted, l[t=1:T_max],lower_bound = 0);
	@constraint(m_restricted, [t=1:T_max], l[t]<=data_demand[t]);
	# Market clearing constraint
	@constraint(m_restricted, balance_constraint[t=1:T_max], sum( sum( sum( z[g,i]*prod for (time,prod) in enumerate(p_schedules[g,i]) if time==t)  for i=1:counter_schedules[g]) for g=1:nb_gen) - l[t] == 0 );
	@constraint(m_restricted, convex_combination[g=1:nb_gen], sum( z[g,i] for i=1:counter_schedules[g]) == 1 );
	@objective(m_restricted, Min, sum( sum( z[g,i]*cost_schedules[g,i] for i=1:counter_schedules[g]) for g=1:nb_gen) - VOLL*sum( l[t] for t=1:T_max) );
	optimize!(m_restricted);
	# return objective_value(m_restricted), dual.(balance_constraint), dual.(convex_combination), MOI.get(m_restricted, MOI.SolveTime());
	return objective_value(m_restricted), dual.(balance_constraint), dual.(convex_combination), 0;
end


function Column_Generation_DW_1(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer=false)
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
		# sub_prob = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
		sub_prob = JuMP.direct_model(Gurobi.Optimizer(GUROBI_ENV)); # OutputFlag=0
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
	

	### Initialisation
	(p_matching, pbar_matching, u_matching, v_matching, w_matching, l_matching, obj_matching, solving_time_matching) = Matching(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
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
	end
	if Printer
		println("[Column_Generation_DW_1] Unit Commitment problem provides the initial schedules considered.");
		println("[Column_Generation_DW_1] p_schedules : $(p_schedules)");
		println("[Column_Generation_DW_1] u_schedules : $(u_schedules)");
		println("[Column_Generation_DW_1] v_schedules : $(v_schedules)");
		println("[Column_Generation_DW_1] w_schedules : $(w_schedules)");
		println("[Column_Generation_DW_1] cost_schedules : $(cost_schedules)");
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
			println("[Column_Generation_DW_1] iter : $(iter)");
			println("[Column_Generation_DW_1] obj_master : $(obj_master)");
			println("[Column_Generation_DW_1] price : $(price)");
		end
		stopping_criteria = 1;
		for g in 1:nb_gen
			sub_problem = tab_subProblems[g];
			@objective(sub_problem.model, Min, sum( NoLoadConsumption[g]*C[g]*sub_problem.u[t] + F[g]*sub_problem.v[t] + C[g]*sub_problem.p[t] for t=1:T_max) - sum( price[t]*sub_problem.p[t] for t=1:T_max));
			optimize!(sub_problem.model);
			# Solving_time_slaves += MOI.get(sub_problem.model, MOI.SolveTime());
			Solving_time_slaves += 0;
			if Printer
				println("[Column_Generation_DW_1] g : $(g) -> obj = $(objective_value(sub_problem.model))");
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
			println("[Column_Generation_DW_1] p_schedules : $(p_schedules)");
			println("[Column_Generation_DW_1] u_schedules : $(u_schedules)");
			println("[Column_Generation_DW_1] v_schedules : $(v_schedules)");
			println("[Column_Generation_DW_1] w_schedules : $(w_schedules)");
			println("[Column_Generation_DW_1] cost_schedules : $(cost_schedules)");
		end
		if stopping_criteria == 1 # no column were added.
			println("[Column_Generation_DW_1] No columns were added, STOP criterion met.");
			println("[Column_Generation_DW_1] iter_stop : $(iter_stop)");
			break;
		end
	end
	return (obj_master, price, iter_stop, Solving_time_master+Solving_time_slaves, Solving_time_master, Solving_time_slaves, obj_vec);
end


function Column_Generation_DW(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, delta_criterion=-10^(-5), Printer=false)
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
	
	### Build Feasible Schedules
	p_schedules = Dict();
	u_schedules = Dict();
	v_schedules = Dict();
	w_schedules = Dict();
	cost_schedules = Dict();
	counter_schedules = Int64[0 for g=1:nb_gen];

	### Initialisation
	(p_matching, pbar_matching, u_matching, v_matching, w_matching, l_matching, obj_matching, solving_time_matching) = Matching(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
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
	end
	if Printer
		println("Initial schedule considered.");
		println(p_schedules);
		println(u_schedules);
		println(v_schedules);
		println(w_schedules);
		println(cost_schedules);
	end

	### Build Restricted Master Program
	m_restricted = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
	JuMP.set_optimizer_attribute(m_restricted,"OutputFlag",0);
	dict_z = Dict();
	z = @variable(m_restricted, [g=1:nb_gen, i=1:counter_schedules[g] ], lower_bound = 0, upper_bound=1, base_name="z");
	for g=1:nb_gen
		for i=1:counter_schedules[g]
			dict_z[g,i] = z[g,i];
		end
	end
	l = @variable(m_restricted, l[t=1:T_max],lower_bound = 0);
	@constraint(m_restricted, [t=1:T_max], l[t]<=data_demand[t]);
	# Market clearing constraint
	balance_constraint = @constraint(m_restricted, [t=1:T_max], sum( sum( sum( dict_z[g,i]*prod for (time,prod) in enumerate(p_schedules[g,i]) if time==t)  for i=1:counter_schedules[g]) for g=1:nb_gen) - l[t] == 0 );
	convex_combination = @constraint(m_restricted, [g=1:nb_gen], sum( dict_z[g,i] for i=1:counter_schedules[g]) == 1 );
	# @objective(m_restricted, Min, sum( sum( dict_z[g,i]*cost_schedules[g,i] for i=1:counter_schedules[g]) for g=1:nb_gen) - VOLL*sum( l[t] for t=1:T_max) );
	RMP = RestrictedMasterProgram(m_restricted, dict_z, l, balance_constraint, convex_combination);

	### Iterations
	iter_max = 500; # 500
	#eps = -10^(-5);
	#eps = -0.001;
	eps = delta_criterion;
	stopping_criteria = 0;
	obj_master = 0;
	iter_stop = 0;
	Solving_time_master = 0;
	Solving_time_slaves = 0;
	obj_vect = [];
	prices_iterates = [];

	nb_var_vec = [];
	nb_cons_vec = [];
	Solving_time_master_vec = [];
	Solving_time_slave_vec = [];

	for iter=1:iter_max
		iter_stop = iter;
		#obj_master, price, pi_dual, sol_time_master = Restricted_Master_Program_Column_Generation(data_demand, counter_schedules, p_schedules, cost_schedules, nb_gen, T_max, VOLL);
		@objective(RMP.model, Min, sum( sum( dict_z[g,i]*cost_schedules[g,i] for i=1:counter_schedules[g]) for g=1:nb_gen) - VOLL*sum( l[t] for t=1:T_max) );
		optimize!(RMP.model);
		obj_master = objective_value(RMP.model);
		push!(obj_vect, obj_master);
		price = dual.(RMP.balance_constraint);
		push!(prices_iterates, price);
		pi_dual = dual.(RMP.convex_combination);
		# Solving_time_master += MOI.get(RMP.model, MOI.SolveTime());
		Solving_time_master += 0;

		push!(nb_var_vec, num_variables(RMP.model));
		push!(nb_cons_vec, sum([num_constraints(RMP.model, i, j) for (i,j) in list_of_constraint_types(RMP.model)]));
		# push!(Solving_time_master_vec, MOI.get(RMP.model, MOI.SolveTime()));
		push!(Solving_time_master_vec, 0);

		if Printer
			println();println();
			println("iter : $(iter)");
			println("obj_master : $(obj_master)");
			println("price : $(price)");
			for i=1:counter_schedules[1]
				println("z[1,$(i)] = $(value(dict_z[1,i]))");
			end
		end
		solving_time_mean_slaves = 0;
		stopping_criteria = 1;
		for g in 1:nb_gen
			sub_problem = tab_subProblems[g];
			@objective(sub_problem.model, Min, sum( NoLoadConsumption[g]*C[g]*sub_problem.u[t] + F[g]*sub_problem.v[t] + C[g]*sub_problem.p[t] for t=1:T_max) - sum( price[t]*sub_problem.p[t] for t=1:T_max));
			optimize!(sub_problem.model);
			if Printer
				println("g : $(g) -> obj = $(objective_value(sub_problem.model))");
			end
			# Solving_time_slaves += MOI.get(sub_problem.model, MOI.SolveTime());
			Solving_time_slaves += 0;
			# solving_time_mean_slaves+= MOI.get(sub_problem.model, MOI.SolveTime());
			solving_time_mean_slaves+= 0;
			#reduced_cost = ComputeReducedCost(objective_value(sub_problem.model), pi_dual[g]);
			reduced_cost = objective_value(sub_problem.model) - pi_dual[g];
			if reduced_cost < eps
				stopping_criteria = 0
				p_sub = value.(sub_problem.p).data;
				u_sub = value.(sub_problem.u).data;
				v_sub = value.(sub_problem.v).data;
				w_sub = value.(sub_problem.w).data;
				#if isInFunction(p_sub[2:T_max+1],p_schedules)==false
				counter_schedules[g]+=1;
				p_schedules[g,counter_schedules[g]] = p_sub[2:T_max+1];
				u_schedules[g,counter_schedules[g]] = u_sub[2:T_max+1];
				v_schedules[g,counter_schedules[g]] = v_sub[2:T_max+1];
				w_schedules[g,counter_schedules[g]] = w_sub[2:T_max+1];
				cost_schedules[g,counter_schedules[g]] = sum( NoLoadConsumption[g]*C[g]*u_sub[t] + F[g]*v_sub[t] + C[g]*p_sub[t] for t=2:T_max+1);
				
				new_var_z = @variable(RMP.model, [g, counter_schedules[g] ], lower_bound = 0, upper_bound=1, base_name="z");
				RMP.z[g, counter_schedules[g]] = new_var_z[g, counter_schedules[g]];
				for t=1:T_max
					#println("Before adding RMP.balance_constraint : $(RMP.balance_constraint[t])");
					set_normalized_coefficient(RMP.balance_constraint[t], RMP.z[g,counter_schedules[g]], p_sub[t+1]);
					#println("After adding RMP.balance_constraint : $(RMP.balance_constraint[t])");println();println();
				end
				set_normalized_coefficient(RMP.convex_combination[g], RMP.z[g,counter_schedules[g]], 1);
			end
		end
		push!(Solving_time_slave_vec, solving_time_mean_slaves);
		if Printer
			println("After adding schedules.")
			println(p_schedules);
			println(u_schedules);
			println(v_schedules);
			println(w_schedules);
			println(cost_schedules);
		end
		if stopping_criteria == 1 # no column were added.
			iter_stop
			break;
		end
	end
	return (obj_master, price, iter_stop, Solving_time_master+Solving_time_slaves, Solving_time_master, Solving_time_slaves, obj_vect, prices_iterates, nb_var_vec ,nb_cons_vec, Solving_time_master_vec, Solving_time_slave_vec)
end