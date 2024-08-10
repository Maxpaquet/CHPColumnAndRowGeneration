###################################################################################################
###################################################################################################
###################################### COLUMN-AND-ROW METHOD ######################################
###################################################################################################
###################################################################################################

"""
Function used in the dynamic column-and-row generation algorithm.
"""
function Pricing_Problem(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, price, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior)
	m_matching = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
	JuMP.set_optimizer_attribute(m_matching,"OutputFlag",0);
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
		#@constraint(m_matching, u_prior_constraint[g=1:nb_gen], u[g,0]==u_prior[g]);@constraint(m_matching, v_prior_constraint[g=1:nb_gen], v[g,0]==v_prior[g]);@constraint(m_matching, w_prior_constraint[g=1:nb_gen], w[g,0]==w_prior[g]);@constraint(m_matching, p_prior_constraint[g=1:nb_gen], p[g,0]==p_prior[g]);@constraint(m_matching, pbar_prior_constraint[g=1:nb_gen], pbar[g,0]==pbar_prior[g]);
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
		#@constraint(m_matching, u_posterior_constraint[g=1:nb_gen], u[g,T_max+1]==u_posterior[g]);@constraint(m_matching, v_posterior_constraint[g=1:nb_gen], v[g,T_max+1]==v_posterior[g]);@constraint(m_matching, w_posterior_constraint[g=1:nb_gen], w[g,T_max+1]==w_posterior[g]);@constraint(m_matching, p_posterior_constraint[g=1:nb_gen], p[g,T_max+1]==p_posterior[g]);@constraint(m_matching, pbar_posterior_constraint[g=1:nb_gen], pbar[g,T_max+1]==pbar_posterior[g]);
		for g=1:nb_gen
			JuMP.fix(u[g,T_max+1],u_posterior[g]; force = true);
			JuMP.fix(v[g,T_max+1],v_posterior[g]; force = true);
			JuMP.fix(w[g,T_max+1],w_posterior[g]; force = true);
			JuMP.fix(p[g,T_max+1],p_posterior[g]; force = true);
			JuMP.fix(pbar[g,T_max+1],pbar_posterior[g]; force = true);
		end
	end
	### Objective function
	@objective(m_matching, Min, sum( sum( NoLoadConsumption[g]*C[g]*u[g,t] + F[g]*v[g,t] + C[g]*p[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum( l[t] for t=1:T_max) - sum( price[t]*(sum( p[g,t] for g=1:nb_gen) - l[t]) for t=1:T_max));
	optimize!(m_matching);
	u_val = value.(u); v_val = value.(v); w_val = value.(w); p_val = value.(p); pbar_val = value.(pbar);l_val = value.(l);
	return (p_val.data, pbar_val.data, u_val.data, v_val.data, w_val.data, l_val, objective_value(m_matching));
end

"""
Function used in the dynamic column-and-row generation.
"""
function Restricted_Extended_Formulation(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, A, B, nb_intervals_gen,u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior)
	# m = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
	m = JuMP.direct_model(Gurobi.Optimizer(GUROBI_ENV)); # OutputFlag=0
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
	### Feasible Dispatch Polytope (equations 13 from source), For interval i of generator g, a = A[g,i] and b=A[g,i]
	@constraint(m, no_production[g=1:nb_gen, i=1:nb_intervals_gen[g], t=0:T_max+1; (t<A[g,i] || t>B[g,i])], p[g,i,t] <= 0);
	@constraint(m, min_output[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]:B[g,i] ], -p[g,i,t] <= -gamma[g,i] * MinRunCapacity[g]);
	@constraint(m, max_output[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]:B[g,i] ], p[g,i,t] - pbar[g,i,t] <= 0);
	@constraint(m, upper_bound_on_power_output_1[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]:B[g,i] ] , pbar[g,i,t] <= gamma[g,i] * MaxRunCapacity[g]);
	@constraint(m, upper_bound_on_power_output_2[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]:B[g,i] ] , pbar[g,i,t] <= gamma[g,i] * (SU[g] + (t - A[g,i])*RU[g] ));
	@constraint(m, upper_bound_on_power_output_3[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]:B[g,i] ] , pbar[g,i,t] <= gamma[g,i] * (SD[g] + (B[g,i] - t)*RD[g] ));
	@constraint(m, limit_power_jumps_up_1[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], pbar[g,i,t] - p[g,i,t-1] <= gamma[g,i] * RU[g]);
	@constraint(m, limit_power_jumps_up_2[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], pbar[g,i,t] - p[g,i,t-1] <= gamma[g,i] * (SD[g] + (B[g,i] - t)*RD[g] - MinRunCapacity[g]));
	@constraint(m, limit_power_jumps_down_1[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], pbar[g,i,t-1] - p[g,i,t] <= gamma[g,i] * RD[g]);
	@constraint(m, limit_power_jumps_down_2[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], pbar[g,i,t-1] - p[g,i,t] <= gamma[g,i] * (SU[g] + (t - A[g,i])*RU[g] - MinRunCapacity[g]));
	### Packing Dispatch Polytopes
	# 15e : to eliminate the impossible combinations (to have the convex hull)
	@constraint(m, minimum_down_time_constraint[g=1:nb_gen, t=1:T_max], sum(gamma[g,i] for i=1:nb_intervals_gen[g] if (A[g,i] <= t && t <= B[g,i]+DT[g])) <= 1);
	
	@constraint(m, production_time[g=1:nb_gen, t=1:T_max], p_time[g,t] == sum(p[g,i,t]  for i=1:nb_intervals_gen[g]));
	@constraint(m, productionbar_time[g=1:nb_gen, t=1:T_max], pbar_time[g,t] == sum(pbar[g,i,t]  for i=1:nb_intervals_gen[g]));
	@constraint(m, commitment_status[g=1:nb_gen,t=0:T_max+1], sum(gamma[g,i] for i=1:nb_intervals_gen[g] if (A[g,i] <= t && t <= B[g,i])) == u[g,t]);
	@constraint(m, startup_status[g=1:nb_gen, t=0:T_max+1], sum(gamma[g,i]  for i=1:nb_intervals_gen[g] if (t == A[g,i])) == v[g,t]);
	@constraint(m, shutdown_status[g=1:nb_gen, t=0:T_max+1], sum(gamma[g,i] for i=1:nb_intervals_gen[g] if (t == B[g,i]+1)) == w[g,t]);
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
	@objective(m, Min, sum( sum( NoLoadConsumption[g]*C[g]*u[g,t] + F[g]*v[g,t] + C[g]*p_time[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum(l[t] for t=1:T_max));
	# @objective(m, Min, sum( sum( NoLoadConsumption[g]*C[g]*sum(gamma[g,i] for i=1:nb_intervals_gen[g] if (A[g,i] <= t && t <= B[g,i])) + F[g]*sum(gamma[g,i]  for i=1:nb_intervals_gen[g] if (t == A[g,i])) + C[g]*sum(p[g,i,t]  for i=1:nb_intervals_gen[g]) for t=1:T_max) for g=1:nb_gen) - VOLL*sum(l[t] for t=1:T_max));
	optimize!(m);
	price = dual.(loads);
	# return (value.(p).data, value.(pbar).data, value.(gamma).data, value.(p_time).data, value.(pbar_time).data, value.(u).data, value.(v).data, value.(w).data, price, value.(l), objective_value(m), MOI.get(m, MOI.SolveTime()));
	return (value.(p).data, value.(pbar).data, value.(gamma).data, value.(p_time).data, value.(pbar_time).data, value.(u).data, value.(v).data, value.(w).data, price, value.(l), objective_value(m), 0);

end

"""
Need to deal case when nb_gen = 1.
"""
function Compute_A_B(u, nb_gen)
	nb_intervals_gen = zeros(Int64, nb_gen);
	if nb_gen>1
		for g=1:nb_gen
			flag = false;
			for elem in u[g,:]
				if elem>0
					if flag==false
						flag = true;
						nb_intervals_gen[g]+=1;
					end
				elseif flag==true
					flag = false;
				end
			end
		end
	else
		flag = false;
		for elem in u
			if elem>0
				if flag==false
					flag = true;
					nb_intervals_gen[1]+=1;
				end
			elseif flag==true
				flag = false;
			end
		end
	end
	A = (-1)*ones(Int64, nb_gen, maximum(nb_intervals_gen));
	B = (-1)*ones(Int64, nb_gen, maximum(nb_intervals_gen));
	a = 0;
	b = 0;
	if nb_gen>1
		for g=1:nb_gen
			flag = false;
			count = 0
			index = 1;
			for elem in u[g,:]
				count+=1;
				#println("elem : ",elem)
				if elem>0
					if flag==false
						flag = true;
						#A[g,index] = count;
						a = count;
					end
				elseif flag==true
					B[g,index] = count-1;
					A[g,index] = a;
					index+=1;
					flag = false;
				end
			end
		end
	else
		flag = false;
		count = 0
		index = 1;
		for elem in u
			count+=1;
			#println("[Compute_A_B] elem : ",elem)
			if elem>0
				if flag==false
					flag = true;
					#A[g,index] = count;
					a = count;
				end
			elseif flag==true
				B[1,index] = count-1;
				A[1,index] = a;
				index+=1;
				flag = false;
			end
		end
	end
	return A,B,nb_intervals_gen
end


function A_B_In_Intervals(a,b,g,A,B,nb_intervals_gen)
	for i=1:nb_intervals_gen[g]
		if a==A[g,i] && b==B[g,i]
			return true
		end
	end
	return false;
end


function Update_A_B(u, nb_gen, A, B, nb_intervals_gen)
	a = 0;
	b = 0;
	nb_intervals_added = zeros(Int64, nb_gen);
	nb_intervals_gen_old = copy(nb_intervals_gen);
	if nb_gen>1
		for g=1:nb_gen
			flag = false;
			count = 0;
			for elem in u[g,:]
				count+=1;
				if elem>0
					if flag==false
						flag = true;
						# nb_intervals_gen[g]+=1;
						a = count;
					end
				elseif flag==true
					b = count-1;
					if A_B_In_Intervals(a,b,g, A, B, nb_intervals_gen_old)==false
						nb_intervals_gen[g]+=1;
						nb_intervals_added[g]+=1;
					end
					flag = false;
				end
			end
		end
	else
		flag = false;
		count = 0;
		for elem in u
			count+=1;
			if elem>0
				if flag==false
					flag = true;
					# nb_intervals_gen[g]+=1;
					a = count;
				end
			elseif flag==true
				b = count-1;
				println("[Update_A_B] A = $(A)");
				println("[Update_A_B] B = $(B)");
				if A_B_In_Intervals(a,b,1, A, B, nb_intervals_gen_old)==false
					nb_intervals_gen[1]+=1;
					nb_intervals_added[1]+=1;
				end
				flag = false;
			end
		end
	end
	a = 0;
	b = 0;
	A_new = (-1)*ones(Int64, nb_gen, maximum(nb_intervals_gen));
	B_new = (-1)*ones(Int64, nb_gen, maximum(nb_intervals_gen));
	A_added = (-1)*ones(Int64, nb_gen, maximum(nb_intervals_added));
	B_added = (-1)*ones(Int64, nb_gen, maximum(nb_intervals_added));
	(i,j) = size(A);
	A_new[1:i,1:j] = A;
	(i,j) = size(B);
	B_new[1:i,1:j] = B;
	if nb_gen>1
		for g=1:nb_gen
			flag = false;
			index = nb_intervals_gen_old[g]+1;
			index_added = 1;
			count = 0;
			for elem in u[g,:]
				count+=1;
				if elem>0
					if flag==false
						flag = true;
						a = count;
					end
				elseif flag==true
					b = count-1;
					if A_B_In_Intervals(a,b,g, A, B, nb_intervals_gen_old)==false
						A_new[g,index] = a;
						B_new[g,index] = b;
						A_added[g,index_added] = a;
						B_added[g,index_added] = b;
						index+=1;index_added+=1;
					end
					flag = false;
				end
			end
		end
	else
		#println("One generator");
		flag = false;
		index = nb_intervals_gen_old[1]+1;
		index_added = 1;
		count = 0;
		for elem in u
			count+=1;
			if elem>0
				if flag==false
					flag = true;
					a = count;
				end
			elseif flag==true
				b = count-1;
				if A_B_In_Intervals(a,b,1, A, B, nb_intervals_gen_old)==false
					A_new[1,index] = a;
					B_new[1,index] = b;
					A_added[1,index_added] = a;
					B_added[1,index_added] = b;
					index+=1;index_added+=1;
				end
				flag = false;
			end
		end
	end
	return A_new,B_new,nb_intervals_gen, A_added, B_added, nb_intervals_added
end


"""
Dynamic column-and-row generation.
No hot-start model implementation.
"""
function Column_And_Row_Generation_1(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer=false)
	# if Printer
	# 	println();
	# 	println();
	# 	println("##############################################################################################################");
	# 	println("Pre-process");
	# 	println("nb_gen : ",nb_gen);
	# end
	### Solve Matching and get "on-intervals" to initialize \bar{S}
	#(p_matching, pbar_matching, u_matching, v_matching, w_matching, l_matching, obj_real_matching, Solving_time) = Matching(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
	
	# m_matching = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
	m_matching = JuMP.direct_model(Gurobi.Optimizer(GUROBI_ENV)); # OutputFlag=0
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
	loads = @constraint(m_matching, [t=1:T_max], sum( p[g,t] for g=1:nb_gen) == l[t]);
	pricing_problem = PricingProblem(m_matching,p,pbar,u,v,w,l,logical_constraint,minimum_up_time,minimum_down_time,generation_limits_1,generation_limits_2,generation_limits_3,ramp_up_constraint,ramp_down_consraint);
	# if Printer
	# 	println("Unit Commitment Problem. It provides the first production intervals for initializing the restricted set of intervals.");
	# 	println(m_matching);
	# end
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
			JuMP.fix(pricing_problem.u[g,T_max+1], u_posterior[g]; force = true);
			JuMP.fix(pricing_problem.v[g,T_max+1], v_posterior[g]; force = true);
			JuMP.fix(pricing_problem.w[g,T_max+1], w_posterior[g]; force = true);
			JuMP.fix(pricing_problem.p[g,T_max+1], p_posterior[g]; force = true);
			JuMP.fix(pricing_problem.pbar[g,T_max+1], pbar_posterior[g]; force = true);
		end
	end
	### Objective function
	@objective(pricing_problem.model, Min, sum( sum( NoLoadConsumption[g]*C[g]*pricing_problem.u[g,t] + F[g]*pricing_problem.v[g,t] + C[g]*pricing_problem.p[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum( pricing_problem.l[t] for t=1:T_max) );
	optimize!(pricing_problem.model);
	u_matching = value.(pricing_problem.u).data;
	delete(pricing_problem.model, loads); # Only need the constraint for initialization. Don't need after because it will be dualize.
	if nb_gen>1
		u_matching = u_matching[:,2:T_max+1];
	else
		u_matching = u_matching[2:T_max+1];
	end	
	(A,B,nb_intervals_gen) = Compute_A_B(u_matching, nb_gen);
	if Printer
		println("[Column_And_Row_Generation_1] u_matching = $(u_matching)");
		println("[Column_And_Row_Generation_1] Initial intervals, A : ",A);
		println("[Column_And_Row_Generation_1] Initial intervals, B : ",B);
		println("[Column_And_Row_Generation_1] nb_intervals_gen : ",nb_intervals_gen);
	end
	beta = 0;
	iter_max = 200; # 5000
	obj_restricted_vector = [];
	obj_pricing_vector = [];
	p_time_restricted = 0;
	pbar_time_restricted = 0;
	u_restricted = 0;
	v_restricted = 0;
	w_restricted = 0;
	price = 0;
	price_opt = 0;
	price_iterates = [];
	l_restricted = 0;
	obj_restricted = 0;
	iter_stop = 0;
	obj_pricing = 0;
	# phi = 1e-3;
	# phi = 10^(-3);
	phi = 0;
	Solving_time = 0;
	for iter=1:iter_max
		if Printer
			println();
			println();
			println("[Column_And_Row_Generation_1] iteration ",iter);
		end
		iter_stop = iter;
		### Solve the restricted problem
		(p_restricted, pbar_restricted, gamma_restricted, p_time_restricted, pbar_time_restricted, u_restricted, v_restricted, w_restricted, price, l_restricted, obj_restricted, time) = Restricted_Extended_Formulation(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, A, B, nb_intervals_gen,u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior);
		if Printer
			println("[Column_And_Row_Generation_1] p_time_restricted = $(p_time_restricted)");
			println("[Column_And_Row_Generation_1] u_restricted = $(u_restricted)");
		end
		if Printer
			println("[Column_And_Row_Generation_1] price = $(price)");
			println("[Column_And_Row_Generation_1] obj_restricted : ",obj_restricted);
			println("[Column_And_Row_Generation_1] obj_pricing    : ",obj_pricing);
		end
		push!(obj_restricted_vector, obj_restricted);
		push!(price_iterates, price);
		price_opt = price;
		Solving_time += time;
		### Solve the pricing problem
		@objective(pricing_problem.model, Min, sum( sum( NoLoadConsumption[g]*C[g]*pricing_problem.u[g,t] + F[g]*pricing_problem.v[g,t] + C[g]*pricing_problem.p[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum( pricing_problem.l[t] for t=1:T_max) - sum( price[t]*(sum( pricing_problem.p[g,t] for g=1:nb_gen) - pricing_problem.l[t]) for t=1:T_max) );
		optimize!(pricing_problem.model);
		# Solving_time += MOI.get(pricing_problem.model, MOI.SolveTime());
		Solving_time += 0;
		obj_pricing = objective_value(pricing_problem.model);
		push!(obj_pricing_vector, obj_pricing);
		u_pricing = value.(pricing_problem.u).data;
		### Compute the Lagrangian dual bound
		if iter==1
			beta = obj_pricing;
		else
			beta = maximum([obj_pricing beta]);
		end
		if Printer
			println("[Column_And_Row_Generation_1] Update dual bound beta = $(beta)");
		end
		if  obj_restricted <= beta + phi
			if Printer
				println();
				println();
				if Printer
					println("[Column_And_Row_Generation_1] OPTIMAL SOLUTION FOUND");
				end
			end
			break; # STOP algorithm.
		end
		### Update current bundle
		if nb_gen>1
			u = u_pricing[:,2:T_max+1];
		else
			u = u_pricing[2:T_max+1];
		end
		A_new,B_new,nb_intervals_gen_new, A_added, B_added, nb_intervals_added = Update_A_B(u, nb_gen, A, B, nb_intervals_gen);
		A = copy(A_new);
		B = copy(B_new);
		nb_intervals_gen = copy(nb_intervals_gen_new);
		if Printer
			println("[Column_And_Row_Generation_1] u : $(u)");
			println("[Column_And_Row_Generation_1] Update starting of intervals A : ",A);
			println("[Column_And_Row_Generation_1] Update ending of intervals   B : ",B);
			println("[Column_And_Row_Generation_1] nb_intervals_gen : ",nb_intervals_gen);
			println("[Column_And_Row_Generation_1] A_added : ",A_added);
			println("[Column_And_Row_Generation_1] B_added : ",B_added);
			println("[Column_And_Row_Generation_1] nb_intervals_added : ",nb_intervals_added);
		end
	end
	if Printer
		println("[Column_And_Row_Generation] Stop at iteration ",iter_stop);
	end
	return (p_time_restricted, pbar_time_restricted, u_restricted, v_restricted, w_restricted, price_opt, l_restricted, obj_restricted, obj_pricing, Solving_time, obj_restricted_vector, obj_pricing_vector, price_iterates);
end



"""
Dynamic column-and-row generation. Implementation of hot-start model.
"""
function Column_And_Row_Generation(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer, phi_criteria=10^(-5))
	# if Printer
	# 	println();
	# 	println("##############################################################################################################");
	# 	println("Pre-process");
	# 	println("nb_gen : ",nb_gen);
	# end
	### Solve Matching and get "on-intervals" to initialize \bar{S}
	m_matching = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
	JuMP.set_optimizer_attribute(m_matching,"OutputFlag",0);
	@variable(m_matching, 0 <= p[g=1:nb_gen, t=0:T_max+1]);		# power output of each generator g at time period t
	@variable(m_matching, 0 <= pbar[g=1:nb_gen, t=0:T_max+1]);	# maximum power available of each generator g at time period t
	@variable(m_matching, u[g=1:nb_gen, t=0:T_max+1], Bin);		# commitment status of geenrator g at time period t
	@variable(m_matching, v[g=1:nb_gen, t=0:T_max+1], Bin); 	# startup status of generator g at time period t
	@variable(m_matching, w[g=1:nb_gen, t=0:T_max+1], Bin);		# shutdown status of generator g at time period t
	@variable(m_matching, 0 <= l[t=1:T_max]);
	@constraint(m_matching, [t=1:T_max], l[t] <= data_demand[t]);			# upper bound for elasticity of the demand
	# The following constraints described the feasible set for technical constraints, the relaxation of 3-bin space
	logical_constraint 	= @constraint(m_matching, [g=1:nb_gen, t=1:T_max+1], u[g,t] - u[g,t-1] == v[g,t] - w[g,t]);					# Constraint 2, logical constraint on the different status (commitment status, startup status, shutdown status)
	minimum_up_time 	= @constraint(m_matching, [g=1:nb_gen, t=UT[g]:T_max+1], sum(v[g,i]  for i=t-UT[g]+1:t) <= u[g,t]);			# Constraint 3, minimum up time constraints
	minimum_down_time 	= @constraint(m_matching, [g=1:nb_gen, t=DT[g]:T_max+1], sum(w[g,i]  for i=t-DT[g]+1:t) <= (1 - u[g,t])); 	# Constraint 4, minimum down time constraints
	generation_limits_1	= @constraint(m_matching, [g=1:nb_gen, t=0:T_max+1], MinRunCapacity[g]*u[g,t] <= p[g,t]);					# Constraint 5_1, generation limits constraints
	generation_limits_2 = @constraint(m_matching, [g=1:nb_gen, t=0:T_max+1], p[g,t] <= pbar[g,t]);									# Constraint 5_2, generation limits constraints
	generation_limits_3 = @constraint(m_matching, [g=1:nb_gen, t=0:T_max+1], pbar[g,t] <= MaxRunCapacity[g]*u[g,t]);				# Constraint 5_3, generation limits constraints
	ramp_up_constraint 	= @constraint(m_matching, [g=1:nb_gen, t=1:T_max+1], pbar[g,t] - p[g,t-1] <= RU[g]*u[g,t-1] + SU[g]*v[g,t]);# Constraint 6, ramp-up constraints and start-up mode
	ramp_down_consraint = @constraint(m_matching, [g=1:nb_gen, t=1:T_max+1], pbar[g,t-1] - p[g,t] <= RD[g]*u[g,t] + SD[g]*w[g,t]);	# Constraint 7, ramp-down constraints and shutdown mode
	loads = @constraint(m_matching, [t=1:T_max], sum( p[g,t] for g=1:nb_gen) == l[t]);
	pricing_problem = PricingProblem(m_matching,p,pbar,u,v,w,l,logical_constraint,minimum_up_time,minimum_down_time,generation_limits_1,generation_limits_2,generation_limits_3,ramp_up_constraint,ramp_down_consraint);
	# Prior data
	if length(u_prior) > 0
		for g=1:nb_gen
			JuMP.fix(pricing_problem.u[g,0], u_prior[g]; force = true);
			JuMP.fix(pricing_problem.v[g,0], v_prior[g]; force = true);
			JuMP.fix(pricing_problem.w[g,0], w_prior[g]; force = true);
			JuMP.fix(pricing_problem.p[g,0], p_prior[g]; force = true);
			JuMP.fix(pricing_problem.pbar[g,0], pbar_prior[g]; force = true);
		end
	end
	# Posterior data
	if length(u_posterior) > 0
		for g=1:nb_gen
			JuMP.fix(pricing_problem.u[g,T_max+1],u_posterior[g]; force = true);
			JuMP.fix(pricing_problem.v[g,T_max+1],v_posterior[g]; force = true);
			JuMP.fix(pricing_problem.w[g,T_max+1],w_posterior[g]; force = true);
			JuMP.fix(pricing_problem.p[g,T_max+1],p_posterior[g]; force = true);
			JuMP.fix(pricing_problem.pbar[g,T_max+1],pbar_posterior[g]; force = true);
		end
	end
	# Objective function
	@objective(pricing_problem.model, Min, sum( sum( NoLoadConsumption[g]*C[g]*pricing_problem.u[g,t] + F[g]*pricing_problem.v[g,t] + C[g]*pricing_problem.p[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum( pricing_problem.l[t] for t=1:T_max) );
	optimize!(pricing_problem.model);
	u_matching = value.(pricing_problem.u).data;
	delete(pricing_problem.model, loads);
	if nb_gen>1
		u_matching = u_matching[:,2:T_max+1];
	else
		u_matching = u_matching[2:T_max+1];
	end	
	(A,B,nb_intervals_gen) = Compute_A_B(u_matching, nb_gen);
	A_added = copy(A); B_added = copy(B); nb_intervals_added = copy(nb_intervals_gen);
	if Printer
		println("[Column_And_Row_Generation_1] u_matching : ",u_matching);
		println("[Column_And_Row_Generation] Initialization - A_added : $(A_added)");
		println("[Column_And_Row_Generation] Initialization - B_added : $(B_added)");
		println("[Column_And_Row_Generation] Initialization - nb_intervals_added : $(nb_intervals_added)");
	end
	###############################################################################################
	### Restricted model
	restricted_model = JuMP.direct_model(Gurobi.Optimizer());
	JuMP.set_optimizer_attribute(restricted_model,"OutputFlag",0);
	# Variables
	@variable(restricted_model, 0 <= p[g=1:nb_gen, i=1:nb_intervals_added[g], t=0:T_max+1], base_name="p");	# Power output vector of the generators		p[generator g, interval(a,b) i, time t]
	@variable(restricted_model, 0 <= pbar[g=1:nb_gen, i=1:nb_intervals_added[g], t=0:T_max+1]);				# Maximum power available from generators	pbar[generator g, interval(a,b) i, time t]
	@variable(restricted_model, 0 <= gamma[g=1:nb_gen, i=1:nb_intervals_added[g]]);							# Variable that indicates whether a vertex (interval) is in the packing or not   gamma[generator g, interval(a,b) i]
	@variable(restricted_model, 0 <= p_time[g=1:nb_gen, t=0:T_max+1]);										# (to simplify computations) production of each generator on each time period
	@variable(restricted_model, 0 <= pbar_time[g=1:nb_gen, t=0:T_max+1]);									# (to simplify computations) maximum power available of each generator on each time period
	@variable(restricted_model, 0 <= u[g=1:nb_gen, t=0:T_max+1]); 											# Commitment status variables
	@variable(restricted_model, 0 <= v[g=1:nb_gen, t=0:T_max+1]); 											# Startup status variables
	@variable(restricted_model, 0 <= w[g=1:nb_gen, t=0:T_max+1]); 											# Shutdown status variables
	@variable(restricted_model, 0 <= l[t=1:T_max]); 														# Modelling an elastic demand
	@constraint(restricted_model, [t=1:T_max], l[t] <= data_demand[t]);
	@constraint(restricted_model, [g=1:nb_gen, i=1:nb_intervals_added[g]], gamma[g,i] <=1); 				# By security
	
	dict_gamma = Dict();
	dict_p = Dict();
	dict_pbar = Dict();
	for g=1:nb_gen
		for i=1:nb_intervals_added[g]
			dict_gamma[g,i] = gamma[g,i];
			for t=0:T_max+1
				dict_p[g,i,t] = p[g,i,t];
				dict_pbar[g,i,t] = pbar[g,i,t];
			end
		end
	end

	count_p = 0;
	count_gamma = 0;
	count_p_time = 0;
	count_u = 0;
	count_l = T_max;
	for g=1:nb_gen
		count_p_time += T_max+2;
		count_u += T_max+2;
		for i=1:nb_intervals_added[g]
			count_p+=T_max+2;
			count_gamma+=1;
		end
	end
	# Feasible Dispatch Polytope (equations 13 from source), For interval i of generator g, a = A[g,i] and b=A[g,i]
	no_production 				  = @constraint(restricted_model, [g=1:nb_gen, i=1:nb_intervals_added[g], t=0:T_max+1; (t<A_added[g,i] || t>B_added[g,i])], dict_p[g,i,t] <= 0);
	min_output 					  = @constraint(restricted_model, [g=1:nb_gen, i=1:nb_intervals_added[g], t=A_added[g,i]:B_added[g,i] ], -dict_p[g,i,t] <= -dict_gamma[g,i] * MinRunCapacity[g]);
	max_output 					  = @constraint(restricted_model, [g=1:nb_gen, i=1:nb_intervals_added[g], t=A_added[g,i]:B_added[g,i] ], dict_p[g,i,t] - dict_pbar[g,i,t] <= 0);
	upper_bound_on_power_output_1 = @constraint(restricted_model, [g=1:nb_gen, i=1:nb_intervals_added[g], t=A_added[g,i]:B_added[g,i] ] , dict_pbar[g,i,t] <= dict_gamma[g,i] * MaxRunCapacity[g]);
	upper_bound_on_power_output_2 = @constraint(restricted_model, [g=1:nb_gen, i=1:nb_intervals_added[g], t=A_added[g,i]:B_added[g,i] ] , dict_pbar[g,i,t] <= dict_gamma[g,i] * (SU[g] + (t - A_added[g,i])*RU[g] ));
	upper_bound_on_power_output_3 = @constraint(restricted_model, [g=1:nb_gen, i=1:nb_intervals_added[g], t=A_added[g,i]:B_added[g,i] ] , dict_pbar[g,i,t] <= dict_gamma[g,i] * (SD[g] + (B_added[g,i] - t)*RD[g] ));
	limit_power_jumps_up_1        = @constraint(restricted_model, [g=1:nb_gen, i=1:nb_intervals_added[g], t=A_added[g,i]+1:B_added[g,i] ], dict_pbar[g,i,t] - dict_p[g,i,t-1] <= dict_gamma[g,i] * RU[g]);
	limit_power_jumps_up_2 		  = @constraint(restricted_model, [g=1:nb_gen, i=1:nb_intervals_added[g], t=A_added[g,i]+1:B_added[g,i] ], dict_pbar[g,i,t] - dict_p[g,i,t-1] <= dict_gamma[g,i] * (SD[g] + (B_added[g,i] - t)*RD[g] - MinRunCapacity[g]));
	limit_power_jumps_down_1 	  = @constraint(restricted_model, [g=1:nb_gen, i=1:nb_intervals_added[g], t=A_added[g,i]+1:B_added[g,i] ], dict_pbar[g,i,t-1] - dict_p[g,i,t] <= dict_gamma[g,i] * RD[g]);
	limit_power_jumps_down_2 	  = @constraint(restricted_model, [g=1:nb_gen, i=1:nb_intervals_added[g], t=A_added[g,i]+1:B_added[g,i] ], dict_pbar[g,i,t-1] - dict_p[g,i,t] <= dict_gamma[g,i] * (SU[g] + (t - A_added[g,i])*RU[g] - MinRunCapacity[g]));
	minimum_down_time_constraint  = @constraint(restricted_model, [g=1:nb_gen, t=1:T_max], sum(dict_gamma[g,i] for i=1:nb_intervals_added[g] if (A_added[g,i] <= t && t <= B_added[g,i]+DT[g])) - 1 <= 0);
	production_time               = @constraint(restricted_model, [g=1:nb_gen, t=1:T_max], p_time[g,t] - sum(dict_p[g,i,t]  for i=1:nb_intervals_added[g]) == 0);
	productionbar_time            = @constraint(restricted_model, [g=1:nb_gen, t=1:T_max], pbar_time[g,t] - sum(dict_pbar[g,i,t]  for i=1:nb_intervals_added[g]) == 0);
	commitment_status             = @constraint(restricted_model, [g=1:nb_gen,t=0:T_max+1], sum(dict_gamma[g,i] for i=1:nb_intervals_added[g] if (A_added[g,i] <= t && t <= B_added[g,i])) - u[g,t] == 0);
	startup_status                = @constraint(restricted_model, [g=1:nb_gen, t=0:T_max+1], sum(dict_gamma[g,i]  for i=1:nb_intervals_added[g] if (t == A_added[g,i])) - v[g,t] == 0);
	shutdown_status               = @constraint(restricted_model, [g=1:nb_gen, t=0:T_max+1], sum(dict_gamma[g,i] for i=1:nb_intervals_added[g] if (t == B_added[g,i]+1)) - w[g,t] == 0);
	loads                         = @constraint(restricted_model, [t=1:T_max], sum( p_time[g,t] for g=1:nb_gen) - l[t] == 0);
	restricted_problem = RestrictedProblem(restricted_model,dict_p,dict_pbar,dict_gamma,p_time,pbar_time,u,v,w,l,no_production,min_output,max_output,upper_bound_on_power_output_1,upper_bound_on_power_output_2,upper_bound_on_power_output_3,limit_power_jumps_up_1,limit_power_jumps_up_2,limit_power_jumps_down_1,limit_power_jumps_down_2,minimum_down_time_constraint,production_time,productionbar_time,commitment_status,startup_status,shutdown_status,loads);
	# Prior data
	if length(u_prior) > 0
		for g=1:nb_gen
			JuMP.fix(restricted_problem.u[g,0], u_prior[g]; force = true);
			JuMP.fix(restricted_problem.v[g,0], v_prior[g]; force = true);
			JuMP.fix(restricted_problem.w[g,0], w_prior[g]; force = true);
			JuMP.fix(restricted_problem.p_time[g,0], p_prior[g]; force = true);
			JuMP.fix(restricted_problem.pbar_time[g,0], pbar_prior[g]; force = true);
		end
	end
	# Posterior data
	if length(u_posterior) > 0
		for g=1:nb_gen
			JuMP.fix(restricted_problem.u[g,T_max+1],u_posterior[g]; force = true);
			JuMP.fix(restricted_problem.v[g,T_max+1],v_posterior[g]; force = true);
			JuMP.fix(restricted_problem.w[g,T_max+1],w_posterior[g]; force = true);
			JuMP.fix(restricted_problem.p_time[g,T_max+1],p_posterior[g]; force = true);
			JuMP.fix(restricted_problem.pbar_time[g,T_max+1],pbar_posterior[g]; force = true);
		end
	end
	# Objective function
	@objective(restricted_problem.model, Min, sum( sum( NoLoadConsumption[g]*C[g]*restricted_problem.u[g,t] + F[g]*restricted_problem.v[g,t] + C[g]*restricted_problem.p_time[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum(restricted_problem.l[t] for t=1:T_max));
	################## Start iterating
	beta = 0;
	iter_max = 500; # 500
	#obj_restricted_vector = zeros(iter_max);
	obj_restricted_vector = [];
	#obj_pricing_vector = zeros(iter_max);
	obj_pricing_vector = [];
	p_time_restricted = 0;
	pbar_time_restricted = 0;
	u_restricted = 0;
	v_restricted = 0;
	w_restricted = 0;
	price = 0;
	l_restricted = 0;
	obj_restricted = 0;
	iter_stop = 0;
	obj_pricing = 0;
	# phi = 10^(-3);
	phi = phi_criteria;
	Solving_time_master = 0;
	Solving_time_slave = 0;
	nb_intervals_prev = copy(nb_intervals_gen);
	
	price_vect = [];
	nb_var_vec = [];
	nb_cons_vec = [];
	Solving_time_master_vec = [];
	Solving_time_slave_vec = [];

	for iter=1:iter_max
		iter_stop = iter;
		# Solve the restricted problem to have the (dual) prices of load constraint.
		optimize!(restricted_problem.model);
		obj_restricted = objective_value(restricted_problem.model);
		price = dual.(restricted_problem.loads);
		push!(price_vect, price);
		push!(obj_restricted_vector,obj_restricted);
		
		push!(nb_var_vec, num_variables(restricted_problem.model));
		push!(nb_cons_vec, sum([num_constraints(restricted_problem.model, i, j) for (i,j) in list_of_constraint_types(restricted_problem.model)]));
		# push!(Solving_time_master_vec, MOI.get(restricted_problem.model, MOI.SolveTime()));
		push!(Solving_time_master_vec, 0);

		# Solving_time_master += MOI.get(restricted_problem.model, MOI.SolveTime());
		Solving_time_master += 0;
		# Solve the pricing problem
		@objective(pricing_problem.model, Min, sum( sum( NoLoadConsumption[g]*C[g]*pricing_problem.u[g,t] + F[g]*pricing_problem.v[g,t] + C[g]*pricing_problem.p[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum( pricing_problem.l[t] for t=1:T_max) - sum( price[t]*(sum( pricing_problem.p[g,t] for g=1:nb_gen) - pricing_problem.l[t]) for t=1:T_max) );
		optimize!(pricing_problem.model);
		# push!(Solving_time_slave_vec, MOI.get(pricing_problem.model, MOI.SolveTime()));
		push!(Solving_time_slave_vec, 0);
		# Solving_time_slave += MOI.get(pricing_problem.model, MOI.SolveTime());
		Solving_time_slave += 0;
		obj_pricing = objective_value(pricing_problem.model);
		if Printer
			println();
			println();
			println("[Column_And_Row_Generation] iteration ",iter);
			# println("[Column_And_Row_Generation] p_opt : $(value.(pricing_problem.p).data)");
			println("[Column_And_Row_Generation] (restricted prob) p_time : $(value.(restricted_problem.p_time).data)");
			println("[Column_And_Row_Generation] (restricted prob) u : $(value.(restricted_problem.u).data)");
			println("[Column_And_Row_Generation] price : $(price)");
			println("[Column_And_Row_Generation] obj_restricted : ",obj_restricted);
			println("[Column_And_Row_Generation] obj_pricing    : ",obj_pricing);
			#println(pricing_problem.model);
		end
		push!(obj_pricing_vector,obj_pricing)
		u_pricing = value.(pricing_problem.u).data;
		# Compute the Lagrangian dual bound
		if iter==1
			beta = obj_pricing;
		else
			beta = maximum([obj_pricing beta]);
		end
		if  obj_restricted <= beta + phi
			if Printer
				println();
				println();
				println("[Column_And_Row_Generation] OPTIMAL SOLUTION FOUND");
			end
			break; # STOP algorithm.
		end
		# Update current bundle
		if nb_gen>1
			u = u_pricing[:,2:T_max+1];
		else
			u = u_pricing[2:T_max+1];
		end
		A_new,B_new,nb_intervals_gen_new,A_added,B_added,nb_intervals_added = Update_A_B(u, nb_gen, A, B, nb_intervals_gen);
		A = copy(A_new);
		B = copy(B_new);
		nb_intervals_gen = copy(nb_intervals_gen_new);
		if  Printer
			println("[Column_And_Row_Generation_1] u : ",u);
			println("[Column_And_Row_Generation] A : ",A);
			println("[Column_And_Row_Generation] B : ",B);
			println("[Column_And_Row_Generation] nb_intervals_gen : ",nb_intervals_gen);
		end

		### Create variable
		for g=1:nb_gen
			for i=nb_intervals_gen[g]-nb_intervals_added[g]+1:nb_intervals_gen[g]
				# println("i : ",i,"-> [A,B] = [$(A[g,i]),$(B[g,i])]");
				new_var_p = @variable(restricted_problem.model, [g,i,0:T_max+1], base_name="p",lower_bound=0);
				new_var_pbar = @variable(restricted_problem.model, [g,i,0:T_max+1], base_name="pbar",lower_bound=0);
				for t=0:T_max+1
					restricted_problem.p[g,i,t] = new_var_p[g,i,t];
					restricted_problem.pbar[g,i,t] = new_var_pbar[g,i,t];
				end
				new_var_gamma = @variable(restricted_problem.model, [g,i],base_name="gamma",lower_bound=0);
				restricted_problem.gamma[g,i] = new_var_gamma[g,i];
			end
		end
		### Create new constraints
		@constraint(restricted_problem.model, [g=1:nb_gen, i=nb_intervals_gen[g]-nb_intervals_added[g]+1:nb_intervals_gen[g], t=0:T_max+1; (t<A[g,i] || t>B[g,i])], restricted_problem.p[g,i,t] <= 0);
		@constraint(restricted_problem.model, [g=1:nb_gen, i=nb_intervals_gen[g]-nb_intervals_added[g]+1:nb_intervals_gen[g], t=A[g,i]:B[g,i] ], -restricted_problem.p[g,i,t] <= -restricted_problem.gamma[g,i] * MinRunCapacity[g]);
		@constraint(restricted_problem.model, [g=1:nb_gen, i=nb_intervals_gen[g]-nb_intervals_added[g]+1:nb_intervals_gen[g], t=A[g,i]:B[g,i] ], restricted_problem.p[g,i,t] - restricted_problem.pbar[g,i,t] <= 0);
		@constraint(restricted_problem.model, [g=1:nb_gen, i=nb_intervals_gen[g]-nb_intervals_added[g]+1:nb_intervals_gen[g], t=A[g,i]:B[g,i] ] , restricted_problem.pbar[g,i,t] <= restricted_problem.gamma[g,i] * MaxRunCapacity[g]);
		@constraint(restricted_problem.model, [g=1:nb_gen, i=nb_intervals_gen[g]-nb_intervals_added[g]+1:nb_intervals_gen[g], t=A[g,i]:B[g,i] ] , restricted_problem.pbar[g,i,t] <= restricted_problem.gamma[g,i] * (SU[g] + (t - A[g,i])*RU[g] ));
		@constraint(restricted_problem.model, [g=1:nb_gen, i=nb_intervals_gen[g]-nb_intervals_added[g]+1:nb_intervals_gen[g], t=A[g,i]:B[g,i] ] , restricted_problem.pbar[g,i,t] <= restricted_problem.gamma[g,i] * (SD[g] + (B[g,i] - t)*RD[g] ));
		@constraint(restricted_problem.model, [g=1:nb_gen, i=nb_intervals_gen[g]-nb_intervals_added[g]+1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], restricted_problem.pbar[g,i,t] - restricted_problem.p[g,i,t-1] <= restricted_problem.gamma[g,i] * RU[g]);
		@constraint(restricted_problem.model, [g=1:nb_gen, i=nb_intervals_gen[g]-nb_intervals_added[g]+1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], restricted_problem.pbar[g,i,t] - restricted_problem.p[g,i,t-1] <= restricted_problem.gamma[g,i] * (SD[g] + (B[g,i] - t)*RD[g] - MinRunCapacity[g]));
		@constraint(restricted_problem.model, [g=1:nb_gen, i=nb_intervals_gen[g]-nb_intervals_added[g]+1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], restricted_problem.pbar[g,i,t-1] - restricted_problem.p[g,i,t] <= restricted_problem.gamma[g,i] * RD[g]);
		@constraint(restricted_problem.model, [g=1:nb_gen, i=nb_intervals_gen[g]-nb_intervals_added[g]+1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], restricted_problem.pbar[g,i,t-1] - restricted_problem.p[g,i,t] <= restricted_problem.gamma[g,i] * (SU[g] + (t - A[g,i])*RU[g] - MinRunCapacity[g]));
		### Update existing constraints
		for g=1:nb_gen
			for t=0:T_max+1
				for i=nb_intervals_gen[g]-nb_intervals_added[g]+1:nb_intervals_gen[g]
					if 1<=t && t<=T_max
						if A[g,i]<=t && t<=B[g,i]+DT[g]
							set_normalized_coefficient(restricted_problem.minimum_down_time_constraint[g,t], restricted_problem.gamma[g,i], 1);
						end
						set_normalized_coefficient(restricted_problem.production_time[g,t], restricted_problem.p[g,i,t], -1);
						set_normalized_coefficient(restricted_problem.productionbar_time[g,t], restricted_problem.pbar[g,i,t], -1);
					end
					if A[g,i]<=t && t<=B[g,i]
						set_normalized_coefficient(restricted_problem.commitment_status[g,t], restricted_problem.gamma[g,i], 1);
					end
					if t==A[g,i]
						set_normalized_coefficient(restricted_problem.startup_status[g,t], restricted_problem.gamma[g,i], 1);
					end
					if t==B[g,i]+1
						set_normalized_coefficient(restricted_problem.shutdown_status[g,t], restricted_problem.gamma[g,i], 1);
					end
				end
			end
		end
	end
	if Printer
		println("[Column_And_Row_Generation] Stop at iteration ",iter_stop);
	end
	p_time_restricted = value.(restricted_problem.p_time).data;
	pbar_time_restricted = value.(restricted_problem.pbar_time).data;
	u_restricted = value.(restricted_problem.u).data;
	v_restricted = value.(restricted_problem.v).data;
	w_restricted = value.(restricted_problem.w).data;
	l_restricted = value.(restricted_problem.l);
	return (p_time_restricted, pbar_time_restricted, u_restricted, v_restricted, w_restricted, price, l_restricted, obj_restricted, obj_pricing, Solving_time_slave+Solving_time_master, obj_restricted_vector, obj_pricing_vector, price_vect, iter_stop, Solving_time_master, Solving_time_slave, nb_var_vec, nb_cons_vec, Solving_time_master_vec, Solving_time_slave_vec);
end
