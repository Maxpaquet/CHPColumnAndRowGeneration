using CSV, DataFrames, JuMP, Gurobi, Plots, GLPK, TickTock, PyPlot, LaTeXStrings, JSON;

# const MOI = MathOptInterface
const GUROBI_ENV = Gurobi.Env();

"""
@in :   - UT (array of size 1 x nb_gen) : minimum Up Time of the considered generator
		- T_max (scalar) : Number of time periods considered
		- nb_gen (scalar) : number of generators that are considered

@out :	- A (array of size 1 x nb_gen) : beginning (a) of all intervals [a,b] such that 1 <= a <= a+UT[g] <= b <= T_max, [0,b], [a,T+1] and [0,T+1]
		- B (array of size 1 x nb_gen) : ending (b) of all intervals [a,b] such that 1 <= a <= a+UT[g] <= b <= T_max, [0,b], [a,T+1] and [0,T+1]

(A,B,index_gen) = set_T(1,5, 1);
println("A : ",A);
println("B : ",B);
println("index_gen : ",index_gen)
"""
function set_T(UT,T_max, nb_gen)
	index_gen = zeros(nb_gen);
	for g=1:nb_gen
		index = 0;
		# Intervals for generator "on" on prior, intervals [0,b]
		# index+=T_max+1; # [0,b] for b=0:T_max
		# Intervals [a,b] such that 1 <= a <= a+(UT-1) <= b <= T_max
		for a=1:T_max
			for b=a+UT[g]-1:T_max
				if(1 <= a && a <= a+UT[g]-1 && a+UT[g]-1 <= b && b <= T_max)
					index += 1;
				end
			end
		end
		# Intervals for generator "on" past time T_max, intervals [a,T_max+1]
		# index+=T_max+1; # [a,T_max+1] for a=1:T_max+1
		# Interval for all period [0,T_max+1]
		# index+=1;
		index_gen[g] = index;
	end

	A = zeros(nb_gen,convert(Int64,maximum(index_gen)));
	B = zeros(nb_gen,convert(Int64,maximum(index_gen)));
	for g=1:nb_gen
		index = 1

		#=a = 0;
		for b=0:T_max
			A[g,index] = a; B[g,index] = b;
			index+=1;
		end=#

		for a=1:T_max
			for b=a+UT[g]-1:T_max
				if(1 <= a && a <= a+UT[g]-1 && a+UT[g]-1 <= b && b <= T_max)
					A[g,index] = a;
					B[g,index] = b;
					index+=1
				end
			end
		end

		#=b = T_max+1;
		for a=1:T_max+1
			A[g,index] = a;B[g,index] = b;
			index+=1;
		end=#
		#A[g,index] = 0; B[g,index] = T_max+1; index+=1;
	end
	A = convert(Matrix{Int64}, A);
	B = convert(Matrix{Int64}, B);
	index_gen = convert(Array{Int64}, index_gen);
	return (A,B,index_gen);
end

"""
Add the vector vect (containing optimal dual values) and the corresponding generator g, to the matrice mat.
"""
function Add_vector(mat, vect, gen)
	(a,b) = size(mat);
	mat_new = zeros(a,b+1);
	mat_new[1:a,1:b] = mat;
	mat_new[2:a,b+1] = vect;
	mat_new[1,b+1] = gen;
	return mat_new;
end

"""
Initialize matrix mat with the vector vect and the generator g.
"""
function Add_vector_first(mat, vect, gen)
	(a,b) = size(mat);
	mat[2:a,1] = vect;
	mat[1,1] = gen;
	return mat;
end

"""
Change values in vectors UT and DT higher than T_max by T_max
"""
function check_UT_DT(vec, T_max)
	new_vec = zeros(length(vec))
	for i=1:length(vec)
		if vec[i]>T_max
			new_vec[i] = T_max;
		else
			new_vec[i] = vec[i];
		end
	end
	return convert(Array{Int64}, new_vec);
end

"""
Look which generators are the same and gather them together.
"""
function Gather_Generators(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, G_c)
	gather_gen = []
    indices_gen = []
    #a = Int64[i for i=1:length(MinRunCapacity)]
    indices_gen = Int64[0 for i=1:length(MinRunCapacity)]
    visited = Int64[];
	count = 1;
	for i in G_c
        if i ∉ visited
            l = Int64[i]
            if i+1<=length(MinRunCapacity)
                for j=i+1:length(MinRunCapacity)
                    if j∉visited && MinRunCapacity[i]==MinRunCapacity[j] && MaxRunCapacity[i]==MaxRunCapacity[j] && RU[i]==RU[j] && RD[i]==RD[j] && UT[i]==UT[j] && DT[i]==DT[j] && SD[i]==SD[j] && SU[i]==SU[j]
                        ### generator i and j have the same parameters
                        append!(l,j);
                        append!(visited,j);
                        indices_gen[j] = count;
                        # filter!(e->e!=j,a)
                    end
                end
            end
			indices_gen[i] = count;
            append!(visited, i);
            append!(gather_gen,[l]);
			count+=1;
        end
	end
    return gather_gen,indices_gen
end

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
	m_matching = JuMP.direct_model(Gurobi.Optimizer(GUROBI_ENV)); # OutputFlag=0
	# m_matching = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
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
	m = JuMP.direct_model(Gurobi.Optimizer(GUROBI_ENV)); # OutputFlag=0
	JuMP.set_optimizer_attribute(m,"OutputFlag",0);

	#JuMP.set_optimizer_attribute(m,"FeasibilityTol",1e-9);
	#JuMP.set_optimizer_attribute(m,"OptimalityTol",1e-9);
	#JuMP.set_optimizer_attribute(m,"BarConvTol",1e-8);
	# Variables
	@variable(m, p[g=1:nb_gen, i=1:nb_intervals_gen[g], t=0:T_max+1], lower_bound = 0);	# Power output vector of the generators		p[generator g, interval(a,b) i, time t]
	@variable(m, pbar[g=1:nb_gen, i=1:nb_intervals_gen[g], t=0:T_max+1], lower_bound = 0);# Maximum power available from generators	pbar[generator g, interval(a,b) i, time t]
	@variable(m, gamma[g=1:nb_gen, i=1:nb_intervals_gen[g]], lower_bound = 0, upper_bound = 1);			# Variable that indicates whether a vertex (interval) is in the packing or not   gamma[generator g, interval(a,b) i]
	@variable(m, p_time[g=1:nb_gen, t=0:T_max+1], lower_bound = 0);						# (to simplify computations) production of each generator on each time period
	@variable(m, pbar_time[g=1:nb_gen, t=0:T_max+1], lower_bound = 0);						# (to simplify computations) maximum power available of each generator on each time period
	@variable(m, u[g=1:nb_gen, t=0:T_max+1], lower_bound = 0, upper_bound = 1); 							# Commitment status variables
	@variable(m, v[g=1:nb_gen, t=0:T_max+1], lower_bound = 0, upper_bound = 1); 							# Startup status variables
	@variable(m, w[g=1:nb_gen, t=0:T_max+1], lower_bound = 0, upper_bound = 1); 							# Shutdown status variables
	@variable(m, l[t=1:T_max], lower_bound = 0); # Modelling an inelastic demand
	@constraint(m, [t=1:T_max], l[t] <= data_demand[t]);
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
	@constraint(m, limit_power_jumps_up_1[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], pbar[g,i,t] - p[g,i,t-1] <= gamma[g,i] * RU[g]);
	@constraint(m, limit_power_jumps_up_2[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], pbar[g,i,t] - p[g,i,t-1] <= gamma[g,i] * (SD[g] + (B[g,i] - t)*RD[g] - MinRunCapacity[g]));
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
	@constraint(m, limit_power_jumps_down_1[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], pbar[g,i,t-1] - p[g,i,t] <= gamma[g,i] * RD[g]);
	@constraint(m, limit_power_jumps_down_2[g=1:nb_gen, i=1:nb_intervals_gen[g], t=A[g,i]+1:B[g,i] ], pbar[g,i,t-1] - p[g,i,t] <= gamma[g,i] * (SU[g] + (t - A[g,i])*RU[g] - MinRunCapacity[g]));
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
	@constraint(m, minimum_down_time_constraint[g=1:nb_gen, t=1:T_max], sum(gamma[g,i] for i=1:nb_intervals_gen[g] if (A[g,i] <= t && t <= B[g,i]+DT[g])) <= 1);
	# 15b : to be able to work with production on time instead of intervals in the objective function
	@constraint(m, production_time[g=1:nb_gen, t=1:T_max], p_time[g,t] == sum(p[g,i,t]  for i=1:nb_intervals_gen[g]));
	# 15c : to be able to work with production on time instead of intervals in the objective function
	@constraint(m, productionbar_time[g=1:nb_gen, t=1:T_max], pbar_time[g,t] == sum(pbar[g,i,t]  for i=1:nb_intervals_gen[g]));
	### Binary variable of UC (commitment status(on/off), startup status, shutdown status)
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
	optimize!(m);
	# Solving_time = MOI.get(m, MOI.SolveTime());
	Solving_time = 0;
	price = dual.(loads);
	u_val = value.(u);
	v_val = value.(v);
	p_time_val = value.(p_time);
	l_val = value.(l);

	if Printer
		nb_con = sum([num_constraints(m, i, j) for (i,j) in list_of_constraint_types(m)]);
		nb_var = num_variables(m);
		println("[EF] nb_con : $(nb_con) and nb_var : $(nb_var)");
	end
	nb_con = sum([num_constraints(m, i, j) for (i,j) in list_of_constraint_types(m)]);
	nb_var = num_variables(m);

	return (value.(p).data, value.(pbar).data, value.(gamma), value.(p_time).data, value.(pbar_time).data, value.(u).data, value.(v).data, value.(w).data, price, l_val, objective_value(m), Solving_time, nb_var, nb_con);
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
					println("[Bender_Decomposition_1] Need to add a cut for generator $(g), z_opt = $(z_opt)");
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


struct ModelCGLP
	model;
	z;
	p;
    pbar;
    p_time;
    pbar_time;
    u;
    v;
    w;
    gamma;
    min_output;
    max_output;
    upper_bound_on_power_output_1;
    upper_bound_on_power_output_2;
    upper_bound_on_power_output_3;
    limit_power_jumps_up_1;
    limit_power_jumps_up_2;
    limit_power_jumps_down_1;
    limit_power_jumps_down_2;
    minimum_down_time_constraint;
    power_output_global;            # rhs to be modified
	power_output_global_2;          # rhs to be modified
    maximum_power_available;        # rhs to be modified
    maximum_power_available_2;      # rhs to be modified
    commitment_status;              # rhs to be modified
    commitment_status_2;            # rhs to be modified
    startup_status;                 # rhs to be modified
    startup_status_2;               # rhs to be modified
    shutdown_status;                # rhs to be modified
    shutdown_status_2;              # rhs to be modified
end


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


function Updata_A_B(u, nb_gen, A, B, nb_intervals_gen)
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


struct PricingProblem
	model;
	p;
	pbar;
	u;
	v;
	w;
	l;
	logical_constraint;
	minimum_up_time;
	minimum_down_time;
	generation_limits_1;
	generation_limits_2;
	generation_limits_3;
	ramp_up_constraint;
	ramp_down_consraint;
end


struct RestrictedProblem
	model;
	p;
	pbar;
	gamma;
	p_time;
	pbar_time;
	u;
	v;
	w;
	l;
	no_production;
	min_output;
	max_output;
	upper_bound_on_power_output_1;
	upper_bound_on_power_output_2;
	upper_bound_on_power_output_3;
	limit_power_jumps_up_1;
	limit_power_jumps_up_2;
	limit_power_jumps_down_1;
	limit_power_jumps_down_2;
	minimum_down_time_constraint;
	production_time;
	productionbar_time;
	commitment_status;
	startup_status;
	shutdown_status;
	loads;
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
			JuMP.fix(pricing_problem.u[g,T_max+1],u_posterior[g]; force = true);
			JuMP.fix(pricing_problem.v[g,T_max+1],v_posterior[g]; force = true);
			JuMP.fix(pricing_problem.w[g,T_max+1],w_posterior[g]; force = true);
			JuMP.fix(pricing_problem.p[g,T_max+1],p_posterior[g]; force = true);
			JuMP.fix(pricing_problem.pbar[g,T_max+1],pbar_posterior[g]; force = true);
		end
	end
	### Objective function
	@objective(pricing_problem.model, Min, sum( sum( NoLoadConsumption[g]*C[g]*pricing_problem.u[g,t] + F[g]*pricing_problem.v[g,t] + C[g]*pricing_problem.p[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum( pricing_problem.l[t] for t=1:T_max) );
	optimize!(pricing_problem.model);
	u_matching = value.(pricing_problem.u).data;
	delete(pricing_problem.model, loads); # Only need the constraint for initialization. Don't need after because it will be dualize.
	if Printer
		println("[Column_And_Row_Generation_1] u_matching : ",u_matching);
	end
	if nb_gen>1
		u_matching = u_matching[:,2:T_max+1];
	else
		u_matching = u_matching[2:T_max+1];
	end	
	(A,B,nb_intervals_gen) = Compute_A_B(u_matching, nb_gen);
	if Printer
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
	price_iterates = [];
	l_restricted = 0;
	obj_restricted = 0;
	iter_stop = 0;
	obj_pricing = 0;
	# phi = 1e-3;
	phi = 10^(-3);
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
		push!(obj_restricted_vector, obj_restricted);
		push!(price_iterates, price);
		Solving_time += time;
		### Solve the pricing problem
		@objective(pricing_problem.model, Min, sum( sum( NoLoadConsumption[g]*C[g]*pricing_problem.u[g,t] + F[g]*pricing_problem.v[g,t] + C[g]*pricing_problem.p[g,t] for t=1:T_max) for g=1:nb_gen) - VOLL*sum( pricing_problem.l[t] for t=1:T_max) - sum( price[t]*(sum( pricing_problem.p[g,t] for g=1:nb_gen) - pricing_problem.l[t]) for t=1:T_max) );
		optimize!(pricing_problem.model);
		# Solving_time += MOI.get(pricing_problem.model, MOI.SolveTime());
		Solving_time += 0;
		obj_pricing = objective_value(pricing_problem.model);
		push!(obj_pricing_vector, obj_pricing);
		u_pricing = value.(pricing_problem.u).data;

		if Printer
			println("[Column_And_Row_Generation_1] price = $(price)")
			println("[Column_And_Row_Generation_1] obj_restricted : ",obj_restricted);
			println("[Column_And_Row_Generation_1] obj_pricing    : ",obj_pricing);
		end
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
		A_new,B_new,nb_intervals_gen_new, A_added, B_added, nb_intervals_added = Updata_A_B(u, nb_gen, A, B, nb_intervals_gen);
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
	return (p_time_restricted, pbar_time_restricted, u_restricted, v_restricted, w_restricted, price, l_restricted, obj_restricted, obj_pricing, Solving_time, obj_restricted_vector, obj_pricing_vector, price_iterates);
end



"""
Dynamic column-and-row generation. Implementation of hot-start model.
"""
function Column_And_Row_Generation(MinRunCapacity, MaxRunCapacity, RU, RD, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, data_demand, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, Printer, phi_criteria=10^(-5))
	if Printer
		println();
		println("##############################################################################################################");
		println("Pre-process");
		println("nb_gen : ",nb_gen);
	end
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
			println("[Column_And_Row_Generation : while loop] obj_restricted : ",obj_restricted);
			println("[Column_And_Row_Generation : while loop] obj_pricing    : ",obj_pricing);
			println("[Column_And_Row_Generation : while loop] price : $(price)");
			println("[Column_And_Row_Generation : while loop] p_opt : $(value.(pricing_problem.p))");
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
		A_new,B_new,nb_intervals_gen_new,A_added,B_added,nb_intervals_added = Updata_A_B(u, nb_gen, A, B, nb_intervals_gen);
		A = copy(A_new);
		B = copy(B_new);
		nb_intervals_gen = copy(nb_intervals_gen_new);
		if  Printer
			println("[Column_And_Row_Generation : while loop] A : ",A);
			println("[Column_And_Row_Generation : while loop] B : ",B);
			println("[Column_And_Row_Generation : while loop] nb_intervals_gen : ",nb_intervals_gen);
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



#######################################################################################################
#######################################################################################################
############################################ DW METHODS ###############################################
#######################################################################################################
#######################################################################################################

struct SubProblem
	model;
	p;
	pbar;
	u;
	v;
	w;
	cost;
	logical_constraint;
	minimum_up_time;
	minimum_down_time;
	generation_limits_1;
	generation_limits_2;
	generation_limits_3;
	ramp_up_constraint;
	ramp_down_consraint;
end

struct RestrictedMasterProgram
	model;
	z;
	l;
	balance_constraint;
	convex_combination;
end


function Restricted_Master_Program_Column_Generation(data_demand, counter_schedules, p_schedules, cost_schedules, nb_gen, T_max, VOLL)
	### Build Master Restricted Problem
	m_restricted = JuMP.direct_model(Gurobi.Optimizer()); # OutputFlag=0
	JuMP.set_optimizer_attribute(m_restricted,"OutputFlag",0);
	JuMP.set_optimizer_attribute(m_restricted,"FeasibilityTol",1e-3);

	@variable(m_restricted, z[g=1:nb_gen, i=1:counter_schedules[g] ], lower_bound = 0, upper_bound=1);
	@variable(m_restricted, l[t=1:T_max],lower_bound = 0);
	@constraint(m_restricted, [t=1:T_max], l[t]<=data_demand[t]);
	# Market clearing constraint
	@constraint(m_restricted, balance_constraint[t=1:T_max], sum( sum( sum( z[g,i]*prod for (time,prod) in enumerate(p_schedules[g,i]) if time==t)  for i=1:counter_schedules[g]) for g=1:nb_gen) - l[t] == 0 );
	# println();println();println(balance_constraint);
	@constraint(m_restricted, convex_combination[g=1:nb_gen], sum( z[g,i] for i=1:counter_schedules[g]) == 1 );
	# println();println();println(convex_combination);
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
			# Solving_time_slaves += MOI.get(sub_problem.model, MOI.SolveTime());
			Solving_time_slaves += 0;
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
