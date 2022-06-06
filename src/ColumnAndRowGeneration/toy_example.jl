using CSV, DataFrames, JuMP, Gurobi, Plots, GLPK, TickTock, LaTeXStrings, JSON;
include("src.jl");

### DATA
# Physical Constraints Generator(s)
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
SU                  = MinRunCapacity;
SD                  = MinRunCapacity;
# Value Of Lost Load
VOLL = 3000; # [€/MWh]
# Demand/Load
L = [6 11 16 11];
# L = [6 11 16 11 11 16 11 16 20];
T_max = length(L);

### Optimization programs
# Prior binary variables
u_prior 	   = zeros(nb_gen);
v_prior 	   = zeros(nb_gen);
w_prior 	   = zeros(nb_gen);
p_prior 	   = zeros(nb_gen);
pbar_prior 	   = zeros(nb_gen);
# Posterior binary variables
u_posterior    = zeros(nb_gen);
v_posterior    = zeros(nb_gen);
w_posterior    = zeros(nb_gen);
p_posterior    = zeros(nb_gen);
pbar_posterior = zeros(nb_gen);

### The set of intervals deals with the minimum up time when building the intervals such that : 1 <= a <= a+(UT-1) <= b <= T_max
(A,B,nb_intervals_gen) = set_T(UT, T_max, nb_gen);
println("A : ",A);
println("B : ",B);
println("nb_intervals_gen : ",nb_intervals_gen);
println();

### Unit Commitment problem
println("---------------------------------------------------------------------------------------");
println("Unit Commitment Problem");
(p_matching, pbar_matching, u_matching, v_matching, w_matching, l_matching, obj_matching, _) = Matching(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, L, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, true);
println("p = $(p_matching)");
println("pbar_UC = $(pbar_matching)");
println("u_UC = $(u_matching)");
println("v_UC = $(v_matching)");
println("w_UC = $(w_matching)");
println("l_UC = $(l_matching)");
println("obj_UC = $(obj_matching)");
println();

### Extended Formulation
println("---------------------------------------------------------------------------------------");
println("Extended Formulation");
(p_EF, pbar_EF, gamma_EF, p_time_EF, pbar_time_EF, u_EF, v_EF, w_EF, price_EF, l_val_EF, obj_EF, _, nb_var_EF, nb_con_EF) = Extended_Formulation(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, L, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, true);
gamma_data = gamma_EF.data;
println("gamma_data = $(gamma_data)");
for (g,int_interval) in keys(gamma_data)
    if nb_gen==1
        println("[a = $(A[int_interval]),b = $(B[int_interval])] : $(gamma_data[(g,int_interval)])");
    else
        println("[a = $(A[g,int_interval]),b = $(B[g,int_interval])] : $(gamma_data[(g,int_interval)])");
    end
end
println("p_time_EF    = $(p_time_EF) ");
println("pbar_time_EF = $(pbar_time_EF)");
println("u_EF = $(u_EF)");
println("v_EF = $(v_EF) ");
println("w_EF = $(w_EF)");
println("price_EF = $(price_EF)");
println("l_val_EF = $(l_val_EF)");
println("obj_EF = $(obj_EF)");

(_, _, _, _, _, _, obj_dual_lagrangian_EF, _) = Dual_Lagrangian(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, L, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price_EF);
println("obj_dual_lagrangian_EF = $(obj_dual_lagrangian_EF)");
println("UC - Dual Lagrangian(price_EF) = $(obj_matching) - $(obj_dual_lagrangian_EF) = $(obj_matching-obj_dual_lagrangian_EF)");

uplift_producer = zeros(nb_gen);
for g=1:nb_gen
    (max_profit, _, _, _, _, _) = MaxProfit_Producer(MinRunCapacity[g], MaxRunCapacity[g], RampUp[g], RampDown[g], UT[g], DT[g], SU[g], SD[g], T_max, F[g], C[g], NoLoadConsumption[g], price_EF, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, g);
    received_profit = sum(price_EF[t]*p_matching[g,t+1] - (NoLoadConsumption[g]*C[g]*u_matching[g,t+1] + F[g]*v_matching[g,t+1] + C[g]*p_matching[g,t+1]) for t=1:T_max);
    uplift_producer[g] = max_profit-received_profit;
end
uplift_consumer = 0;
(max_profit_consumer,_) = MaxProfit_Consumer(L, T_max, price_EF, VOLL);
received_profit_consumer = sum( (VOLL - price_EF[t])*l_matching[t] for t=1:T_max);
uplift_consumer = max_profit_consumer - received_profit_consumer;
println("Max Profit = sum($(uplift_producer)) + $(uplift_consumer) = $(sum(uplift_producer[i] for i=1:length(uplift_producer)) + uplift_consumer)");
println();


### Row Generation : Toy Example
println("---------------------------------------------------------------------------------------");
println("Row Generation (Bender decomposition)");
(iter_max, price_RG, obj_vec_RG, p_time_RG, pbar_time_RG, u_RG, v_RG, w_RG, l_RG, obj_RG, _) = Bender_Decomposition_1(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, L, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, 0, true);
println("price_RG = $(price_RG)");
# price = price_RG

println();
(_, _, _, _, _, _, obj_dual_lagrangian_RG, _) = Dual_Lagrangian(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, L, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price_RG);
println("\nobj_dual_lagrangian_RG = $(obj_dual_lagrangian_RG)");
println("UC - Dual Lagrangian(price_RG) = $(obj_matching) - $(obj_dual_lagrangian_RG) = $(obj_matching-obj_dual_lagrangian_RG)");

uplift_producer_RG = zeros(nb_gen);
for g=1:nb_gen
    (max_profit, _, _, _, _, _) = MaxProfit_Producer(MinRunCapacity[g], MaxRunCapacity[g], RampUp[g], RampDown[g], UT[g], DT[g], SU[g], SD[g], T_max, F[g], C[g], NoLoadConsumption[g], price_RG, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, g);
    received_profit = sum(price_RG[t]*p_matching[g,t+1] - (NoLoadConsumption[g]*C[g]*u_matching[g,t+1] + F[g]*v_matching[g,t+1] + C[g]*p_matching[g,t+1]) for t=1:T_max);
    uplift_producer_RG[g] = max_profit-received_profit;
end
uplift_consumer_RG = 0;
(max_profit_consumer,_) = MaxProfit_Consumer(L, T_max, price_RG, VOLL);
received_profit_consumer = sum( (VOLL - price_RG[t])*l_matching[t] for t=1:T_max);
uplift_consumer_RG = max_profit_consumer - received_profit_consumer;
println("Max Profit = sum($(uplift_producer_RG)) + $(uplift_consumer_RG) = $(sum(uplift_producer_RG[i] for i=1:length(uplift_producer_RG)) + uplift_consumer_RG)");
println();

### Column Generation - Dantzig Wolfe : Toy Example
println("---------------------------------------------------------------------------------------");
println("Column Generation (Dantzig-Wolfe)");
(obj_CG, price_CG, iter_stop_CG, _, _, _, obj_vec_CG) = Column_Generation_DW_1(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, L, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, true);
println("price_CG = $(price_CG)");
# price = price_CG;
println();
(_, _, _, _, _, _, obj_dual_lagrangian_CG, _) = Dual_Lagrangian(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, L, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price_CG);
println("\nobj_dual_lagrangian_CG = $(obj_dual_lagrangian_CG)");
println("UC - Dual Lagrangian(price_CG) = $(obj_matching) - $(obj_dual_lagrangian_CG) = $(obj_matching-obj_dual_lagrangian_CG)");

uplift_producer_CG = zeros(nb_gen);
for g=1:nb_gen
    (max_profit, _, _, _, _, _) = MaxProfit_Producer(MinRunCapacity[g], MaxRunCapacity[g], RampUp[g], RampDown[g], UT[g], DT[g], SU[g], SD[g], T_max, F[g], C[g], NoLoadConsumption[g], price_CG, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, g);
    received_profit = sum(price_CG[t]*p_matching[g,t+1] - (NoLoadConsumption[g]*C[g]*u_matching[g,t+1] + F[g]*v_matching[g,t+1] + C[g]*p_matching[g,t+1]) for t=1:T_max);
    uplift_producer_CG[g] = max_profit-received_profit;
end
uplift_consumer_CG = 0;
(max_profit_consumer,_) = MaxProfit_Consumer(L, T_max, price_CG, VOLL);
received_profit_consumer = sum( (VOLL - price_CG[t])*l_matching[t] for t=1:T_max);
uplift_consumer_CG = max_profit_consumer - received_profit_consumer;
println("Max Profit = sum($(uplift_producer_CG)) + $(uplift_consumer_CG) = $(sum(uplift_producer_CG[i] for i=1:length(uplift_producer_CG)) + uplift_consumer_CG)");
println();


# ### Column And Row Generation : Toy Example
# println("---------------------------------------------------------------------------------------");
# println("Column And Row Generation - 1");
# (p_time_CRG, pbar_time_CRG, u_CRG, v_CRG, w_CRG, price_CRG, l_CRG, obj_CRG, obj_CRG, _, obj_restricted_vector, obj_pricing_vector, price_CRG_iterates) = Column_And_Row_Generation_1(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, L, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, true);
# println("Price CRG = $(price_CRG)");
# println("price_CRG_iterates = $(price_CRG_iterates)");
# println();
# (_, _, _, _, _, _, obj_dual_lagrangian_CRG, _) = Dual_Lagrangian(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, L, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price_CRG);
# println("\nobj_dual_lagrangian_CRG = $(obj_dual_lagrangian_CRG)");
# println("UC - Dual Lagrangian(price_CRG) = $(obj_matching) - $(obj_dual_lagrangian_CRG) = $(obj_matching-obj_dual_lagrangian_CRG)");

# uplift_producer_CRG = zeros(nb_gen);
# for g=1:nb_gen
#     (max_profit, _, _, _, _, _) = MaxProfit_Producer(MinRunCapacity[g], MaxRunCapacity[g], RampUp[g], RampDown[g], UT[g], DT[g], SU[g], SD[g], T_max, F[g], C[g], NoLoadConsumption[g], price_CRG, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, g);
#     received_profit = sum(price_CRG[t]*p_matching[g,t+1] - (NoLoadConsumption[g]*C[g]*u_matching[g,t+1] + F[g]*v_matching[g,t+1] + C[g]*p_matching[g,t+1]) for t=1:T_max);
#     uplift_producer_CRG[g] = max_profit-received_profit;
# end
# uplift_consumer_CRG = 0;
# (max_profit_consumer,_) = MaxProfit_Consumer(L, T_max, price_CRG, VOLL);
# received_profit_consumer = sum( (VOLL - price_CRG[t])*l_matching[t] for t=1:T_max);
# uplift_consumer_CRG = max_profit_consumer - received_profit_consumer;
# println("Max Profit = sum($(uplift_producer_CRG)) + $(uplift_consumer_CRG) = $(sum(uplift_producer_CRG[i] for i=1:length(uplift_producer_CRG)) + uplift_consumer_CRG)");
# println();


println("---------------------------------------------------------------------------------------");
println("Column And Row Generation");
(_, _, _, _, _, price_CRG, _, obj_restricted, obj_pricing, _, _, _, price_vect, _, _, _, _, _, _, _) = Column_And_Row_Generation(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, L, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, true)
println("price_CRG = $(price_CRG)");
# println("price_vect");
# println(price_vect);
# println("obj_restricted = $(obj_restricted)");
# println("obj_pricing = $(obj_pricing)");

println();
(_, _, _, _, _, _, obj_dual_lagrangian_CRG, _) = Dual_Lagrangian(MinRunCapacity, MaxRunCapacity, RampUp, RampDown, UT, DT, SU, SD, T_max, nb_gen, F, C, NoLoadConsumption, L, VOLL, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, price_CRG);
println("\nobj_dual_lagrangian_CRG = $(obj_dual_lagrangian_CRG)");
println("UC - Dual Lagrangian(price_CRG) = $(obj_matching) - $(obj_dual_lagrangian_CRG) = $(obj_matching-obj_dual_lagrangian_CRG)");

uplift_producer_CRG = zeros(nb_gen);
for g=1:nb_gen
    (max_profit, _, _, _, _, _) = MaxProfit_Producer(MinRunCapacity[g], MaxRunCapacity[g], RampUp[g], RampDown[g], UT[g], DT[g], SU[g], SD[g], T_max, F[g], C[g], NoLoadConsumption[g], price_CRG, u_prior, v_prior, w_prior, p_prior, pbar_prior, u_posterior, v_posterior, w_posterior, p_posterior, pbar_posterior, g);
    received_profit = sum(price_CRG[t]*p_matching[g,t+1] - (NoLoadConsumption[g]*C[g]*u_matching[g,t+1] + F[g]*v_matching[g,t+1] + C[g]*p_matching[g,t+1]) for t=1:T_max);
    uplift_producer_CRG[g] = max_profit-received_profit;
end
uplift_consumer_CRG = 0;
(max_profit_consumer,_) = MaxProfit_Consumer(L, T_max, price_CRG, VOLL);
received_profit_consumer = sum( (VOLL - price_CRG[t])*l_matching[t] for t=1:T_max);
uplift_consumer_CRG = max_profit_consumer - received_profit_consumer;
println("Max Profit = sum($(uplift_producer_CRG)) + $(uplift_consumer_CRG) = $(sum(uplift_producer_CRG[i] for i=1:length(uplift_producer_CRG)) + uplift_consumer_CRG)");


