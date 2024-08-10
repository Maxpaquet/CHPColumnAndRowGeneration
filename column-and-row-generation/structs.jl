

function hello_world()
    print("Test hello world")
end

# Column generation structs

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


# DW structs

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


# Column and row structs

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