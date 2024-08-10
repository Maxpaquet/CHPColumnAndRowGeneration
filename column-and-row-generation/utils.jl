


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