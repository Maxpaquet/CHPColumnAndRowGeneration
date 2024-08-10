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

### STEP1 : Initialize A, B, nb_intervals_gen
u_matching = [1.0, 1.0, 1.0, 0.0];
nb_gen = 1;
(A, B, nb_intervals_gen) = Compute_A_B(u_matching, nb_gen);


u = [1.0, 0.0, 1.0, 0.0];
# nb_intervals_gen = [2];
(A_new, B_new, nb_intervals_gen, A_added, B_added, nb_intervals_added) = Update_A_B(u, nb_gen, A, B, nb_intervals_gen);
println("A_new = $(A_new)");
println("B_new = $(B_new)");
println("nb_intervals_gen = $(nb_intervals_gen)");
println("A_added = $(A_added)");
println("B_added = $(B_added)");
println("nb_intervals_added = $(nb_intervals_added)");
