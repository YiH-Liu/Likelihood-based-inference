using Random
using Measures
using CSV
using DataFrames
using StatsBase
using FilePathsBase

# In this script we genrate the data required to plot Figure 2(b).

# Get the path of the current script
path = dirname(@__FILE__)

# Init(): Set up the Initial condition for stochstic model involve two subpopulation.
function init(W,H,Δ,x1_0,d1_0,x2_0,d2_0)
    #------------------------------------------------------------------------------------------------------------------
    # Input:
    # W: (number) The width of the lattice.
    # H: (number) The height of the lattice.
    # Δ: (number) The lattice size.
    # x1_0: (list) List of tuple, each tuple contain the location of agent from subpopulation 1 initially spread from.
    # d1_0: (vector) Lsit of number, each number represent the initial occupancy of agent from subpopulation 1 at associate loaction in x1_0.
    # x2_0: (list) List of tuple, each tuple contain the location of agent from subpopulation 2 initially spread from.
    # d2_0: (vector) Lsit of number, each number represent the initial occupancy of agent from subpopulation 2 at associate loaction in x2_0.
    # Output:
    # A: (matrix) Lattice in matrix form, A[j,i] is associate with site (J+1-j,i) of the lattice.
    # N: (number) Totol number of agents on the lattice.
    #------------------------------------------------------------------------------------------------------------------
    # Number of agent in subpopulation 1 and 2.
    N1 = 0; N2 = 0; 
    # Genrate the empty lattice matrix based on W and H.
    I = Int(W/Δ + 1); J = Int(H/Δ + 1); 
    A = zeros(J,I);
    # Place agent in subpopulation 1 on the lattice based on x1_0 and d1_0.
    index = 1
    for location in x1_0
        # Transform location in to lattice site index.
        i_low = Int(location[1]/Δ + 1); i_up = Int(location[2]/Δ + 1);
        # Get the number of agents in each colomn from colomn "i_low" to "i_up" from d1_0.
        C1_i = Int(d1_0[index] * J);index = index + 1
        # Place agent on the ith colomn of lattice based on d1_0.
        for i = i_low:i_up
            # Place C1_i agents randomly on the ith colomn.
            js = sample(1:J, C1_i, replace=false)
            for j in js
                A[j,i] = 1
                N1 = N1 + 1
            end
        end
    end
    # Place agent in subpopulation 2 on the lattice based on x2_0 and d2_0.
    index = 1
    for location in x2_0
        # Transform location in to lattice site index.
        i_low = Int(location[1]/Δ + 1); i_up = Int(location[2]/Δ + 1);
        # Get the number of agents in each colomn from colomn "i_low" to "i_up" from d2_0.
        C2_i = Int(d2_0[index] * J); index = index + 1
        # Place agent on the ith colomn of lattice based on d2_0.
        for i = i_low:i_up
            # find all empty site in ith colomn 
            A_i = A[:,i]
            empty_site = findall(j -> j == 0, A_i)
            # Place C2_i agents randomly on the empty site in ith colomn.
            js = sample(empty_site, C2_i, replace=false)
            for j in js
                A[j,i] = 2
                N2 = N2 + 1
            end
        end
    end
    N = N1 + N2
    return A,N
end 

# ctrfer(): Transform jth row and ith colomn of the lattice matrix A into lattice site (i,j).
function ctrfer(i,j,J)
    #------------------------------------------------------------------------------------------------------------------
    # Input:
    # i: (number) i-th colomn of the lattice matrix A.
    # j: (number) j-th row of the lattice matrix A.
    # J: (number) Number of rows in the lattice.
    # Output:
    # (site_i,site_j): (tuple) Lattice site (site_i,site_j) correspond to jth row and ith colomn of the lattice matrix A.
    #------------------------------------------------------------------------------------------------------------------
    site_i = i
    site_j = J+1 - j
    return (site_i,site_j)
end

# Indices(): Extract the Indices of each agent in lattice matrix A and store them in dataframe.
function Indices(A) 
    #------------------------------------------------------------------------------------------------------------------
    # Input:
    # A: (matrix) lattice in matrix form, A[j,i] is associate with site (J+1-j,i) of the lattice.
    # Output:
    # index1: (dataframe) index1.i stores all the i-indices of agents in subpopulation 1 on the lattice and the corresponding j-indices are stored in index1.j.
    # index2: (dataframe) index2.i stores all the i-indices of agents in subpopulation 2 on the lattice and the corresponding j-indices are stored in index2.j.
    #------------------------------------------------------------------------------------------------------------------ 
    (J,I) = size(A)
    #Indices for subpopulation 1 and 2
    indices_1 = []
    indices_2 = []
    # Extract the indices of each agent on lattice and store them in as list of tuple.
    for j in 1:J
        for i in 1:I
            if A[j,i] == 1
                site = ctrfer(i,j,J)
                push!(indices_1, site)
            elseif A[j,i] == 2
                site = ctrfer(i,j,J)
                push!(indices_2, site)
            end
        end
    end
    # i j index for subpopulation 1
    i_indices_1 = [index[1] for index in indices_1]
    j_indices_1 = [index[2] for index in indices_1]
    # i j index for subpopulation 2
    i_indices_2 = [index[1] for index in indices_2]
    j_indices_2 = [index[2] for index in indices_2]
    # store in a dataframe
    index1 = DataFrame(i = vec(i_indices_1),j = vec(j_indices_1))
    index2 = DataFrame(i = vec(i_indices_2),j = vec(j_indices_2))
    return index1,index2
end

# Count(): Simulate the count data for subpopulation 1 and 2.
function Count(index1,index2,W,Δ)
    #------------------------------------------------------------------------------------------------------------------
    # Input:
    # index1: (dataframe) index1.i stores all the i-indices of agents in subpopulation 1 on the lattice and the corresponding j-indices are stored in index1.j.
    # index2: (dataframe) index2.i stores all the i-indices of agents in subpopulation 2 on the lattice and the corresponding j-indices are stored in index2.j.
    # Δ: (number) The lattice size.
    # Output:
    # C_1: (vector) Colomn count data for subpopulation 1.
    # C_2: (vector) Colomn count data for subpopulation 2.
    #------------------------------------------------------------------------------------------------------------------
    # Genrate empty C1 and C2.
    I = Int(W/Δ + 1);
    C_1 = zeros(I,1); C_2 = zeros(I,1)
    # Count number of agent in subpopulation 1 in each colomn from the the i_indices.
    for i in index1.i
        C_1[i] = C_1[i]+1
    end
    # Count number of agent in subpopulation 2 in each colomn from the the i_indices.
    for i in index2.i
        C_2[i] = C_2[i]+1
    end
    return C_1,C_2
end

# realisation(): update the lattice matrix A from time t to time t + τ.
function realisation(A,P1,P2,ρ1,ρ2,I,J,N)
    #------------------------------------------------------------------------------------------------------------------
    # Input:
    # A: (matrix) Lattice in matrix form, A[j,i] is associate with site (J+1-j,i) of the lattice.
    # P1: (number) Probability of move for agnets in subpopulation 1.
    # ρ1: (number) Bias parametefor agnets in subpopulation 1.
    # P2: (number) Probability of move for agnets in subpopulation 2.
    # ρ2: (number) Bias parametefor agnets in subpopulation 2.
    # I: (number) Number of colomns in the lattice.
    # J: (number) Number of rows in the lattice.
    # N: (number) Totol number of agents on the lattice.
    # Output:
    # A: (matrix) Updated lattice in matrix form, A[j,i] is associate with site (J+1-j,i) of the lattice.
    #------------------------------------------------------------------------------------------------------------------
    # Probability to move to right for subpopulation 1 and 2.
    prob_right1 = (1 + ρ1) / 4; prob_right2 = (1 + ρ2) / 4
    # Number of select agent from subpopulation 1 and 2.
    n_1 = 0; n_2 = 0;
    # Random sequential update- select until we have total N agent.
    while n_1 + n_2 < N
        # Select a random ith colomn and jth row of the lattice matrix.
        i = rand(1:I); j = rand(1:J)
        if A[j,i] == 1 
            # An agent in subpopulation 1 is selected.
            n_1 = n_1 + 1
            if rand() <= P1
                # If the agent move decide the direction of movement.
                prob = rand()
                if prob <= 1/4   # move up
                    if j==1 # at boundary-periodic boundary condition
                        if A[J,i] == 0 
                            A[J,i] = 1
                            A[j,i] = 0
                        end
                    else # move as normal
                        if A[j-1,i] == 0 
                            A[j-1,i] = 1
                            A[j,i] = 0
                        end
                    end
                elseif prob <= 1/2 # move down
                    if j==J # at boundary-periodic boundary condition
                        if A[1,i] == 0 
                            A[1,i] = 1
                            A[j,i] = 0
                        end
                    else # move as normal
                        if A[j+1,i] == 0 
                            A[j+1,i] = 1
                            A[j,i] = 0
                        end
                    end
                elseif prob <= 1/2 + prob_right1 # move right
                    if i != I # not at boundary
                        if A[j,i+1] == 0
                            A[j,i+1] = 1
                            A[j,i] = 0
                        end
                    end
                else # move left
                    if i != 1 # not at boundary
                        if A[j,i-1] == 0
                            A[j,i-1] = 1
                            A[j,i] = 0
                        end
                    end
                end
            end
        elseif A[j,i] == 2
            # An agent in subpopulation 2 is selected. 
            n_2 = n_2 + 1
            if rand() <= P2
            # If the agent move decide the direction of movement.
                prob = rand()
                if prob <= 1/4   # move up 
                    if j == 1 # at boundary-periodic boundary condition
                        if A[J,i] == 0
                            A[J,i] = 2
                            A[j,i] = 0
                        end
                    else # move as normal
                        if A[j-1,i] == 0
                            A[j-1,i] = 2
                            A[j,i] = 0
                        end
                    end
                elseif  prob <= 1/2 # move down
                    if j == J # at boundary-periodic boundary condition
                        if A[1,i] == 0
                            A[1,i] = 2
                            A[j,i] = 0
                        end
                    else # move as normal
                        if A[j+1,i] == 0
                            A[j+1,i] = 2
                            A[j,i] = 0
                        end
                    end
                elseif  prob <= 1/2 + prob_right2 # move right
                    if i != I # not at boundary
                        if A[j,i+1] == 0
                            A[j,i+1] = 2
                            A[j,i] = 0
                        end
                    end
                else #move left
                    if i != 1 # not at boundary
                        if A[j,i-1] == 0
                            A[j,i-1] = 2
                            A[j,i] = 0
                        end
                    end
                end
            end
        end
    end
    return A
end

# discrete_simulation(): update the lattice matrix A from time 0 to time t.
function discrete_simulation(A,N,P1,P2,ρ1,ρ2,W,H,Δ,τ,t)
    #------------------------------------------------------------------------------------------------------------------
    # Input:
    # A: (matrix) Lattice in matrix form at initial condition, A[j,i] is associate with site (J+1-j,i) of the lattice.
    # N: (number) Totol number of agents on the lattice.
    # P1: (number) Probability of move for agnets in subpopulation 1.
    # ρ1: (number) Bias parametefor agnets in subpopulation 1.
    # W: (number) The width of the lattice.
    # H: (number) The height of the lattice.
    # Δ: (number) The lattice size.
    # τ: (number) The discrete time step duration.
    # Output:
    # A: (matrix) Updated lattice in matrix form at time t, A[j,i] is associate with site (J+1-j,i) of the lattice.
    #------------------------------------------------------------------------------------------------------------------
    I = Int(W/Δ + 1); J = Int(H/Δ + 1); 
    for i = τ:τ:t
        A = realisation(A,P1,P2,ρ1,ρ2,I,J,N)
    end 
    return A
end

# Lattice size.
W = 199; H = 19; Δ = 1
# Time step duration.
τ=1
# Discrete model parameters.
P1 = 1; P2 = 1
ρ1 = 0;ρ2 = 0
# Simuation time.
t = 5000
# Inital condition
#---------------------------
# Agent in subpopulation 1 are initially occpied 50% of site at 1<=i<=55 and 146<=i<=200, 
# which correspond to Cartesian coordinate 0<=x<=54 and 145<=x<=199.
# The remaining sites in these columns are occupied by agents from subpopulation 2 
x1_0 = [(0,54),(145,199)]; x2_0 = [(0,54),(145,199)]
d1_0 = [0.5,0.5]; d2_0 = [0.5,0.5]


# genrate the initial condition
A0,N= init(W,H,Δ,x1_0,d1_0,x2_0,d2_0)

# genrate index for subpopulation 1 and 2 at initial condition
index1_0,index2_0 = Indices(A0) 

# save the index for subpopulation 1 and 2 at initial condition
CSV.write("$path\\index1_0.csv", index1_0)
CSV.write("$path\\index2_0.csv", index2_0)


# genrate count data for subpopulation 1 and 2 at initial condition
C_1_0,C_2_0 = Count(index1_0,index2_0,W,Δ)

# save the count data for subpopulation 1 and 2 at initial condition.
datas_0 = DataFrame(a = vec(C_1_0),b = vec(C_2_0))
CSV.write("$path\\data_0.csv", datas_0)

# simulate the discrete model at time t
@time A=discrete_simulation(A0,N,P1,P2,ρ1,ρ2,W,H,Δ,τ,t)

# genrate index for subpopulation 1 and 2 at time t
index1,index2 = Indices(A) 

# save the index for subpopulation 1 and 2 at time t
CSV.write("$path\\index1.csv", index1)
CSV.write("$path\\index2.csv", index2)

# genrate count data 
C_1,C_2 = Count(index1,index2,W,Δ)

# save the observation data: count data for subpopulation 1 and 2.
datas = DataFrame(a = vec(C_1),b = vec(C_2))
CSV.write("$path\\data.csv", datas)

