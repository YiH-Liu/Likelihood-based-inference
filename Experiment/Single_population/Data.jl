using Random
using CSV
using DataFrames
using StatsBase
using FilePathsBase

# In this script we genrate the data required to plot Figure 2(a).

# Get the path of the current script.
path = dirname(@__FILE__)

# Init(): Set up the Initial condition for stochstic model involve single population.
function init(W,H,Δ,x_0,d_0)
    #------------------------------------------------------------------------------------------------------------------
    # Input:
    # W: (number) The width of the lattice.
    # H: (number) The height of the lattice.
    # Δ: (number) The lattice size.
    # x_0: (list) List of tuple, each tuple contain the location of agent initially spread from.
    # d_0: (vector) Lsit of number, each number represent the initial occupancy of agent at associate loaction in x_0.
    # Output:
    # A: (matrix) Lattice in matrix form, A[j,i] is associate with site (J+1-j,i) of the lattice.
    # N: (number) Totol number of agents on the lattice.
    #------------------------------------------------------------------------------------------------------------------
    # Number of agent in subpopulation 1.
    N = 0; 
    # Genrate the empty lattice based on W and H.
    I = Int(W/Δ + 1); J = Int(H/Δ + 1); 
    A = zeros(J,I);
    # Place agent on the lattice based on x_0 and d_0.
    index = 1
    for location in x_0
        # Transform location in to lattice site index.
        i_low = Int(location[1]/Δ + 1); i_up = Int(location[2]/Δ + 1);
        # Get the number of agents in each colomn from colomn "i_low" to "i_up" from d_0.
        C_i = Int(d_0[index] * J); index = index + 1;
        # Place agent on the ith colomn of lattice based on d_0.
        for i = i_low:i_up
            # Place the C_i agents randomly on the ith colomn.
            js = sample(1:J, C_i, replace=false)
            for j in js
                A[j,i] = 1
                N = N + 1
            end
        end
    end
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
    # index: (dataframe) index.i stores all the i-indices of agents on the lattice and the corresponding j-indices are stored in index.j.
    #------------------------------------------------------------------------------------------------------------------
    (J,I) = size(A)
    #Indices for subpopulation 1.
    indices = []
    # Extract the indices of each agent on lattice and store them in as list of tuple.
    for j in 1:J
        for i in 1:I
            if A[j,i] == 1
                site = ctrfer(i,j,J)
                push!(indices, site)
            end
        end
    end
    # i j index for subpopulation 1.
    i_indices = [index[1] for index in indices]
    j_indices = [index[2] for index in indices]
    # Store in a dataframe
    index = DataFrame(i = vec(i_indices),j = vec(j_indices))
    return index
end

# Count(): Simulate the count data for subpopulation 1.
function Count(index,W,Δ)
    #------------------------------------------------------------------------------------------------------------------
    # Input:
    # index: (dataframe) index.i stores all the i-indices of agents on the lattice and the corresponding j-indices are stored in index.j.
    # W: (number) The width of the lattice.
    # Δ: (number) The lattice size.
    # Output:
    # C_1: (vector) Colomn count data for subpopulation 1.
    #------------------------------------------------------------------------------------------------------------------
    # Genrate empty C_1.
    I = Int(W/Δ + 1);
    C_1 = zeros(I,1)
    # Count number of agent in each colomn from the the i_indices.
    for i in index.i
        C_1[i] = C_1[i]+1
    end
    return C_1
end

# realisation(): update the lattice matrix A from time t to time t + τ.
function realisation(A,P1,ρ1,I,J,N)
    #------------------------------------------------------------------------------------------------------------------
    # Input:
    # A: (matrix) Lattice in matrix form, A[j,i] is associate with site (J+1-j,i) of the lattice.
    # P1: (number) Probability of move for agnets in subpopulation 1.
    # ρ1: (number) Bias parametefor agnets in subpopulation 1.
    # I: (number) Number of colomns in the lattice.
    # J: (number) Number of rows in the lattice.
    # N: (number) Totol number of agents on the lattice.
    # Output:
    # A: (matrix) Updated lattice in matrix form, A[j,i] is associate with site (J+1-j,i) of the lattice.
    #------------------------------------------------------------------------------------------------------------------
    # Probability to move to right for subpopulation 1.
    prob_right = (1 + ρ1) / 4
    # Number of select agent from subpopulation 1
    n1 = 0;
    # Random sequential update- select until we have total N agent.
    while n1 < N
        # Select a random ith colomn and jth row of the lattice matrix.
        i = rand(1:I); j = rand(1:J)
        if A[j,i] == 1 
            # An agent in subpopulation 1 is selected.
            n1 = n1 + 1
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
                elseif prob <= 1/2 + prob_right # move right
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
        end
    end
    return A
end

function discrete_simulation(A,N,P1,ρ1,W,H,Δ,τ,t)
    #------------------------------------------------------------------------------------------------------------------
    # Input:
    # A: (matrix) Lattice in matrix form at initial condition, A[j,i] is associate with site (J+1-j,i) of the lattice.
    # N: (number) Totol number of agents on the lattice.
    # P1: (number) Probability of move for agnets in subpopulation 1.
    # ρ1: (number) Bias parametefor agnets in subpopulation 1.
    # W: (number) The width of the lattice.
    # H: (number) The height of the lattice.
    # Δ: (number) The lattice site size.
    # τ: (number) The discrete time step duration.
    # Output:
    # A: (matrix) Updated lattice in matrix form at time t, A[j,i] is associate with site (J+1-j,i) of the lattice.
    #------------------------------------------------------------------------------------------------------------------
    I = Int(W/Δ + 1); J = Int(H/Δ + 1); 
    for i = τ:τ:t
        A = realisation(A,P1,ρ1,I,J,N)
    end 
    return A
end

# Lattice size.
W = 199; H = 19; Δ = 1
# Time step duration.
τ=1
# Discrete model parameters.
P1 = 1; ρ1 = 0
# Simuation time.
t = 5000

# Inital condition
#---------------------------
# Agent in subpopulation 1 are initial fully occpied site at 1<=i<=55 and 146<=i<=200, 
# which correspond to Cartesian coordinate 0<=x<=54 and 145<=x<=199. 
x_0 = [(0,54),(145,199)]
d_0 = [1,1]

# Genrate the initial condition.
A0,N= init(W,H,Δ,x_0,d_0)

# Extract indices for subpopulation 1 at initial condition.
index_0 = Indices(A0) 

# save the indices for subpopulation 1 at initial condition for snapshot.
CSV.write("$path\\index_0.csv", index_0)


# Genrate count data for subpopulation 1 at initial condition.
C_0 = Count(index_0,W,Δ)

# save the count data for subpopulation 1 at initial condition.
datas_0 = DataFrame(a = vec(C_0))
CSV.write("$path\\data_0.csv", datas_0)

# Simulate the discrete model at time t.
@time A=discrete_simulation(A0,N,P1,ρ1,W,H,Δ,τ,t)

# Extract indices for subpopulation 1 at time t.
index = Indices(A) 

# Save the indices for subpopulation 1 at time t.
CSV.write("$path\\index.csv", index)


# Genrate count for subpopulation 1 at t.
C = Count(index,W,Δ)

# Save the count data for subpopulation 1 at t.
datas = DataFrame(a = vec(C))
CSV.write("$path\\data.csv", datas) 

