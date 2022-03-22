from gurobipy import *
import math
import numpy as np
import time
import networkx as nx
import random

def createGraph(input_file):
    global G, n, m, a0, start_time, apsp_time, n_nodes, n_edges, conn_comp, av_degree, density, card_comp
    
    G = nx.Graph()
    for j in range(0,n):
        G.add_node(j)
   
    f = open(input_file, "r")
    string = f.readline()
    for i in range(0, m):
        string = f.readline()
        string = string.split()
        j = int(string[0])-1
        k = int(string[1])-1
        G.add_edge(j, k)
    f.close()
    
    for i in range(0, n):
        G.add_edge(i, i)
    
    print("graph created")
                
    conn_comp = nx.number_connected_components(G)
    n_nodes   = len(nx.nodes(G))
    n_edges   = len(nx.edges(G))
    av_degree = 0
    for e in nx.degree(G):
        av_degree = av_degree + e[1]
    av_degree = av_degree / n_nodes
    density   = nx.density(G)
    print("conn. comp.: " + str(conn_comp))
    print("num vertices: " + str(n_nodes))
    print("num edges: " + str(n_edges))
    print("density: " + str(density))
    print("av degree: " + str(av_degree))
    print("---Compute all pairs shortest path - running time: %s seconds ---" % (time.time() - start_time))

def run():
    global total_runtime, runtime, n, m, a0, feasible, best_sequence, T_input, F
    try:
        m = Model("mip1")
        m.Params.outputFlag = 1  # 0 - Off  //  1 - On
        m.setParam("MIPGap", 0.0);
        m.setParam("Presolve", 2); # -1 - Automatic // 0 - Off // 1 - Conservative // 2 - Aggresive 
        #m.setParam("PreQLinearize", -1); # -1 - Automatic // 0 - Off // 1 - Strong LP relaxation // 2 - Compact relaxation
        #m.params.BestObjStop = k
        
        # ------------------------------INPUT --------------------------------
        T = T_input
        #---------------------------- VARIABLES -----------------------------------------------------------

        z = []
        for i in range(n):
            temp = []
            for j in range(T):
                temp.append(0)
            z.append(temp)
        for i in range(n):  #-------------------------------------------------- (A.1)
            for j in range(T):
                z[i][j] = m.addVar(vtype=GRB.BINARY, name="z,%s" % str(i+1) + "," + str(j+1))

        b = []
        for i in range(n):
            temp = []
            for j in range(T):
                temp.append(0)
            b.append(temp)
        for i in range(n):  #-------------------------------------------------- (A.1)
            for j in range(T):
                b[i][j] = m.addVar(vtype=GRB.BINARY, name="b,%s" % str(i+1) + "," + str(j+1))
        
        b_prime = [] # -------------------------------------------------------- (B)
        for i in range(T):
            b_prime.append(0)
            b_prime[i] = m.addVar(vtype=GRB.BINARY, name="b_prime,%s" % str(i+1))

        s = [] #--------------------------------------------------------------- (A.2)
        for i in range(n):
            temp = []
            for j in range(T):
                temp.append(0)
            s.append(temp)
        for i in range(n):
            for j in range(T):
                s[i][j] = m.addVar(vtype=GRB.BINARY, name="s,%s" % str(i+1) + "," + str(j+1))
        
        #---------------------------- CONSTRAINTS ---------------------------------------------------------
        
        for i in range(n): #---------------------------------------------------(2)--
            for j in range(1,T):
                m.addConstr(b[i][j] >= b[i][j-1])
        
        for i in range(n): #---------------------------------------------------(3)--
            for j in range(1,T):
                for k in G.neighbors(i):
                    m.addConstr(b[i][j] >= b[k][j-1])
                        
        for i in range(n): #---------------------------------------------------(4)--
            for j in range(1,T):
                for k in G.neighbors(i):
                    m.addConstr(b[i][j] >= s[k][j-1])

        sum_ = 0
        for i in range(n): # -------------------------------(5)--
            sum_ = sum_ + b[i][0]
        m.addConstr(sum_ == 0)

        s_transpose = np.array(s).T.tolist()# --------------(6 y 7)--
        for i in range(T):
            if i == 0:
                m.addConstr(sum(s_transpose[i]) <= F)
            else:
                m.addConstr(sum(s_transpose[i]) <= sum(s_transpose[i-1]) + F)

        for i in range(n): #-------------------------------(8)--
            for j in range(1,T):
                m.addConstr(s[i][j] >= s[i][j-1])
                
        for i in range(n): #--------------------------------(9, 10 y 11)--
            for j in range(T):
                m.addConstr(z[i][j] >= b[i][j])
                m.addConstr(z[i][j] >= s[i][j])
                m.addConstr(z[i][j] <= b[i][j] + s[i][j])
                
        epsilon = 0.0001 #------------------------------------------(12)--
        for j in range(T): 
            sum_ = 0
            for i in range(n):
                sum_ = sum_ + z[i][j]
            m.addConstr(b_prime[j] <= (sum_ / n))
            m.addConstr(b_prime[j] + 1 >= (sum_ / n) + epsilon)
            
        m.addConstr(b_prime[0] == 0) #--------------------------------(13)--
 
        for j in range(1,T): # --------------------------------------- (14)--
            m.addConstr(b_prime[j] >= b_prime[j-1])               
                        
        for i in range(n): #--------------------------------------------(15)
            for j in range(1,T):
                sum_ = 0
                for k in G.neighbors(i):
                    sum_ = sum_ + b[k][j-1]
                for k in G.neighbors(i):
                    sum_ = sum_ + s[k][j-1]
                m.addConstr(b[i][j] <= sum_)

        #---------------------------- OBJECTIVE FUNCTION --------------------------------------------------
        
        b_transpose = np.array(b).T.tolist()
        m.setObjective(T - sum(b_prime), GRB.MINIMIZE)#------------------------(1)--
                
        #---------------------------- OPTIMIZATION -------------------------------------------------------
        
        #m.computeIIS()
        #m.write("model_or.lp")
        
        m.optimize()
        runtime = m.Runtime
        #print("The run time is %f" % runtime)
        if m.status == GRB.INFEASIBLE:
            #print('Model is infeasible')
            feasible = False
        else:
            print("Obj:", m.objVal)
            feasible = True
            b_out       = []
            s_out       = []
            b_prime_out = []
            for v in m.getVars():
                varName = v.varName
                varNameSplit = varName.split(',')
                if varNameSplit[0] == 'b':
                    b_out.append(v.x)
                    #print(str(varName) + ": " + str(v.x))
                if varNameSplit[0] == 's':
                    #print(str(varName) + ": " + str(v.x))
                    s_out.append(v.x)
                if varNameSplit[0] == 'b_prime':
                    b_prime_out.append(v.x)
            #print(b_out)
            #print(b_prime_out)
            
            b_prime = []
            for i in range(len(b_prime_out)):
                if b_prime_out[i] > 0.9:
                    b_prime.append(1)
                else:
                    b_prime.append(0)
            
            sol_size = T - sum(b_prime) + 1
            
            s = []
            for i in range(n):
                temp = []
                for j in range(T):
                    temp.append(0)
                s.append(temp)
            j = 0
            k = 0
            for i in range(len(s_out)):
                if s_out[i] > 0.9:
                    s[j][k] = 1
                else:
                    s[j][k] = 0
                if (i+1)%T == 0:
                    j = j + 1
                    k = -1
                k =  k + 1
            sol = []
            for j in range(sol_size):
                temp = []
                for i in range(n):
                    if j == 0:
                        if s[i][j] == 1:
                            temp.append(i)
                    else:
                        if s[i][j] == 1 and s[i][j-1] == 0:
                            temp.append(i)
                sol.append(temp)
            print("OPT: " + str(sol))
            print("b(G): " + str(sol_size))
            print("The run time is %f" % runtime)
    except GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

def main(n_in, m_in, input_file_in, F_input):
    global n, m, start_time, T_input, F
    start_time = time.time()
    input_file  = input_file_in
    n = n_in
    m = m_in
    F = F_input
    createGraph(input_file)
    T_input = math.ceil(math.sqrt(n/F)) + 0
    run()
    
if __name__ == "__main__":
    global instance, T_input
    folder_dataset = 'C:/Users/jgd/Documents/FireFighter/'
    dataset = [
        #['test.mtx' , 10, 9, [0], 1], # File name, n, m, fire sources, number of firefighters
        #['test2.mtx', 20, 19, [0,19], 1]
        #['ba300_1.mtx' , 300, 299]
        #['ba100_1.mtx' , 100, 99]
        #['polbooks.mtx' , 105, 441]
        #['reds200_10.mtx' , 200, 380]
        #['lattice2d.mtx', 1089, 2112, 1]
        #['ca-netscience.mtx', 379, 914, 1]
        #['web-polblogs.mtx', 643, 2280, 1]
        #['dolphins.mtx', 62, 159, [0], 1]
        #['karate.mtx', 34, 78, [0], 1]
        #['lattice10x10.mtx', 100, 180]
        #["path4.mtx", 4, 3]
        #["path4.mtx", 4, 3, 1] # file, n, m, number of fire sources
        #["path9.mtx", 9, 8, 1] # file, n, m, number of fire sources
        #["path25.mtx", 25, 24, 2] # file, n, m, number of fire sources
        #["path16.mtx", 16, 15, 4] # file, n, m, number of fire sources
        #["path36.mtx", 36, 35, 9] # file, n, m, number of fire sources
        ["line49nodes.mtx", 49, 48, 1]
        ]
    for i in range(len(dataset)):
        print("--------------------------------------------------------------")
        instance = dataset[i][0]
        print("instance: " + instance)
        main(dataset[i][1], dataset[i][2], folder_dataset + dataset[i][0], dataset[i][3])
        