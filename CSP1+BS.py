
from gurobipy import *
import math
import numpy as np
import time
import networkx as nx
import random

def createGraph(input_file):
    global G, n, m, start_time, n_nodes, n_edges, conn_comp, av_degree, density, card_comp
    
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

def run(k, max_time):
    global G, runtime, n, m, feasible, best_sequence, sol_size
    try:
        m = Model("mip1")
        m.Params.outputFlag = 0  # 0 - Off  //  1 - On
        m.setParam("MIPGap", 0.0);
        m.setParam("Presolve", 1); # -1 - Automatic // 0 - Off // 1 - Conservative // 2 - Aggresive 
        #m.setParam("PreQLinearize", -1); # -1 - Automatic // 0 - Off // 1 - Strong LP relaxation // 2 - Compact relaxation
        #m.params.BestObjStop = k
        #m.setParam("TimeLimit", 10*3600)
        m.setParam("TimeLimit", max_time)

        T = k
        #---------------------------- VARIABLES -----------------------------------------------------------

        b = []
        for i in range(n):
            temp = []
            for j in range(T):
                temp.append(0)
            b.append(temp)
        for i in range(n):
            for j in range(T):
                b[i][j] = m.addVar(vtype=GRB.BINARY, name="b,%s" % str(i+1) + "," + str(j+1))

        s = []
        for i in range(n):
            temp = []
            for j in range(T):
                temp.append(0)
            s.append(temp)
        for i in range(n):
            for j in range(T):
                s[i][j] = m.addVar(vtype=GRB.BINARY, name="s,%s" % str(i+1) + "," + str(j+1))

        z = []
        for i in range(n):
            z.append(0)
        for i in range(n):
            z[i] = m.addVar(vtype=GRB.BINARY, name="z,%s" % str(i+1))
        
        #---------------------------- CONSTRAINTS ---------------------------------------------------------

        b0 = [] #-------------------------------------------------------------(1.11)--
        s0 = [] #-------------------------------------------------------------(1.12)--
        for i in range(n):
            b0.append(0)
            s0.append(0)
        
        for i in range(n): #--------------------------------------------------(1.2)--
            for j in range(T):
                if j == 0:
                    m.addConstr(s[i][j] >= s0[i])
                else:
                    m.addConstr(s[i][j] >= s[i][j-1])

        for i in range(n): #--------------------------------------------------(1.3)--
            for j in range(T):
                for k in G.neighbors(i):
                    if j == 0:
                        m.addConstr(b[i][j] >= b0[k])
                    else:
                        m.addConstr(b[i][j] >= b[k][j-1])
                        
        for i in range(n): #---------------------------------------------------(1.4)--
            for j in range(T):
                for k in G.neighbors(i):
                    if j == 0:
                        m.addConstr(b[i][j] >= s0[k])
                    else:
                        m.addConstr(b[i][j] >= s[k][j-1])
                    
        for i in range(n): #--------------------------------------------------(1.5)
            for j in range(T):
                sum_ = 0
                for k in G.neighbors(i):
                    if j == 0:
                        sum_ = sum_ + b0[k] + s0[k]
                    else:
                        sum_ = sum_ + b[k][j-1] + s[k][j-1]
                m.addConstr(b[i][j] <= sum_)

        s_transpose = np.array(s).T.tolist()# --------------------------------(1.6)--
        for i in range(T):
            if i == 0:
                sum_ = 0
                for j in range(n):
                    sum_ = sum_ + s[j][i] - s0[j]
                m.addConstr(sum_ == 1)
            else:
                sum_ = 0
                for j in range(n):
                    sum_ = sum_ + s[j][i] - s[j][i-1]
                m.addConstr(sum_ == 1)

        # most important constraint                        
        for i in range(n): #--------------------------------------------------( )--
            m.addConstr(z[i] >= b[i][T-1])
            m.addConstr(z[i] >= s[i][T-1])
            m.addConstr(z[i] <= b[i][T-1] + s[i][T-1])
        m.addConstr(sum(z) == n)

        #---------------------------- OBJECTIVE FUNCTION --------------------------------------------------
        
        m.setObjective(0, GRB.MINIMIZE)#------------------------(0)--
                
        #---------------------------- OPTIMIZATION -------------------------------------------------------
        
        m.optimize()
        runtime = m.Runtime
        if m.status == GRB.INFEASIBLE:
            feasible = False
        else:
            try:
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
                    if varNameSplit[0] == 's':
                        s_out.append(v.x)
                
                b = []
                for i in range(n):
                    temp = []
                    for j in range(T):
                        temp.append(0)
                    b.append(temp)
                j = 0
                k = 0
                for i in range(len(b_out)):
                    if b_out[i] > 0.9:
                        b[j][k] = 1
                    else:
                        b[j][k] = 0
                    if (i+1)%T == 0:
                        j = j + 1
                        k = -1
                    k =  k + 1
                
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
                    
                # Compute the solution size from s and b
                sol_size = 0
                for j in range(T):
                    num_burned = 0
                    for i in range(n):
                        if s[i][j] or b[i][j]:
                            num_burned += 1
                    if num_burned == n:
                        sol_size = j+1
                        break                        
                    
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
                print("MIP gap value: %f" % m.MIPGap)
                best_sequence = sol
            except:
                print("Something happened")
    except GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

    
if __name__ == "__main__":
    global instance, T_input, best_sequence, sol_size, start_time, time_out
    folder_dataset = 'C:/Users/jgd/Documents/GBP/dataset/paths/'
    dataset = [
        ['path16.mtx' ,16 ,15 ,2 ,6] # instance, n, m, L, U
        ]
    for i in range(len(dataset)):
        print("--------------------------------------------------------------")
        #max_time = 10*3600
        max_time = 2000
        input_file = folder_dataset + dataset[i][0]
        instance = dataset[i][0]
        n = dataset[i][1]
        m = dataset[i][2]
        L = dataset[i][3]
        U = dataset[i][4]
        start_time = time.time()
        createGraph(input_file)
        print("instance: " + instance)
        best_sequence = []
        total_runtime = 0
        upper = U
        lower = L
        time_out = False
        while lower <= upper and not time_out:
            mid = math.floor((upper + lower) /2)
            print("MID: " + str(mid))
            run(mid, max_time)
            max_time = max_time - runtime
            if feasible == True:
                upper = sol_size - 1
            else:
                lower = mid + 1
            if max_time <= 0:
                time_out = True
        print("burning number: " + str(len(best_sequence)))
        print("optimal burning sequence: " + str(best_sequence))
        total_runtime = time.time() - start_time
        print("---Total running time: %s seconds ---" % total_runtime)
