from gurobipy import *
import math
import numpy as np
import time
import networkx as nx
import random

def callback_incumbent_logger(model, where):
    global t_incumbent
    if where == GRB.Callback.MIPSOL:
        t_incumbent = model.cbGet(GRB.Callback.RUNTIME)

def createGraph(input_file):
    global t_incumbent, G, n, m, a0, start_time, apsp_time, n_nodes, n_edges, conn_comp, av_degree, density, card_comp
    
    t_incumbent = math.inf
    
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

def run():
    global G, total_runtime, runtime, n, m, a0, feasible, best_sequence, L_input, T_input, density, av_degree, n_nodes, n_edges
    try:
        m = Model("mip1")
        m.Params.outputFlag = 1  # 0 - Off  //  1 - On
        m.setParam("MIPGap", 0.0);
        m.setParam("Presolve", 2); # -1 - Automatic // 0 - Off // 1 - Conservative // 2 - Aggresive 
        #m.setParam("PreQLinearize", -1); # -1 - Automatic // 0 - Off // 1 - Strong LP relaxation // 2 - Compact relaxation
        #m.params.BestObjStop = k
        m.setParam("TimeLimit", 2000)
        
        # ------------------------------INPUT --------------------------------
        T = T_input
        #---------------------------- VARIABLES -----------------------------------------------------------

        z = []
        for i in range(n):
            temp = []
            for j in range(T):
                temp.append(0)
            z.append(temp)
        for i in range(n):
            for j in range(T):
                z[i][j] = m.addVar(vtype=GRB.BINARY, name="z,%s" % str(i+1) + "," + str(j+1))

        b = []
        for i in range(n):
            temp = []
            for j in range(T):
                temp.append(0)
            b.append(temp)
        for i in range(n):
            for j in range(T):
                b[i][j] = m.addVar(vtype=GRB.BINARY, name="b,%s" % str(i+1) + "," + str(j+1))
        
        b_prime = []
        for i in range(T):
            b_prime.append(0)
            b_prime[i] = m.addVar(vtype=GRB.BINARY, name="b_prime,%s" % str(i+1))

        s = []
        for i in range(n):
            temp = []
            for j in range(T):
                temp.append(0)
            s.append(temp)
        for i in range(n):
            for j in range(T):
                s[i][j] = m.addVar(vtype=GRB.BINARY, name="s,%s" % str(i+1) + "," + str(j+1))
        
        #---------------------------- CONSTRAINTS ---------------------------------------------------------

        b0 = [] 
        s0 = [] 
        for i in range(n):
            b0.append(0)
            s0.append(0)

        for i in range(n): #--------------------------------------------------(1.3)--
            for j in range(T):
                m.addConstr(b[i][j] >= s[i][j])
        
        for i in range(n): #--------------------------------------------------(1.2)--
            for j in range(T):
                if j == 0:
                    m.addConstr(s[i][j] >= s0[i])
                else:
                    m.addConstr(s[i][j] >= s[i][j-1])

        for i in range(n): #--------------------------------------------------(1.4)--
            for j in range(T):
                for k in G.neighbors(i):
                    if j == 0:
                        m.addConstr(b[i][j] >= b0[k])
                    else:
                        m.addConstr(b[i][j] >= b[k][j-1])
                        
        #for i in range(n): #---------------------------------------------------(1.4)--
        #    for j in range(T):
        #        closed_neighbors = [n for n in G.neighbors(i)]
        #        closed_neighbors.append(i)
                #for k in G.neighbors(i):
        #        for k in closed_neighbors:
        #            if j == 0:
        #                m.addConstr(b[i][j] >= s0[k])
        #            else:
        #                m.addConstr(b[i][j] >= s[k][j-1])
        
        for i in range(n): #--------------------------------------------------(1.5)
            for j in range(T):
                sum_ = 0
                closed_neighbors = [n for n in G.neighbors(i)]
                closed_neighbors.append(i)
                for k in closed_neighbors:
                    if j == 0:
                        sum_ = sum_ + b0[k]
                    else:
                        sum_ = sum_ + b[k][j-1]
                sum_ += s[i][j]
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
        
        for j in range(T): #--------------------------------------------------(1.7)--
            sum_ = 0
            for i in range(n):
                sum_ = sum_ + b[i][j]
                m.addConstr(b_prime[j] <= b[i][j])
            m.addConstr(b_prime[j] >= sum_ - (n-1))
                        

        #---------------------------- OBJECTIVE FUNCTION --------------------------------------------------
        
        b_transpose = np.array(b).T.tolist()
        m.setObjective(T - sum(b_prime), GRB.MINIMIZE)#------------------------(1)--
        
        # Set bounds
        m.addConstr(T - sum(b_prime) <= T_input - 1)
        m.addConstr(T - sum(b_prime) >= L_input - 1)
                
        #---------------------------- OPTIMIZATION -------------------------------------------------------
        
        m.optimize(callback_incumbent_logger)
        runtime = m.Runtime
        if m.status == GRB.INFEASIBLE:
            feasible = False
        else:
            print("Obj:", m.objVal)
            print("t_incumbent:" + str(t_incumbent))
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
                if varNameSplit[0] == 'b_prime':
                    b_prime_out.append(v.x)
            
            b_prime = []
            for i in range(len(b_prime_out)):
                if b_prime_out[i] > 0.9:
                    b_prime.append(1)
                else:
                    b_prime.append(0)
            
            sol_size = T - sum(b_prime) + 1
            
            if T - sum(b_prime) == len(b_prime):
                print("No feasible solution")
            else:
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
                print("Final MIP gap value: %f" % m.MIPGap)
    except GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

def main(n_in, m_in, input_file_in, L, U):
    global n, m, start_time, T_input, L_input
    start_time = time.time()
    input_file  = input_file_in
    n = n_in
    m = m_in
    createGraph(input_file)
    T_input = U
    L_input = L
    run()
    
if __name__ == "__main__":
    global instance, T_input
    folder_dataset = 'C:/Users/perro/Documents/GBP/paperILP/'
    #folder_dataset = 'C:/Users/perro/Documents/GBP/erdos/'
    dataset = [
        ['ca-netscience.mtx',379,914,2,10], # instance, n, m, L, U
        ['karate.mtx',34,78,1,6], # instance, n, m, L, U
        ]
    for i in range(len(dataset)):
        print("--------------------------------------------------------------")
        instance = dataset[i][0]
        print("instance: " + instance)
        n = dataset[i][1]
        m = dataset[i][2]
        L = dataset[i][3]
        U = dataset[i][4]
        main(n, m, folder_dataset + dataset[i][0], L, U)
