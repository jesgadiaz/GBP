from gurobipy import *
import math
import numpy as np
import time
import networkx as nx
import random

def createGraph(input_file):
    global G, n, m, d, start_time, apsp_time, n_nodes, n_edges, conn_comp, av_degree, density
    
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
    
    
    if nx.is_connected(G):
        sp = dict(nx.all_pairs_shortest_path_length(G))
        d = []
        for i in range(0,n):
            list = []
            for j in range(0,n):
                list.append(sp[i][j])
            d.append(list)
        for i in range(0, n):
            for j in range(0, n):
                if d[i][j] == float("inf"):
                    d[i][j] = (n+1);
        del sp
    else:
        d = []
        for i in range(0,n):
            list = []
            for j in range(0,n):
                list.append(float("inf"))
            d.append(list)
        connected_components = nx.connected_components(G)
        for component in connected_components:
            # create subgraph    
            G_sub = nx.Graph()
            for v in component:
                G_sub.add_node(v)
                for e in G.edges:
                    if e[0]==v or e[1]==v:
                        G_sub.add_edge(e[0],e[1])
            
            sp = dict(nx.all_pairs_shortest_path_length(G_sub))
            for item1 in sp.items():
                for item2 in item1[1].items():
                    d[item1[0]][item2[0]] = item2[1]
                    d[item2[0]][item1[0]] = item2[1]
            del G_sub
            del sp
            
        for i in range(0, n):
            for j in range(0, n):
                if d[i][j] == float("inf"):
                    d[i][j] = (n+1)
                
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
    
def run(mid):
    global total_runtime, runtime, n, m, d, feasible, best_sequence
    try:
        m = Model("mip1")
        m.Params.outputFlag = 0  # 0 - Off  //  1 - On
        m.setParam("MIPGap", 0.0);
        m.setParam("Presolve", 2); # -1 - Automatic // 0 - Off // 1 - Conservative // 2 - Aggresive 
        #m.setParam("PreQLinearize", -1); # -1 - Automatic // 0 - Off // 1 - Strong LP relaxation // 2 - Compact relaxation
        #m.params.BestObjStop = k
        #m.setParam("TimeLimit", 2*3600)
        k = mid
        M       =  n + 2 #(3*n)
        epsilon =  0.001
        
        #---------------------------- VARIABLES -----------------------------------------------------------
        
        x = [] 
        for i in range(n):
            x.append(0)
            x[i] = m.addVar(vtype=GRB.BINARY, name="x%s" % str(i+1))
            
        z = []
        for i in range(n):
            temp = []
            for j in range(k+1):
                temp.append(0)
            z.append(temp)
        for i in range(n):
            for j in range(k+1):
                z[i][j] = m.addVar(vtype=GRB.BINARY, name="z%s" % "_" + str(i+1) + "," + str(j+1))
                
        p = []
        for i in range(n):
            p.append(0)
            p[i] = m.addVar(vtype=GRB.INTEGER, name="p%s" % str(i+1))
                
        y = []
        for i in range(n):
            temp = []
            for j in range(n):
                temp.append(0)
            y.append(temp)
        for i in range(n):
            for j in range(n):
                y[i][j] = m.addVar(vtype=GRB.BINARY, name="y%s" % "_" + str(i+1) + "," + str(j+1))                
        
        #---------------------------- OBJECTIVE FUNCTION --------------------------------------------------
        
        m.setObjective(0, GRB.MINIMIZE)#--------------------------------------( )
        
        #---------------------------- CONSTRAINTS ---------------------------------------------------------
        
        m.addConstr(sum(x) == k) #--------------------------------------------(2.2)
        
        for i in range(n): #--------------------------------------------------(2.3)
            m.addConstr(x[i] <= p[i])
            
        for i in range(n): #--------------------------------------------------(2.4)
            m.addConstr(p[i] <= x[i] * k)
            
        z_transpose = np.array(z).T.tolist()
        for j in range(k): #--------------------------------------------------(2.5)
            m.addConstr(sum(z_transpose[j]) == 1)     
        
        for i in range(n): #--------------------------------------------------(2.6)
            m.addConstr(sum(z[i]) == 1)
        
        q = []
        for i in range(k):
            q.append(i+1)
        q_np = np.array(q)
        for i in range(n): #--------------------------------------------------(2.7)
            z_np = np.array(z[i])
            m.addConstr(p[i] == sum(q_np * z_np[:-1]))
            
        for i in range(n): #--------------------------------------------------(2.8)
            for j in range(n):
                m.addConstr(y[i][j] <= x[i])        
                
        for i in range(n): #--------------------------------------------------(2.9)
            for j in range(n):
                m.addConstr((x[i] * d[i][j]) + p[i] <= x[i] * (k) + (1-y[i][j]) * M )
        
        for i in range(n): #--------------------------------------------------(2.10)
            for j in range(n):
                m.addConstr((x[i] * d[i][j]) + p[i] >= x[i] * (k + epsilon) - (y[i][j] * M))
                
        y_transpose = np.array(y).T.tolist()
        for j in range(n): #--------------------------------------------------(2.11)
            m.addConstr(sum(y_transpose[j]) >= 1)    
                
        #---------------------------- OPTIMIZATION -------------------------------------------------------
        
        m.optimize()
        runtime = m.Runtime
        if m.status == GRB.INFEASIBLE:
            feasible = False
        else:
            feasible = True
            p_out = []
            x_out = []
            y_out = []
            for v in m.getVars():
                varName = v.varName
                if varName[0] == 'p':
                    p_out.append(v.x)
                if varName[0] == 'x':
                    x_out.append(v.x)
                if varName[0] == 'y':
                    y_out.append(v.x)
            best_sequence = []
            for j in range(k):
                for i in range(len(p_out)):
                    if p_out[i] <= (j+1) + 0.1 and p_out[i] >= (j+1) - 0.1:
                        best_sequence.append(i)
                        break
    except GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

def binarySearch(L, B):
    global total_runtime, n, runtime, feasible, mid, best_sequence, G, n_nodes, n_edges, conn_comp, instance, av_degree, density, plus
    best_sequence = []
    total_runtime = 0
    
    if B == -1:
        if nx.is_connected(G):
            upper = math.ceil(math.sqrt(n))
            lower = 0
        else:
            upper = n
            lower = nx.number_connected_components(G) - 1
    else:
        upper = B
        lower = L        
    
    checked_feasible_k = []
    while lower <= upper:
        mid = math.floor((upper + lower) /2)
        run(mid)
        print("MID: " + str(mid))
        total_runtime = total_runtime + runtime
        if feasible == True:
            upper = mid - 1
        else:
            lower = mid + 1
    print("burning number: " + str(len(best_sequence)))
    print("optimal burning sequence: " + str(best_sequence))
    total_runtime = time.time() - start_time
    print("---Total running time: %s seconds ---" % total_runtime)
    
def main(n_in, m_in, input_file_in, L, U):
    global n, m, start_time
    start_time = time.time()
    input_file  = input_file_in
    n = n_in
    m = m_in
    createGraph(input_file)
    binarySearch(L,U)
    
if __name__ == "__main__":
    global instance
    folder_dataset = 'C:/Users/jgd/Documents/GBP/dataset/paths/'
    dataset = [
        ['path16.mtx',16,15,2,6], # instance, n, m, L, U
        ['path25.mtx',25,24,3,7], # instance, n, m, L, U  
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
