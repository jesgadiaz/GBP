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
        
        x = []   #-----------------------------------------------------------------------------------------(13)
        for i in range(n):
            x.append(0)
            x[i] = m.addVar(vtype=GRB.BINARY, name="x%s" % str(i+1))
            
        z = []   #-----------------------------------------------------------------------------------------(14)
        for i in range(n):
            temp = []
            for j in range(k+1):
                temp.append(0)
            z.append(temp)
        for i in range(n):
            for j in range(k+1):
                z[i][j] = m.addVar(vtype=GRB.BINARY, name="z%s" % "_" + str(i+1) + "," + str(j+1))
                
        p = []   #-----------------------------------------------------------------------------------------(15)
        for i in range(n):
            p.append(0)
            p[i] = m.addVar(vtype=GRB.INTEGER, name="p%s" % str(i+1))
                
        y = []   #-----------------------------------------------------------------------------------------(17)
        for i in range(n):
            temp = []
            for j in range(n):
                temp.append(0)
            y.append(temp)
        for i in range(n):
            for j in range(n):
                y[i][j] = m.addVar(vtype=GRB.BINARY, name="y%s" % "_" + str(i+1) + "," + str(j+1))                
        
        #---------------------------- OBJECTIVE FUNCTION --------------------------------------------------
        
        m.setObjective(0, GRB.MINIMIZE)#--------------------------------------------------------------(1)
        
        #---------------------------- CONSTRAINTS ---------------------------------------------------------
        
        m.addConstr(sum(x) == k) #-------------------------------------------------------------------------(25)
        
        for i in range(n): #-------------------------------------------------------------------------------(26)
            m.addConstr(x[i] <= p[i])
            
        for i in range(n): #-------------------------------------------------------------------------------(C4)
            m.addConstr(p[i] <= x[i] * k)
            
        z_transpose = np.array(z).T.tolist()
        for j in range(k): #-------------------------------------------------------------------------------(27)
            m.addConstr(sum(z_transpose[j]) == 1)     
        
        for i in range(n): #-------------------------------------------------------------------------------(28)
            m.addConstr(sum(z[i]) == 1)
            
        #--------------------------------------------------------------------------------------------------(C7)
        #m.addConstr(sum(z_transpose[k]) == n - k)
        
        q = []
        for i in range(k):
            q.append(i+1)
        q_np = np.array(q)
        for i in range(n): #-------------------------------------------------------------------------------(29)
            z_np = np.array(z[i])
            m.addConstr(p[i] == sum(q_np * z_np[:-1]))
            
        #for i in range(n): #-------------------------------------------------------------------------------(C9)
        #    for j in range(n):
        #        m.addConstr((x[i] * d[i][j]) + p[i] <= x[i] * (k + (1-y[i][j]) * M) )
        
        for i in range(n): #-------------------------------------------------------------------------------(31)
            for j in range(n):
                m.addConstr((x[i] * d[i][j]) + p[i] <= x[i] * (k) + (1-y[i][j]) * M )
                
        #for i in range(n): #------------------------------------------------------------------------------(C10)
        #    for j in range(n):
        #        m.addConstr((x[i] * d[i][j]) + p[i] >= x[i] * (k + epsilon - (y[i][j] * M) ) )
        
        for i in range(n): #------------------------------------------------------------------------------(32)
            for j in range(n):
                m.addConstr((x[i] * d[i][j]) + p[i] >= x[i] * (k + epsilon) - (y[i][j] * M))
            
        for i in range(n): #------------------------------------------------------------------------------(30)
            for j in range(n):
                m.addConstr(y[i][j] <= x[i])
                
        y_transpose = np.array(y).T.tolist()
        for j in range(n): #------------------------------------------------------------------------------(33)
            m.addConstr(sum(y_transpose[j]) >= 1)    
                
        #---------------------------- OPTIMIZATION -------------------------------------------------------
        
        m.optimize()
        runtime = m.Runtime
        #print("The run time is %f" % runtime)
        if m.status == GRB.INFEASIBLE:
            #print('Model is infeasible')
            feasible = False
        else:
            #print("Obj:", m.objVal)
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
            #print("p_out: " + str(p_out))
            best_sequence = []
            for j in range(k):
                for i in range(len(p_out)):
                    if p_out[i] <= (j+1) + 0.1 and p_out[i] >= (j+1) - 0.1:
                        best_sequence.append(i+1)
                        break
            #print("c_out: " + str(c_out))
            #print("t_out: " + str(t_out))
            #print("t_out: " + str(t))
            #print("r_out: " + str(r))
       
        
    #except GurobiError:
    #    print("Error reported")
    #    print(str(GurobiError.message))
    #    print(str(GurobiError.errno))
        #presolved_model = m.presolve()
        #print("Presolve model - IsQCP?: " + str(presolved_model.IsQCP))
        #print("Presolve model - IsQP?: " + str(presolved_model.IsQP))
        #presolved_model.printStats()
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
    
    #print("lower: " + str(lower))
    #print("upper: " + str(upper))
    #while lower <= upper:
    checked_feasible_k = []
    #while lower != upper:
    while lower <= upper:
        mid = math.floor((upper + lower) /2)
        #mid = math.ceil((upper + lower) /2)
        #if mid not in checked_feasible_k:
        run(mid)
        print("MID: " + str(mid))
        total_runtime = total_runtime + runtime
        if feasible == True:
            #checked_feasible_k.append(mid)
            upper = mid - 1
            #upper = mid - 1
        else:
            lower = mid + 1
        #else:
        #    lower = upper
    print("burning number: " + str(len(best_sequence)))
    print("optimal burning sequence: " + str(best_sequence))
    total_runtime = time.time() - start_time
    print("---Total running time: %s seconds ---" % total_runtime)
    file_log = open("log.txt","a")
    
    file_log.write(instance[:-4] + " & "    + str(n_nodes) + " & " + str(n_edges) 
                   + " & " + '%.3f'%(density)   + " & " + '%.3f'%(av_degree)
                   #+ " & " + str(conn_comp) + " & " + str(len(best_sequence))
                   + " & " + str(len(best_sequence))
                   + " & " + str(math.ceil(total_runtime)) + "\\\\" + '\n')

    file_dispersion = open("disper.txt","a")
    #file_dispersion.write(str(n_nodes) + "," + str(total_runtime) + '\n')
    file_dispersion.write(str(total_runtime) + '\n')
    
    file_sols = open("sols.txt","a")
    file_sols.write(instance[:-4] + ", " + str(total_runtime) + ", " + str(best_sequence) + '\n')
    
    file_log.close()
    file_sols.close()
    file_dispersion.close()
    
def main(n_in, m_in, input_file_in, L, B):
    global n, m, start_time
    start_time = time.time()
    input_file  = input_file_in
    n = n_in
    m = m_in
    createGraph(input_file)
    binarySearch(L,B)
    
if __name__ == "__main__":
    global instance
    #folder_dataset = 'C:/Users/jgd/Documents/GBP/dataset/nr/'
    #folder_dataset = 'C:/Users/jgd/Documents/GBP/dataset/paths/'
    folder_dataset = 'C:/Users/jgd/Documents/GBP/dataset/paperILP/'
    dataset = [
        #['er100_1.mtx',  100, 233],
        #['er100_2.mtx',  100, 245],
        #['er100_3.mtx',  100, 231],
        #['er100_4.mtx',  100, 241],
        #['er100_5.mtx',  100, 236],
        #['er100_6.mtx',  100, 233],
        #['er100_7.mtx',  100, 242],
        #['er100_8.mtx',  100, 231],
        #['er100_9.mtx',  100, 255],
        #['er100_10.mtx', 100, 215],
        
        #['er100_1_m.mtx', 100, 352],
        #['er100_2_m.mtx', 100, 354],
        #['er100_3_m.mtx', 100, 366],
        #['er100_4_m.mtx', 100, 374],
        #['er100_5_m.mtx', 100, 382],
        #['er100_6_m.mtx', 100, 383],
        #['er100_7_m.mtx', 100, 365],
        #['er100_8_m.mtx', 100, 386],
        #['er100_9_m.mtx', 100, 363],
        #['er100_10_m.mtx', 100, 383],
        
        #['er100_1_ss.mtx', 100, 151],
        #['er100_2_ss.mtx', 100, 148],
        #['er100_3_ss.mtx', 100, 144],
        #['er100_4_ss.mtx', 100, 163],
        #['er100_5_ss.mtx', 100, 161],
        #['er100_6_ss.mtx', 100, 162],
        #['er100_7_ss.mtx', 100, 149],
        #['er100_8_ss.mtx', 100, 152],
        #['er100_9_ss.mtx', 100, 161],
        #['er100_10_ss.mtx', 100, 162],   
        
        ['er100_1_007.mtx', 100, 340],   
        ['er100_2_007.mtx', 100, 328],   
        ['er100_3_007.mtx', 100, 362],   
        ['er100_4_007.mtx', 100, 379],   
        ['er100_5_007.mtx', 100, 332],   
        ['er100_6_007.mtx', 100, 340],   
        ['er100_7_007.mtx', 100, 336],   
        ['er100_8_007.mtx', 100, 315],   
        ['er100_9_007.mtx', 100, 316],   
        ['er100_10_007.mtx', 100, 365],   
        
        
        #['reds100_1.mtx', 100, 183],
        #['reds100_2.mtx', 100, 199],
        #['reds100_3.mtx', 100, 187],
        #['reds100_4.mtx', 100, 186],
        #['reds100_5.mtx', 100, 179],
        #['reds100_6.mtx', 100, 201],
        #['reds100_7.mtx', 100, 200],
        #['reds100_8.mtx', 100, 185],
        #['reds100_9.mtx', 100, 184],
        #['reds100_10.mtx', 100, 195],
        
        
        #['path16.mtx',16,15],
        #['path25.mtx',25,24],
        #['path36.mtx',36,35],
        #['path49.mtx',49,48],
        #['path64.mtx',64,63],
        #['path81.mtx',81,80],
        #['path100.mtx',100,99],
        #['path121.mtx',121,120],
        #['path144.mtx',144,143],
        #['path169.mtx',169,168],
        #['path196.mtx',196,195],
        #['path225.mtx',225,224],
        #['path256.mtx',256,255],
        #['path289.mtx',289,288],
        #['path324.mtx',324,323],
        #['path361.mtx',361,360],        
        
        #['karate.mtx'            , 34  , 78],
        #['chesapeake.mtx'        , 39  , 170],
        #['line49nodes.mtx'       , 49  , 48],
        #['dolphins.mtx'          , 62  , 159],
        #['rt-retweet.mtx'        , 96  , 117],
        #['polbooks.mtx'          , 105 , 441],
        #['adjnoun.mtx'           , 112 , 425],
        #['ia-infect-hyper.mtx'   , 113 , 2196],
        #['ia-enron-only.mtx'     , 143 , 623],
        #['c-fat200-1.mtx'        , 200 , 1534],
        #['c-fat200-2.mtx'        , 200 , 3235],
        #['c-fat200-5.mtx'        , 200 , 8473],
        #['squaredIdealBurn7.mtx' , 231 , 418],
        #['DD244.mtx'             , 291 , 822],
        #['ca-netscience.mtx'     , 379 , 914],
        #['sanr400-0-5.mtx'       , 400 , 39984],
        #['infect-dublin.mtx'     , 410 , 2765],
        #['frb30-15-4.mtx'        , 450 , 83194],
        #['c-fat500-1.mtx'        , 500 , 4459],
        #['c-fat500-2.mtx'        , 500 , 9139],
        #['c-fat500-5.mtx'        , 500 , 23191],
        #['c-fat500-10.mtx'       , 500 , 46627],
        #['bio-diseasome.mtx'     , 516 , 1188],
        #['web-polblogs.mtx'      , 643 , 2280],
        #['p-hat700-1.mtx'        , 700 , 60999],
        #['p-hat700-2.mtx'        , 700 , 121728], 
        #['DD68.mtx'              , 775 , 2093],
        #['BA-1_10_60-L5.mtx'     , 804 , 46410],
        #['ia-crime-moreno.mtx'   , 829 , 1475],
        #['soc-wiki-Vote.mtx'     , 889 , 2914],
        #['socfb-Reed98.mtx'      , 962 , 18812],
        #['ba-1k-2k.mtx'          , 1000, 1996],
        #['lattice3D.mtx'         , 1000, 2700],
        #['san1000.mtx'           , 1000, 250500],
        #['bal_bin_tree_9.mtx'    , 1023, 1022],
        #['delaunay_n10.mtx'      , 1024, 3056],
        #['stufe.mtx'             , 1036, 1868],
        #['lattice2D.mtx'         , 1089, 2112],
        #['bal_ter_tree_6.mtx'    , 1093, 1092],
        #['email-univ.mtx'        , 1133, 5451],
        #['econ-mahindas.mtx'     , 1258, 7513],
        #['ia-fb-messages.mtx'    , 1266, 6451],
        #['bio-yeast.mtx'         , 1458, 1948]
        ]
    for i in range(len(dataset)):
        print("--------------------------------------------------------------")
        instance = dataset[i][0]
        print("instance: " + instance)
        n = dataset[i][1]
        L = 4
        #main(n, dataset[i][2], folder_dataset + dataset[i][0], 1, math.ceil(math.sqrt(n))) # con B = -1 se desactiva
        #for j in range(math.ceil(math.sqrt(n)), n+1):
        #for j in range(6, n+1):
        for j in range(4, n+1):
        #for j in [10]:
            main(n, dataset[i][2], folder_dataset + dataset[i][0], L, j) # con B = -1 se desactiva
            
#            if B == -1:
#                if nx.is_connected(G):
#                   upper = math.ceil(math.sqrt(n))
#                    lower = 0
#                else:
#                    upper = n
#                    lower = nx.number_connected_components(G) - 1
#            else:
#                upper = B
#                if nx.is_connected(G):
#                    lower = 0
#                else:
#                    lower = nx.number_connected_components(G) - 1            
