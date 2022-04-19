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
    global G, total_runtime, runtime, n, m, a0, feasible, best_sequence, L_input, T_input, density, av_degree, n_nodes, n_edges
    try:
        m = Model("mip1")
        m.Params.outputFlag = 1  # 0 - Off  //  1 - On
        m.setParam("MIPGap", 0.0);
        m.setParam("Presolve", 2); # -1 - Automatic // 0 - Off // 1 - Conservative // 2 - Aggresive 
        #m.setParam("PreQLinearize", -1); # -1 - Automatic // 0 - Off // 1 - Strong LP relaxation // 2 - Compact relaxation
        #m.params.BestObjStop = k
        m.setParam("TimeLimit", 2*3600)
        
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

        b0 = []
        s0 = []
        for i in range(n):
            b0.append(0)
            s0.append(0)
        
        for i in range(n): #-------------------------------(6)--
            for j in range(T):
                if j == 0:
                    m.addConstr(s[i][j] >= s0[i])
                else:
                    m.addConstr(s[i][j] >= s[i][j-1])                    

        for i in range(n): #---------------------------------------------------(7)--
            for j in range(T):
                for k in G.neighbors(i):
                    if j == 0:
                        m.addConstr(b[i][j] >= b0[k])
                    else:
                        m.addConstr(b[i][j] >= b[k][j-1])
                        
        for i in range(n): #---------------------------------------------------(8)--
            for j in range(T):
                for k in G.neighbors(i):
                    if j == 0:
                        m.addConstr(b[i][j] >= s0[k])
                    else:
                        m.addConstr(b[i][j] >= s[k][j-1])
                    
        for i in range(n): #--------------------------------------------(9)
            for j in range(T):
                sum_ = 0
                for k in G.neighbors(i):
                    if j == 0:
                        sum_ = sum_ + b0[k]
                    else:
                        sum_ = sum_ + b[k][j-1]                        
                for k in G.neighbors(i):
                    if j == 0:
                        sum_ = sum_ + s0[k]
                    else:
                        sum_ = sum_ + s[k][j-1]
                m.addConstr(b[i][j] <= sum_)

        #sum_ = 0
        #for i in range(n): # -------------------------------(12)--
        #    sum_ = sum_ + b[i][0]
        #m.addConstr(sum_ == 0)

        s_transpose = np.array(s).T.tolist()# --------------(12)--
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
                
        for i in range(n): #--------------------------------(13)--
            for j in range(T):
                m.addConstr(z[i][j] >= b[i][j])
                m.addConstr(z[i][j] >= s[i][j])
                m.addConstr(z[i][j] <= b[i][j] + s[i][j])
        
        for j in range(T): #--------------------------------(13)--
            sum_ = 0
            for i in range(n):
                sum_ = sum_ + z[i][j]
                m.addConstr(b_prime[j] <= z[i][j])
            m.addConstr(b_prime[j] >= sum_ - (n-1))
                        

        #---------------------------- OBJECTIVE FUNCTION --------------------------------------------------
        
        b_transpose = np.array(b).T.tolist()
        m.setObjective(T - sum(b_prime), GRB.MINIMIZE)#------------------------(1)--
        
        # Set bounds
        m.addConstr(T - sum(b_prime) <= T_input - 1)
        m.addConstr(T - sum(b_prime) >= L_input - 1)
                
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
    
                
                file_log = open("log.txt","a")
                
                file_log.write(instance[:-4] + " & "    + str(n_nodes) + " & " + str(n_edges - n)
                               + " & " + '%.3f'%(density)   + " & " + '%.3f'%(av_degree)
                               + " & " + str(round(runtime,3))
                               + " & " + str(m.MIPGap)
                               + " & " + str(sol_size)
                               + "\\\\" + '\n')
            
                file_dispersion = open("disper.txt","a")
                #file_dispersion.write(str(n) + "," + str(runtime) + '\n')
                file_dispersion.write(str(runtime) + '\n')
                
                file_sols = open("sols.txt","a")
                file_sols.write(instance[:-4] + ", " + str(runtime) + ", " + str(sol) + '\n')
                
                file_log.close()
                file_sols.close()
                file_dispersion.close()
    except GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

def main(n_in, m_in, input_file_in, L, B):
    global n, m, start_time, T_input, L_input
    start_time = time.time()
    input_file  = input_file_in
    n = n_in
    m = m_in
    createGraph(input_file)
    T_input = B #math.ceil(math.sqrt(n)) + B
    L_input = L
    run()
    
if __name__ == "__main__":
    global instance, T_input
    #folder_dataset = 'C:/Users/jgd/Documents/FireFighter/nr/'
    folder_dataset = 'C:/Users/jgd/Documents/GBP/dataset/paths/'
    #folder_dataset = 'C:/Users/jgd/Documents/GBP/dataset/paperILP/'
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
        
        #['er100_1_007.mtx', 100, 340],   
        #['er100_2_007.mtx', 100, 328],   
        #['er100_3_007.mtx', 100, 362],   
        #['er100_4_007.mtx', 100, 379],   
        #['er100_5_007.mtx', 100, 332],   
        #['er100_6_007.mtx', 100, 340],   
        #['er100_7_007.mtx', 100, 336],   
        #['er100_8_007.mtx', 100, 315],   
        #['er100_9_007.mtx', 100, 316],   
        #['er100_10_007.mtx', 100, 365],           
        
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
        
        #['ba100_1.mtx', 100, 99],
        #['ba100_2.mtx', 100, 99],
        #['ba100_3.mtx', 100, 99],
        #['ba100_4.mtx', 100, 99],
        #['ba100_5.mtx', 100, 99],
        #['ba100_6.mtx', 100, 99],
        #['ba100_7.mtx', 100, 99],
        #['ba100_8.mtx', 100, 99],
        #['ba100_9.mtx', 100, 99],
        #['ba100_10.mtx', 100, 99],        
        
        ['path16.mtx',16,15],
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
        
        #['path16',16,15],
        #['path25',25,24],
        #['path36',36,35],
        #['path49',49,48],
        #['path64',64,63],
        #['path81',81,80],
        #['path100',100,99],
        #['path121',121,120],
        #['path144',144,143],
        #['path169',169,168],
        #['path196',196,195],
        #['path225',225,224],
        #['path256',256,255],
        #['path289',289,288],
        #['path324',324,323],
        #['path361',361,360],
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
        #['frb30-15-4.mtx'        , 450 , 83194], # No feasible after 2 hours
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
        #['karate.mtx', 34, 78, 1]
        #['lattice10x10.mtx', 100, 180]
        #["path4.mtx", 4, 3]
        #["path4.mtx", 4, 3, 1] # file, n, m, number of fire sources
        #["path9.mtx", 9, 8, 1] # file, n, m, number of fire sources
        #["path25.mtx", 25, 24, 2] # file, n, m, number of fire sources
        #["path16.mtx", 16, 15, 4] # file, n, m, number of fire sources
        #["path36.mtx", 36, 35, 9] # file, n, m, number of fire sources
        #["line49nodes.mtx", 49, 48, 1]
        ]
    for i in range(len(dataset)):
        print("--------------------------------------------------------------")
        instance = dataset[i][0]
        print("instance: " + instance)
        n = dataset[i][1]
        L = 1
        #main(n, dataset[i][2], folder_dataset + dataset[i][0], math.ceil(math.sqrt(n)))
        #for j in range(math.ceil(math.sqrt(n)), n+1):
        #for j in range(4, n+1):
        for j in [16]:
            main(n, dataset[i][2], folder_dataset + dataset[i][0], L, j)
        
