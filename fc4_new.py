import pandas as pd
import numpy as np
import networkx as nx
import os
import sys
from itertools import permutations, repeat
import itertools
import math
from io import StringIO
from GraphRicciCurvature.FormanRicci import FormanRicci

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', None)
pd.set_option('display.width', None)
np.set_printoptions(threshold=sys.maxsize)
N=1000
thres = 0.4
#print("threshold:")
#print(thres)
sum1 = np.empty((200, 200))
count = 1
intersect = []

#reading the matrices and getting [0,1] matrices from them

for filename in os.listdir('/Users/dasgupta/Desktop/brain2/disease'):
    if filename.endswith('txt'):
        mat1 = np.genfromtxt(open('/Users/dasgupta/Desktop/brain2/disease/'+filename, 'rb'), delimiter=",", dtype='float',
                             filling_values="0.0", autostrip=True, usecols=range(1,201))
        mat1 = np.nan_to_num(mat1)
        sum1 = mat1
        sum1 = np.nan_to_num(sum1)
        #count = count+1

#sum1 = sum1/(count-1)
data_main = np.where(sum1 > thres,sum1, 0)
#print(data_main)
network_matrix1 = nx.from_numpy_matrix(data_main, create_using=None)
edge_list = list(network_matrix1.edges)
#print(edge_list)
data=pd.DataFrame(data_main)

#finding connected components, largest connected component, diameters of connected components, diameter of largest connected component

#largest_cc1 = max(nx.connected_components(network_matrix1), key=len)
#component_size1 = len(largest_cc1)

#list_cc1 = nx.connected_components(network_matrix1)
#diams1= []
#S1 = [network_matrix1.subgraph(c).copy() for c in list_cc1]
#subGraph1 = []
#for c in S1:
 #   subGraph1.append(c.nodes(data=True))
  #  diams1.append(nx.diameter(c))

#function to find common value between 2 lists

def intersection(a:list,b:list):
    return(set(a).intersection(b))

#reading the regions matrices

for filename in os.listdir('/Users/dasgupta/Desktop/brain2/regions'):
    if filename.endswith('txt'):
        control1 = np.genfromtxt(open('/Users/dasgupta/Desktop/brain2/regions/'+filename, 'rb'),delimiter=","
                            ,dtype=str)

data_rc=pd.DataFrame(data = control1,index = None,columns = None)
rc_np = data_rc.to_numpy()
rc_list = list(rc_np)

#function to find area of triangles

def area_c(triangles_nodes:list):
    sum=0.0
    edges=[]
    combo = list(permutations(triangles_nodes, 2))
    for i in range(1,4) : 
           #print(data.iloc[list(combo[i])[0]][list(combo[i])[1]])
           p = (data.iloc[list(combo[i])[0]][list(combo[i])[1]])
           sum = sum + p
           edges.append(data.iloc[list(combo[i])[0]][list(combo[i])[1]])

    t_area = (sum/3)
    #summ = (sum/2)
    #t_area = math.sqrt(abs(summ*(summ-edges[0])*(summ-edges[1])*(summ-edges[2])))
    return t_area    

#function to find edge weights

def weight_edge_c(node:int):
    wt_e = 0
    #for i in range(0,len(S1[0])):
    for i in range(0,len(network_matrix1)):
        wt_e = wt_e + data.iloc[node][i]
    return wt_e


#for future : setting triangles for higher order faces
def triangles(nodes:list):
    combo = list(permutations(nodes, 3))
    return combo



def list_diff(list1, list2):
    out = [item for item in list1 if not item in list2]
    return out

def intersection(lst1, lst2): 
    lst3 = [value for value in lst1 if value in lst2] 
    return lst3    



#driver function and forman curvature for disease




print("DISEASE")
print()
print()

#ls = nx.cycle_basis(S1[0],0)
ls = nx.cycle_basis(network_matrix1,0)
edg=[]
edg_ls=[]
final=[]
#print(ls)
n_A = []
n_A4 = []
n_x = []
n_fcy = []
n_ficy = []
fi_cycle = []
f_cycle = []
i_not_j = []
par_edg = []
for x in ls:
    if len(x)==3:
       	#set_triangles = triangles(x)
       	#for t in set_triangles:
            #S = area(t)
    #else:
       	#print(x)
       	n_x.append(x)
        A = area_c(x)
        #A = math.sqrt(3/16)
        n_A.append(A)
        for items in x:
            comb=list(permutations(x,2)) 
            #n_comb.append(comb)   
        for q in range(1,4):
            edg.append(comb[q])   
        for i in range(1,len(edg)):
            edg_ls.append(edg[i])     
        #print(edg_ls) 
        for num in edg_ls:
            if num not in final:
               final.append(num)
#print(final)
#Finding 4-cycles:
for i in range(0,(len(n_x)-2)): 
    for j in range(i+1,(len(n_x)-1)):
        common = len(intersection(n_x[i],n_x[j]))
        #print(intersection(n_x[i],n_x[j]))
        if common == 2:
            #print("4- cycle found:")
            i_not_j = list_diff(n_x[i],n_x[j])
            j_not_i = list_diff(n_x[j],n_x[i])
            f_cycle = i_not_j + n_x[j]
            #print(f_cycle)
            n_fcy.append(f_cycle)
            A = area_c(n_x[i]) + area_c(n_x[j])
            n_A4.append(A)
#print(n_fcy) 




#finding C(G) for 3 cycles:
#C_G3 = network_matrix1.number_of_nodes() - network_matrix1.number_of_edges() + len(n_x) 
#print(C_G3)

#finding the C(G) for 4 cycles
#C_G4 = network_matrix1.number_of_nodes() - network_matrix1.number_of_edges() + len(n_x) + len(n_fcy)       
#print(C_G4)

#Finding Forman Curvature till 4 cycles:
#terms needed are term1 term1A term2 term3 term4(containing term4A) term5A term5B:
#finding term1A (under for i in range(0,len(final))):

#n_term1 = [] 
n_FC = []

#el = list(nx.generate_edgelist(network_matrix1, data=False))
#NG = nx.parse_edgelist(el, nodetype = int)
#edge_list = NG.edges()
#print(edge_list)
for i in range(0,len(edge_list)):
    #print(edge_list[i])
    wtt_e=data.iloc[edge_list[i][0]][edge_list[i][1]]
    #print(wtt_e)
    term1 = 0
    #term2 = weight_edge_c(final[i][0])/(S1[0].degree(final[i][0]))
    term2 = weight_edge_c(edge_list[i][0])/(network_matrix1.degree(edge_list[i][0]))
    #term3 = weight_edge_c(final[i][1])/(S1[0].degree(final[i][1]))
    term3 = weight_edge_c(edge_list[i][1])/(network_matrix1.degree(edge_list[i][1])) 
    for j in range(0,len(n_x)):
        if edge_list[i][0] in n_x[j] and edge_list[i][1] in n_x[j]:
            term1 = term1 + (wtt_e/n_A[j]) 
            #print (term1)

          
    term1A = 0
    term4 = 0
    term4A = 0
    term5A = 0
    term5B = 0
    term5 = 0 
    for k in range (0,len(n_fcy)):
         

        if edge_list[i][0] in n_fcy[k] and edge_list[i][1] in n_fcy[k]:
            term1A = term1A + (wtt_e/(n_A4[k]))
            #finding term4:
            par_edg = list_diff(n_fcy[k],edge_list[i])
            #print(par_edg)
            wtt_par_edg = data.iloc[par_edg[0]][par_edg[1]]
            #print(wtt_par_edg)
            term4A = math.sqrt((wtt_e)*(wtt_par_edg)/(n_A4[k]))
            term4 = term4 + term4A
            #print(term4)
 
            #finding term5 = ((term2+term3)+(term5A + term5B))/term4A:
            #term5A = term5A + weight_edge_c(par_edg[0])/(S1[0].degree(par_edg[0]))
            #term5B = term5B + weight_edge_c(par_edg[1])/(S1[0].degree(par_edg[1]))
            #term5A = term5A + weight_edge_c(par_edg[0])/(network_matrix1.degree(par_edg[0]))
            #term5B = term5B + weight_edge_c(par_edg[1])/(network_matrix1.degree(par_edg[1]))
            #if term4 == 0:
             #   term5 = 0
            #else:
            #    term5 = term5 + ((term2+term3)+(term5A + term5B))/term4

  #Final FC:
    #FC = wtt_e * (term1 + term1A + ((term2 + term3)/wtt_e)-abs(term4 - term5))
    FC = wtt_e * (term1 + term1A + ((term2 + term3)/wtt_e)-abs(term4))
    
    #print(FC)
    n_FC.append(FC)
    #print(n_FC)

#print(final)

final_node = []

for i in range(0,len(edge_list)):
    a = rc_np[edge_list[i][0]]
    b = rc_np[edge_list[i][1]]
    d = n_FC[i]
    c = [a,b,d]
    #c = [a,b]
    print(c)
    final_node.append(c)
print()
#print(final_node)
print()
print()

