import preprocessing as pre
from input import create_circuit
import numpy as np
import sys
from test2 import *
import time


if __name__ == "__main__":
    file_name=sys.argv[1]
    circuit_info = open(file_name, "r").read()
    circuit=create_circuit("Test",circuit_info)
    two_qubit_gates=pre.get_two_qubit_gates(circuit)
    n=8
    physical_architecture=pre.get_nnGrid(n)
#    physical_architecture=[(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,0),(3,8),(4,15),(8,9),(9,10),(10,11),(11,12),(12,13),(13,14),(14,15),(15,8)]
#    physical_architecture=[(0, 6), (1, 6), (1, 7), (2, 7), (2, 8), (3, 8), (3, 9), (4, 9), (4, 10), (5, 10), (5, 11),
#                                    (6, 12), (6, 13), (7, 13), (7, 14), (8, 14), (8, 15), (9, 15), (9, 16), (10, 16), (10, 17), (11, 17),
#                                    (12, 18), (13, 18), (13, 19), (14, 19), (14, 20), (15, 20), (15, 21), (16, 21), (16, 22), (17, 22), (17, 23),
#                                    (18, 24), (18, 25), (19, 25), (19, 26), (20, 26), (20, 27), (21, 27), (21, 28), (22, 28), (22, 29), (23, 29),
#                                    (24, 30), (25, 30), (25, 31), (26, 31), (26, 32), (27, 32), (27, 33), (28, 33), (28, 34), (29, 34), (29, 35),
#                                    (30, 36), (30, 37), (31, 37), (31, 38), (32, 38), (32, 39), (33, 39), (33, 40), (34, 40), (34, 41), (35, 41),
#                                    (36, 42), (37, 42), (37, 43), (38, 43), (38, 44), (39, 44), (39, 45), (40, 45), (40, 46), (41, 46), (41, 47),
#                                    (42, 48), (42, 49), (43, 49), (43, 50), (44, 50), (44, 51), (45, 51), (45, 52), (46, 52), (46, 53), (47, 53)]
#    number_of_physical_qubit=54
#    number_of_physical_qubit=np.array(physical_architecture).max()+1
#    neighbours=pre.get_list_of_neighbour(physical_architecture,number_of_physical_qubit)
#    distance_matrix=pre.construct_device_distance_matrix(neighbours,number_of_physical_qubit)
    
#    test_optimization(physical_architecture,two_qubit_gates,dependency,10,600,60)
#    print("Problem solved")
    n,gates,mappings,swaps=find_mappings(physical_architecture,two_qubit_gates)
    print("Valid solution:",verifier(physical_architecture,two_qubit_gates,gates,mappings,swaps))
    print("Number of SWAPs:",n)
