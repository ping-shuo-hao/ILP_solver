from collections import deque

def get_two_qubit_gates(circuit):
  gates=[]
  for gate in circuit.list_gate_program_qubit:
    if len(gate)==2:
      gates.append(gate)
  return gates

def get_nnGrid(n: int):
    my_coupling = []
    for i in range(n):
        for j in range(n):
            if j < n-1:
                my_coupling.append((i*n +j,i*n +j+1))
            if i < n-1:
                my_coupling.append((i*n +j,i*n +j+n))
    return my_coupling


def construct_device_distance_matrix(list_qubit_neighbors,n):
  qubit_distance_matrix = [[n for j in range(n)] for i in range(n)]
  for i in range(n):
    qubit_distance_matrix[i][i] = 0
    traverse_set = set()
    queue = deque()
    queue.append(i)
    while len(queue) > 0:
      q = queue.popleft()
      traverse_set.add(q)
      for j in list_qubit_neighbors[q]:
        if not (j in traverse_set):
          traverse_set.add(j)
          qubit_distance_matrix[i][j] = qubit_distance_matrix[i][q] + 1
          queue.append(j)
  return qubit_distance_matrix

def get_list_of_neighbour(qubit_pairs,n):
  list_qubit_neighbors = [[] for i in range(n)]
  for qubit_pair in qubit_pairs:
    list_qubit_neighbors[qubit_pair[0]].append(qubit_pair[1])
    list_qubit_neighbors[qubit_pair[1]].append(qubit_pair[0])
  return list_qubit_neighbors

