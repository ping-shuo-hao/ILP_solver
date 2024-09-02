from gurobipy import *
import numpy as np
from collections import deque
import time
import copy
import sys


def get_dependecy_list(gate_connectivity_list,number_of_physical_qubit,non_execute_list):
  num_of_gates=len(gate_connectivity_list)
  qubit_occupation=[-1 for _ in range(number_of_physical_qubit)]
  previous_gates=[]

  for i in non_execute_list:
    q0,q1=gate_connectivity_list[i]
    if qubit_occupation[q0]!=-1:
      if not (qubit_occupation[q0],i) in previous_gates:
        previous_gates.append((qubit_occupation[q0],i))
    qubit_occupation[q0]=i

    if qubit_occupation[q1]!=-1:
      if not (qubit_occupation[q1],i) in previous_gates:
        previous_gates.append((qubit_occupation[q1],i))
    qubit_occupation[q1]=i

  return previous_gates

def get_edge_list(model:Model()):
  edge_list=[[] for i in range(model._number_of_physical_qubit)]
  for i,j in model._qubit_connectivity_list:
    edge_list[i].append(j)
    edge_list[j].append(i)
  model._edge_list=edge_list

def construct_device_distance_matrix(model:Model()):
  n=model._number_of_physical_qubit
  list_qubit_neighbors=model._edge_list
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
  model._distance_matrix=qubit_distance_matrix

def get_related_qubits(model:Model(),gate_layers):
  list_of_related_qubit=[[] for _ in range(model._number_of_physical_qubit)]
  related_pairs=[]
  layer_index=0
  for layer in gate_layers:
    for q0,q1 in layer:
      list_of_related_qubit[q0].append((q1,layer_index))
      list_of_related_qubit[q1].append((q0,layer_index))
    layer_index+=1


  for related_qubits in list_of_related_qubit:
    for i in range(len(related_qubits)):
      for j in range(i+1,len(related_qubits)):
        i_qubit,i_layer=related_qubits[i]
        j_qubit,j_layer=related_qubits[j]
        if i_qubit!=j_qubit:
          related_pairs.append((i_qubit,j_qubit,abs(i_layer-j_layer),(i_layer+j_layer)/2))

  return related_pairs

def get_next_executable_gate(distance_matrix,gate_layers,gate_connectivity_list,prev_mapping=[]):
  print("Current gate layer:",gate_layers)
  if len(gate_layers)==0:
    return None

  if len(gate_layers[0])==0:
    return None

  if len(prev_mapping)==0:
    return gate_layers[0][0]

  candidate_list=[]
  for gate in gate_layers[0]:
    q0,q1=gate_connectivity_list[gate]
    p0,p1=prev_mapping[q0],prev_mapping[q1]
    candidate_list.append((gate,distance_matrix[p0][p1]))

  candidate_list.sort(key=lambda x:x[1])

  return candidate_list[0][0]

def get_gate_layers(gate_dependecy_list,gate_connectivity_list,execute_gate,non_executed_list):
  if execute_gate!=-1:
    new_dependecy_list=[]
    for g1,g2 in gate_dependecy_list:
      if g1 != execute_gate:
        new_dependecy_list.append((g1,g2))
    gate_dependecy_list=new_dependecy_list
  else:
    new_dependecy_list=gate_dependecy_list

  record={}
  max_layer=0
  for gate in non_executed_list:
    prev_list=[]
    for prev,later in new_dependecy_list:
      if later==gate:
        prev_list.append(prev)
    current_layer=0
    for prev in prev_list:
      if current_layer<record[prev]+1:
        current_layer=record[prev]+1
    if current_layer>max_layer:
      max_layer=current_layer
    record[gate]=current_layer

  gate_layer=[[] for i in range(max_layer+1)]
  gate_layer_connectivity=[[] for i in range(max_layer+1)]

  for i in non_executed_list:
    gate_layer[record[i]].append(i)
    gate_layer_connectivity[record[i]].append(gate_connectivity_list[i])

  return gate_layer,new_dependecy_list

def create_qubit_mapping_variables_and_constraints(model:Model()):
  model._qubit_mapping_variable_list=model.addVars(np.arange(model._time),
                                                   np.arange(model._number_of_physical_qubit),
                                                   np.arange(model._number_of_physical_qubit),name='x',vtype=GRB.BINARY)
  model.addConstrs(model._qubit_mapping_variable_list.sum(t,i,'*')==1 for t in range(model._time) for i in range(model._number_of_physical_qubit))
  model.addConstrs(model._qubit_mapping_variable_list.sum(t,'*',i)==1 for t in range(model._time) for i in range(model._number_of_physical_qubit))
  model.update()

def find_program_to_physical_mapping(model:Model(),qubit_mapping,t,final_result=False):
  if final_result:
    current_mapping=[]
    for program_q in range(model._number_of_physical_qubit):
      for physical_q in range(model._number_of_physical_qubit):
        if qubit_mapping[t,program_q,physical_q].x>0.9:
          current_mapping.append(physical_q)
          break
    return current_mapping
  else:
    current_mapping=[]
    for program_q in range(model._number_of_physical_qubit):
      for physical_q in range(model._number_of_physical_qubit):
        if qubit_mapping[t,program_q,physical_q]>0.9:
          current_mapping.append(physical_q)
    return current_mapping

  return locations

def find_physical_to_program_mapping(model:Model(),qubit_mapping,t,final_result=False):
  if final_result:
    current_mapping=[]
    for physical_q in range(model._number_of_physical_qubit):
      for program_q in range(model._number_of_physical_qubit):
        if qubit_mapping[t,program_q,physical_q].x>0.5:
          current_mapping.append(program_q)
          break
    return current_mapping
  else:
    current_mapping=[]
    for physical_q in range(model._number_of_physical_qubit):
      for program_q in range(model._number_of_physical_qubit):
        if qubit_mapping[t,program_q,physical_q]>0.5:
          current_mapping.append(program_q)
    return current_mapping

  return locations

def mapper_initilization(qubit_connectivity_list:list,gate_connectivity_list):
  model=Model()
  model._qubit_connectivity_list=qubit_connectivity_list
  model._gate_connectivity_list=gate_connectivity_list
  model._number_of_physical_qubit=np.max(qubit_connectivity_list)+1
  model._time=1
  create_qubit_mapping_variables_and_constraints(model)
  get_edge_list(model)
  construct_device_distance_matrix(model)
  model._non_connect_list=[]
  model._is_initial_mapping=False
  model._mapping_candidate=[]


  for i in range(model._number_of_physical_qubit):
    non_connecting_qubit=[]
    for j in range(model._number_of_physical_qubit):
     if i!=j:
        if not (j in model._edge_list[i]):
          non_connecting_qubit.append(j)
    model._non_connect_list.append(non_connecting_qubit)

  model._previous_distance_constraints=None
  model._current_gate_distance_constraint=None
  model._max_runtime=0
  model._start_time=None
  model._previous_solution=None
  model.setParam('OutputFlag', 0)

  return model

def create_gate_distance_constraints(model:Model(),excution_connectivity):
  model._current_gate_distance_constraint=model.addConstrs(model._qubit_mapping_variable_list[0,q0,p0]
                                                           +model._qubit_mapping_variable_list[0,q1,p1]<=1
                                                           for q0,q1 in excution_connectivity
                                                           for p0 in range(model._number_of_physical_qubit)
                                                           for p1 in model._non_connect_list[p0])
  model.update()


def set_objective_for_mapper(model:Model(),execution_list,gate_layers,prev_location=[],a=10,b=0.9,c=0.9,d=0.5):
  objective=0
  a=a*((len(gate_layers)+1)+model._number_of_physical_qubit)
  related_pairs=get_related_qubits(model,[execution_list]+gate_layers)
  if len(prev_location)>0:
#    print(prev_location)
    objective=quicksum(model._distance_matrix[prev_location[program_qubit]][physical_qubit]*model._qubit_mapping_variable_list[0,program_qubit,physical_qubit]
                       for program_qubit in range(model._number_of_physical_qubit) for physical_qubit in range(model._number_of_physical_qubit))
#    objective=objective-0.5*quicksum(1-model._qubit_mapping_variable_list[0,p,prev_location[p]] for p in range(model._number_of_physical_qubit))
    objective=a*objective

  for t in range(model._time):
    current_constant=len(gate_layers)
    for layer in gate_layers:
      for q0,q1 in layer:
        objective=objective+b*current_constant*quicksum((model._distance_matrix[p0][p1]-1) *
                                                        model._qubit_mapping_variable_list[t,q0,p0]*model._qubit_mapping_variable_list[t,q1,p1]
                                                        for p0 in range(model._number_of_physical_qubit) for p1 in model._non_connect_list[p0])
        current_constant*=b

  constant=1
  objective=objective+quicksum(constant*c**diff*d**aver*(model._distance_matrix[p0][p1]-2)*model._qubit_mapping_variable_list[t,q0,p0]*model._qubit_mapping_variable_list[t,q1,p1]
                                 for q0,q1,diff,aver in related_pairs for p0 in range(model._number_of_physical_qubit)
                                 for p1 in model._non_connect_list[p0] if model._distance_matrix[p0][p1]>1)

  model.setObjective(objective,GRB.MINIMIZE)
  model.update()

def mapping_callback(model,where):
  if model._start_time is None:
    model._start_time=time.time()
  elif time.time()-model._start_time>model._max_runtime and not(model._previous_solution is None):
      print("Mapping callback: Reach time limit")
      model._start_time=None
      model._previous_solution=None
      model.terminate()

  if where==GRB.Callback.MIPSOL:
    obj_value = model.cbGet(GRB.Callback.MIPSOL_OBJ)
    solution = model.cbGetSolution(model._qubit_mapping_variable_list)
    if model._previous_solution is None:
      print("Mapping callback: Find initial mapping with objective value",obj_value)
      model._previous_solution=obj_value
    elif abs(obj_value-model._previous_solution)>model._improvement_limit:
      print("Mapping callback: Find better mapping with objective value",obj_value," and improvement",abs(obj_value-model._previous_solution))
      model._start_time=time.time()
      model._previous_solution=obj_value
      candidate_mapping=find_program_to_physical_mapping(model,solution,0)
      model._mapping_candidate.append((obj_value,candidate_mapping))


def mapper_optimization(mapper:Model(),execution_list,future_list,prev_location=[],run_time_limit=60,
                        improvement_limit=1,number_of_candidate=3):
  mapper._max_runtime=run_time_limit
  mapper._improvement_limit=improvement_limit
  execution_list=[mapper._gate_connectivity_list[i] for i in execution_list]
  future_list=[[mapper._gate_connectivity_list[i] for i in layer]for layer in future_list]
  create_gate_distance_constraints(mapper,execution_list)
  set_objective_for_mapper(mapper,execution_list,future_list,prev_location)

  if len(prev_location)>0:
    for i in range(len(prev_location)):
      mapper._qubit_mapping_variable_list[0,i,prev_location[i]].start=1 #program_to_physical

  mapper.optimize(mapping_callback)
  mappings=mapper._mapping_candidate
  last_obj_value = mapper.ObjVal
  last_program_to_physical=find_program_to_physical_mapping(mapper,mapper._qubit_mapping_variable_list,0,True)
  mappings.append((last_obj_value,last_program_to_physical))
  mappings.sort(key=lambda x:x[0])
  number_of_candidate=min(number_of_candidate,len(mappings))
  mappings=mappings[:number_of_candidate]
  mapper.remove(mapper._current_gate_distance_constraint)
  mapper._mapping_candidate=[]
  mapper._start_time=None
  mapper._previous_solution=None
  mapper.update()
  mapper.reset()

  return mappings

def find_initial_mapping(mapper:Model(),qubit_connectivity_list,gate_connectivity_list,gate_dependecy_list,max_runtime=90):
  number_of_gates=len(gate_connectivity_list)
  non_execute_list=[i for i in range(number_of_gates)]
  gate_layer,gate_dependecy_list=get_gate_layers(gate_dependecy_list,gate_connectivity_list,-1,non_execute_list)
#  mapper.setParam('OutputFlag', 0)
  execution_gate=get_next_executable_gate(mapper._distance_matrix,gate_layer,gate_connectivity_list)
  non_execute_list.remove(execution_gate)
  gate_layer,gate_dependecy_list=get_gate_layers(gate_dependecy_list,gate_connectivity_list,execution_gate,non_execute_list)
  initial_mappings=mapper_optimization(mapper,[execution_gate],gate_layer,run_time_limit=max_runtime)

  return initial_mappings,non_execute_list,[execution_gate],gate_layer,gate_dependecy_list

def check_mapping(qubit_connectivity_list,gate_connectivity_list,gate_dependecy_list,program_to_physical,
                  execution_list,gate_layer,non_execute_list):
  i=0
  terminate=False
  execution_list=copy.deepcopy(execution_list)
  non_execute_list=copy.deepcopy(non_execute_list)
  while len(gate_layer)>0 and (not terminate):
    terminate=True
    if len(gate_layer[0])==0:
      gate_layer.pop(0)
      continue

    for gate in gate_layer[0]:
      q0,q1= gate_connectivity_list[gate]
      q0,q1=program_to_physical[q0],program_to_physical[q1]
      if (q0,q1) in qubit_connectivity_list or (q1,q0) in qubit_connectivity_list:
        terminate=False
        non_execute_list.remove(gate)
        if gate not in execution_list:
          execution_list.append(gate)


    gate_dependecy_list=[(g0,g1) for g0,g1 in gate_dependecy_list if g0 in non_execute_list]
    gate_layer,gate_dependecy_list=get_gate_layers(gate_dependecy_list,gate_connectivity_list,-1,non_execute_list)



  return execution_list,non_execute_list,gate_layer,gate_dependecy_list

def fix_location(router:Model(),previous_mapping:list,next_mapping:list):
#  print(next_mapping)
  previous=router.addConstrs(router._qubit_mapping_variable_list[0,program_qubit,previous_mapping[program_qubit]]==1
                             for program_qubit in range(router._number_of_physical_qubit))
  next=router.addConstrs(router._qubit_mapping_variable_list[1,program_qubit,next_mapping[program_qubit]]==1
                             for program_qubit in range(router._number_of_physical_qubit))

  return previous,next

def create_swap_variables_and_constraints(model:Model()):
  model._routing_variable_list=model.addVars(np.arange(model._time-1),np.arange(model._number_of_physical_qubit),model._expand_connectivity_list,
                                             name="w",vtype=GRB.BINARY)

  model.addConstrs(model._qubit_mapping_variable_list[t+1,p,q]==model._qubit_mapping_variable_list[t,p,q]+
                   model._routing_variable_list.sum(t,p,"*",q)-model._routing_variable_list.sum(t,p,q,"*")
                   for t in range(model._time-1) for p in range(model._number_of_physical_qubit)
                   for q in range(model._number_of_physical_qubit))

  model.addConstrs(model._routing_variable_list.sum(t,"*",physical1,physical2)==model._routing_variable_list.sum(t,"*",physical2,physical1)
  for t in range(model._time-1) for physical1 in range(model._number_of_physical_qubit) for physical2 in model._edge_list[physical1])

  model.addConstrs(model._routing_variable_list.sum(t,program_q,physical_q,"*")<=1 for t in range(model._time-1)
  for program_q in range(model._number_of_physical_qubit) for physical_q in range(model._number_of_physical_qubit))

  model.update()

def create_ordering_variable_and_constraints(model:Model()):
  n=model._number_of_physical_qubit
  model._ordering_variable_list=model.addVars(np.arange(model._time-1),np.arange(model._number_of_physical_qubit),
                                              np.arange(model._number_of_physical_qubit),name="o",vtype=GRB.CONTINUOUS)

  model.addConstrs(model._ordering_variable_list[t,p,j]>=model._ordering_variable_list[t,p,i]+1-n*(1-model._routing_variable_list[t,p,i,j])
  for t in range(model._time-1) for p in range(model._number_of_physical_qubit) for i,j in model._expand_connectivity_list)

  model.update()

def set_objective_for_router(model:Model()):
  model.setObjective(model._routing_variable_list.sum(),GRB.MINIMIZE)
  model.update()

def router_initilization(qubit_connectivity_list:list):
  model=Model()
  model._qubit_connectivity_list=qubit_connectivity_list
  model._number_of_physical_qubit=np.max(qubit_connectivity_list)+1
  model._time=2
  model._expand_connectivity_list=qubit_connectivity_list+[(q1,q0) for q0,q1 in qubit_connectivity_list]
  get_edge_list(model)
  create_qubit_mapping_variables_and_constraints(model)
  create_swap_variables_and_constraints(model)
  create_ordering_variable_and_constraints(model)
  set_objective_for_router(model)
  model._max_runtime=0
  model._start_time=None
  model._find_feasible_solution=False
  model._number_of_swaps=float("inf")
  model._swap_list=None
  model.setParam('OutputFlag', 0)
  model.update()

  return model

def get_routes(model:Model(),swap_variable_list):
  routes=[]
  for program_q in range(model._number_of_physical_qubit):
    program_qubit_layer=[]
    for physical_q in range(model._number_of_physical_qubit):
      for neighbor in model._edge_list[physical_q]:
        if swap_variable_list[0,program_q,physical_q,neighbor]>=1:
            program_qubit_layer.append((physical_q,neighbor))
    routes.append(program_qubit_layer)

  return routes

def sort_route(list_of_routes):
  if len(list_of_routes)==0:
    return [],-1,-1
  ordered_routes=[list_of_routes.pop(0)]
  start_position=ordered_routes[0][0]
  end_position=ordered_routes[0][1]
  while len(list_of_routes)>0:
    for current_start,current_end in list_of_routes:
      if start_position==current_end:
        ordered_routes=[(current_start,current_end)]+ordered_routes
        list_of_routes.remove((current_start,current_end))
        start_position=current_start
      elif end_position==current_start:
        ordered_routes.append((current_start,current_end))
        list_of_routes.remove((current_start,current_end))
        end_position=current_end

  return ordered_routes,start_position,end_position

def order_routes(model:Model(),program_to_physical_mapping,unordered_routes):
  reorder_routes=[]
  unorder_copy=[]
  for program_q in range(model._number_of_physical_qubit):
    program_list=[]
    for operation in unordered_routes[program_q]:
      program_list.append(operation)
    unorder_copy.append(program_list)

  for program_q in range(model._number_of_physical_qubit):
    ordered_route,start,end=sort_route(unordered_routes[program_q])
    if start!=-1:
      if (start!=program_to_physical_mapping[0][program_q]) or (end!=program_to_physical_mapping[1][program_q]):
        return False,None
    reorder_routes.append(ordered_route)

  return True,reorder_routes

def find_next_target(model:Model(),current_routes):
  target_list=[]
  length_list=[]
  for program_q1 in range(model._number_of_physical_qubit):
    length_list.append(len(current_routes[program_q1]))
    if len(current_routes[program_q1])>0:
      target_list.append(current_routes[program_q1][0])
  max_length=max(length_list)

  for i in range(1,max_length):
    for program_q1 in range(model._number_of_physical_qubit):
      if length_list[program_q1]>i:
        start,end=current_routes[program_q1][i]
        if (end,start) in target_list:
          return current_routes[program_q1][:i],program_q1

def find_pairs(model:Model(),current_routes,program_to_physical,physical_to_program):
  for program_q1 in range(model._number_of_physical_qubit):
      if len(current_routes[program_q1])>0:
        current_start,current_end=current_routes[program_q1][0]
      else:
        continue

      for program_q2 in range(program_q1+1,model._number_of_physical_qubit):
        if len(current_routes[program_q2])>0:
          matching_start,matching_end=current_routes[program_q2][0]
          if (current_end==matching_start) and (current_start==matching_end):
            program_to_physical[program_q1]=current_end
            program_to_physical[program_q2]=matching_end
            physical_to_program[current_start]=program_q2
            physical_to_program[current_end]=program_q1
            current_routes[program_q1].pop(0)
            current_routes[program_q2].pop(0)
#              print(current_routes)
            return [(current_start,current_end)],current_routes,program_to_physical,physical_to_program

  operations,program_q1=find_next_target(model,current_routes)
  for operation in operations:
    start,end=operation
    program_q2=physical_to_program[end]
    current_routes[program_q2].insert(0,(start,end))
    current_routes[program_q1].pop(0)
    program_to_physical[program_q1]=end
    program_to_physical[program_q2]=start
    physical_to_program[start]=program_q2
    physical_to_program[end]=program_q1

  return operations,current_routes,program_to_physical,physical_to_program

def get_swap(model:Model(),routing_variable_list,qubit_mapping_variable_list):
  swap_list=[]
  current_routes=get_routes(model,routing_variable_list)
  routes_copy=get_routes(model,routing_variable_list)
  program_to_physical_mappings=[find_program_to_physical_mapping(model,qubit_mapping_variable_list,0),
                      find_program_to_physical_mapping(model,qubit_mapping_variable_list,1)]
  physical_to_program_mappings=[find_physical_to_program_mapping(model,qubit_mapping_variable_list,0),
                       find_physical_to_program_mapping(model,qubit_mapping_variable_list,1)]
  program_to_physical=program_to_physical_mappings[0].copy()
  physical_to_program=physical_to_program_mappings[0].copy()
#  print("unordered routes:",current_routes)
  success,current_routes=order_routes(model,program_to_physical_mappings,current_routes)
  _,routes_copy=order_routes(model,program_to_physical_mappings,routes_copy)

  prev=sum([len(r) for r in current_routes])
 # print("Initial pairs:",prev)
  repeate=0
  while sum([len(r) for r in current_routes])>0:
#    print(sum([len(r) for r in current_routes]),"route pairss left")
    operation,current_routes,program_to_physical,physical_to_program=find_pairs(model,current_routes,program_to_physical,physical_to_program)
    swap_list=swap_list+operation

#    if prev==sum([len(r) for r in current_routes]):
#      repeate+=1
#    else:
#      repeate=0
#    if repeate>3:
#      print("Error happens in get swap")
#      print("Routes list:",current_routes)
#      print("Program to physical:",program_to_physical)
#      print("Orginal:",program_to_physical_mappings[0])
#      sys.exit(1)

  number_of_swaps=len(swap_list)
 # print("Number of SWAP Found:",number_of_swaps)

  if number_of_swaps<prev/2:
    print("SWAP Error!")
    print("Ordered routes:",routes_copy)
    print("SWAP details:",swap_list)
    sys.exit(1)
  return number_of_swaps,swap_list

def router_callback(model,where):
  if model._start_time is None:
    model._start_time=time.time()
  elif time.time()-model._start_time>model._max_runtime and model._find_feasible_solution:
    model._start_time=None
    model._find_feasible_solution=False
    model._current_iteration=0
    model.terminate()

  if where==GRB.Callback.MIPSOL:
    model._find_feasible_solution=True
    routing=model.cbGetSolution(model._routing_variable_list)
    mapping=model.cbGetSolution(model._qubit_mapping_variable_list)
    number_of_swaps,swap_list=get_swap(model,routing,mapping)
    if number_of_swaps<model._number_of_swaps:
      model._number_of_swaps=number_of_swaps
      model._swap_list=swap_list

def router_optimization(router,previous_mapping,next_mapping,run_time_limit=10):
  router._max_runtime=run_time_limit
  previous,next=fix_location(router,previous_mapping,next_mapping)
  router.optimize(router_callback)
  number_of_swaps=router._number_of_swaps
  swap_list=router._swap_list
  router._number_of_swaps=float('inf')
  router._swap_list=None
  router.remove(previous)
  router.remove(next)
  router._start_time=None
  router._find_feasible_solution=False
  router.reset()
  router.update()

  return number_of_swaps,swap_list

def find_mappings(qubit_connectivity_list,gate_connectivity_list):
  non_execute_list=[i for i in range(len(gate_connectivity_list))]
  gate_dependecy_list=get_dependecy_list(gate_connectivity_list,np.max(qubit_connectivity_list)+1,non_execute_list)
  mapper=mapper_initilization(qubit_connectivity_list,gate_connectivity_list)
  router=router_initilization(qubit_connectivity_list)
  print("Looking for initial solutions")
  list_of_initial_mapping,non_execute_list_orginal,execution_list_orginal,gate_layer_orginal,gate_dependecy_list_orginal=find_initial_mapping(mapper,qubit_connectivity_list,
                                                                                                              gate_connectivity_list,gate_dependecy_list)
  print("Found initial solutions")

  answer_candidate=[]
  k=0
  for _,initial_mapping in list_of_initial_mapping:
    print("Start processing initial solution:",k)
    gate_execution_list=[]
    qubit_mapping_list=[]
    total_swap_list=[]
    total_number_of_swaps=0
    gate_layer=copy.deepcopy(gate_layer_orginal)
    current_mapping=initial_mapping
    non_execute_list=copy.deepcopy(non_execute_list_orginal)
    execution_list=copy.deepcopy(execution_list_orginal)
    gate_dependecy_list=copy.deepcopy(gate_dependecy_list_orginal)
    execution_list,non_execute_list,gate_layer,gate_dependecy_list=check_mapping(qubit_connectivity_list,gate_connectivity_list,gate_dependecy_list,
                                                                                 current_mapping,execution_list,gate_layer,non_execute_list)
    print("Initial mapping:",execution_list)
    gate_execution_list.append(execution_list)
    qubit_mapping_list.append(current_mapping)
 
    while len(non_execute_list)>0:
      print("Move to next gate")
      print(len(non_execute_list),"gates left")
#      print("Non execute list:",non_execute_list)
#      print("Gate dependecy list:",gate_dependecy_list)
#      print("Before checking executable gate:",)
      execution_gate=get_next_executable_gate(mapper._distance_matrix,gate_layer,gate_connectivity_list,current_mapping)
      execution_list=[execution_gate]
#      print("After checking execute gate:",execution_gate)
      next_mappings=mapper_optimization(mapper,[execution_gate],gate_layer,current_mapping,30)
#      print("Found next mappings")
      number_of_swaps,swap_list=router_optimization(router,current_mapping,next_mappings[0][1],10)
      print(number_of_swaps,"SWAPS Found")
      current_mapping=next_mappings[0][1]
      execution_list,non_execute_list,gate_layer,gate_dependecy_list=check_mapping(qubit_connectivity_list,gate_connectivity_list,gate_dependecy_list,
                                                                                   current_mapping,execution_list,gate_layer,non_execute_list)
      gate_execution_list.append(execution_list)
      qubit_mapping_list.append(current_mapping)
      total_swap_list.append(swap_list)
      total_number_of_swaps+=number_of_swaps

    print("Total number of SWAPS:",total_number_of_swaps)
    answer_candidate.append((total_number_of_swaps,gate_execution_list,qubit_mapping_list,total_swap_list))
    print("Finish processing initial solution:",k)
#    print((total_number_of_swaps,gate_execution_list,qubit_mapping_list,total_swap_list))
    k+=1
  answer_candidate.sort(key=lambda x:x[0])

  return answer_candidate[0]

def verifier(qubit_connectivity_list,gate_connectivity_list,gates,mappings,swaps):
  for i in range(len(gates)):
    for gate in gates[i]:
      q0,q1=gate_connectivity_list[gate]
      q0,q1=mappings[i][q0],mappings[i][q1]
      if not ((q0,q1) in qubit_connectivity_list or (q1,q0) in qubit_connectivity_list):
        print("Gate Constraint Not Satisfied")
        return False
    if i>0:
      previous_mapping=mappings[i-1].copy()
      for q0,q1 in swaps[i-1]:
        q0_index=previous_mapping.index(q0)
        q1_index=previous_mapping.index(q1)
        previous_mapping[q0_index]=q1
        previous_mapping[q1_index]=q0

      if previous_mapping!=mappings[i]:
        print("SWAP not Match Mapping")
        return False

  return True
