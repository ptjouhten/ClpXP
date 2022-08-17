import cobra
import cplex
import docplex.mp.model as cpx
from collections import Counter
from itertools import combinations, chain
from math import isinf
from math import floor, ceil
from cobra.medium import minimal_medium
import pandas as pd
import numpy as np
import json


def get_flux_distances(sbml_model, source_flux, outfile_name): 

    model = sbml_model.copy()
    #solver = cpx.Model(name="fba")

    #let's read which metabolites should not be considered to create adjacency between reactions
    with open("yeast_GEM_86_mets_remove.txt") as file:
        lines = file.readlines()
        lines = [line.rstrip() for line in lines]
    remove_mets_list = lines
    #file_content.split(",").strip()
    file.close()
    print(len(remove_mets_list))
    #print("The list is: ", remove_mets_list)
    #variables = {}

    #create variables
    #for rxn in model.reactions:
        #print(rxn.id)
    #    variables[rxn.id] = solver.continuous_var(name=rxn.id, lb=rxn.lower_bound, ub=rxn.upper_bound)

    #create the stoichiometric constraints for the flux variables
    A = cobra.util.array.create_stoichiometric_matrix(model, array_type="DataFrame")
    print(A.shape)
    #remove mets to be removed
    A=A.drop(labels=remove_mets_list, axis=0)
    print(A.shape)
    #A=A.T[A.any()].T
    #A = A.loc[:, (A!= 0).any(axis=1)]
    A = A.loc[(A).any(1), (A!=0).any(0)]
    #print(list(A.index))
    print(A.shape)
    Amat = A.to_numpy()
    #get the metabolite names
    #names = A.index
    fluxes = list(A.columns)
    #print(A['r_1110'])
    #for v in A['r_1110']:
    #    print(v)
    A = A.transpose()
    Amat_t=A.to_numpy()
    #create a constraint for each metabolite
    #for name in names:
    #    r = A[A[name] != 0]
        #get the fluxes
    #    index = list(r.index)
    #    var = None
        #get the variables corresponding to the fluxes
    #    c = [variables[rxn_id] for rxn_id in index]
    #    solver.add_constraint(solver.sum(c * r[name]) == 0)


    #create undirected reaction adjacency matrix
    Ad = np.matmul(Amat_t,Amat)
    #print(Ad.size)
    #calculate distances from a source flux to all other fluxes
    solver_dist = cpx.Model(name="flux_dist")
    variables_dist = {}
    for rxn in fluxes:
        #if rxn == 'r_1110':
        #    break
        variables_dist[rxn] = solver_dist.continuous_var(name= rxn, lb=0, ub=10000)
        if rxn == source_flux:
            solver_dist.add_constraint(variables_dist[rxn] == 0)

    for j, rxn_j in enumerate(fluxes):
        #print(j)
    #for rxn_j in model.reactions:
    	#rxn_j_ind = fluxes.index(rxn_j)
        #for rxn_i in model.reactions:
        i = j+1
        while i < len(fluxes):
            #print(i)
            rxn_i = fluxes[i]
        	#rxn_i_ind = fluxes.index(rxn_i)
        	#add constraint if the reactions are adjacent to each other
            #print(Ad[i,j])
            if Ad[i,j] != 0:
                #print("Adjacent: "+rxn_i+ " and "+rxn_j)
                solver_dist.add_constraint(variables_dist[rxn_j] <= variables_dist[rxn_i] + 1)
                solver_dist.add_constraint(variables_dist[rxn_i] <= variables_dist[rxn_j] + 1)
            i = i+1

    solver_dist.maximize(solver_dist.sum(variables_dist))
    sol_dist = solver_dist.solve()

    flux_distance_to_source = {}
    for rxn_id in fluxes:
        flux_distance_to_source[rxn_id] = sol_dist.get_value(variables_dist[rxn_id])

    with open('flux_distances_'+outfile_name+'.txt', 'w') as f:
        for key,value in flux_distance_to_source.items():
            f.write(key + "\t" + str(value) + "\n")
    f.close()

    return flux_distance_to_source

        
        

