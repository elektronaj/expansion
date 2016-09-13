"""
Computing edge expansion and vertex expansion of a given graph.

Usage: 'python expansion.py filename',

where filename is the name of a .cg file representing a graph;
in a .cg file, the vertices are named 1..n,
the first line shall be 'n',
and for each edge {u,v} in the graph (u,v \in 1..n)
there shall be one line of the format 'u v'.

The program assumes Python 2.7 and the existence of IBM's CPLEX (including the Python API).
"""


import os
import sys
import cplex
import random
import itertools
from cplex.exceptions import CplexError


def cg_to_linked_list(filename):
    """
    Takes a .cg file and returns
    a representation of the graph as a linked list.
    """
    lines = [line.rstrip('\n') for line in open(filename)]
    n = lines[0]
    linked_list = [[] for _ in range(int(n))]    
    for line in lines[1:]:
        i, j = line.split(' ')
        i = int(i) - 1
        j = int(j) - 1
        linked_list[i].append(j)
        linked_list[j].append(i)    
    return linked_list


def linked_list_to_n_and_m_and_edgemap(linked_list):
    """
    Takes a linked list representation of a graph and returns
    the number n of vertices,
    the number m of edges,
    and the list of edges.
    """
    n = len(linked_list)
    m = 0
    edgemap = []
    for i in range(n):
        for j in linked_list[i]:
            if j > i:
                m += 1
                edgemap.append([i, j])
    return n, m, edgemap


def linked_list_to_edge_expansion_lp(linked_list, set_size):
    """
    Takes a linked list representation of a graph and returns
    a string containing the contents of a corresponding .lp file
    for computing its edge expansion.
    """
    n, m, edgemap = linked_list_to_n_and_m_and_edgemap(linked_list)
    lp = ''

    # Objective
    lp += 'Minimize\n'
    lp += ' obj: '
    for j in range(m):
        lp += 'z' + str(j)
        if j < m - 1:
            lp += ' + '
        else:
            lp += '\n'

    # Constraints
    lp += 'Subject to\n'

    # Constraint: size of the set is bounded
    lp += ' c1: '
    for i in range(n):
        lp += 'x' + str(i)
        if i < n - 1:
            lp += ' + '
        else:
            lp += ' = ' + str(set_size) + '\n'

    # Constraint: z's shall be 1 if exactly one endpoint is in
    for j in range(m):
        lp += ' c' + str(2 + j) + ': ' + 'z' + str(j) + ' + ' + 'x' + str(edgemap[j][0])  + ' - ' + 'x' + str(edgemap[j][1]) + ' >= 0 ' + '\n'
        lp += ' c' + str(2 + j + m) + ': ' + 'z' + str(j) + ' + ' + 'x' + str(edgemap[j][1])  + ' - ' + 'x' + str(edgemap[j][0]) + ' >= 0 ' + '\n'
            
    # Variable declarations
    lp += 'Binary\n'
    for i in range(n):
       lp += ' x' + str(i) + '\n'
    for j in range(m):
       lp += ' z' + str(j) + '\n'
    
    # End
    lp += 'End\n'

    # Return
    return lp


def linked_list_to_vertex_expansion_lp(linked_list, set_size):
    """
    Takes a linked list representation of a graph and returns
    a string containing the contents of a corresponding .lp file
    for computing its edge expansion.
    """
    n = len(linked_list)
    lp = ''

    # Objective
    lp += 'Minimize\n'
    lp += ' obj: '
    for i in range(n):
        lp += 'y' + str(i)
        if i < n - 1:
            lp += ' + '
        else:
            lp += '\n'

    # Constraints
    lp += 'Subject to\n'

    # Constraint: size of the set is bounded
    lp += ' c1: '
    for i in range(n):
        lp += 'x' + str(i)
        if i < n - 1:
            lp += ' + '
        else:
            lp += ' = ' + str(set_size) + '\n'

    # Constraint: considering closed neighborhoods
    for i in range(n):
        lp += ' c' + str(i + 2) + ': x' + str(i) + ' + ' + 'y' + str(i) + ' <= 1\n'

    # Constraint: y's shall be x's neighbors
    for j in range(n):
        lp += ' c' + str(n + 2 + j) + ': '
        for l in range(len(linked_list[j])):
            lp += 'x' + str(linked_list[j][l])
            if l < len(linked_list[j]) - 1:
                lp += ' + '
            else:
                lp += ' - y' + str(j) + ' >= 0\n'

    c_counter = 2 * n + 1
    for j in range(n):
        for i in linked_list[j]:
            c_counter += 1
            lp += ' c' + str(c_counter) + ': y' + str(j) + ' - ' + 'x' + str(i) + ' + ' + 'x' + str(j) + ' >= 0\n'
            
    # Variable declarations
    lp += 'Binary\n'
    for i in range(n):
       lp += ' x' + str(i) + '\n'
    for i in range(n):
       lp += ' y' + str(i) + '\n'
    
    # End
    lp += 'End\n'

    # Return
    return lp


def run_cplex(lp, linked_list):
    """
    Takes a string representating an .lp file
    and a linked list representing the graph,
    creates a temporary .lp file,
    runs CPLEX on it,
    removes the temporary .lp file,
    and returns a CPLEX solution object's objective and values.
    
    """
    n, m, edgemap = linked_list_to_n_and_m_and_edgemap(linked_list)
    fname = 'tmp.lp'
    f = open(fname, 'w')
    f.write(lp)
    f.close()

    cpx = cplex.Cplex(fname)

    cpx.parameters.threads.set(1)
    cpx.set_log_stream(None)
    cpx.set_error_stream(None)
    cpx.set_results_stream(None)

    try:
        cpx.solve()
    except CplexError, exc:
        print 'problem from within CPLEX ERROR:\n'
        print exc
        print 'problem ended from CPLEX.\n'
        os.remove(f)
        return

#    os.remove(f)

    return cpx.solution.get_objective_value(), cpx.solution.get_values()


def cplex_edge_expansion_solution_to_expansion_and_vertex_set(solution_objective, solution_values, linked_list):
    """
    Takes a CPLEX solution object's objective and values for edge expansion
    and a linked list representation of a graph and returns
    a tuple (expansion, vertex_set) where
    expansion is the edge expansion computed
    and vertex_set is the vertex set defining the expansion (the best cut).
    """
    expansion = solution_objective
    vertex_set = []
    n, m, edgemap = linked_list_to_n_and_m_and_edgemap(linked_list)
    for i in range(m, n + m):
        if solution_values[i] > 0.1:
            vertex_set.append(i - m)
   
    return expansion, vertex_set


def cplex_vertex_expansion_solution_to_expansion_and_vertex_set(solution_objective, solution_values, linked_list):
    """
    Takes a CPLEX solution object's objective and values for edge expansion
    and a linked list representation of a graph and returns
    a tuple (expansion, vertex_set) where
    expansion is the vertex expansion computed
    and vertex_set is the vertex set defining the expansion (the best cut).
    """
    expansion = solution_objective
    vertex_set = []
    n, m, edgemap = linked_list_to_n_and_m_and_edgemap(linked_list)
    for i in range(n, 2 * n):
        if solution_values[i] > 0.1:
            vertex_set.append(i - n)

    return expansion, vertex_set


def compute_expansion(linked_list):
    """
    Takes a linked list representation of a graph,
    computes its edge expansion and edge expansion,
    and returns a tuple (ex_expansion, vx_expansion, ex_vertex_set, vx_vertex_set)
    where ex_expansion (vx_expansion) is the edge expansion (vertex expansion)
    and ex_vertex_set (vx_vertex_set) is the corresponding vertex set (the best cut).
    """
    ex_expansion = 666
    vx_expansion = 666
    ex_vertex_set = ['inf']
    ve_vertex_set = ['inf']
    # Looping over various set sizes
    for i in range(1, len(linked_list) / 2 +  1):
        ex_lp = linked_list_to_edge_expansion_lp(linked_list, i)
        vx_lp = linked_list_to_vertex_expansion_lp(linked_list, i)
        ex_solution_objective, ex_solution_values = run_cplex(ex_lp, linked_list)
        vx_solution_objective, vx_solution_values = run_cplex(vx_lp, linked_list)
        ex_i_expansion, ex_i_vertex_set = cplex_edge_expansion_solution_to_expansion_and_vertex_set(ex_solution_objective, ex_solution_values, linked_list)
        ex_i_expansion = float(ex_i_expansion) / float(i)
        vx_i_expansion, vx_i_vertex_set = cplex_vertex_expansion_solution_to_expansion_and_vertex_set(vx_solution_objective, vx_solution_values, linked_list)
        vx_i_expansion = float(vx_i_expansion) / float(i)
        if ex_i_expansion < ex_expansion:
            ex_expansion = ex_i_expansion
            ex_vertex_set = ex_i_vertex_set
        if vx_i_expansion < vx_expansion:
            vx_expansion = vx_i_expansion
            vx_vertex_set = vx_i_vertex_set
    return (ex_expansion, vx_expansion, ex_vertex_set, vx_vertex_set)


def main():
    """
    Takes a filename as a command-line argument,
    computes its edge expansion and vertex expansion,
    and print the corresponding values as well as the corresponding sets (the best cuts).
    """
    filename = sys.argv[1]
    linked_list = cg_to_linked_list(filename)
    print '\nComputing edge expansion and vertex expansion of %s\n'%(filename)
    ex_expansion, vx_expansion, ex_vertex_set, vx_vertex_set = compute_expansion(linked_list)
    FORMAT = '%-27s of %s is %s'
    print FORMAT%('edge expansion', filename, ex_expansion)
    print FORMAT%('vertex expansion', filename, vx_expansion)
    print FORMAT%('edge expansion vertex set', filename, ex_vertex_set)
    print FORMAT%('vertex expansion vertex set', filename, vx_vertex_set)
    print '\n%s;%s;%s;%s\n'%(ex_expansion, vx_expansion, ex_vertex_set, vx_vertex_set)
main()
