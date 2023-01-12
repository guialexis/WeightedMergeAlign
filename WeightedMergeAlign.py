# -*- coding: utf-8 -*-
"""
Weighted MergeAlign for MSA Consensus

Authors : Antoine, Renwen, Guillaume

Based on the MergeAlign Algorithm
"""

import sys
import pymsa
import os
import pandas as pd

class Node:
    """ A position in alignment space with edges to the previous nodes in all alignments """
    
    def __init__(self):
        self.previous_nodes = {}
        self.path_score = 0
        self.path_length = 0
        self.path_average = 0
        self.best_previous_node = None

def weightedMA(filenames, score_type = "PAM_SOP", weights = None):
    """ Main function : Return an alignment that is a mix of the entry alignements, based on their scores with the asked scoring method"""
    lst_alignments = [getFASTA(filename) for filename in filenames]
    
    original = dict([(name, seq.replace('-', '')) for name, seq in lst_alignments[0].items()])
    sequence_names = list(original.keys())
    
    indices = [[seqtoind(alignment[seq]) for seq in sequence_names] for alignment in lst_alignments]
    
    if weights == None:
        w = [score_align(alignment, score_type) for alignment in lst_alignments]
    else:
        w = weights
    coordinates = [[a for a in zip(*alignment)] for alignment in indices]
    
    graph = create_graph(coordinates, w)
    path = best_path(graph) #Liste de tuples correspondant aux noeuds
    align_ind = [[t[i] for t in path] for i in range(len(path[0]))]
    
    final_align = {}
    
    for i, seq in enumerate(align_ind):
        seq_name = sequence_names[i]
        final_align[seq_name] = indtoseq(seq, original[seq_name])
        
    full_score = total_score(final_align)
    return final_align, full_score
    
def getFASTA(filename):
    """ Read file in FASTA format.
        Return a dictionary with key=sequence name, value=sequence. """

    try:
        f = open(filename, 'r')
    except IOError:
        print("Unable to open file" + str(filename))
        sys.exit(2)

    sequences = {}
    seq_id = None

    for line in f.readlines():
        if line.startswith('>'):
            seq_id = line.rstrip()[1:]
            sequences[seq_id] = ''
        elif seq_id != None:
            sequences[seq_id] += line.rstrip()

    if len(sequences) == 0:
        print(str(filename)+ "contains no sequences")
        sys.exit(2)
    
    
    return sequences

def seqtoind(sequence):
    """ Convert a sequence to a list of amino acid indices.
        Gap are converted to the index of the preceeding amino acid. """
    
    indices = []
    i = 0
    
    for aa in sequence:
        if aa != '-':
            i += 1
        indices.append(i)
        
    return indices

def indtoseq(seq_ind, original):
    """ Convert list of indices back to original sequence plus gaps """
    
    sequence = ''

    previous_i = 0
    for i in seq_ind:
        if i != previous_i:
            sequence += original[i-1]
        else:
            sequence += '-'
        previous_i = i
        
    return sequence

def score_align(alignment, score_type):
    """ Different scoring method implement by the library pymsa"""
    
    seq_id = list(alignment.keys())
    sequences = [alignment[ids] for ids in seq_id]
    
    msa = pymsa.MSA(sequences, seq_id)
    
    if score_type == "Non_Gaps":
        non_gaps = pymsa.PercentageOfNonGaps(msa)
        score = non_gaps.compute()
    
    elif score_type == "Totally_Conserved_Columns":
        tot_cons = pymsa.PercentageOfTotallyConservedColumns(msa)
        score = tot_cons.compute()
    
    elif score_type == "Entropy":
        ent = pymsa.Entropy(msa)
        score = ent.compute()
    
    elif score_type == "PAM_SOP":
        score = pymsa.SumOfPairs(msa, pymsa.PAM250()).compute()
    
    elif score_type == "BLOSUM_SOP":
        score = pymsa.SumOfPairs(msa, pymsa.Blosum62()).compute()
    
    elif score_type == "PAM_Star":
        score = pymsa.Star(msa, pymsa.PAM250()).compute()
    
    elif score_type == "BLOSUM_Star":
        score = pymsa.Star(msa, pymsa.Blosum62()).compute()
    
    else:
        raise Exception("Invalid_Score_Type")
        
    return score

def total_score(alignment):
    """Compute all the different scores of the alignment"""
    return [score_align(alignment, "Non_Gaps"), score_align(alignment, "Totally_Conserved_Columns"), score_align(alignment, "Entropy"), 
                  score_align(alignment, "PAM_SOP"), score_align(alignment, "BLOSUM_SOP"), score_align(alignment, "PAM_Star"), score_align(alignment, "BLOSUM_Star")]

def create_graph(alignments, w):
    """ Traverse each alignment creating nodes at each point found in alignment space. 
        For each node record which nodes preceded it and how often. """
    
    dimensions = len(alignments[0][0])
    first_node = tuple([0]*dimensions)
    nodes = {first_node: Node()}
    
    for alignment, weight in zip(alignments, w):
        previous_node = first_node
        
        for point in alignment:
            if not point in nodes:
                nodes[(point)] = Node()
            node_dict = nodes[(point)].previous_nodes
            # Add previous node to current nodes count of previous nodes
            node_dict[previous_node] = node_dict.get(previous_node, 0) + weight
            previous_node = point

    return nodes

def best_path(graph, num_paths = 100):
    """ Travserse nodes, finding the best route to each one with a dynamic programming approach. 
        Return the best path to the final node. """
   
    dimensions = len(list(graph.keys())[0])
    first_node = tuple([0]*dimensions)
    
    for coord, current_node in sorted(graph.items())[1:]:
        previous_nodes = current_node.previous_nodes.items()
        best_node = max(previous_nodes, key=lambda n: n[1] + graph[n[0]].path_average)
        #best_node = max(previous_nodes, key=lambda n: (n[1] + graph[n[0]].path_score)/(graph[n[0]].path_length+1))
        
        # Update this node with the best path so far
        current_node.path_length = graph[best_node[0]].path_length + 1
        current_node.path_score = best_node[1] + graph[best_node[0]].path_score
        current_node.path_average = current_node.path_score* 1.0 / current_node.path_length
        current_node.best_previous_node = best_node[0]
    
    
    
    # Start at the final node and work backwards, moving from each node to its best previous node
    path = []
    scores = []
    while coord != first_node:
        path.append(coord)
        scores.append(graph[coord].previous_nodes[graph[coord].best_previous_node]*1.0/num_paths)
        coord = graph[coord].best_previous_node

    path.reverse()
    scores.reverse()
    
    return path

def outputAlignmentAsFASTA(filename, final_alignment, threshold=None, scores=None):
    """ Output alignment in form of a dictionary as a FASTA file. """
    
    with open(filename, 'w') as output:
        
        for name, sequence in final_alignment.items():
            output.write(">%s\n" % name)
            if threshold:
                sequence = ''.join([aa for (aa, score) in zip (sequence, scores) if score>threshold])
            output.write("%s\n" % sequence)


##Example of use

#data dictionnary enables to output results as an excel file
data = {}

#Folder with the alignments as FASTA Files
r = 'C:/Users/guill/Documents/M1_MPRI/Bioinfo'

#Computation of the consensus alignment with each possible type of scoring
for score_type in ["Non_Gaps", "Totally_Conserved_Columns", "Entropy", "PAM_SOP", "BLOSUM_SOP", "PAM_Star", "BLOSUM_Star"]:
    
    align, score = weightedMA([r+"/Tests/BB3/"+file for file in os.listdir(r+"/Tests/BB3")], score_type)
    outputAlignmentAsFASTA("BB3_result"+score_type+".fasta", align)
    data[score_type] = score

#Computation of the consensus alignment with equal weights
align, score = weightedMA([r+"/Tests/BB3/"+file for file in os.listdir(r+"/Tests/BB3")], weights = [1, 1, 1, 1, 1])
outputAlignmentAsFASTA("BB3_result_eqweights.fasta", align)
data["EqWeights"] = score

#Computation of the scores of the base alignments
for file in os.listdir(r+"/Tests/BB3"):
    
    alignment = getFASTA(r+"/Tests/BB3/"+file)
    data[file] = total_score(alignment)


#Output results as an xlsx file
df = pd.DataFrame(data)
df.to_excel('Results_BB3.xlsx', sheet_name='sheet1', index=False)
