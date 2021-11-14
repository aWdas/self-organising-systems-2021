import array
import random
import numpy as np

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from deap import algorithms
from deap import base
from deap import creator
from deap import tools

from .util import generalize, seq_str, dist, build_pairwise_alignments, fix_odd_sequence_count


def mwm_ga(sequences: [SeqRecord], alignment_result):
    # build a map of the pairwise distances between the sequences,
    # and keep track of the pair with the minimum distance
    if alignment_result is None:
        print("No pre-computed pairwise alignments were provided -> calculating them now")
        alignment_result = build_pairwise_alignments(sequences)
    pairwise_alignments = alignment_result["pairwise_alignments"]
    min_dist_pair = alignment_result["min_dist_pair"]

    # apply the fix suggested in the paper to get to an even number of sequences
    pairwise_alignments = fix_odd_sequence_count(sequences, pairwise_alignments, min_dist_pair)


    # Build the configuration required for DEAP
    # This example was used as a guideline: https://github.com/DEAP/deap/blob/master/examples/ga/tsp.py
    # Parameters have then been tuned to improve the results

    # we use one fitness metric, which is minimized
    creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
    # an individual is represented by an int-array
    creator.create("Individual", array.array, typecode='i', fitness=creator.FitnessMin)

    toolbox = base.Toolbox()

    # to generate a new individual, the indices present in pairwise_alignments are sampled in random order
    indices = pairwise_alignments.keys()
    ind_size = len(indices)
    toolbox.register("indices", random.sample, indices, ind_size)
    toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.indices)

    # the population is a list of individuals
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    # use Partially Mapped Crossover (PMX) for crossover
    toolbox.register("mate", tools.cxPartialyMatched)
    # use index shuffling as mutation operator; indpb is the probability with which each slot value is moved
    toolbox.register("mutate", tools.mutShuffleIndexes, indpb=0.05)
    # use roulette wheel selection
    toolbox.register("select", tools.selTournament, tournsize=3)

    # use the sum of the distances between the paired sequences as fitness value
    def eval_pairing(individual):
        dist_sum = 0
        for k in range(0, len(individual), 2):
            i, j = individual[k], individual[k+1]
            if i > j:
                i, j = j, i
            dist_sum += pairwise_alignments[i][j][2]
        return dist_sum,
    toolbox.register("evaluate", eval_pairing)

    population = toolbox.population(n=300)
    hof = tools.HallOfFame(1)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("std", np.std)
    stats.register("min", np.min)
    stats.register("max", np.max)

    # Execute the GA algorithm
    algorithms.eaSimple(population, toolbox, 0.7, 0.2, 1000, stats=stats, halloffame=hof, verbose=True)

    best_pairing = hof[0]

    generalized_sequences = []
    generalized_distance = 0
    for k in range(0, len(best_pairing), 2):
        i, j = best_pairing[k], best_pairing[k+1]
        if i > j:
            i, j = j, i
        pairwise_alignment = pairwise_alignments[i][j]

        # if f is part of the sequence match, replace it by its sources a and b and
        # include their distance d in the distance result
        if j == len(sequences):
            a, b, d = min_dist_pair
            generalized_sequences.append({
                "group": [pairwise_alignment[0].id, sequences[a].id, sequences[b].id],
                "generalizedSeq": generalize(seq_str(pairwise_alignment[0]), seq_str(pairwise_alignment[1])),
                "distance": pairwise_alignment[2] + d
            })
            generalized_distance = generalized_distance + pairwise_alignment[2] + d
        else:
            generalized_sequences.append({
                "group": [pairwise_alignment[0].id, pairwise_alignment[1].id],
                "generalizedSeq": generalize(seq_str(pairwise_alignment[0]), seq_str(pairwise_alignment[1])),
                "distance": pairwise_alignment[2]
            })
            generalized_distance = generalized_distance + pairwise_alignment[2]

    return generalized_sequences, generalized_distance
