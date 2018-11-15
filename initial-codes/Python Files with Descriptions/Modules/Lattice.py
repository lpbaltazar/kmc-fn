import random as rd
import numpy as np


def count(lattice, length):
    """ Counts the population of the lattice
    
    Parameter:
    -------------------
    lattice: 2-dimensional numpy array
        lattice to be counted
    length: integer
        length of the lattice
    """
    n_count = 0
    g_count = 0
    s_count = 0
    t_count = 0
    b_count = 0

    for i in range(length):
        for j in range(length):
            if lattice[j][i] == 0:
                s_count += 1
            elif 1 >= lattice[j][i] >= 0.5:
                t_count += 1
            elif 0.5 > lattice[j][i] > 0:
                g_count += 1
            elif lattice[j][i] > 1:
                n_count += 1
            else:
                b_count += 1  # negative numbers

    m_sum = s_count + t_count + g_count + n_count + b_count
    print("Source:", s_count)
    print("Tolerant: ", t_count)
    print("Gullible: ", g_count)
    # print("Neutral: ", n_count)
    # print("b_count", b_count)
    print("TOTAL", m_sum)
    return


def lattice_initialize(lattice, length, coverage):
    """Initializes the lattice with the preferred coverage of gullible
    
    Parameters:
    ---------------
    lattice: 2-dimensional numpy array
        lattice to be updated
    length: integer
        length of the lattice
    coverage: float
        preferred coverage percentage
        
    Returns:
    -------------
    lattice: 2-dimensional numpy array
        updated lattice
    """
    
    counter = 0  # main counter
    gullible_count = 0
    tolerant_count = 0
    maximum = length ** 2 - 1

    max_gullible = maximum * coverage
    max_tolerant = maximum - 1 - max_gullible

    while counter < maximum:
        # guess random position
        xs = rd.choice(range(length))
        ys = rd.choice(range(length))
        v = rd.triangular()
        # print(v, counter, tolerant_count, gullible_count, xs, ys)

        if lattice[ys][xs] == 10000:
            if v > 0.5 and tolerant_count < max_tolerant:
                lattice[ys][xs] = v
                tolerant_count += 1
                counter += 1
            elif v < 0.5 and gullible_count < max_gullible:
                lattice[ys][xs] = v
                gullible_count += 1
                counter += 1

    count(lattice, length)
    np.savetxt("Initial Population " + str(coverage) + ".csv", lattice, delimiter = ",")  # Save init Population to csv
    return lattice


def gullible_counter(lattice, length):
    """ Counts the population of the lattice

    Parameter:
    -------------------
    lattice: 2-dimensional numpy array
        lattice to be counted
    length: integer
        length of the lattice
    """

    counter = 0
    for i in range(length):
        for j in range(length):
            if lattice[j][i] < 0.5 and lattice[j][i] != 0:
                counter += 1
    return counter


def count_source(lattice, length):
    """ Counts the population of the lattice

    Parameter:
    -------------------
    lattice: 2-dimensional numpy array
        lattice to be counted
    length: integer
        length of the lattice
    """

    counter = 0
    for i in range(length):
        for j in range(length):
            if lattice[j][i] == 0:
                counter += 1
    return counter
