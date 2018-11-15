import random as rd
import math as mt


def main_equation(value, x, y, i, j, multiplier):
    """ Implements the original reduction equation

    new_value = old_value + multiplier * (-alpha * ( 1 / distance))

    Parameters
    -----------
    value: float
        value of the current lattice position
    i: integer
        x coordinate of the source
    j: integer
        y coordinate of the source
    x: integer
        x coordinate of the target cell
    y: integer
        y coordinate of the target cell
    multiplier: float
        exponential decay term

    Returns
    -----------
        out: float
            new value of the cell
    """
    alpha = rd.triangular()  # random number
    out = value + multiplier * (-alpha * (1/(mt.sqrt(((x - i)**2 + (y - j)**2)))))
    return out


def update_lattice_equation(lattice, length, initial_n, time, decay_constant):
    """ Updates the lattice cells which are neighbors to the source

       Parameters
       -----------
       lattice: numpy 2-dimensional array
           lattice to be updated
       length: integer
           length of the side of a lattice
       initial_n: integer
           initial number of sources
       time: integer
           number of time scans already done
       decay_constant: float
           exponential decay term

       Returns
       -----------
           lattice: numpy 2-dimensional lattice
               updated lattice
       """

    for i in range(length):
        for j in range(length):
            if lattice[j][i] == 0:

                multi = initial_n * (mt.exp(-time / decay_constant))
                
                try:  # below the source
                    lattice[j+1][i] = main_equation(lattice[j + 1][i], i, j, i, j+1, multi)
                    # print("1")
                except IndexError:
                    pass

                try:  # above the source
                    lattice[j-1][i] = main_equation(lattice[j + 1][i], i, j, i, j-1, multi)
                    # print("2")
                except IndexError:
                    pass

                try:  # right to the source
                    lattice[j][i+1] = main_equation(lattice[j + 1][i], i, j, i + 1, j, multi)
                    # print("3")
                except IndexError:
                    pass

                try:  # left to the source
                    lattice[j][i-1] = main_equation(lattice[j + 1][i], i, j, i - 1, j, multi)
                    # print("4")
                except IndexError:
                    pass

                try:  # diagonal down and right
                    lattice[j+1][i + 1] = main_equation(lattice[j + 1][i], i, j, i+1, j+1, multi)
                    # print(5)
                except IndexError:
                    pass

                try:  # diagonal down and left
                    lattice[j+1][i-1] = main_equation(lattice[j + 1][i], i, j, i-1, j+1, multi)
                    # print("6")
                except IndexError:
                    pass

                try:  # diagonal up and right
                    lattice[j-1][i+1] = main_equation(lattice[j + 1][i], i, j, i+1, j-1, multi)
                    # print("7")
                except IndexError:
                    pass

                try:  # diagonal up and left
                    lattice[j-1][i-1] = main_equation(lattice[j + 1][i], i, j, i-1, j-1, multi)
                    # print("9")
                except IndexError:
                    pass

    for i in range(length):
        for j in range(length):
            if lattice[j][i] < 0:
                lattice[j][i] = 0

    return lattice


def update_matrix_equation(lattice, decay_constant, initial_n, length, time, x0, y0):  # Equation Based Approach
    """ Updates all lattice cells depending on the location of the chosen source

    Parameters:
    -----------------
    lattice: 2-dimensional numpy array
        lattice to be updated
    decay_constant: float
        exponential decay multiplier
    initial_n: integer
        initial number of sources
    length: integer
        length of the lattice
    time: integer
        number of time scans
    x0: integer
        x coordinate of source
    y0: integer
        y coordinate of source

    Returns:
    ------------------
    lattice: 2-dimensional numpy array
        updated lattice
    """

    multiplier = initial_n * (mt.exp(-time / decay_constant))

    for i in range(length):
        for j in range(length):
            if i != x0 and j != y0:
                alpha = rd.triangular()
                lattice[j][i] = lattice[j][i] + multiplier * (-alpha * 1 / (mt.sqrt(((x0 - i) ** 2 + (y0 - j) ** 2))))

    for i in range(length):
        for j in range(length):
            if lattice[j][i] < 0:
                lattice[j][i] = 0

    return lattice


def new_model(value):
    """ New model for reducing lattice values

    Parameters
    -------------------
    value: float
        value of the lattice cell to be updated

    Returns
    ------------------
    value: float
        updated value of the lattice cell
    """
    return value


def gullible_selection(lattice, x0, y0):
    """ Updating lattice using gullible selection

    Parameters:
    -----------------
    lattice: 2-dimensional numpy array
        lattice to be updated
    x0: integer
        x coordinate of the source selected
    y0: integer
        y coordinate of the source selected

    Returns:
    -----------------
    lattice: 2-dimensional numpy array
        updated lattice
    """

    test = {1: [x0 - 1, y0 - 1], 2: [x0, y0 - 1], 3: [x0 + 1, y0 - 1], 4: [x0 - 1, y0], 5: [x0 + 1, y0],
            6: [x0 - 1, y0 + 1], 7: [x0, y0 + 1], 8: [x0 + 1, y0 + 1]}  # near neighbors of the source
    gul = []  # list for the locations for each gullible
    tol = []  # list for the locations for each tolerant
    for i in range(1, 9):
        if test[i][1] >= 0 and test[i][0] >= 0:  # checks if the location is valid
            try:
                # print(lattice[[test[i][1]], [test[i][0]]])
                if 0 < lattice[[test[i][1]], [test[i][0]]] <= 0.5:  # checks if value is gullible
                    gul.append(i)
                else:
                    tol.append(i)
            except IndexError:
                pass
    if gul == 0:
        total = tol
    else:
        total = gul
    select = rd.choice(total)

    lattice[[test[select][1]], [test[select][0]]] = new_model(lattice[[test[select][1]], [test[select][0]]])
    # changes value of the lattice

    return lattice


def monte_carlo_selection(lattice, x0, y0):
    """ Updating lattice using monte carlo method

    Parameters:
    -----------------
    lattice: 2-dimensional numpy array
        lattice to be updated
    x0: integer
        x coordinate of the source selected
    y0: integer
        y coordinate of the source selected

    Returns:
    -----------------
    lattice: 2-dimensional numpy array
        updated lattice
    """

    test = {1: [x0 - 1, y0 - 1], 2: [x0, y0 - 1], 3: [x0 + 1, y0 - 1], 4: [x0 - 1, y0], 5: [x0 + 1, y0],
            6: [x0 - 1, y0 + 1], 7: [x0, y0 + 1], 8: [x0 + 1, y0 + 1]}  # near neighbors of the source
    total_gul = 0  # total for gullible
    total_tol = 0  # total for tolerant
    mode = 0

    for i in range(1, 9):
        if test[i][1] >= 0 and test[i][0] >= 0:  # checks if the location is valid
            try:
                # print(lattice[[test[i][1]], [test[i][0]]])
                if 0 < lattice[[test[i][1]], [test[i][0]]] <= 0.5:  # checks if value is gullible
                    total_gul += lattice[[test[i][1]], [test[i][0]]]
                else:
                    total_tol += lattice[[test[i][1]], [test[i][0]]]
            except IndexError:
                pass

    if total_gul == 0:
        total = total_tol
        mode = 1
    else:
        total = total_gul

    select = rd.triangular(total)

    mov = 0
    change = 0

    for i in range(1, 9):
        if test[i][1] >= 0 and test[i][0] >= 0:  # checks if the location is valid
            if mode == 0:
                try:
                    # print(lattice[[test[i][1]], [test[i][0]]])
                    if 0 < lattice[[test[i][1]], [test[i][0]]] <= 0.5:  # checks if value is gullible
                        mov += lattice[[test[i][1]], [test[i][0]]]
                except IndexError:
                    pass
            else:
                try:
                    # print(lattice[[test[i][1]], [test[i][0]]])
                    if 0.5 < lattice[[test[i][1]], [test[i][0]]] < 1.0:  # checks if value is gullible
                        mov += lattice[[test[i][1]], [test[i][0]]]
                except IndexError:
                    pass

        if mov >= select:
            change = i
            break

    lattice[[test[change][1]], [test[change][0]]] = new_model(lattice[[test[change][1]], [test[change][0]]])

    return lattice
