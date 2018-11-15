import Modules.Lattice as lat
import Modules.Lattice_Update as upd
import numpy as np

decay_const: int = 1
n_0 = 1
step = 0

n = 100    # length of the matrix
N = np.full(shape=(n, n), fill_value=10000, dtype=np.float32)  # creates matrix will all values having 10000
tm = 0  # time

source_counter = 0
max_source = (n**2) * 0.9  # stopping criteria

xs = 50
ys = 50

# xs = rd.choice(range(n))  # initial x and y for the source, (value == 0)
# ys = rd.choice(range(n))

N[ys][xs] = 0
print(xs, ys)

print("------------------------------A")
lat.lattice_initialize(N, n, 0.1)  # Initializes the population to have 1/3 of each type of persons
print("------------------------------B")

lat.count(N, n)
while source_counter < max_source:
    print("-----------------------------------BEFORE")
    upd.update_lattice_equation(N, n, n_0, 1, decay_const)  # updates lattice
    print("-----------------------------------AFTER")
    source_counter = lat.count_source(N, n)  # counts the source population of the lattice
    print("Step: ", tm)
    # lat.count(N, n)
    print(source_counter, max_source, tm)
    tm += 1
np.savetxt("Final Population.csv", N, delimiter = ",")
print("END----------------------------")
print()
lat.count(N, n)
print()
print("Step END: ", tm)
