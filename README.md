# Ising-Model-CUDA
The evolution of an Ising model, in CUDA, in two dimensions for a given number of steps 𝑘. The Ising model (/ˈaɪsɪŋ/; German: [ˈiːzɪŋ]), named after the physicist Ernst Ising, is a mathematical model of ferromagnetism in statistical mechanics. The model consists of discrete magnetic dipole moments of atomic “spins” that can be in one of two states (+1 or −1). The spins are arranged in a square lattice with periodic boundary conditions, allowing each spin to interact with its neighbors. The dipole moments update in discrete time steps according to the majority of the spins within the 3 × 3 window centered to each lattice point. Windows centered on the edge lattice points wrap around to the other side (known as toroidal or periodic boundary conditions).

## Validate Results

To validate the results of all versions run:

```sh
make all
```

To validate a specific version(e.g. "V1") run:

```sh
make V1
```

For cleaning:
```sh
make clean
```
