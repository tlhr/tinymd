# tinymd
A very simple and primitive molecular dynamics simulation in rust.

## Installation and Usage
To build:
```sh
cargo build
```

Configuration is currently done by editing the source files. The simulation takes plain files as input and only works on monoatomic rare-gas clusters, using the Lennard-Jones potential. The arrays used are fixed size and have to be changed to the appropriate system size. Correctness of the output has not been verified!
