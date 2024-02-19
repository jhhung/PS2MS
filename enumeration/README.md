# Enumeration
- PS<sup>2</sup>MS build the synthetic database by enumerating possible derivatives.

## Requirement

- GNU [g++-10](https://gcc.gnu.org/gcc-10/) or higher
- [CMake 3.16.0](https://cmake.org/download/) or higher to build the enumeration step and the detection step
- [rdkit](https://www.rdkit.org/docs/Install.html)
  - build the c++ code from the source and install python package from conda

## Usage:

### Install conda env 
```bash
conda env create --name drug_detection -f drug-detection/environment.yml
```

### Build project
```
cd drug-detection
conda activate drug_detection

mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j
```

### Run permutation
```bash
cd drug-detection
conda activate drug_detection

./build/build_database <split_count> <split_idx> <output_dir>

# for example: ./build/build_database 200 0 /tmp
```

