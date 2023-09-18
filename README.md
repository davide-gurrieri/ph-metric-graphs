# ph-metric-graphs

Computation of persistent homology on metric graphs, with particular focus on vascular networks.

This repository is developed as part of the “Advanced Programming for Scientific Computing” course at Polytechnic University of Milan.

The aim is to use topological data analysis to study geometrical features of vascular networks. This allows the detection of anomalies like devascularized regions, loops and chaotic patterns which can be used for disease diagnosis purposes.

The project is developed in C++ and Python, using the [GUDHI](https://gudhi.inria.fr/) library for the computation of persistent homology.

A detailed description of the work can be found in the [report](https://github.com/davide-gurrieri/ph-metric-graphs/blob/main/report/report.pdf).

To use the repository `cd` to the folder you wish to install it, and clone it:

```shell
git clone https://github.com/davide-gurrieri/ph-metric-graphs.git
```

## Structure

```
├── README.md
├── data
│   ├── complex_network
│   └── synthetic_network
├── docs
│   └── Doxyfile
├── output
│   ├── README
│   ├── plot
│   └── plot_python
├── python
│   ├── network_analysis.py
│   ├── plot_ph_from_cpp_file.py
│   ├── real_networks_analysis.ipynb
│   ├── requirements.txt
│   └── synthetic_network_analysis.ipynb
├── report
│   └── report.pdf
└── src
    ├── CMakeLists.txt
    ├── GetPot
    ├── chrono.hpp
    ├── geometry.h
    ├── main.cpp
    ├── network_analysis.cpp
    ├── network_analysis.h
    ├── options1.txt
    └── options2.txt
```

## C++ version

### Prerequisites

To build `GUDHI` you will need `cmake`, `Boost`, `CGAL`, `EIGEN3`, `OpenGL` and `TBB`.

On Linux machines, it is sufficient to run:

```shell
sudo snap install cmake --classic # cmake
sudo apt-get install libboost-all-dev # Boost
sudo apt-get install libcgal-dev # CGAL
sudo apt-get install libeigen3-dev # Eigen
sudo apt install freeglut3-dev # OpenGL
sudo apt install libtbb-dev # TBB
```

To build and install `GUDHI`:

```shell
git clone https://github.com/GUDHI/gudhi-devel.git
cd gudhi-devel
mkdir build
cd build/
cmake -DCMAKE_BUILD_TYPE=Release ..
make
make install
```

### Compilation and execution

To compile the C++ source code `cd` to the repository folder and run:

```bash
cd src/
mkdir build && cd build/
cmake ..
make
```

You can run the analysis on both the real (`options1.txt`) and the synthetic networks (`options2.txt`) by running:

```bash
./network_analysis file=options1.txt
```

To select a different synthetic network, change `options2.txt`. The available networks are in the `data` folder.

To obtain the documentation run:

```shell
sudo apt install doxygen
cd docs/
doxygen
```

Open the `index.html` file in the `html` folder to view the documentation.

## Python version

### Prerequisites

You need pip to get te prerequisites for the Python version.

```shell
sudo apt install python3-pip
```

To install the prerequisites run from the repository folder:

```shell
cd python
pip install -U -r requirements.txt
```

### Execution

You can execute jupyter notebooks in a web browser by running:

```shell
jupyter notebook
```

## Results

All the plots will be in the `output` folder.
