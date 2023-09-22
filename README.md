# Spiking RNN simulator

## Introduction
This package provides an implementation of a spiking recurrent neural network simulator with c++.
The network can have multiple neural populations, different connectivity profiles (all to all, sparse, tuned, ...).
For more info look at the config files in ./conf/.

## Installation
This package requires
A C++ Compiler supporting C++17, CMake (optional) and yaml-cpp

To install with CMake:
```bash
cd bin/
cmake ..
make
```

Otherwise you can just use the *Makefile*
```bash
make
```


## Usage
Assuming the dependencies are installed, here is how to run the model (see notebooks folder or org folder for more doc)

from python
```python
from run_model import run_cpp

run_cpp('session_name', 'bin_path', 'config_file_path')
```

## Contributing
Feel free to contribute.
```
MIT License
Copyright (c) [2023] [A. Mahrach]
```
