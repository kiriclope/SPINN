# Spiking RNN simulator

## Introduction
This package provides an implementation of a spiking recurrent neural network simulator with c++.
The network can have multiple neural populations, different connectivity profiles (all to all, sparse, tuned, ...).
For more info look at the config files in ./conf/.

## Installation
Provide clear instructions on how to get your development environment running.
```bash
cd bin/
cmake ..
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
