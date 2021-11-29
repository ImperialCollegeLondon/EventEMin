# EventEMin
[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg?style=flat-square)](https://creativecommons.org/licenses/by-nc-sa/4.0/)

Code for the following papers:
- [Entropy Minimisation Framework for Event-based Vision Model Estimation](http://www.ecva.net/papers/eccv_2020/papers_ECCV/papers/123500154.pdf)
- [Robust Event-based Vision Model Estimation by Dispersion Minimisation](https://ieeexplore.ieee.org/document/9625712)

The authors provide this code in the hope it will be useful for understanding the proposed method, as well as for reproducibility of the results.

For more information and more open-source software please visit the Personal Robotic Lab's website: <https://www.imperial.ac.uk/personal-robotics/>.

## Requirements
This code was tested on Ubuntu 16.04, 18.04 and 20.04 distros.

### Dependencies
CMake: <https://cmake.org/download/>
  ```
  sudo apt install cmake cmake-curses-gui
  ```

OpenCV: <https://docs.opencv.org/trunk/d7/d9f/tutorial_linux_install.html>

Eigen3: <http://eigen.tuxfamily.org/index.php?title=Main_Page#Download>
  ```
  git clone https://gitlab.com/libeigen/eigen.git
  cd eigen
  mkdir build && cd build
  cmake ..
  sudo make install
  ```
or
  ```
  sudo apt install libeigen3-dev
  ```

GSL - GNU (only used for the batch mode): <https://www.gnu.org/software/gsl/>
  ```
  git clone https://github.com/ampl/gsl.git
  cd gsl
  ./configure
  sudo make install
  ```

OpenMP (optional): <https://www.openmp.org/>
  ```
  sudo apt install libomp-dev
  ```

### General
  ```bash
  git clone https://github.com/ImperialCollegeLondon/EventEMin.git
  cd EventEMin
  mkdir build && cd build
  cmake .. <cmake arguments>
  ```
The cmake arguments can be set as follows:
  ```
  -DEventEMin_BATCH_MODE=ON/OFF
                              Build batch mode.
                              (default: OFF)
  -DEventEMin_INCREMENTAL_MODE=ON/OFF
                              Build incremental mode.
                              (default: OFF)
  -DEventEMin_FAST_EXP=ON/OFF
                              Use fast exponentiation.
                              (default: ON)
  -DEventEMin_USE_OPENMP=ON/OFF
                              Uses the OpenMP library for parallelization.
                              (default: ON)
  ```
Lastly, ensure all environment path variables are well set, and compile everything:
  ```bash
  make
  ```

## Datasets
We provide samples in the [dataset](./dataset) folder, corresponding to some results obtained in the paper.
Most of the examples provided work with these samples, so you are not required to download any dataset.

## Examples
The source files are located in the [test](./test) directory and the binary files will be located in the [bin](./bin) directory.
We provide estimation examples per model, and the dispersion measure to be used can be chosen on the corresponding source file.
Please note that the exact entropy-based measures have quadratic complexity with the number of events and the respective examples are expected to take longer (especially if you do not use OpenMP).

### Batch Mode
For each example, two images should be seen, corresponding to the original events and the transformed events, according to the estimated parameters, accumulated on the image plane.
The status of the optimisation procedure should be displayed in the following format:
  ```
  iteration, restart iteration: score, gradient magnitude
  gradient
  parameters
  ```
In the end, the estimated parameters are also displayed.

#### 2D Translation Estimation
To run the example, on a terminal type:
  ```bash
  ./example_translation2d
  ```

#### Rotation Estimation
To run the example, on a terminal type:
  ```bash
  ./example_rotation
  ```

#### Motion Estimation in Planar Scenes
To run the example, on a terminal type:
  ```bash
  ./example_homography
  ```

#### 6-DOF in 3D
To run the example, on a terminal type:
  ```bash
  ./example_6dof
  ```

### Incremental Mode
For each example, two images should be seen, corresponding to the original events and the transformed events, according to the estimated parameters, accumulated on the image plane.
The timestamp and corresponding estimates should be displayed in the following format:
  ```
  ts: timestamp, v: motion parameter estimates
  ```
In the end, the estimated parameters are also displayed.

#### 2D Translation Estimation
To run the example, on a terminal type:
  ```bash
  ./example_incremental_translation2d
  ```

#### Rotation Estimation
To run the example, on a terminal type:
  ```bash
  ./example_incremental_rotation
  ```

#### 6-DOF in 3D
To run the example, on a terminal type:
  ```bash
  ./example_incremental_6dof
  ```

## Test Sequences
To run the sequences test, you need to download at least one sequence of the dataset provided in <http://rpg.ifi.uzh.ch/davis_data.html>.

### Batch Mode
To run the test, on a terminal type:
  ```bash
  ./example_test_sequence <path-to-events-dir> <batch-size> <path-to-estimates-saving-dir> <estimates-file-name>
  ```
The executable arguments are as follows:

- path-to-events-dir:
Path to the events' directory, following the format proposed in <http://rpg.ifi.uzh.ch/davis_data.html>, e.g. `/poster_rotation`. The events' directory must contain two files, namely, `events.txt` (list of events) and `calib.txt` (camera parameters). Please check the folders under [dataset](./dataset) for examples.
- batch-size:
Number of events of each batch, e.g. `20000`.
- path-to-estimates-saving-dir:
Path to where the estimates are to be stored, e.g. `/poster_rotation/estimates`.
- estimates-file-name:
File name of the estimates, e.g. `approx_tsallis_2`.

For example, if you downloaded the `poster_rotation` sequence and stored it under `foo` directory, by running
  ```bash
  ./example_test_sequence /foo/poster_rotation 20000 /foo/poster_rotation/estimates approx_tsallis_2
  
  ```
a file containig the estimates using the *Approx. Tsallis* measure should be created under the `/foo/poster_rotation/estimates` directory (`/estimates` directory should be created before running the command).

### Incremental Mode
To run the test, on a terminal type:
  ```bash
  ./example_incremental_test_sequence <path-to-events-dir> <number-events> <path-to-estimates-saving-dir> <estimates-file-name>
  ```
The executable arguments are as follows:

- path-to-events-dir:
Path to the events' directory, following the format proposed in <http://rpg.ifi.uzh.ch/davis_data.html>, e.g. `/poster_rotation`. The events' directory must contain two files, namely, `events.txt` (list of events) and `calib.txt` (camera parameters). Please check the folders under [dataset](./dataset) for examples.
- number-events:
Number of the most recent events to maintain, e.g. `10000`.
- path-to-estimates-saving-dir:
Path to where the estimates are to be stored, e.g. `/poster_rotation/estimates`.
- estimates-file-name:
File name of the estimates, e.g. `incremental_potential`.

For example, if you downloaded the `poster_rotation` sequence and stored it under `foo` directory, by running
  ```bash
  ./example_incremental_test_sequence /foo/poster_rotation 10000 /foo/poster_rotation/estimates incremental_potential
  
  ```
a file containig the estimates using the *Incremental Potential* measure should be created under the `/foo/poster_rotation/estimates` directory (`/estimates` directory should be created before running the command).

### Compute Errors
To compute the errors for rotational motion estimation, run the MATLAB script [sequence_error.m](./dataset/poster_rotation/sequence_error.m).

## License
The EventEMin code is licensed under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/). Commercial usage is not permitted.
