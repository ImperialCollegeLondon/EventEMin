# EventEMin
[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg?style=flat-square)](https://creativecommons.org/licenses/by-nc-sa/4.0/)

Code for the ECCV 2020 paper entitled [Entropy Minimisation Framework for Event-based Vision Model Estimation](http://www.ecva.net/papers/eccv_2020/papers_ECCV/papers/123500154.pdf).
The authors provide this code in the hope it will be useful for understanding the proposed method, as well as for reproducibility of the results.

For more information and more open-source software please visit the Personal Robotic Lab's website: <https://www.imperial.ac.uk/personal-robotics/>.

## Requirements
This code was tested on Ubuntu 16.04 and 18.04 distros.

### Dependencies
CMake: <https://cmake.org/download/>

	$ sudo apt install cmake cmake-curses-gui

OpenCV: <https://docs.opencv.org/trunk/d7/d9f/tutorial_linux_install.html>

Eigen3: <http://eigen.tuxfamily.org/index.php?title=Main_Page#Download>

	$ sudo apt install libeigen3-dev

GSL - GNU: <https://www.gnu.org/software/gsl/>

	$ git clone https://github.com/ampl/gsl.git
	$ cd gsl
	$ ./configure
	$ make
	$ sudo make install

OpenMP (optional): <https://www.openmp.org/>

	$ sudo apt install libomp-dev

## General

	$ git clone https://github.com/ImperialCollegeLondon/EventEMin.git
	$ cd EventEMin
	$ mkdir build && cd build
	$ cmake .. && make

Please ensure all environment path variables are well set.

## Datasets
Two datasets were used in the experiments: <http://rpg.ifi.uzh.ch/davis_data.html> and <https://daniilidis-group.github.io/mvsec/>.
We provide samples in the [dataset](./dataset) folder, corresponding to some results obtained in the paper.
All the examples provided work with these samples, so you are not required to download any dataset.

## Examples
The source files are located in the [test](./test) directory and the binary files will be located in the [bin](./bin) directory.
We provide two estimation examples per model, one corresponding to the exact entropy measures (Eqs. 7 and 9-12) and the other corresponding to the proposed efficient approximate entropy measures (Eqs. 15-19).
Please note that the exact entropy measures have quadratic complexity with the number of events and the respective examples are expected to take longer (specially if you do not allow OpenMP).

### Motion Estimation in the Image Plane
For each example, two images should be seen, corresponding to the original events and the modelled events, according to the estimated parameters, accummulated on the image plane.
Also, the status of the optimisation procedure should be displayed in the following format:

	- iteration, restart iteration: score, gradient magnitude
	- gradient
	- parameters

In the end, the estimated parameters are also displayed.

#### Optical Flow Estimation
To run examples, on a terminal type:

	$ ./example_entropy_optical_flow
	$ ./example_approx_entropy_optical_flow

#### Rotation Estimation
To run examples, on a terminal type:

	$ ./example_entropy_rotation
	$ ./example_approx_entropy_rotation

#### Motion Estimation in Planar Scenes
To run examples, on a terminal type:

	$ ./example_entropy_homography
	$ ./example_approx_entropy_homography

### Motion Estimation in 3D
No images are displayed.
Instead, a file will be generated containing the 3D modelled events according to the estimated parameters.
These events can then be visualised running the MATLAB script [visualisation.m](./dataset/indoor_flying1/visualisation.m).

#### 6-DOF in 3D
To run examples, on a terminal type:

	$ ./example_entropy_6dof
	$ ./example_approx_entropy_6dof

### Test Sequences
To run the sequences tests, you need to download at least one sequence of the dataset provided in <http://rpg.ifi.uzh.ch/davis_data.html> and then on a terminal type:

        $ ./example_test_sequence_distorted <path-to-events-dir> <batch-size> <path-to-estimates-saving-dir> <estimates-file-name>

- path-to-events-dir:
Path to the events' directory, following the format proposed in <http://rpg.ifi.uzh.ch/davis_data.html> (e.g. `/poster_rotation`).
The events' directory must contain two files, namely, events.txt (list of events) and calib.txt (camera parameters).
Please check the folders under [dataset](./dataset) for examples.
- batch-size:
Number of events of each batch (e.g. `20000`).
- path-to-estimates-saving-dir:
Path to where the estimates are to be stored (e.g. `/poster_rotation_estimates`)
- estimates-file-name:
Estimates file name (e.g. `approx_tsallis_2`).

For example, if you downloaded the `poster_rotation` sequence and stored it under `foo` directory, by running

        $ ./example_test_sequence_distorted /foo/poster_rotation 20000 /foo/poster_rotation/estimates approx_tsallis_2

a file containig the estimates using the *Approx. Tsallis* function should be created under `/foo/poster_rotation/estimates` directory (`/estimates` directory should be created before running the test).

To compute the errors, by running the MATLAB script [sequence_error.m](./dataset/poster_rotation/sequence_error.m), you should get the following in deg/s (approximately), using the *Approx. Tsallis* function:

| Sequence         | mean  | std   | RMS   |
| ---------------- | ----- | ----- | ----- |
| boxes_rotation   | 8.24  | 10.74 | 10.85 |
| dynamic_rotation | 8.07  | 11.37 | 11.48 |
| poster_rotation  | 11.19 | 14.82 | 14.91 |

We provide these results as a reference to check whether the code is running as intended.
If you have troubles reproducing these results, please do not hesitate and let me know.

## License
The EventEMin code is licensed under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/). Commercial usage is not permitted.
