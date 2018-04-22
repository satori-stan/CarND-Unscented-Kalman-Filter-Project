# Unscented Kalman Filter Project Starter Code
Self-Driving Car Engineer Nanodegree Program

In this project utilize an Unscented Kalman Filter to estimate the state of a moving object of interest with noisy lidar and radar measurements. Passing the project requires obtaining RMSE values that are lower that the tolerance outlined in the project rubric. 

## Objective

The goal of this project is to track a simulated object using the CTRV (Constant Turn Rate and Velocity Magnitude) model (which is non-linear) and an Unscented Kalman Filter, which uses a statistical approximation for non-linear equations instead of the Jacobian used by the Extended Kalman Filter.

A successful implementation requires:

* Compilation on any platform (Linux, Mac, Windows)
* Follow the Prediction -> Update cycle using the CTRV motion model equations and observations for a 2d Kalman Filter.
* Handling the first measurement appropriately to account for the fact that there is no prior prediction.
* Handling both Lidar and Radar measurements correctly.
* Obtaining RMSE values that are lower than [.09, .10, .40, .30] for the x, y positions and x and y velocities respectively.
* Performing all of the above in an efficient manner.

[image1]: ./SimulationRMSE.gif

### Here I will consider each point individually and explain how I addressed them in my implementation

#### 1. Compilation on any platform

The starter code provided for the project was already meant to compile correctly in Linux, Mac and Windows; although compilation in Windows required installing Windows 10 Bash on Ubuntu or using Docker to create a virtual environment.

Since my current platform is Windows but I have low resources, I decided to go for a native approach to avoid the installation of such resource-intensive environments. This decision was originally taken during the implementation of the first project of the term [Extended-KalmanFilter](https://github.com/satori-stan/CarND-Extended-Kalman-Filter-Project). Then, it had a non-trivial time cost, which I was able to capitalize on for this project. There is little in the repository that will give my choice away, except for a couple of lines at the end of the CMakeLists.txt file where the names of the libraries are redefined when a Windows environment is found.

The main program file (main.cpp) was modified slightly to accommodate api changes in uWS 0.14.4. This file receives a measurement package message from the simulator and passes it to the UKF class, where the information is processed.

#### 2. Follow the process of a 2d Kalman filter using the CTRV motion model

In this case, the process is implemented in the UKF::ProcessMeasurement method in ukf.cpp. In general, it is a Predict -> Update process as defined by a Kalman Filter. The first measurement processed, however, doesn't follow this flow but rather initializes the belief of where the object is. This goes on from line 81 to 110. The only processing done is transforming from polar to cartesian coordinates if the first reading is from a Lidar sensor.

In every call, the time elapsed is measured. This is a simple subtraction of the new reading's timestamp minus the last.

From the second reading on, the prediction and update steps are executed, once we have enough data to make the first prediction.

The prediction step in this case is more involved than in the Extended Kalman Filter. This is because we are using a non-linear motion model and the Unscented Kalman Filter uses a stochastic transformation (using the mean and noise). The steps for the unscented transformation should be easily followed in the code. The Prediction function starts in line ukf.cpp:142 and each step is detailed in the comments.

After the prediction step comes the update step. This is (as with the Extended Kalman Filter) split between Lidar and Radar measurements. Here, in both methods (UpdateLidar in line 238 and UpdateRadar in line 301) the steps to complete are detailed in the comments. The Kalman Filter steps can be easily followed: transform the prediction to the measurement space, calculate the difference between the prediction and the measurement and update our belief. The unscented transformation can also easily be spotted in the creation of matrices for our sigma points (lines 244 and 307) and the loops used to process them (261-270, 277-283, 325-340 and 347-357).

Finally, the Normalized Innovation Squared (NIS) is calculated and displayed to be used as a measure of the consistency of the noise of the system. This is done because the noise values of the model (not the measurements) are not provided for the project but rather assumed by the implementer (me).

#### 3. Handling of the first measurement

As was already explained in the previous section, the first measurement is never used for the prediction / update cycle since there would be no prior belief of the position of the object. It is only used to initialize said belief and then the system waits for the next measurement.

#### 4. Handling Lidar and Radar measurements correctly

The nature of Lidar and Radar measurement differences is in display in the initialization of the ProcessMeasurement function and the split Update functions for Lidar and Radar. The key difference here is the transformation from polar to cartesian coordinates (or vice versa) to make sure our calculations are consistent.

#### 5. Obtain low RMSE values

Using the simulator, RMSE values obtained after a complete run were below the maximums allowed ([.09, .10, .40, .30]).

![final frame of the simulator run displaying the RMSE values][image1]

#### 6. Efficient code

As much as possible, the code performs each operation once. For example, in FusionEKF.cpp:113-119 where the time delta is used power 4 and divided by 4 in several operations so the individual expression is calculated once and the result used later.

Other examples are the Update functions in kalman_filter.cpp that perform transpose and inverse operations only once and reuse the result as needed.

## This repository

This project involves the Term 2 Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases)

This repository includes two files that can be used to set up and install [uWebSocketIO](https://github.com/uWebSockets/uWebSockets) for either Linux or Mac systems. For windows you can use either Docker, VMware, or even [Windows 10 Bash on Ubuntu](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/) to install uWebSocketIO. Please see [this concept in the classroom](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/f758c44c-5e40-4e01-93b5-1a82aa4e044f/concepts/16cf4a78-4fc7-49e1-8621-3450ca938b77) for the required version and installation scripts.

Once the install for uWebSocketIO is complete, the main program can be built and ran by doing the following from the project top directory.

1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./UnscentedKF

Tips for setting up your environment can be found [here](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/f758c44c-5e40-4e01-93b5-1a82aa4e044f/concepts/23d376c7-0195-4276-bdf0-e02f1f3c665d)

Here is the main protocol that main.cpp uses for uWebSocketIO in communicating with the simulator.

INPUT: values provided by the simulator to the c++ program

["sensor_measurement"] => the measurement that the simulator observed (either lidar or radar)


OUTPUT: values provided by the c++ program to the simulator

["estimate_x"] <= kalman filter estimated position x
["estimate_y"] <= kalman filter estimated position y
["rmse_x"]
["rmse_y"]
["rmse_vx"]
["rmse_vy"]

---

## Other Important Dependencies

* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF` Previous versions use i/o from text files.  The current state uses i/o
from the simulator.

