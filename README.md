# 3D Lower Envelope Algorithms: C++ Implementation, Parallelization and Benchmarking 
## About

This repository contains the developments and results I achieved during the small GMT project (INFOMSPGMT) run from September 5th to December 1st 2023, supervised by F. Staals at Utrecht University. In his (not yet published) notes, F. Staals introduces an efficient parallel algorithm for computing 3d lower envelopes. The logic of such an algorithm is assembled with multiple sub-routines solving specific problems in a parallel fashion. This project aims to investigate the practical feasibility of such an efficient algorithm and its components and provide a working implementation in C++. In addition, we want to investigate the relationship between theoretical complexities and running times of such implementations as the input size and number of processors is increased. For more information, you can refer to the project's final report and presentation available [here](https://github.com/gianmarcopicarella/ga-spgmt-uu/blob/main/FinalProjectReport.pdf) and [here](https://github.com/gianmarcopicarella/ga-spgmt-uu/blob/main/slides_presentation.pdf).

## Getting Started

The building process has been tested on Windows 11 and Visual Studio 17 2022. It should work also for MacOS and Linux but additional steps are required. For further information, please refer to [CGAL's documentation](https://github.com/gianmarcopicarella/ga-spgmt-uu/blob/main/FinalProjectReport.pdf).

Clone the project

```bash
  $ git clone --recurse-submodules https://github.com/gianmarcopicarella/ga-spgmt-uu.git
```

Go to the project directory

```bash
  $ cd ga-spgmt-uu
```

Run the install script

```bash
  $ .\script\cgal-install-windows.ps1
```

Run the CMake script

```bash
  $ .\script\src-build-windows.ps1
```

## Tech Stack
CMake, C++, CGAL, Intel OneTBB, Catch2 and SFML
