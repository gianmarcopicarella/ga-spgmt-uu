# 3D Lower Envelope Algorithms: C++ Implementation, Parallelization and Benchmarking 
## About
This repository contains the work done for the small GMT project (**INFOMSPGMT**) I worked on from September 5th to December 1st 2023, supervised by Frank Staals at Utrecht University. The project's main goal involves designing, implementing and benchmarking parallel algorithms required for computing 3D envelopes based on Staals' theoretical notes (not yet published). For more information, you can refer to the project's final report available [here](https://github.com/gianmarcopicarella/ga-spgmt-uu/blob/main/FinalProjectReport.pdf).

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
