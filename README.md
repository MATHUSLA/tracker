# MATHUSLA Particle Tracker

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/bbbe2cb1269e4de68a9780534652a3d2)](https://app.codacy.com/app/MATHUSLA/tracker?utm_source=github.com&utm_medium=referral&utm_content=MATHUSLA/tracker&utm_campaign=badger)

_a time-based tracking library for particle detection experiments_

## Projects

| Name        | Description               | Location (`demo/`) | Status           |
|:-----------:|:-------------------------:|:------------------:|:----------------:|
| _Prototype_ | MATHUSLA Test Stand at P1 | `prototype/`       | **`INCOMPLETE`** |

The tracking library comes with a project directory `demo` which holds the MATHUSLA Test Stand Prototype tracking code (see `demo/prototype`). More projects to be added soon.

## Tracking Script

The tracking script is a configuration file which allows the user to specify details of the tracking algorithm. Here is a list of allowed keys and their associated value types.

| Key                           | Value Type                    | Description                                 |
|:-----------------------------:|:-----------------------------:|:-------------------------------------------:|
| `geometry-file`               | file path                     | path to `GDML` geometry file                |
| `geometry-map-file`           | file path                     | path to geometry map                        |
| `geometry-default-time-error` | positive `real`               | default time resolution of detector volumes |
| `root-data`                   | directory                     | path to ROOT data to be processed           |
| `root-keys`                   | `T, X, Y, Z` OR `T, Detector` | ROOT keys for reading in `root-data`        |
| `collapse-size`               | `dt, dx, dy, dz`              | spacetime interval for point collapse       |
| `layer-depth`                 | positive `real`               | layer size for geometry approximation       |
| `line-width`                  | positive `real`               | tolerance for line approximation            |
| `seed-size`                   | positive `integer`            | number of points per seed                   |

## Building

The tracker comes with a simple build script called `install` which allows for build customization of the tracking library.

Here is a list of useful command line options:

| Action             | Options for `./install` |
|:------------------:|:-----------------------:|
| Build Only         | none                    |
| Build and Auto Run | `--run`  (see below)    |
| Clean CMake Build  | `--cmake --clean`       |
| More Options       | `--help`                |

After a successful run of `./install`, the entire tracking library is built and all of the `demo` projects are installed. The `build` directory is created which stores all of the binaries and libraries associated with the tracker.

At the moment the tracking library only comes with one "project tracker" so running `./install --run` will automatically run the `prototype` executable from `build/demo/`. Here is a list of command line options that the _Prototype_ uses:

| Action                 | Short Options    | Long Options            |
|:----------------------:|:----------------:|:-----------------------:|
| Geometry File          | `-g <filepath>`  | `--geometry=<filepath>` |
| Detector Map           | `-m <filepath>`  | `--map=<filepath>`      |
| ROOT Data Directory    | `-d <directory>` | `--data=<directory>`    |
| Custom Tracking Script | `-s <file>`      | `--script=<file>`       |
| Quiet Mode             | `-q`             | `--quiet`               |
| Help                   | `-h`             | `--help`                |

The tracking script can provide the all of the detail above to the tracker and will eventually replace all command line arguments.
