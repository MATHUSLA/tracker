# MATHUSLA Particle Tracker

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/bbbe2cb1269e4de68a9780534652a3d2)](https://app.codacy.com/app/MATHUSLA/tracker?utm_source=github.com&utm_medium=referral&utm_content=MATHUSLA/tracker&utm_campaign=badger)

_tracking algorithm for particle detector_

## Build and Run

The tracker comes with a simple build script called `install` which allows for build customization and execution of the tracking algorithm.

Here is a list of useful commands:

| Action             | Options for `./install` |
|:------------------:|:-----------------------:|
| Build Only         | none                    |
| Build and Auto Run | `--run`                 |
| Clean CMake Build  | `--cmake --clean`       |
| More Options       | `--help`                |

The tracking executable also comes with a set of command line arguments to specify geometry files, data files, ... etc. Here is a list of options:

| Action                 | Short Options    | Long Options            |
|:----------------------:|:----------------:|:-----------------------:|
| Geometry File          | `-g <filepath>`  | `--geometry=<filepath>` |
| Detector Map           | `-m <filepath>`  | `--map=<filepath>`      |
| ROOT Data Directory    | `-d <directory>` | `--data=<directory>`    |
| Custom Tracking Script | `-s <file>`      | `--script=<file>`       |
| Quiet Mode             | `-q`             | `--quiet`               |
| Help                   | `-h`             | `--help`                |

## Tracking Script

The tracking script is a configuration file which allows the user to specify details of the tracking algorithm. Here is a list of allowed keys and their associated value types.

| Key             | Value Type                              | Description                                     |
|:---------------:|:---------------------------------------:|:-----------------------------------------------:|
| `geometry-file` | file path                               | path to `GDML` geometry file                    |
| `geometry-map`  | file path                               | path to geometry map                            |
| `root-data`     | directory                               | path to ROOT data to be processed               |
| `root-keys`     | `T, X, Y, Z` keys OR `T, Detector` keys | ROOT keys for reading from files in `root-data` |
| `collapse-size` | `dt, dx, dy, dz`                        | spacetime interval for point collapse           |
| `layer-depth`   | positive `real`                         | layer size for geometry approximation           |
| `line-width`    | positive `real`                         | tolerance for line approximation                |
| `seed-size`     | positive `integer`                      | number of points per seed                       |

The tracking script also includes the arguments for file paths so it replaces the  `-gmd` options defined in [__Build and Run__](#build-and-run).
