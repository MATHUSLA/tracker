# MATHUSLA Particle Tracker

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/bbbe2cb1269e4de68a9780534652a3d2)](https://app.codacy.com/app/MATHUSLA/tracker?utm_source=github.com&utm_medium=referral&utm_content=MATHUSLA/tracker&utm_campaign=badger)

_a time-based tracking library for particle physics experiments_

## Projects

| Name        | Description                    | Location (`demo/`) | Status           |
|:-----------:|:------------------------------:|:------------------:|:----------------:|
| _Prototype_ | MATHUSLA Test Stand at P1      | `prototype/`       | **`FUNCTIONAL`** |
| _Box_       | MATHUSLA Full Box              | `box/`             | **`INCOMPLETE`** |
| _ModuleBox_ | MATHUSLA Full Box with Modules | `module_box/`      | **`INCOMPLETE`** |

The tracking library comes with a project directory `demo` which holds the MATHUSLA Test Stand Prototype tracking code (see `demo/prototype`). More projects to be added soon.

## Tracking Script

The tracking script is a configuration file which allows the user to specify details of the tracking algorithm. Here is a list of allowed keys and their associated value types.

| Key                           | Value Type                    | Description                                 |
|:-----------------------------:|:-----------------------------:|:-------------------------------------------:|
| `verbose-output`              | boolean                       | type of output to terminal                  |
| `draw-events`                 | boolean                       | type of event display mode                  |
| `geometry-file`               | file path                     | path to `GDML` geometry file                |
| `geometry-map-file`           | file path                     | path to geometry map                        |
| `geometry-default-time-error` | positive `real`               | default time resolution of detector volumes |
| `statistics-directory`        | directory path                | path to output statistics                   |
| `statistics-file-prefix`      | prefix string                 | prefix name for output file                 |
| `statistics-file-extension`   | extension string              | extension to file (or filetype)             |
| `data-directory`              | directory path                | path to ROOT data to be processed           |
| `data-file-extension`         | extension string              | extension to file (or filetype)             |
| `data-position-keys`          | string for each of R4         | data key read-in for position               |
| `data-position-error-keys`    | string for each of R4         | data key read-in for position error         |
| `data-detector-key`           | string                        | data key read-in for detector               |
| `data-track-id-key`           | string                        | data key read-in for track id               |
| `data-parent-id-key`          | string                        | data key read-in for parent id              |
| `data-momentum-keys`          | string for each of R4         | data key read-in for momentum               |
| `time-smearing`               | boolean                       | mode for smearing time input on resolution  |
| `layer-axis`                  | R3 coordinate                 | direction for track parameterization        |
| `layer-depth`                 | positive `real`               | layer size for geometry approximation       |
| `line-width`                  | positive `real`               | tolerance for line approximation            |
| `seed-size`                   | positive `integer`            | number of points per seed                   |
| `event-density-limit`         | positive `real < 1`           | density limit before event is dropped       |

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

At the moment the tracking library only comes with one complete "project tracker" so running `./install --run` will automatically run the `prototype` executable from `build/demo/`. Here is a list of command line options that the _Prototype_ uses:

| Action                 | Short Options    | Long Options            |
|:----------------------:|:----------------:|:-----------------------:|
| Geometry File          | `-g <filepath>`  | `--geometry=<filepath>` |
| Detector Map           | `-m <filepath>`  | `--map=<filepath>`      |
| ROOT Data Directory    | `-d <directory>` | `--data=<directory>`    |
| Custom Tracking Script | `-s <file>`      | `--script=<file>`       |
| Quiet Mode             | `-q`             | `--quiet`               |
| Help                   | `-h`             | `--help`                |

The tracking script can provide the all of the detail above to the tracker and will eventually replace all command line arguments.
