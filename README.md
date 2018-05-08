# MATHUSLA Particle Tracker
_tracking algorithm for particle detector_

## Build & Run

The tracker comes with a simple build script called `install` which allows for build customization and execution of the tracking algorithm.

Here is a list of useful commands:

| Action | Options after `./install` |
|:-:|:-:|
| Build Only | none  |
| Build and Auto Run | `--run`  |
| Clean CMake Build | `--cmake --clean` |
| More Options | `--help` |

The tracking executable also comes with a set of command line arguments to specify geometry files, data files, ... etc. Here is a list of options:

| Action | Short Options after `./tracker` | Long Options after `./tracker` |
|:-:|:-:|:-:|
| Geometry File | `-g <filepath>` | `--geometry=<filepath>`  |
| Detector Map | `-m <filepath>` | `--map=<filepath>` |
| ROOT Data Directory | `-d <directory>` | `--data=<directory>` |
| Custom Tracking Script  | `-s <file>`  | `--script=<file>`  |
| Quiet Mode  | `-q` | `--quiet` |
| Help | `-h` | `--help` |

## Tracking Script

The tracking script is a configuration file which allows the user to specify details of the tracking algorithm. Here is a list of allowed keys and their associated value types.

| Key | Value Type | Description |
|:---:|:----------:|:-----------:|
| `geometry-file` | file path | path to `GDML` geometry file |
| `geometry-map`  | file path | path to geometry map |
| `root-data`     | directory | path to ROOT data to be processed |
| `root-keys`     | `T-key, X-key, Y-key, Z-key` OR `T-key, Detector-key` | ROOT keys for reading from files in `root-data` |
| `collapse-size` | `dt, dx, dy, dz` | spacetime interval for point collapse |
| `layer-depth`   | positive `real` | layer size for geometry approximation |
| `line-width`    | positive `real` | tolerance for line approximation |
| `seed-size`     | positive `integer` | number of points per seed |

The tracking script also includes the arguments for file paths so it replaces the  `-gmd` options defined in the __Build & Run__ section.