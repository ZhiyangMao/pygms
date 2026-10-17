# PyGMS: Python Ground Motion Selection

## Overview

PyGMS is a Python package for selecting and scaling ground motion records to match target response spectra. The package implements various target spectra including Conditional Spectrum (CS), Conditional Mean Spectrum (CMS), and CS with duration constraints (CS+duration). 

The implementation provides a user-friendly interface through JSON configuration files.

## Features

- Support for multiple target intensity measure types:
  - Conditional Spectrum (CS)
  - Conditional Mean Spectrum (CMS)
  - Conditional Spectrum plus duration constraints (CS+duration)
- Integration with USGS hazard API for ground motion prediction equations
- Built-in significant duration models
- Includes dataset of NGA-West ground motions (Sa from 0.01s to 10s, Ds5-75, and Ds5-95)
- Optional per-earthquake penalty to limit the number of records selected from the same event

## Structure

The package consists of:

- **pygms.py**: Main script that interfaces with user-defined JSON configuration files
- **Classes**:
  - `TargetIntensityMeasure`: Processes target spectra configuration
  - `GroundMotionSelection`: Handles the ground motion selection algorithm
- **Modules**:
  - `USGSHazardGMM`: Fetches USGS API for different ground motion prediction equations
  - `SignificantDurationModel`: Implements significant duration models
  - `CorrelationModel`: Implements intensity-measure correlation models

## Usage

### Basic Usage

1. Create a JSON configuration file (see examples below)
2. Run the selection algorithm using the provided Jupyter notebook [examples.ipynb](https://github.com/bakerjw/pygms/blob/main/examples.ipynb)  or directly with the pygms.py Python script:

```python
import pygms

# Load configuration from JSON file
pygms.run_selection("example1.json")

# Optionally specify the output file (or directory)
pygms.run_selection("example1.json", "results/example1.csv")
```

### Limiting Records per Earthquake

Large, well-recorded events (e.g., Chi-Chi) can dominate a selected set. An optional `EventPenalty` section in `GroundMotionSelection` limits the number of records taken from the same earthquake:

```json
"EventPenalty": {
  "MaxPerEvent": 5,
  "Weight": 1.0,
  "Exponent": 2,
  "Mode": "soft+hard"
}
```

- `MaxPerEvent`: maximum number of records per earthquake (EQID). The penalty is disabled when `EventPenalty` is absent.
- `SoftCap`: number of records per earthquake above which the soft penalty starts (default `MaxPerEvent`). With the default, up to `MaxPerEvent` records per earthquake are allowed without penalty. Set `SoftCap` below `MaxPerEvent` to further discourage earthquakes from approaching the cap.
- `Weight`: penalty weight. The penalty is scaled by the initial fitting error, so values of order 1 work across sites and targets.
- `Exponent`: the penalty for each earthquake is `max(0, n - SoftCap) ** Exponent`, where `n` is its number of selected records.
- `Mode`: `"soft+hard"` (default), `"soft"` (penalty only), or `"hard"` (cap only).

Earthquakes are identified by the NGA-West2 `EQID` from `gmdb/nga_west_event.csv` (default database) or `gmdb/nga_west2/event.csv` (`"NGA-West2"`), built from the NGA-West2 flatfile with `tools/build_event_table.py`. When the penalty is enabled, the output file includes `EQID` and `EarthquakeName` columns and a summary of records per earthquake is printed.

## Configuration Options

Configuration is handled through a JSON file with the following main sections:

1. **Target Intensity Measure**: Define the target spectrum type, periods, and conditional intensity measure
2. **Ground Motion Models**: Select models and parameters (magnitude, distance, Vs30)
3. **Selection Parameters**: Define scaling approach, number of ground motions, and error weights

## Examples

The package includes three example configurations:

1. **example1.json**: CS+duration target
2. **example2.json**: CS only target
3. **example3.json**: CMS only target
4. **example_event_penalty.json**: CS+duration target with at most 5 records per earthquake

## Selection Algorithm

The ground motion selection algorithm follows these steps:

1. Sample a set of pseudo spectra from the target distribution
2. Perform initial selection to find the best match to individual pseudo spectra (minimizing weighted total error)
3. Apply a greedy algorithm to optimize the weighted total error for:
   - Mean and standard deviation of spectral acceleration (Sa)
   - Mean and standard deviation of significant duration (for CS+duration targets)

## Output

The default output is a CSV file containing:
1. Ground motion filenames (RSN identifiers)
2. Unscaled spectra values
3. Significant duration values
4. Scaling factors

## References

Baker, J. W., and Lee, C. (2018). “[An Improved Algorithm for Selecting Ground Motions to Match a Conditional Spectrum](https://www.jackwbaker.com/Publications/Baker_Lee_(2018)_GM_Selection,_JEE.pdf).” Journal of Earthquake Engineering, 22(4), 708–723. This algorithm is also implemented in [Matlab](https://github.com/bakerjw/CS_Selection).


## License

[Add license information]

## Contributors

- Kuanshi Zhong (Implementation in Python)
- Jack Baker (Development of a previous implementation in Matlab)
- Marguerite Dubertret (Performance and functional improvements)
- Zhiyang Mao (Per-earthquake record penalty)

## Contact

- Kuanshi Zhong [zhongki@ucmail.uc.edu]
- Jack Baker [bakerjw@stanford.edu]

## TODO
- Integrating site-specific seismic disaggregation info
- Multi-threading/processing function
