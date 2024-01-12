# Visual Plumes Model Code

## Code Repositories for Visual Plumes:

- https://github.com/sfei/Visual-Plumes-Back-End
- https://github.com/sfei/Visual-Plumes-Front-End
- https://github.com/sfei/Visual-Plumes-Models

## Introduction

This repository represents the core model code that runs the Visual Plumes simulation. It handles both the plume (currently only UM3) and far-field simulations.

## Folder and file description

Note however that this code is meant to be imported as a module into a larger project, wherein all the exposed classes and objects are provided in a flat organization from `\__init__.py`.

### Root

- `\.`
  - `\__init__.py` - Module init defining classes and objects exposed in the main module.
  - `\all_releastic_outputs.py` - For example only. Use as reference for adding parameters to output handler.
  - `\Element.py` - Simple class for holding a set a values related to control-volume element.
  - `\globals.py` - Commonly used, generic variables and classes.
  - `\GraphOutputs.py` - Class for tracking coordinates to be used as data for graph outputs.
  - `\helpers.py` - Useful helper functions.
  - `\Middleware.py` - Middleware and entry point function for running the entirety of the model.
  - `\Output.py` - Class for storing data representing tabular outputs.
  - `\UMUnit.py` - Class representing the UM3 model which simulates the jet/plume.
  - `\units.py` - Classes representing supported units of measure and conversations between.
  - `\vectors.py` - Commonly used vector functions.

### Ambient

- `\ambient\`
  - `\ambient\Ambient.py` - Simple class for holding a set of ambient values.
  - `\ambient\AmbientHandler.py` - Handler for filling and interpolating a stack of ambient values representing a depth profile.
  - `\ambient\AmbientStore.py` - Metadata for ambient values (e.g. units, timeseries info, extrapolation).
  - `\ambient\calcuations.py` - Ambient calculation functions.

### Farfield

- `\farfield\`
  - `\farfield\BrooksFarField.py` - Classes and functions for handling the Brooks Far Field model.
  - `\farfield\TidalPollutionBuildup.py` - Classes and functions for handling the Tidal Pollution Buildup model.

### Parameters

- `\params\`
  - `\params\DiffuserParameters.py` - Simple class for holding a set of diffuser values.
  - `\params\DiffuserStore.py` - Metadata for ambient values (e.g. units, timeseries info).
  - `\params\ModelParameters.py` - Simple class for holding model parameter values.

### Timeseries

- `\timeseries\`
  - `\timeseries\AbstractTimeseries.py` - Abstract class with commonly used functions to inherit for timeseries handlers.
  - `\timeseries\AmbientTimeseries.py` - Timeseries handler for ambient timeseries data.
  - `\timeseries\DiffuserTimeseries.py` - Timeseries handler for diffuser timeseries data.
  - `\timeseries\TimeseriesHandler.py` - Handler class for syncing and incrementing a set of timeseries data.
  - `\timeseries\TimeseriesStore.py` - Start time, end time, and time increment metadata store for timeseries data.

## Middleware and entry-point

The main entry point for running the model is in `Middleware.py` as the function `run()`. This function handles all setup, coordination of objects, and transitions to far-field models. Currently runs UM3 plume model and optionally a far-field model. Far-field may be either Brooks far-field or tidal pollution buildup, but not both.

### Arguments

- `model_params`: The model parameters, as set in an object instance of the class `ModelParameters`.
- `diffuser_params`: The diffuser parameters, as set in an object instance of the class `DiffuserParameters`.
- `diffuser_store`: The metadata for the diffuser parameters, as set in an object instance of the class `DiffuserStore`.
- `timeseries_handler`: Optional. Set as `None` or empty for no timeseries data. An instance of the class `TimeseriesHandler` to signify the presence of timeseries data to use in the model run.
- `ambient_stack`: A stack of ambient conditions in order of increasing depth. Defined as a list of `Ambient` objects.
- `ambient_store`: The metadata for the diffuser parameters, as set in an object instance of the class `AmbientStore`.
- `output_handler`: Optional. Set as `None` or empty to use default output configuration. If set (as a instance of the class `OutputHandler`), use presupplied handler for tracking tabular outputs.

## Setup and example

At minimum, initialize the objects of classes `ModelParameters`, `DiffuserParameters`, `DiffuserStore`, and the necessary ambient conditions (detailed in next subsection). Fill in and modify values as appropriate (or leave in defaults values).

```python
# adjust as needed based on how module was setup in your project
from visualplumes import units, Middleware, OutputUM3, Ambient, AmbientStore, ModelParameters, BacteriaModel, \
                         SimilarityProfile, FarfieldDiffusivity, MaxVerticalReversals, DiffuserParameters, \
                         DiffuserStore, AmbientTimeseries, DiffuserTimeseries, TimeseriesHandler


model_params = ModelParameters()
diff_params  = DiffuserParameters()
diff_store   = DiffuserStore()

# modify some model parameters
model_params.write_step_freq = 100
model_params.max_reversals   = MaxVerticalReversals.SECOND_MAX_RISE_OR_FALL
model_params.max_dilution    = 10000

# setup tidal pollution build model to run
model_params.tidal_pollution_buildup  = True
model_params.tpb_channel_width        = 39.6  # all units assumed meters or degrees
model_params.tpb_segment_length       = 80.46
model_params.tpb_upstream_dir         = 180
model_params.tpb_coast_bin            = 65
model_params.tpb_coast_concentration  = 9e-6
model_params.tpb_mixing_zone_ceil     = 3.657

# set up diffuser parameters and define units
diff_params.diameter                 = 0.0254
diff_store.diameter.units            = units.Length.METERS
diff_params.offset_x                 = 0.0
diff_store.offset_x.units            = units.Length.FEET
diff_params.offset_y                 = 0.0
diff_store.offset_y.units            = units.Length.FEET
diff_params.vertical_angle           = 22.5
diff_store.vertical_angle.units      = units.Angle.DEGREES
diff_params.horizontal_angle         = 270.0
diff_store.horizontal_angle.units    = units.Angle.DEGREES
diff_params.num_ports                = 1
diff_store.num_ports.units           = units.Unitless.UNITLESS
diff_params.acute_mixing_zone        = 30
diff_store.acute_mixing_zone.units   = units.Length.FEET
diff_params.isopleth                 = 60.0
diff_store.isopleth.units            = units.Isopleth.CONCENTRATION
diff_params.depth                    = 5.0
diff_store.depth.units               = units.Length.FEET
diff_params.salinity                 = 15.0
diff_store.salinity.units            = units.Salinity.PRACTICAL_SALINITY_UNITS
diff_params.temperature              = 25.0
diff_store.temperature.units         = units.Temperature.CELSIUS
diff_params.effluent_flow            = 0.021087
diff_store.effluent_flow.units       = units.FlowRate.MEGAGALLONS_PER_DAY
diff_params.concentration            = 8300.0
diff_store.concentration.units       = units.Concentration.PARTS_PER_MILLION
```

### Ambient

Ambient data is supplied as a stack of `Ambient[]` objects, in order of increasing depth. Not all values have to be filled in the stack (though at least one for each value/parameter type), and the unfilled values will be interpolated and extrapolated as needed.

```python
# start by defining units
ambient_store = AmbientStore()
ambient_store.z.units             = units.Length.FEET
ambient_store.current_speed.units = units.Speed.FEET_PER_SECOND
ambient_store.current_dir.units   = units.Angle.DEGREES
ambient_store.salinity.units      = units.Salinity.PRACTICAL_SALINITY_UNITS
ambient_store.temperature.units   = units.Temperature.CELSIUS
ambient_store.bg_conc.units       = units.Concentration.PARTS_PER_MILLION
ambient_store.ff_diff_coeff.units = units.EddyDiffusivity.DIFFUSIVITY

# first layer, here representing surface conditions
ambient_depth0 = Ambient()
ambient_depth0.z             = 0
ambient_depth0.current_speed = 0
ambient_depth0.current_dir   = 0
ambient_depth0.salinity      = 0.2
ambient_depth0.temperature   = 25
ambient_depth0.bg_conc       = 9
ambient_depth0.decay_rate    = 0
ambient_depth0.ff_diff_coeff = 0.0003

# second and only other defined layer (at 12ft depth)
# here, only the temperature is defined as changing by depth
# unfilled values we be interpolated (using default param of Interpolation.Constant) from surface layer
ambient_depth12 = Ambient()
ambient_depth12.z            = 12
ambient_depth12.temperature  = 22.2

# the ambient stack
ambient_stack = (ambient_depth0, ambient_depth12)
```

### Timeseries

Timeseries data must be defined using a created `TimeseriesHandler`.

```python
timeseries = TimeseriesHandler()
timeseries.start_time           = 0  # inclusive
timeseries.end_time             = 3  # inclusive
timeseries.time_increment       = 1  # will result in 4 cases to run from 0-3 inclusive
timeseries.units.start_time     = units.Time.HOURS
timeseries.units.end_time       = units.Time.HOURS
timeseries.units.time_increment = units.Time.HOURS
```

Timeseries data from CSV files are fed into the handler using the appropriate timeseries data handling class.

```python
diff_store.effluent_flow.units        = units.FlowRate.MEGAGALLONS_PER_DAY
diff_store.effluent_flow.ts_increment = 1.0
timeseries.diffuser.effluent_flow     = DiffuserTimeseries("./TRwtp_tsfiles/flowrate.csv", diff_store.effluent_flow)

ambient_store.current_dir.z_is_depth       = True  # timeseries may be depth/height layers indepedently
ambient_store.current_dir.ts_depth_units   = units.Length.FEET
ambient_store.current_dir.units            = units.Angle.DEGREES
ambient_store.current_dir.ts_increment     = 1.0  # this is assumed as hours
timeseries.ambient.current_dir             = AmbientTimeseries("./TRwtp_tsfiles/current_dir.csv", ambient_store.current_dir)

ambient_store.current_speed.z_is_depth     = True
ambient_store.current_speed.ts_depth_units = units.Length.FEET
ambient_store.current_speed.units          = units.Speed.FEET_PER_SECOND
ambient_store.current_speed.ts_increment   = 1.0
timeseries.ambient.current_speed           = AmbientTimeseries("./TRwtp_tsfiles/current_speed.csv", ambient_store.current_speed)
```

### Output handler

The output handler is optional, but if wanting to manually configure, you may create your own handler and feed into the middleware function. The handler simply needs to be created and then parameters added for tracking. Knowing the proper regime and parameter names may not be obvious, so use the file `all_releastic_outputs.py` for examples of most likely parameters.

If not supplied, a default list of parameters to track is set up.

```python
output_handler = OutputUM3()
# add parameters by defining the regime, parameter name, label, unit class, and unit type
output_handler.add_parameter('element', 'depth',          'Depth',     units.Length,        diff_store.depth.units)
output_handler.add_parameter('element', 'diameter',       'Width',     units.Length,        diff_store.diameter.units)
output_handler.add_parameter('element', 'vertical_angle', 'V-angle',   units.Angle,         diff_store.vertical_angle.units)
output_handler.add_parameter('element', 'concentration',  'Pollutant', units.Concentration, diff_store.concentration.units)
output_handler.add_parameter('element', 'density',        'Density',   units.Density,       units.Density.SIGMA_T)
output_handler.add_parameter('ambient', 'density',        'Amb-den',   units.Density,       units.Density.SIGMA_T)
output_handler.add_parameter('ambient', 'current_speed',  'Amb-cur',   units.Speed,         ambient_store.current_speed.units)
output_handler.add_parameter('element', 'speed',          'Velocity',  units.Speed,         units.Speed.METERS_PER_SECOND)
output_handler.add_parameter('element', 'dilution',       'Dilution',  units.Unitless,      units.Unitless.UNITLESS)
```

### Running the model

Pass your parameters to the middleware function.

```python
output_dict = Middleware.run(
    model_params=model_params,
    diffuser_params=diff_params,
    diffuser_store=diff_store,
    timeseries_handler=timeseries,
    ambient_stack=ambient_stack,
    ambient_store=ambient_store,
    output_handler=output_handler
)
```

### Output dictionary

The output dictionary (in the above example assigned to `output_dict`) is a large Python dictionary containing all the output information. Many output variables are duplicated by time-slice case run and thus wrapped a list or tuple. Even if no time-slices are defined by the timeseries data, the output may be wrapped a in list or tuple of length one.

#### General outputs

These parameters provided basic information on the model run itself and whether it ran successfully or errored.

- `success` (`bool`): Whether the run successfully completed.
- `error` (`Exception`): The `AssertionError` or `UserInputError` thrown or None if run was successful. Note that other types of errors will not be caught and instead thrown.
- `errortraceback` (`str`): The error traceback message or None if run was successful.
- `cases` (`int`): Number of cases run.
- `casetime` (`float[]`): List of model time (relevant for timeseries data, if supplied), in seconds, per case.

#### Parameters and inputs

These parameters show the input parameters. For timeseries data, can show how diffuser and ambient conditions have changed for each time-slice case run.

- `modelparams` (`dict`): Simple dictionary for model parameters information.
  - `memos` (`str[]`): List of print memos detailing model parameters for output text.
- `timeseries` (`dict`): Simple dictionary for timeseries indices by case.
  - `memos` (`str[]`): List of print memos detailing timeseries indices by case.
- `diffuser` (`dict`): Diffuser values.
  - `headers` (`dict[]`): List of header information dicts.
    - `name` (`str`) : The variable name.
    - `label` (`str`) : The formatted/pretty name for this parameter.
    - `units` (`Unit`) : Subclass of `Unit` (from units.py) with unit type (e.g. `units.Length`).
    - `in_units` (`int`) : The specific unit (a numeric from Unit type member variable) (e.g. `units.Length.METERS`).
    - `units_label` (`str`) : The formatted label for the specific unit given by in_units.
  - `outputs` (`dict[][]`) : List of diffuser values for each sequential, time-slice case. Provided as 2-dimension list by time-slice then values in same order as headers.
- `ambient` (`dict`): Ambient values.
  - `headers` (`dict[]`): List of header information dicts.
    - `name` (`str`) : The variable name.
    - `label` (`str`) : The formatted/pretty name for this parameter.
    - `units` (`Unit`) : Subclass of `Unit` (from units.py) with unit type (e.g. `units.Length`).
    - `in_units` (`int`) : The specific unit (a numeric from Unit type member variable) (e.g. `units.Length.METERS`).
    - `units_label` (`str`) : The formatted label for the specific unit given by in_units.
  - `outputs` (`dict[][][]`) : List of ambient values for each sequential, time-slice case. Provided as 3-dimension list by time-slice then ambient stack (in order of increasing depth) then values in same order as headers.

#### Plume model outputs

Outputs are a list of values to be tracked per timestep as defined by the `OutputHandler` supplied or a set of default parameters if not. The `headers` values defined the parameters and the `outputs` are a list of these values in the same order, per time-slice case run, per timestep increment where values were tracked.

E.g. `output_dict['plume']['outputs'][1][3]` calls the 4th row of tracked values (using index=3) of the 2nd time-slice case run (using index=1). This returns a dictionary with keys `step`, `values`, and `status` with data from this snapshot.

- `plume` (`dict`): UM3 model outputs.
  - `headers` (`dict[]`): List of header information dicts.
    - `regime` (`str`) : The regime name (`ambient`, `element`, `diffuser`, or `model`) for the variable to be tracked.
    - `name` (`str`) : The variable name.
    - `label` (`str`) : The formatted/pretty name for this parameter.
    - `units` (`Unit`) : Subclass of Unit (from units.py) with unit type (e.g. `units.Length`).
    - `in_units` (`int`) : The specific unit (a numeric from Unit type member variable) (e.g. `units.Length.METERS`).
    - `units_label` (`str`) : The formatted label for the specific unit given by in_units.
  - `outputs` (`dict[]`): List of outputs for each sequential, time-slice case.
    - `step` (`int`) : Timestep number.
    - `values` (`float[]`) : List of values for tracked params at that step (in order params were added to outputs).
    - `status` (`str`): : Any notable status update/changes (given as comma separated string).
  - `memos` (`str[][]`): List of print memos that occurred during model run for each case. (`memo()` in o.c.)

#### Far-field model outputs

Far-field data returns different by the time of far-field model returned. The Brooks Far Field is run per time-slice case run individually, while the Tidal Pollution Buildup model runs only once, using each time-slice case run to simulate the next step and update its status.

- `farfield` (`dict`): Farfield model outputs (only filled for Brooks far-field model).
  - `was_run` (`bool`): True if Brooks far-field model was run.
  - `headers` (`dict[]`): List of header information dicts.
    - `name` (`str`) : The variable name.
    - `label` (`str`) : The formatted/pretty name for this parameter.
    - `units` (`Unit`) : Subclass of Unit (from units.py) with unit type (e.g. `units.Length`).
    - `in_units` (`int`) : The specific unit (a numeric from Unit type member variable) (e.g. `units.Length.METERS`).
    - `units_label` (`str`) : The formatted label for the specific unit given by in_units.
  - `outputs` (`dict[]`): List of outputs for each sequential, time-slice case.
    - `step` (`int`) : Timestep number.
    - `values` (`float[]`) : List of values for tracked params at that step (in order params were added to outputs).
  - `memos` (`str[][]`): List of print memos that occurred during model run for each case. (`memo()` in o.c.)
- `tpb` (`dict`): Tidal pollution buildup model outputs (only filled for tidal pollution buildup model).
  - `was_run` (`bool`): True if tidal pollution buildup model was run.
  - `memos` (`str[][]`): List of print memos that occurred during model run. (`memo()` in o.c.)

#### Graph outputs

The graph data is provided per graph series, which are predefined and cannot currently be customized. Each series is provided as a list of x/y coordinates and meta data about the units for each axis. A null coordinate of `[None, None]` is used to signify a break and new line in the same series.

- `graphs` (`dict`): Dictionary of series names to series data.
  - `coords` (`float[][]`): List of graph coordinates for this series. (All `NaN`s converted to `None`).
  - `units` (`dict`): Dictionary of units information for axes of this series.
    - `x` (`dict`): Dictionary of units information for the x-axis of this series.
      - `units` (`Unit`) : Subclass of Unit (from units.py) with unit type (e.g. `units.Length`).
      - `in_units` (`int`) : The specific unit (a numeric from Unit type member variable) (e.g. `units.Length.METERS`).
      - `units_label` (`str`) : The formatted label for the specific unit given by `in_units`.
      - `axis_label` (`str`) : The axis label.
    - `y` (`dict`): Dictionary of units information for the x-axis of this series.
      - `units` (`Unit`) : Subclass of Unit (from units.py) with unit type (e.g. `units.Length`).
      - `in_units` (`int`) : The specific unit (a numeric from Unit type member variable) (e.g. `units.Length.METERS`).
      - `units_label` (`str`) : The formatted label for the specific unit given by `in_units`.
      - `axis_label` (`str`) : The axis label.
