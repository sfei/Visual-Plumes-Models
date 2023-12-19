from visualplumes import units, OutputUM3, AmbientStore, DiffuserStore, ModelParameters, FarfieldDiffusivity


# assumed objects created in preparing scenario run
model_params   = ModelParameters()
diff_store     = DiffuserStore()
ambient_store  = AmbientStore()
output_handler = OutputUM3()

# Note that eddy diffusivity name depends on model parameters -- wonder if better to just say "Eddy Diffusivity" in the
# setup UI though and change the column name as needed in outputs later.
if model_params.farfield_diffusivity == FarfieldDiffusivity.POWER_4_3:
    diff_col_name = "4/3 eddy diffusivity"
else:
    diff_col_name = "Eddy diffusivity"

# Diffuser/case related -- TODO: seems kind of redundant considering outputs also include the diffuser table which is static per case
#--------------------------- regime ---- parameter name ---- pretty/formatted label --- unit class ------------ unit type -----------
output_handler.add_parameter('diffuser', 'diameter',         'Port diameter',           units.Length,           diff_store.diameter.units)
output_handler.add_parameter('diffuser', 'num_ports',        'Num. ports',              units.Unitless,         units.Unitless.UNITLESS)
output_handler.add_parameter('diffuser', 'port_spacing',     'Port spacing ',           units.Length,           diff_store.port_spacing.units)
output_handler.add_parameter('diffuser', 'start_time',       'Start time',              units.Time,             diff_store.start_time.units)
output_handler.add_parameter('diffuser', 'end_time',         'End time',                units.Time,             diff_store.end_time.units)
output_handler.add_parameter('diffuser', 'time_increment',   'Time increment',          units.Time,             diff_store.time_increment.units)
output_handler.add_parameter('diffuser', 'acute_mixing_zone','Mixing zone distance',    units.Length,           diff_store.acute_mixing_zone.units)
output_handler.add_parameter('diffuser', 'effluent_flow',    'Total flow',              units.FlowRate,         diff_store.effluent_flow.units)

# Ambient conditions related
#--------------------------- regime ---- parameter name ---- pretty/formatted label --- unit class ------------ unit type -----------
output_handler.add_parameter('ambient',  'density',          'Ambient density',         units.Density,          units.Density.SIGMA_T)
output_handler.add_parameter('ambient',  'current_speed',    'Current speed',           units.Speed,            ambient_store.current_speed.units)
output_handler.add_parameter('ambient',  'current_dir',      'Current direction',       units.Angle,            ambient_store.current_dir.units)
output_handler.add_parameter('ambient',  'salinity',         'Ambient salinity',        units.Salinity,         ambient_store.salinity.units)
output_handler.add_parameter('ambient',  'temperature',      'Ambient temperature',     units.Temperature,      ambient_store.temperature.units)
output_handler.add_parameter('ambient',  'bg_conc',          'Background concentration',units.Concentration,    ambient_store.bg_conc.units)
output_handler.add_parameter('ambient',  'kt',               'Decay rate',              units.DecayRate,        units.DecayRate.PER_SECOND)  # decay rate output always converted
output_handler.add_parameter('ambient',  'ff_velocity',      'Far-field speed',         units.Speed,            ambient_store.ff_velocity.units)
output_handler.add_parameter('ambient',  'ff_dir',           'Far-field direction',     units.Angle,            ambient_store.ff_dir.units)
output_handler.add_parameter('ambient',  'ff_diff_coeff',    'Dispersion',              units.EddyDiffusivity,  units.EddyDiffusivity.DIFFUSIVITY)
# TODO: below doesn't seem to work in VisualPlumes (selecting 4/3Eddy in special settings doesn't actually seem to output anything new in table)
# output_handler.add_parameter('ambient',  'diffusivity',      diff_col_name,             units.EddyDiffusivity,  units.EddyDiffusivity.DIFFUSIVITY)

# Plume/element or general model related
#--------------------------- regime ---- parameter name ---- pretty/formatted label --- unit class ------------ unit type -----------
output_handler.add_parameter('element',  'total_time',       'Time',                    units.Time,             units.Time.SECONDS)
output_handler.add_parameter('element',  'depth',            'Depth',                   units.Length,           diff_store.depth.units)
output_handler.add_parameter('element',  'diameter',         'Diameter',                units.Length,           diff_store.diameter.units)
output_handler.add_parameter('element',  'height',           'Height',                  units.Length,           diff_store.diameter.units)
output_handler.add_parameter('element',  'vertical_angle',   'Vertical angle',          units.Angle,            diff_store.vertical_angle.units)
output_handler.add_parameter('element',  'horizontal_angle', 'Horizontal angle',        units.Angle,            diff_store.horizontal_angle.units)
output_handler.add_parameter('element',  'salinity',         'Salinity',                units.Salinity,         ambient_store.salinity.units)
output_handler.add_parameter('element',  'temperature',      'Temperature',             units.Temperature,      ambient_store.temperature.units)
output_handler.add_parameter('element',  'concentration',    'Pollutant',               units.Concentration,    diff_store.concentration.units)
output_handler.add_parameter('element',  'density',          'Density',                 units.Density,          units.Density.SIGMA_T)
output_handler.add_parameter('element',  'dilution',         'Dilution',                units.Unitless,         units.Unitless.UNITLESS)
output_handler.add_parameter('element',  'cl_dilution',      'CL-dilution',             units.Unitless,         units.Unitless.UNITLESS)
output_handler.add_parameter('model',    'net_dilution',     'Net dilution',            units.Unitless,         units.Unitless.UNITLESS)
output_handler.add_parameter('element',  'speed',            'Velocity',                units.Speed,            ambient_store.current_speed.units)
# Total distance traveled (e.g. if it loops back, it would still be positive value)
output_handler.add_parameter('model',    'travel',           'Distance travel',         units.Length,           diff_store.offset_x.units)
# Displacement from source (e.g. if it loops back to source, would be zero) -- technically surface distance only, depth not factored
# TODO: adjust label? remove this option?
output_handler.add_parameter('element',  'total_surf_dsp',   'Distance from source',    units.Length,           diff_store.offset_x.units)
output_handler.add_parameter('element',  'x_displacement',   'X-position',              units.Length,           diff_store.offset_x.units)
output_handler.add_parameter('element',  'y_displacement',   'Y-position',              units.Length,           diff_store.offset_y.units)
output_handler.add_parameter('model',    'iso_diameter',     'Iso diameter',            units.Length,           diff_store.diameter.units)

# Useful for debugging, but maybe shouldn't be included as arguably unnecessary and also some of these could be
# confusing. I.e., what they means as a var in a control-volume element of the plume model and what it means in the
# real-life plume it's simulating is subject to misinterpretation. E.g. the mass pollutant is just for one singular
# slice of the plume representing a single control-volume in that position, not the mass pollutant of the entire plume.
#--------------------------- regime ---- parameter name ---- pretty/formatted label --- unit class ------------ unit type -----------
output_handler.add_parameter('model',    'density_diff',     'Delta density',           units.Density,          units.Density.KILOGRAMS_PER_CUBIC_METER)
output_handler.add_parameter('model',    'um3isoplet',       'Isopleth',                units.Unitless,         units.Unitless.UNITLESS)
output_handler.add_parameter('element',  'mass',             'Mass',                    units.Mass,             units.Mass.KILOGRAMS)
output_handler.add_parameter('model',    'mass_pollutnant',  'Mass pollutant',          units.Mass,             units.Mass.KILOGRAMS)
output_handler.add_parameter('element',  'd_mass',           'Step mass entrained',     units.Mass,             units.Mass.KILOGRAMS)

# TODO: unrecognized vars available in Visual Plumes currently
#     - Bottom -> bottom depth? this is constant per case so why would it be needed to output in the table?
#     - ExtraVar -> ?
