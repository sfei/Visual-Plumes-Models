import math
import numpy as np
from . import units
from .GraphOutput import GraphOutput
from .Output import OutputUM3
from .globals import missing, COS_20_DEGREES, maxsteps, UserInputError
from .vectors import DEGREES_TO_RADIAN, RADIANS_TO_DEGREE, \
    ZERO_VECTOR, I_UNIT_VECTOR, K_UNIT_VECTOR, GRAVITY_VECTOR, MAGNITUDE_GRAVITY, \
    make_vector, change_vector, magnitude, rescale_vector, unit_vector, \
    angle, vector_cosine, project_vector
from .params.ModelParameters import ModelParameters, Model, MaxVerticalReversals, SimilarityProfile
from .params.DiffuserStore import DiffuserStore
from .params.DiffuserParameters import DiffuserParameters
from .ambient.Ambient import Ambient
from .ambient.AmbientHandler import AmbientHandler
from .ambient.AmbientStore import AmbientStore
from .ambient.calculations import seawater_density
from .Element import Element


class UMUnit:

    def __init__(self, model_parameters, diffuser_parameters, diffuser_store, ambient_stack, ambient_store,
                 ambient_ts_stacks=None, output_handler=None, graph_handler=None):
        assert isinstance(model_parameters, ModelParameters)
        assert isinstance(diffuser_parameters, DiffuserParameters)
        assert isinstance(diffuser_store, DiffuserStore)
        assert isinstance(ambient_store, AmbientStore)
        assert isinstance(ambient_stack, (list, tuple))
        for layer in ambient_stack:
            assert isinstance(layer, Ambient)
        if ambient_ts_stacks:
            for ts_stack in ambient_ts_stacks.values():
                for layer in ts_stack:
                    assert isinstance(layer, Ambient)
        # model parameters
        self.model_params       = model_parameters
        # diffuser stuff
        self.diff_params        = None                 # replaces `c`    -> param holding object (diffuser)
        self.og_diff_params     = diffuser_parameters
        self.diffuser_store     = diffuser_store
        # ambient conditions stuff
        self.ambient_handler    = None
        self.ambient_store      = ambient_store
        self.ambient_stack      = ambient_stack
        self.ambient_ts_stacks  = ambient_ts_stacks if ambient_ts_stacks else {}
        # variable/condition tracking objects
        self.ambient            = Ambient()            # replaces ambient in `c`
        self.orig_ambient       = Ambient()            # replaces ambient in `r`
        self.element            = Element()            # replaces `n`    -> control volume at current timestep
        self.last_element       = Element()            # replaces `o`    -> control volume at last timestep
        self.orig_element       = Element()            # replaces `orig` -> control volume conditions at origin
        # general variables
        self.dt                 = 0.003  # timestep
        self.v_source           = K_UNIT_VECTOR  # formerly `Source`
        self.v_ambient          = ZERO_VECTOR    # replaces `Uam`  -> ambient current vector
        self.v_align            = ZERO_VECTOR
        self.density_diff       = 0.0  # (formerly drho) difference in density from ambient to control volume
        self.last_density_diff  = 0.0  # (formerly drho0) last density difference
        self.mass_pollutant     = 0.0  # (formerly polm) mass of pollutant
        self.ini_mass_pollutant = 0.0  # (formerly polm0)
        self.net_dilution       = 0.0  # (formerly netdil)
        # stuff related to isopleth something
        self.um3isoplet         = 0.0  # explicitly store to reduce redundant calls to calculating function
        self.iso_diameter       = 0
        # some persistent variables related to entrainment calculation
        self.cth                = 0.0
        self.entrainment_spco2  = 0.0
        self.entrainment_limspc = 0.0
        self.entrainment_a      = [1.0] * 6  # only uses indices 2-5, but for ease of transcription leaving unused 0-1
        self.entrainment_phi    = 0
        self.entrainment_veave  = 0
        self.entrainment_Nv     = None
        self.last_diameters     = [0, 0]  # (formerly w1, w2) indexed like n-1, n-2
        self.last_dz            = [0, 0]  # (formerly z1, z2) indexed like n-1, n-2
        # related to automatic adjusting of timestep (TODO: do these have to persist? can just be local vars in loop)
        self.bint               = 0.02    # something to do with comparing ratio change in mass
        self.thint              = 0.005   # something to do with comparing ratio change in velocity direction
        self.gamma              = 1.0     # greater of ratio of mass change or velocity direction change (ratio-ed to bint and thint)
        # unidentified general variables
        self.px                 = 0.0
        self.mx                 = 0.0
        self.N_                 = None
        self.Nopp_              = None
        self.cdilprev           = 1.0     # (formerly cdilprev) previous timestep value of self.params.vcld (???)
        self.dist_separation    = 0.0     # (formerly tesmer) appears to be related to spacing before plumes merge
        # self.ini_port_elevation = 0.0     # (formerly horzo) initial port elevation
        self.org_dir            = make_vector()  # (formerly OrgDir)
        # things related to output conditions
        self.last_diameter      = 0
        self.last_graf_step     = 0
        # things related to some post memo messages in middleware
        self.dil100mx           = 0.0
        self.dil100mp           = 0.0
        self.dis100mx           = 0.0
        self.dis100mp           = 0.0
        # status and stop condition trackers and limits
        self.travel             = 0.0     # total traveled distance
        self.step               = 0       # count number of timesteps run
        self.reversals          = 0       # number of times the plume as reversed vertical direction
        self.traps              = 0       # number of times the plume reached a trap level
        # current status tracker dict
        self.statuses_changed   = []
        self.status = {
            'done':               False,
            'surfacing':          False,
            'has_surfaced':       False,
            'has_merged':         False,
            'merging':            False,
            'reversal':           False,
            'stream_limit':       False,
            'not_stream_limited': True,
            'bottom':             False,
            'bottomed':           False,
            'was_stopped':        False,
            'is_lap':             False,
            'is_lap_old':         False,
            'nogo':               False,
            'reached_cig':        False,
            'atvmxz':             False,
            'close_isoplet':      False,
            #
            'hitzid':             False,
            'ambconctoobig':      False
        }
        if output_handler:
            assert isinstance(output_handler, OutputUM3)
            self.outputit = output_handler
        else:
            # default set of outputs
            xy_units = (
                units.Length.FEET
                if self.ambient_store.current_speed.units in (
                    units.Speed.FEET_PER_SECOND,
                    units.Speed.KNOTS,
                    units.Speed.MILES_PER_HOUR,
                    units.Speed.FATHOMS_PER_SECOND
                )
                else units.Length.METERS
            )
            self.outputit = OutputUM3()
            self.outputit.add_parameter('element', 'depth',          'Depth',     units.Length,        self.diffuser_store.depth.units)
            self.outputit.add_parameter('element', 'diameter',       'Width',     units.Length,        self.diffuser_store.diameter.units)
            self.outputit.add_parameter('ambient', 'current_speed',  'Amb-vel',   units.Speed,         self.ambient_store.current_speed.units)
            self.outputit.add_parameter('element', 'speed',          'Velocity',  units.Speed,         self.ambient_store.current_speed.units)
            self.outputit.add_parameter('element', 'x_displacement', 'X-pos',     units.Length,        xy_units)
            self.outputit.add_parameter('element', 'y_displacement', 'Y-pos',     units.Length,        xy_units)
            self.outputit.add_parameter('element', 'density',        'Density',   units.Density,       units.Density.KILOGRAMS_PER_CUBIC_METER)
            self.outputit.add_parameter('element', 'concentration',  'Pollutant', units.Concentration, self.diffuser_store.concentration.units)
            self.outputit.add_parameter('element', 'dilution',       'Dilution',  units.Unitless)
        if graph_handler:
            assert isinstance(graph_handler, GraphOutput)
            self.graphit = graph_handler
        else:
            self.graphit = GraphOutput(self.model_params, self.ambient_store, self.diffuser_store)
        self.model_params.max_reversals_limit = 5
        self.model_params.max_traps_limit = 5
        match self.model_params.max_reversals:
            case MaxVerticalReversals.INITIAL_TRAP_LEVEL:
                self.model_params.max_traps_limit = 1
            case MaxVerticalReversals.SECOND_TRAP_LEVEL:
                self.model_params.max_traps_limit = 2
            case MaxVerticalReversals.MAX_RISE_OR_FALL:
                self.model_params.max_reversals_limit = 1
            case MaxVerticalReversals.SECOND_MAX_RISE_OR_FALL:
                self.model_params.max_reversals_limit = 2
            case _:
                raise UserInputError("Error reading max reversals argument")

    def reset_output_booleans(self):
        for key in self.status.keys():
            self.status[key] = key == 'not_stream_limited'

    def initialize(self):
        # standarize diffuser params units (copy from original as values will get modified)
        self.diff_params = self.og_diff_params.copy()
        # depth/temperature must be converted first as some conversions depend on
        self.diff_params.depth = units.convert(
            self.diff_params.depth,
            units.Length,
            self.diffuser_store.depth.units,
            units.Length.METERS
        )
        self.diff_params.temperature = units.convert(
            self.diff_params.temperature,
            units.Temperature,
            self.diffuser_store.temperature.units,
            units.Temperature.CELSIUS
        )
        # convert rest of vars
        for vkey in self.diff_params._vars_:
            if vkey in ('depth', 'temperature'):
                continue
            unit_type = units.from_var_name(vkey)
            from_units = self.diffuser_store.get(vkey).units
            if vkey == 'isopleth':
                # Calculate isopleth parameter
                if self.model_params.model == Model.PDS:
                    self.diff_params.isopleth = 1e9
                    continue
                # isopleth is matched to units from another field, depending on the type
                match from_units:
                    case units.Isopleth.CONCENTRATION:
                        unit_type = units.Concentration
                        from_units = self.diffuser_store.concentration.units
                    case units.Isopleth.SALINITY:
                        unit_type = units.Salinity
                        from_units = self.diffuser_store.salinity.units
                    case units.Isopleth.TEMPERATURE:
                        unit_type = units.Temperature
                        from_units = self.diffuser_store.temperature.units
                    case units.Isopleth.SPEED:
                        unit_type = units.Speed
                        from_units = self.ambient_store.current_speed.units
                    case _:
                        raise UserInputError("Unknown/unsupported isopleth type")
            self.diff_params.set(
                vkey,
                units.convert(
                    value=self.diff_params.get(vkey),
                    units_or_var_name=unit_type,
                    from_units=from_units,
                    model_params=self.model_params,
                    celsius=self.diff_params.temperature,
                    # PSU parameter not necessary, only used in decay rate conversions
                    depth=self.diff_params.depth
                )
            )

        # diffuser input validations (not complete but the most important ones)
        if self.diff_params.diameter <= 0:
            self.outputit.memo("To complete simulation, diffuser diameter must be greater than zero")
            raise UserInputError("Invalid diffuser diameter. Check inputs.")
        if self.diff_params.num_ports <= 0:
            self.outputit.memo("To complete simulation, number of ports must be greater than zero")
            raise UserInputError("Invalid number of ports. Check inputs.")
        if self.diff_params.num_ports > 1 and self.diff_params.port_spacing <= 0:
            self.outputit.memo("To complete simulation, port spacing must be greater than zero")
            raise UserInputError("Invalid port spacing. Check inputs.")
        if self.diff_params.effluent_flow <= 0:
            self.outputit.memo("To complete simulation, diffuser effluent flow must be greater than zero")
            raise UserInputError("Invalid diffuser effluent flow. Check inputs.")

        # from Procedure Initialize in vectors.py
        # some stuff outside of loop in merge() pulled into here as well

        # unit vector perpendicular to the shore
        #uv_perp_shore = unit_vector(np.cross(K_UNIT_VECTOR, self.model_params.v_shore))

        # graph vector here moved into top of merge() code

        self.ambient_handler = AmbientHandler(self.model_params, self.ambient_store)

        # set bottom depth (this was taken from initializeUM in main.pas)
        # before z is converted to depth in ambient_handler.fill()
        self.model_params.bottom_depth = self.diff_params.depth
        if self.ambient_store.z.z_is_depth:
            max_z = max((amb.z for amb in self.ambient_stack))
            max_z = self.ambient_handler.get_z(max_z)
            if max_z > self.model_params.bottom_depth:
                self.model_params.bottom_depth = max_z
        for vkey in self.ambient_store._input_vars_:
            ambient_store_at = self.ambient_store.get(vkey)
            if not ambient_store_at.from_time_series or not ambient_store_at.z_is_depth:
                continue
            max_z = max((amb.z for amb in self.ambient_ts_stacks[vkey]))
            max_z = self.ambient_handler.get_z(max_z)
            if max_z > self.model_params.bottom_depth:
                self.model_params.bottom_depth = max_z

        self.ambient_handler.fill(
            model_params=self.model_params,
            ambient_stack=self.ambient_stack,
            ambient_ts_stacks=self.ambient_ts_stacks,
            diff_params=self.diff_params,
            orig_ambient=self.orig_ambient
        )
        self.ambient_stack = self.ambient_handler.ambient_stack
        # TODO: if casecount != 1: statusbar.panels[1].text:= ' Run '+inttostr(casecount)+' of '+inttostr(nruns)+' or more runs'

        # get ambient conditions (set z-level, then call ambient_level() to reinterpolate into params.ambient)
        self.ambient.depth = self.diff_params.depth
        self.v_ambient = self.ambient_handler.ambient_level(
            umunit=self,
            ambient_cond=self.ambient,
            depth=self.ambient.depth,
            bottom_depth=self.model_params.bottom_depth
        )
        self.orig_ambient.density = self.ambient.density
        self.orig_ambient = self.ambient.copy()

        self.diff_params.num_ports = round(self.diff_params.num_ports)
        if self.diff_params.num_ports == 1 and self.model_params.model != Model.DOS3:
            self.diff_params.port_spacing = missing

        self.element.diameter         = self.diff_params.diameter
        self.element.depth            = self.diff_params.depth
        self.element.vertical_angle   = self.diff_params.vertical_angle
        self.element.horizontal_angle = self.diff_params.horizontal_angle
        self.element.salinity         = self.diff_params.salinity
        self.element.temperature      = self.diff_params.temperature
        self.element.concentration    = self.diff_params.concentration

        # calculate initial speed and diameter
        radius = self.element.diameter*0.5
        if self.model_params.model == Model.PDS:
            self.element.speed = self.diff_params.effluent_flow / self.element.diameter
            # TODO: Divide by zero safeguard, but does this make sense? Speed explodes exponentially as angle gets
            # vertical. Putting threshold at ~5 degrees (~11x multiplier)
            self.element.speed /= (self.element.vertical_angle if abs(self.element.vertical_angle) > 0.09 else 0.09)
        else:
            self.element.diameter *= math.sqrt(self.model_params.contraction_coeff)
            radius = self.element.diameter*0.5
            self.element.speed = self.diff_params.effluent_flow / self.diff_params.num_ports / math.pi / radius**2

        self.last_element = self.element.copy()

        # create velocity vector (for last!)
        vhor_radians = DEGREES_TO_RADIAN*self.diff_params.horizontal_angle
        cos_vhor   = math.cos(vhor_radians)
        sin_vhor   = math.sin(vhor_radians)
        vang_radians = DEGREES_TO_RADIAN*self.diff_params.vertical_angle
        cos_vang   = math.cos(vang_radians)
        sin_vang   = math.sin(vang_radians)
        self.last_element.v_velocity = make_vector(
            x=self.element.speed * cos_vhor * cos_vang,
            y=self.element.speed * sin_vhor * cos_vang,
            z=self.element.speed * sin_vang
        )

        # seawater density
        self.element.density = seawater_density(
            at_equilibrium=self.model_params.at_equilibrium,
            salinity=self.diff_params.salinity,
            temperature=self.diff_params.temperature,
            depth=self.element.depth,
            in_sigma=False
        )
        self.orig_element.density = self.element.density
        self.density_diff = self.ambient.density - self.element.density

        # initial height (10% of diameter)
        self.element.height = self.diff_params.diameter * 0.1

        # initial mass
        self.element.mass = self.element.density*math.pi*(radius**2)*self.element.height
        #self.element.d_mass = self.bint*self.element.mass
        # newmass = self.element.d_mass  # never used

        # source velocity vector in horizontal plane only
        self.v_source = make_vector(
            x=self.diff_params.offset_x,
            y=self.diff_params.offset_y,
            z=0
        )
        # total displacement vector from x/y of source but on water surface
        self.element.v_surface_tdsp = make_vector(
            x=self.diff_params.offset_x,
            y=self.diff_params.offset_y,
            z=-self.ambient.depth
        )
        # displacement vector during timestep only
        self.element.v_displace = make_vector()  #self.dt*self.last_element.v_velocity

        # velocity vector (at first timestep) as mass-averaged momentum of current plus entrained ambient
        self.element.v_velocity = (  # (oc) estimate bend
            (self.element.mass*self.last_element.v_velocity + self.element.d_mass*self.v_ambient) /
            (self.element.mass + self.element.d_mass)
        )

        self.N_, self.Nopp_ = self.element.body_calc(self.last_element, self.dt, in_init=True)
        self.element.Cc = rescale_vector(0.0005*self.element.diameter, self.element.Rc)
        # fill this 37 length list of radii, which is all a cross product of the last against the next velocity vector
        # looks like it just rotates the radius vector around..
        # TODO is any element besides the 0th ever even used?
        b_scale = radius*math.tan(math.pi/18.0)
        for i in range(1, 37):
            last_B = self.element.v_radii[i-1]
            self.element.v_radii[i] = rescale_vector(
                radius,
                last_B + rescale_vector(
                    b_scale,
                    np.cross(self.element.v_velocity, last_B)
                )
            )

        # pollutant mass initialization
        self.mass_pollutant     = self.element.concentration * self.element.mass
        self.ini_mass_pollutant = self.mass_pollutant
        # self.element.v4o3       = self.element.concentration

        # edge to edge port distance is initial plume spacing
        self.dist_separation = self.diff_params.port_spacing - self.element.diameter

        if self.ambient.current_speed <= 0:
            self.ambient.current_speed = 1e-5

        # reset initial velocity
        self.element.v_velocity = self.last_element.v_velocity
        self.orig_element = self.last_element = self.element.copy()
        self.convert_vectors()
        # self.ini_port_elevation = self.diff_params.port_elevation

        self.define_net_dilution()
        self.calc_isoplet()
        self.update_status()

    def define_net_dilution(self):
        # this comes from main.pas
        if self.outputit.has_parameter('element', 'vcld'):  # (oc) bug 2013 enhancement  3.79
            self.net_dilution = self.element.dilution / (2.22 if self.element.diameter > 2 * self.diff_params.port_spacing else 3.79)
        elif not self.model_params.report_effective_dillution or self.ambient.bg_conc == 0:
            self.net_dilution = self.element.dilution
        elif self.element.concentration > 0:
            self.net_dilution = self.ini_mass_pollutant / self.element.concentration
        else:
            self.net_dilution = 0  # shouldn't happen by divide by zero safeguard

    def convert_vectors(self):
        # update coordinates
        self.element.x_displacement = self.element.v_surface_tdsp[0]
        self.element.y_displacement = self.element.v_surface_tdsp[1]
        self.element.total_surf_dsp = magnitude(self.element.v_surface_tdsp)
        # are we randomly using ambient-z to hold something different?
        # TODO: isn't the z-coordinate always 0 for surface displacement, so we're just resetting it to surface?
        self.ambient.depth = self.element.depth = -self.element.v_surface_tdsp[2]
        # update speed
        self.element.speed = magnitude(self.element.v_velocity)
        # update angles
        v_horz = make_vector(x=self.element.v_velocity[0], y=self.element.v_velocity[1])
        self.element.vertical_angle = RADIANS_TO_DEGREE * angle(v_horz, self.element.v_velocity)
        if self.element.v_velocity[2] < 0:
            self.element.vertical_angle *= -1
        self.element.horizontal_angle = RADIANS_TO_DEGREE * angle(v_horz, make_vector(x=1))
        if self.element.v_velocity[2] < 0 and abs(self.element.v_velocity[2]) > 1e-10:
            self.element.horizontal_angle *= -1  # (oc) debug 1/2/00

    def kpro(self):
        match self.model_params.similarity_profile:
            case SimilarityProfile.DEFAULT:
                if self.model_params.model == Model.DKH:
                    poly = 1.2
                    sngl = 1.79
                else:
                    poly = 1.5
                    sngl = 2.0
            case SimilarityProfile.POWER_3_2:
                poly = 2.222
                sngl = 3.89
            case SimilarityProfile.GAUSSIAN:
                poly = 2.1466
                sngl = 3.67
            case _:
                raise UserInputError("Invalid/unrecognized similarity profile")
        if self.status['has_merged'] and self.diff_params.num_ports > 1:
            # TODO: what makes sense from element, what makes sense from diffuser?
            if self.element.diameter > 2.0*self.diff_params.port_spacing:
                arg = poly
            else:
                arg = sngl - (sngl - poly)*(self.element.diameter - self.diff_params.port_spacing)/self.diff_params.port_spacing
        else:
            arg = sngl
        self.element.cl_dilution = self.element.dilution / arg
        if self.step == 1:
            self.cdilprev = 1.0
        if self.element.cl_dilution < self.cdilprev:
            self.element.cl_dilution = self.cdilprev
        else:
            self.cdilprev = self.element.cl_dilution
        return arg

    def calc_isoplet(self):
        kappa  = 1.888
        isoval = 1.0
        kpro   = self.kpro()
        if self.model_params.model != Model.PDS:
            match self.diffuser_store.isopleth.units:
                case units.Isopleth.CONCENTRATION:
                    camb = self.ambient.bg_conc
                    sorc = self.diff_params.concentration
                case units.Isopleth.SALINITY:
                    camb = self.ambient.salinity
                    sorc = self.diff_params.salinity
                case units.Isopleth.TEMPERATURE:
                    camb = self.ambient.temperature
                    sorc = self.diff_params.temperature
                case units.Isopleth.SPEED:
                    camb = self.ambient.current_speed
                    sorc = self.diff_params.effluent_flow / (math.pi*0.25*self.orig_element.diameter**2)
                case _:
                    raise UserInputError("Unknown/unsupported isopleth type")
            if self.status['stream_limit'] and self.ambient.current_speed != 0:
                # this makes no sense, can simplify gf to 2/3
                # gf = 1-c[vdia]/2*c[vdia]/2*2*2/3/c[vdia]/c[vdia]
                #      1 - c/2*c/2*2*2/3/c/c
                #      1 - (2*2*c*c) / (2*2*3*c*c)
                #      1 - 1/3
                # gf = 2/3
                cCL = camb + 1.5*(self.px/self.mx - camb)
                arg = 1.0 - (self.diff_params.isopleth - camb) / (cCL - camb)
                isoval = 0 if arg < 0 else math.sqrt(arg)
                max_isoval = self.model_params.tpb_channel_width / self.element.diameter
                if isoval > max_isoval:
                    isoval = max_isoval
            else:
                this_thing = self.element.dilution / kpro * (self.diff_params.isopleth - camb) / (sorc - camb)
                if this_thing > 1.0:
                    isoval = 0.0
                else:
                    match self.model_params.similarity_profile:
                        case SimilarityProfile.DEFAULT:
                            isoval = math.sqrt(1.0 - this_thing)
                        case SimilarityProfile.POWER_3_2:
                            isoval = (1.0 - math.sqrt(this_thing))**(2.0/3.0)
                        case SimilarityProfile.GAUSSIAN:
                            isoval = 1.0/kappa*math.sqrt(-math.log(this_thing))
                        case _:
                            raise UserInputError("Invalid/unrecognized similarity profile")
                    if self.element.dilution < kpro:
                        a = self.orig_element.diameter*(kpro - self.element.dilution)/kpro
                        isoval = (isoval*(self.element.diameter - a) + a)/self.element.diameter
                    if isoval > 1.0 and not self.status['ambconctoobig']:
                        self.status['ambconctoobig'] = True
                        self.outputit.memo('Ambient species greater than plume isopleth value, physical boundary graphed')
                if isoval > 1.0:
                    isoval = 1.0  # (o.c.) limits max diameter
                # this was weird in o.c. -- was converted to concentration (w/o checking iso units), this is what I
                # think is actually intended
                if (
                    self.diffuser_store.isopleth.units == units.Isopleth.CONCENTRATION
                    and self.diff_params.isopleth > sorc
                ):
                    isoval = 0.0
                    self.status['done'] = True
                    self.outputit.memo('Isopleth value > base value, simulation stopped.')
        self.um3isoplet = isoval
        self.iso_diameter = self.um3isoplet*self.element.diameter

    def calc_entrainment(self, radius, d_radius, displacement):
        # calc new vectors
        v_horiz_velocity = change_vector(self.element.v_velocity, z=0)
        if angle(GRAVITY_VECTOR, self.element.v_velocity) != 0:
            self.entrainment_Nv = unit_vector(np.cross(self.element.v_velocity, GRAVITY_VECTOR))
        Uamv = project_vector(self.v_ambient, self.entrainment_Nv)
        Uamn = self.v_ambient - Uamv
        mag_velocity = magnitude(self.element.v_velocity)
        if mag_velocity == 0:
            mag_velocity = 1e-12  # avoid divide by zero error
        sth = self.element.v_velocity[2] / mag_velocity
        cth = magnitude(v_horiz_velocity) / mag_velocity

        veavgin = 0
        if self.model_params.allow_induced_current and self.diff_params.num_ports > 1:
            veavgin = self.entrainment_veave
            # (oc) induced-current effect
            if self.element.diameter > self.diff_params.port_spacing:
                veavgin *= self.element.diameter / self.diff_params.port_spacing

        # this really only needs to be calculated once
        if self.diff_params.num_ports <= 1:
            self.entrainment_spco2 = missing
        else:
            cos_from_horiz = vector_cosine(self.v_align, v_horiz_velocity)
            if abs(cos_from_horiz) < COS_20_DEGREES:
                angle_from_horiz = angle(self.v_align, v_horiz_velocity)
                self.entrainment_spco2 = abs(0.5*self.diff_params.port_spacing*math.sin(angle_from_horiz))
            else:
                self.entrainment_spco2 = 0.5*self.diff_params.port_spacing*math.sin(math.pi/9.0)
            half_limspc = self.entrainment_limspc * 0.5
            if self.entrainment_spco2 > half_limspc and self.entrainment_limspc != 0:
                self.entrainment_spco2 = half_limspc
        # espc = self.entrainment_spco2

        # adjustment of some other variables (only hits in multiple port scenarios)
        if radius >= self.entrainment_spco2:
            self.entrainment_phi = math.atan(math.sqrt(abs(  # (oc) arctan bug
                (radius**2 - self.entrainment_spco2**2) / self.entrainment_spco2
            )))
            self.entrainment_a[2] = 1.0 - 2.0*self.entrainment_phi/math.pi
            self.entrainment_a[3] = self.entrainment_spco2 / radius
            self.entrainment_a[4] = self.entrainment_a[2] + math.sin(2.0*self.entrainment_phi)/math.pi
            self.entrainment_a[5] = 1.0 / self.diff_params.num_ports

        # calc some other random parameters
        cvaden_dt_rb_h = self.ambient.density * self.dt * radius * self.element.height
        cvaden_dt_rb_h_2 = 2.0 * cvaden_dt_rb_h
        cyl = cvaden_dt_rb_h_2 * abs(sth) * self.entrainment_a[3] * (magnitude(Uamv) + veavgin)
        cyln = cvaden_dt_rb_h_2 * self.entrainment_a[5] * magnitude(Uamn)
        if displacement != 0:
            if (np.dot(Uamv, self.element.v_velocity) > 0) == (d_radius > 0):
                rgro = abs(math.pi / displacement * cth * d_radius * self.entrainment_a[2] * magnitude(Uamv) * cvaden_dt_rb_h)
            else:
                rgro = 0
            curv = (
                math.pi*0.5
                * (radius**2)
                * angle(self.last_element.v_velocity, self.element.v_velocity) / displacement
                * self.element.height * self.ambient.density * self.entrainment_a[4]
                * (np.dot(self.element.v_velocity, unit_vector(self.element.Rc)) + veavgin)
            )
            if np.dot(Uamv, self.element.Rc) < 0:
                curv *= -1
        else:
            # if no movement, values move towards infinity. instead assume static (no growth, no curve)
            rgro = 0
            curv = 0

        self.ambient.depth = -self.element.v_surface_tdsp[2]

        # adjust upstream diameter based on ratio of total vertical distance traveled from last, last timestep
        upstream_diameter = 0
        if (
            (self.last_dz[-2] > self.element.v_surface_tdsp[2] > self.last_dz[-1])
            or (self.last_dz[-2] < self.element.v_surface_tdsp[2] < self.last_dz[-1])
        ):
            upstream_diameter = self.last_diameters[-2] + (
                (self.last_diameters[-1] - self.last_diameters[-2])
                * (self.element.v_surface_tdsp[2] - self.last_dz[-2])
                / (self.last_dz[-1] - self.last_dz[-2])
            )
            if upstream_diameter < 0:
                upstream_diameter = 0
            elif upstream_diameter > self.element.diameter:
                upstream_diameter = self.element.diameter

        # curv and cyl adjusted by inverse ratio of diameter difference on upstream-side
        one_minus_ratio_diameter = 1.0 - upstream_diameter / self.element.diameter
        curv *= one_minus_ratio_diameter
        cyl *= one_minus_ratio_diameter

        # first entrainment factor is some addition of all this
        eins = rgro + cyln
        if cyl + curv > 0:
            eins += cyl + curv

        # second entrainment factor is...
        uwper = abs(magnitude(Uamv) * sth)
        ven = self.model_params.aspiration_coeff * magnitude(
            self.element.v_velocity - rescale_vector(
                np.dot(Uamv, unit_vector(self.element.v_velocity)),
                self.element.v_velocity
            )
        )
        angl = math.atan((uwper/ven - 1.0)**0.5) if uwper > ven else 0
        self.entrainment_veave = ven - angl/math.pi*ven - uwper/math.pi*(1.0 - math.sin(angl))
        zwei = cvaden_dt_rb_h_2 * math.pi * abs(self.entrainment_veave) * self.entrainment_a[2]

        # entrainment increases total mass by this much
        return (eins + zwei), cth

    def merge(self):
        # initial conditions and reset
        self.initialize()

        # graph source to shore vector
        # for now, since only one diffuser which cannot change source x/y, set so only drawn once in first case
        if self.model_params.use_shore_vector and self.model_params.casecount == 1:
            v_shore = make_vector(
                x=(self.model_params.dist_to_shore*math.cos(DEGREES_TO_RADIAN*self.model_params.dir_to_shore)),
                y=(self.model_params.dist_to_shore*math.sin(DEGREES_TO_RADIAN*self.model_params.dir_to_shore))
            )
            self.graphit.graph_vector(self.v_source, v_shore)

        # TODO: this seems related to output graphs, can be pulled out of here I think
        # self.set_mark_size(casecount)

        if self.element.concentration <= 0:
            self.outputit.memo("To complete simulation, effluent concentration must be greater than zero")
            return
        if self.ambient.kt > 0.01:
            self.outputit.memo("Decay rate very fast, may overflow error")
        if self.status['nogo']:
            return
        z_vel = self.element.v_velocity[2]
        if ((self.element.density < 0 < z_vel) or (self.element.density > 0 > z_vel)) \
                and self.model_params.max_traps_limit < 2:
            self.outputit.memo("Potential for more diluition")

        # setup random vars
        denom1 = travel0 = 0

        orig_density = self.element.density
        last_radius = radius = magnitude(self.element.v_radii[0])  # formerly `rsave` and `rb`
        self.last_dz[-2] = self.last_dz[-1] = self.element.v_surface_tdsp[2]    # formerly `z1` and `z2`
        self.last_diameters[-2] = self.last_diameters[-1] = self.element.diameter  # formerly `w1` and `w2`

        # something related to entrainment
        vert_angle = angle(GRAVITY_VECTOR, self.element.v_velocity)
        if abs(vert_angle) <= 1e-10 or abs(math.pi - vert_angle) < 1e-10:
            if magnitude(self.v_ambient) == 0:
                self.entrainment_Nv = make_vector(y=1)
            else:
                self.entrainment_Nv = unit_vector(np.cross(GRAVITY_VECTOR, self.v_ambient))
        else:
            self.entrainment_Nv = unit_vector(np.cross(self.element.v_velocity, GRAVITY_VECTOR))
        self.v_align = self.entrainment_Nv

        Ao_ = self.diff_params.acute_mixing_zone*unit_vector(self.element.v_velocity)

        message = self.graphit.set_org_dir(self.element, self.v_ambient)
        if message:
            self.outputit.memo(message)

        # TODO: isn't this just 0?
        lapstep = self.step

        # self.graphit.add_nulls(self.step, self.model_params.casecount)

        ############################################################
        # ENTER LOOP
        ############################################################
        while not self.status['done'] and not self.status['atvmxz']:
            # add to timestep count
            self.step += 1

            # reset vars
            self.dil100mx = self.dil100mp = self.dis100mx = self.dis100mp = 0

            # entrainment calculation happens inside here (few things only used by it can be calc'd right before)
            mass_entrained, self.cth = self.calc_entrainment(
                radius=radius,
                d_radius=(radius - last_radius),
                displacement=(magnitude(self.element.v_displace) if self.step > 1 else self.element.speed*self.dt)
            )
            self.element.d_mass = mass_entrained

            # automatic timestep adjustment
            # (oc) unrealistic large dt correction if used at first step: Kenwyn Nov 2012
            if self.step > 1:
                self.gamma = max(
                    mass_entrained/(self.bint*self.element.mass),
                    angle(self.element.v_velocity, self.last_element.v_velocity) / self.thint
                )
                self.dt /= self.gamma
                self.element.d_mass /= self.gamma

            # not sure exactly what's going on, but end goal is to set `vab`
            # this was cfentrain function, but since it's only called once
            # d_4_3 = self.element.diameter**(4.0/3.0)
            # d_2   = self.element.diameter**2
            # dmfcf = self.element.mass / math.erf(
            #     # TODO: erf() was called from global, I assume this is error function, equivalent to python math.erf
            #     1.5 / (1.0 + 8.0*self.ambient.ff_diff_coeff*d_4_3*self.dt / d_2)**3
            # )
            # TODO: what is this used for?
            # self.element.vab = dmfcf/(self.element.d_mass + 1e-12)  # (oc) bug 2013, dmcf can be zero

            # something about checking for lap?
            mag_nRc = magnitude(self.element.Rc)
            if mag_nRc < self.element.diameter:
                self.element.d_mass *= mag_nRc / self.element.diameter
                self.status['is_lap'] = True
                lapstep = self.step
            elif self.step > lapstep+1:
                self.status['is_lap'] = False

            # copy element state to last, set 'last' variables
            self.last_element = self.element.copy()
            self.last_density_diff = self.density_diff
            last_radius = radius
            last_dist_separation = self.dist_separation

            # adjust mass and mass-averaged properties
            new_mass = self.element.mass + self.element.d_mass
            self.element.salinity = (
               (self.element.mass * self.element.salinity + self.element.d_mass * self.ambient.salinity)
               / new_mass  # (oc) not Sonntag enthalpy approach
            )
            self.element.temperature = (
                (self.element.mass * self.element.temperature + self.element.d_mass * self.ambient.temperature)
                / new_mass
            )

            # adjust density
            self.element.density = seawater_density(
                at_equilibrium=self.model_params.at_equilibrium,
                salinity=self.element.salinity,
                temperature=self.element.temperature,
                depth=self.element.depth,
                in_sigma=False
            )
            self.density_diff = self.ambient.density - self.element.density
            # TODO: used anywhere?
            #equilibrium = self.ini_density_diff * self.density_diff < 0

            # move time elapsed
            self.element.total_time += self.dt

            # new velocity as mass-averaged momentum of previous plus entrained
            self.element.v_velocity = (
                (self.element.mass*self.element.v_velocity + self.element.d_mass*self.v_ambient)
                / (self.element.mass + self.element.d_mass)
                # plus buoyancy factor
                - self.density_diff/self.element.density*self.dt*GRAVITY_VECTOR
            )

            # adjust movement and travel
            self.element.v_displace = self.element.v_velocity*self.dt
            travel0 = self.travel
            self.travel += magnitude(self.element.v_displace)
            self.element.v_surface_tdsp += self.element.v_displace

            # update An_ vector.. whatever that is
            An_ = unit_vector(project_vector(self.element.v_surface_tdsp, K_UNIT_VECTOR))
            cos_Ao_An = abs(vector_cosine(Ao_, An_))
            const_a = 0.5*self.diff_params.num_ports*self.diff_params.port_spacing + magnitude(Ao_)
            if cos_Ao_An < 1e-6:
                An_ = const_a*unit_vector(An_)
            else:
                An_ = magnitude(Ao_)/cos_Ao_An*An_
                if magnitude(An_) > const_a:
                    An_ = const_a*unit_vector(An_)
            if self.diff_params.num_ports < 2:
                An_ = rescale_vector(self.diff_params.acute_mixing_zone, I_UNIT_VECTOR)
            if magnitude(project_vector(self.element.v_surface_tdsp, K_UNIT_VECTOR)) >= magnitude(An_):
                self.status['hitzid'] = True

            # increase mass and update radius (radius will get adjusted again, diameter calc must wait)
            self.element.mass += self.element.d_mass
            radius = math.sqrt(self.element.mass/self.element.density/math.pi/self.element.height)

            # update separation, check for merging
            self.dist_separation = self.entrainment_spco2 - radius
            self.status['merging'] = (
                self.dist_separation*last_dist_separation < 0.0
                or (self.diff_params.num_ports > 1 and 0 < self.entrainment_limspc < self.element.diameter)
            )

            # adjust radius after separation calculation
            if radius > self.entrainment_spco2:
                phi = self.entrainment_phi
                radius *= math.sqrt(math.pi / (math.pi - 2.0*phi + 2.0*math.cos(phi)*math.sin(phi)))

            # do this calc
            self.calc_isoplet()

            # body calculations
            self.N_, self.Nopp_ = self.element.body_calc(self.last_element, self.dt, um3iso=self.um3isoplet)
            # check bottomed
            self.status['bottom'] = (
                self.N_[2] < -self.model_params.bottom_depth
                or self.Nopp_[2] < -self.model_params.bottom_depth
            )

            # update height as function of change in velocity, with flattening/lengthening proportional to decrease/increase
            mag_v_curr = magnitude(self.element.v_velocity)
            mag_v_last = magnitude(self.last_element.v_velocity)
            if mag_v_curr != 0 and mag_v_last != 0:
                # for complete stop or acceleration from stopped, simply ignore change in height (edge cases 0/infinity)
                self.element.height *= mag_v_curr/mag_v_last

            # update dilution as function of proportional change in mass and density
            self.element.dilution = self.element.mass*orig_density / (self.element.density*self.orig_element.mass)

            # adjust pollutant mass and concentration
            self.mass_pollutant += (
                # loss from decay rate of pollution greater than mass from background concentration
                -self.ambient.kt*(self.mass_pollutant - self.element.mass*self.ambient.bg_conc)*self.dt
                # added from background concentration of entrained
                + self.element.d_mass*self.ambient.bg_conc
                # something related to tidal pollution buildup
                + self.element.d_mass*self.model_params.tpb_bincon
            )
            self.element.concentration = (self.mass_pollutant / self.element.mass) * (self.element.density / orig_density)  # (oc) Paul Taylor bug, email 20051003; was polm/m befor
            # self.element.v4o3 = self.element.concentration  # ? <- original comment

            # update diameter now
            self.element.diameter = 2.0*radius
            self.iso_diameter = self.um3isoplet * self.element.diameter

            v_arg = unit_vector(project_vector(self.orig_element.v_velocity, K_UNIT_VECTOR))
            if np.dot(self.element.v_surface_tdsp, v_arg) > self.diff_params.acute_mixing_zone and self.dil100mx < 1.0:
                self.dil100mx = self.element.dilution
                self.dis100mx = self.diff_params.acute_mixing_zone
            v_arg = unit_vector(np.cross(GRAVITY_VECTOR, v_arg))
            if np.dot(self.element.v_surface_tdsp, v_arg) > self.diff_params.acute_mixing_zone and self.dil100mp < 1.0:
                self.dil100mp = self.element.dilution
                self.dis100mp = self.diff_params.acute_mixing_zone

            self.last_dz[-2] = self.last_dz[-1]
            self.last_dz[-1] = self.element.v_surface_tdsp[2]
            self.last_diameters[-2] = self.last_diameters[-1]
            self.last_diameters[-1] = self.element.diameter

            # update ambient conditions (set new z-level, then call ambient_level() to reinterpolate into params.ambient)
            self.v_ambient = self.ambient_handler.ambient_level(
                umunit=self,
                ambient_cond=self.ambient,
                depth=self.ambient.depth,  # TODO: ambient.depth gets reset in convert vectors, which is called after this, probably right but check later
                bottom_depth=self.model_params.bottom_depth
            )
            self.convert_vectors()

            self.define_net_dilution()

            # check surfacing
            if self.N_[2] > 0 or self.Nopp_[2] > 0:
                self.status['surfacing'] = True

            # status update check (with stop, update, write checks)
            last_denom = denom1
            speed = magnitude(self.element.v_velocity)
            denom1 = (
                MAGNITUDE_GRAVITY*(self.ambient.density - self.element.density)/self.ambient.density*self.ambient.depth
                + 0.5*(0.5*self.element.diameter*speed/self.ambient.depth)**2
                - speed**2
            )
            denomproduct = last_denom*denom1
            self.update_status(magnitude_An=magnitude(An_), denomproduct=denomproduct)

            # last states that must be set after stop/write condition check
            self.status['is_lap_old'] = self.status['is_lap']

        ############################################################
        # EXIT LOOP
        ############################################################

        # this was moved out of loop, because it's tied to break conditions anyway
        # memo stuff post loop moved to middleware
        self.graphit.graph(self, self.model_params.casecount, self.net_dilution)  # (oc) bug 2015 25 Mar added to close isopleth

        if self.ambient.current_speed == 0:
            self.ambient.current_speed = 0.001  # (oc) Salas div by zero bug 2018

        # graph39 call includes what was in graph4panel
        self.graphit.graph_end(self, self.model_params.casecount, self.net_dilution)
        if math.floor(self.travel/2.5) != math.floor(travel0/2.5):
            self.graphit.graph(self, self.model_params.casecount, self.net_dilution)

        # tidal pollution buildup updated
        self.element.buildup = self.model_params.tpb_bincon

        # Brooks Far-Field model moving to middleware

        if self.status['bottom'] and self.step < 3:
            self.outputit.memo('Bottom geometry consistent?  Try increasing port elev and/or ambient depth')

    def reset_status_changes(self):
        self.statuses_changed = []

    def check_status_changed(self, name, value=True):
        return self.status[name] == value and name in self.statuses_changed

    def _update_status(self, name, set_to):
        self.status[name] = set_to
        if name not in self.statuses_changed:
            self.statuses_changed.append(name)

    def update_status(self, magnitude_An=0, denomproduct=0):
        # oc from stop condition
        if self.step > 0:
            if self.last_element.v_velocity[2]*self.element.v_velocity[2] < 0:
                if not self.status['reversal'] and 'reversal' not in self.statuses_changed:
                    self.statuses_changed.append('reversal')
                self.status['reversal'] = True
                self.reversals += 1
            elif self.status['reversal']:
                self._update_status('reversal', False)
            density_reversal = self.last_density_diff * self.density_diff < 0
            if density_reversal:
                self.traps += 1

            # reverse direction (reversal check was during vm_stop_condition)
            if self.status['reversal']:
                self.last_dz[-2] = self.last_dz[-1]
                self.last_dz[-1] = self.element.v_surface_tdsp[2]
                self.last_diameters[-2] = self.last_diameters[-1]
                self.last_diameters[-1] = self.element.diameter*0.5
        else:
            density_reversal = 0

        # check stop conditions
        self.vm_stop_condition(density_reversal)

        # oc from write condition
        if (
            not self.status['reached_cig']
            and self.model_params.tpb_mixing_zone_ceil > 0
            and (-self.N_[2] < self.model_params.tpb_mixing_zone_ceil or -self.Nopp_[2] < self.model_params.tpb_mixing_zone_ceil)
        ):
            self._update_status('reached_cig', True)
        if self.status['surfacing'] and not self.status['has_surfaced']:
            self._update_status('has_surfaced', True)
        if self.status['stream_limit'] and self.status['not_stream_limited']:
            self._update_status('not_stream_limited', False)
        # from outputit
        if self.status['bottom'] and not self.status['bottomed']:
            self._update_status('bottomed', True)

        if self.step > 0:
            # from output it
            if not self.status['atvmxz']:
                kmag_surface_tdsp = magnitude(project_vector(self.element.v_surface_tdsp, K_UNIT_VECTOR))
                if kmag_surface_tdsp >= magnitude_An:
                    self._update_status('atvmxz', True)
                    if self.model_params.model != Model.PDS:
                        self.graphit.graph_mixing_zone(self, self.model_params.casecount, self.net_dilution)
            if self.status['bottom'] and not self.status['bottomed']:
                self._update_status('bottomed', True)
            if self.entrainment_limspc == 0 and self.element.diameter > 2*self.diff_params.depth:
                self.entrainment_limspc = self.diff_params.depth
                self.status['merging'] = True
                if not self.status['has_merged']:
                    self._update_status('has_merged', True)
            if (
                not self.status['stream_limit']
                and magnitude(self.v_ambient) > 1e-5
                and self.element.dilution > (self.ambient_handler.Qstream + self.diff_params.effluent_flow)/self.diff_params.effluent_flow
            ):
                # TODO: because Qstream is multiplier of `model_params.tpb_channel_width` which is 0 if tidal pollution
                # buildup is not set, this always triggers at first timestep in such cases (should it be ignored if TPB
                # is unset?)
                self._update_status('stream_limit', True)
                if self.entrainment_limspc == 0:
                    self.entrainment_limspc = self.diff_params.depth
                self.px = self.mass_pollutant
                self.mx = self.element.mass
            if self.status['merging'] and not self.status['has_merged']:
                self._update_status('has_merged', True)
            self.status['merging'] = False

        # check write conditions
        do_outputs, do_graphs = self.vm_write_condition(denomproduct=denomproduct)
        if do_outputs:
            self.output_it(denomproduct)
        if do_graphs:
            self.graphit.graph(self, self.model_params.casecount, self.net_dilution)

    def vm_write_condition(self, denomproduct=0):
        do_outputs = False  # goes to outputit()
        do_graphs  = False  # goes to grafit()

        if (
            self.step == 0
            or self.element.diameter == 0
            or (self.step > 0 and (self.density_diff == 0 or self.last_density_diff == 0))
            or self.status['done']
            or self.step > maxsteps
            or self.element.v_displace[2] * self.last_element.v_displace[2] < 0
            or (self.step >= self.model_params.write_step_freq and self.step % round(self.model_params.write_step_freq) == 0)
            or self.check_status_changed('has_merged')
            or self.last_density_diff*self.density_diff < 0
            or self.status['reversal']
            or (
                self.element.v_displace[2] == 0 or self.last_element.v_displace[2] == 0
                and (self.reversals == self.model_params.max_reversals_limit or self.traps == self.model_params.max_traps_limit)
            )
            or self.check_status_changed('has_surfaced')
            or self.check_status_changed('bottomed')
            or (self.status['is_lap'] != self.status['is_lap_old'])
            or self.element.dilution > self.model_params.max_dilution
            or self.check_status_changed('not_stream_limited')
            or (self.status['hitzid'] and not self.status['atvmxz'])  # (oc) ((mag(_pp(n.R_, k_)) > mag(An_)) and (not atvmxz))
            # (oc) or ((mag(_pp(n.R_, k_)) > c[vmxz]) and (not atvmxz))
            or (self.check_status_changed('reached_cig') and self.model_params.tpb_mixing_zone_ceil > 0)
            or (denomproduct < 0 < self.step and self.element.diameter > self.ambient.depth)
            # (oc) or ((magRc - magCc) / c[vdia] < 3.6)
        ):
            self.element.dilution = (self.element.mass / self.orig_element.mass) * self.orig_element.density / self.element.density
            self.element.density = self.element.density
            self.convert_vectors()
            if self.model_params.model != Model.DOS3:
                do_outputs = True

        if (
            self.iso_diameter == 0
            or (self.step%10 == 0 and self.step - self.last_graf_step > 5)
            or abs(self.last_diameter / self.iso_diameter - 1 > 0.1)
        ):
            # if step == 0:
            #     liftpens()
            self.last_diameter = self.iso_diameter
            self.last_graf_step = self.step
            do_graphs = True
        if self.step%15 == 0:
            do_graphs = True

        return do_outputs, do_graphs

    def vm_stop_condition(self, density_reversal):
        if (
            self.step > maxsteps
            or (self.step > 0 and density_reversal and self.model_params.max_traps_limit == 0)
            or self.traps >= self.model_params.max_traps_limit
            or self.reversals >= self.model_params.max_reversals_limit
            or (self.status['surfacing'] and not self.model_params.dont_stop_on_surface_hit)
            or (self.status['bottomed'] and self.model_params.stop_on_bottom_hit)
            or self.element.dilution > self.model_params.max_dilution
        ):
            self.status['done'] = True
        if self.step > 0 and self.um3isoplet <= 0:
            self._update_status('close_isoplet', True)
            self.status['done'] = True

    # semi-implementation of outputit in vmunit.pas
    # actual output functions are not transcribed here, but a few status checks are done during this output procedure
    def output_it(self, denomproduct=0):
        self.element.buildup = self.model_params.tpb_bincon

        if (
            self.step == 0
            and self.model_params.model == Model.NRFIELD
            and self.outputit.has_parameter('element', 'total_time')
        ):
            if self.step == 0:
                density_1 = seawater_density(
                    ambient_cond=self.ambient,
                    at_equilibrium=self.model_params.at_equilibrium,
                    in_sigma=False
                )
                density_2 = seawater_density(
                    ambient_cond=self.ambient_stack[0],
                    at_equilibrium=self.model_params.at_equilibrium,
                    in_sigma=False
                )
                dumie = (
                    (density_1 - density_2)
                    *math.sqrt(MAGNITUDE_GRAVITY / self.diff_params.depth / self.element.density)
                )
                self.element.total_time = 1e6 if dumie == 0 else math.pi/dumie
            elif self.step == 1 and self.ambient.current_speed != 0:
                self.element.total_time = magnitude(self.element.v_surface_tdsp) / math.sqrt(self.ambient.current_speed)

        # mixing zone graphing moved to update_status
        # if self.check_status_changed('atvmxz') and self.model_params.model != Model.PDS:
        #     # mixgraf()
        #     pass

        self.outputit.output(self, denomproduct)
        self.reset_status_changes()
