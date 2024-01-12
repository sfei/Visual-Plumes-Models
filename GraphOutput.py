import math
import numpy as np
from enum import Enum
from . import units
from .Element import Element
from .ambient.Ambient import Ambient
from .ambient.AmbientStore import AmbientStore
from .ambient.calculations import seawater_density, mancini
from .globals import UserInputError
from .params.DiffuserParameters import DiffuserParameters
from .params.DiffuserStore import DiffuserStore
from .params.ModelParameters import Model, ModelParameters
from .vectors import GRAVITY_VECTOR, project_vector, magnitude, make_vector, DEGREES_TO_RADIAN, unit_vector, \
    I_UNIT_VECTOR


class ProjectionPlane(Enum):
    EFFLUENT = "Effluent"
    CURRENT  = "Current"
    SET      = "Set"  # basically, custom


_TEMPLATE_AMBIENT  = Ambient()
_TEMPLATE_ELEMENT  = Element()
_TEMPLATE_DIFFUSER = DiffuserParameters()
_TEMPLATE_MODEL    = [
    'density_diff',
    'mass_pollutant',
    'um3isoplet',
    'total_time',
    'travel',
    'steps',
    'reverals'
]


class GraphOutput:
    """ Class for tracking coordinates from graph outputs. """

    def __init__(self, model_params, ambient_store, diffuser_store):
        """
        Args:
            model_params: instance of ModelParameters
            ambient_store: instance of AmbientStore
            diffuser_store: instance of DiffuserStore
        """
        assert isinstance(model_params, ModelParameters)
        assert isinstance(ambient_store, AmbientStore)
        assert isinstance(diffuser_store, DiffuserStore)

        # this controls whether to stop graphing plan/profile view things (the by-case stuff still graphs)
        self.stop_planprofile_graphs = False
        # Persistent variables
        self.ambient_store           = ambient_store
        self.diffuser_store          = diffuser_store
        self.depth_units             = self.ambient_store.z.units
        self.length_units            = self.diffuser_store.offset_x.units
        self.concentration_units     = self.diffuser_store.concentration.units
        # Whether elevation is given in depth or rise (formerly `exdepth`)
        self.graph_by_depth          = True   # TODO: where set from in UI
        self.negative_depth_axis     = False  # TODO: custom var, no equivalent in current UI
        self.plot_centerline         = False  # TODO: this gets prompted if custom_y is depth
        # Which plane in which graph perspective is angles from (o.c. as `projplane`)
        self.projplane               = ProjectionPlane.EFFLUENT # TODO: taken from mainform.cutplanerg
        self.org_dir                 = 0
        # If custom plane ('all' in o.c.), what angle this plane is at
        self.set_plane_deg           = 0      # TODO: taken from mainform.AllPlaneed.text
        self.set_plan_rad            = DEGREES_TO_RADIAN * self.set_plane_deg
        # Message settigns
        self.show_vert_dis_uam_msg   = True   # TODO: from main.pas (formerly shoVertDisUamMsg)
        self.show_vert_dis_const_msg = True   # TODO: from main.pas (formerly shoVertDisConstMsg)
        self.show_no_cur_const_msg   = True   # TODO: from main.pas (formerly shoNoCurConstMsg)
        # Whether to include other graphs
        self.custom_graph            = False  # (formerly `ifivx39`) custom graph
        # Custom axis settings
        self.custom_x                = None   # (formerly `ivx`) "Abscissa" in special settings
        self.custom_x_regime         = "element"
        self.custom_y1               = None   # (formerly `ivy`) "Ordinate-1" in special settings
        self.custom_y1_regime        = "element"
        self.custom_y2               = None   # (formerly `ivy1`) "Ordinate-2" in special settings
        self.custom_y2_regime        = "element"
        # Boundary stuff??
        # self.bdysize                 = 20000
        # self.bdycount                = -1
        # self.bdy                     = np.zeros((self.bdysize, 2))
        # Some sort of persistent vars
        # self._cl                     = [0, 0]
        # self._up                     = [0, 0]
        # self._ark                    = [0, 0]
        # self.gtmax                   = -round((math.pi-1)*1e8-1)
        # Prepare units and common axis metadata
        length_unit_meta1  = self._get_units_meta(units.Length, self.length_units, "Distance from Origin")
        length_unit_meta2  = self._get_units_meta(units.Length, self.length_units, "Distance from Source")
        z_axis_label       = "Depth" if self.graph_by_depth else "Rise"
        z_unit_meta        = self._get_units_meta(units.Length, self.depth_units, z_axis_label)
        dilution_unit_meta = self._get_units_meta(units.Unitless, units.Unitless.UNITLESS, "Dilution")
        case_unit_meta = {
            'units':       units.Unitless,
            'in_units':    units.Unitless.UNITLESS,
            'units_label': "",
            'axis_label':  "Case sequence"
        }
        xz_profile_view = {
            "x": length_unit_meta1,
            "y": z_unit_meta
        }
        xy_plan_view = {
            "x": self._get_units_meta(units.Length, self.length_units, "West-East"),
            "y": self._get_units_meta(units.Length, self.length_units, "South-North")
        }
        xy_density = {
            "x": self._get_units_meta(units.Density, units.Density.SIGMA_T, "Density"),
            "y": z_unit_meta
        }
        xy_dilution = {
            "x": length_unit_meta2,
            "y": dilution_unit_meta
        }
        xy_case_dilution = {
            "x": case_unit_meta,
            "y": dilution_unit_meta
        }
        xy_case_concentration = {
            "x": case_unit_meta,
            "y": self._get_units_meta(units.Concentration, self.concentration_units, "Concentration")
        }
        # Dictionary of series coordinates by name
        self.series = {
            # In "Graphical Output" > "Style" > "4pnl" > top-left
            "trajectory":       [],  # (formerly `serside`)  profile-view centerline along projection plane
            "boundary1":        [],  # (formerly `serbdy`)   profile-view boundary
            "boundary2":        [],  # (formerly ???)        profile-view boundary, other side
            # In "Graphical Output" > "Style" > "4pnl" > top-right
            "ambdensity":       [],  # (formerly `serden`)   ambient density profile by depth
            "density":          [],  # (formerly `serfden`)  plume trajectory density by depth
            # In "Graphical Output" > "Style" > "4pnl" > bottom-left
            "path":             [],  # (formerly `serpath`)  plan-view centerline (x/y displacement)
            "out1":             [],  # (formerly `serout`)   plan-view boundary (displacement + radius)
            "out2":             [],  # (formerly `bdy`)      plan-view boundary, other side (displacement - radius)
            "pathv":            [],  # (formerly `serpathv`) plan-view vector pointing towards shore
            # In "Graphical Output" > "Style" > "4pnl" > bottom-right
            "dilution":         [],  # (formerly `serdil`)   net dilution by distance
            "cldilution":       [],  # (formerly `sercldil`) element cl dilution by distance
            # In "Graphical Output" > "Style" > "diln"
            "enddilution":      [],  # (formerly `seridil`)  net dilution by case at end by case
            "mzdilution":       [],  # (formerly `sermixd`)  concentration at mixing zone depth by case
            # In "Graphical Output" > "Style" > "con"
            "concentration":    [],  # (formerly `sercon`)   concentration by case at end by case
            "mzconcentration":  [],  # (formerly `sermixc`)  net dilution at mixing zone depth by case
            # In "Graphical Output" > "Style" > "cus"
            "custom1":          [],  # (formerly `sergen`)   custom graph 1
            # In "Graphical Output" > "Style" > "cus"
            "custom2":          []   # (formerly `sergenr`)  custom graph 2
        }
        # Dictionary of series units by name
        self.units = {
            "trajectory":       xz_profile_view,
            "boundary1":        xz_profile_view,
            "boundary2":        xz_profile_view,
            "ambdensity":       xy_density,
            "density":          xy_density,
            "path":             xy_plan_view,
            "out1":             xy_plan_view,
            "out2":             xy_plan_view,
            "pathv":            xy_plan_view,
            "dilution":         xy_dilution,
            "cldilution":       xy_dilution,
            "enddilution":      xy_case_dilution,
            "mzdilution":       xy_case_dilution,
            "concentration":    xy_case_concentration,
            "mzconcentration":  xy_case_concentration,
            "custom1":          {"x": None, "y": None},
            "custom2":          {"x": None, "y": None},
        }

    @staticmethod
    def _get_units_meta(unit_type, in_units, axis_label):
        """ Get formatted units metadata for outputs dict.
        Args:
            unit_type: unit subclass
            in_units: unit value for unit type
            axis_label: axis label
        Returns: dict of
            units: unit subclass
            in_units: unit value for unit type
            units_label: units label
            axis_label: axis label
        """
        return {
            'units':       unit_type,
            'in_units':    in_units,
            'units_label': (
                "" if unit_type == units.Unitless else unit_type.label(
                    in_units if in_units is not None else 1
                )
            ),
            'axis_label':  axis_label
        }

    @staticmethod
    def _get_regime(regime):
        """ Get the regime template.
        Args:
            regime: regime name
        Returns: tuple of formatted regime name and template (a list of param names or an object)
        """
        regime = regime.lower().strip()
        match regime:
            case 'ambient':
                template = _TEMPLATE_AMBIENT
            case 'element':
                template = _TEMPLATE_ELEMENT
            case 'diffuser':
                template = _TEMPLATE_DIFFUSER
            case _:
                raise UserInputError(f"Invalid parameter regime ({regime})")
        return regime, template

    def _verify_custom_var(self, regime, varname, errorname):
        """ Verify custom variable for graphing is valid.
        Args:
            regime: regime name
            varname: variable name
            errorname: formatted/pretty var name for if UserInputError is invalid
        """
        assert regime and isinstance(regime, str)
        assert varname and isinstance(varname, str)
        regime, template = self._get_regime(regime)
        if isinstance(template, (list, tuple)):
            if varname not in template:
                raise UserInputError(f"Could not determine regime for custom {errorname} graph variable.")
        else:
            try:
                getattr(template, varname)
            except:
                raise UserInputError(f"Could not find custom {errorname} graph variable in regime.")
        real_template, units = self._get_custom_units(regime, varname)
        if not units:
            raise UserInputError(f"Could not determine units for custom {errorname} graph variable.")
        return units

    def set_custom_graph(self, x_regime=None, x_var=None, y1_regime=None, y1_var=None, y2_regime=None, y2_var=None):
        """ Set the custom graph variables.
        Args:
            x_regime: regime name for x-axis
            x_var: var name for x-axis
            y1_regime: regime name for y-axis
            y1_var: var name for y-axis
            y2_regime: regime name for y2-axis
            y2_var: var name for y2-axis
        """
        if x_var:
            x_units = self._verify_custom_var(x_regime, x_var, 'x')
        else:
            x_var = x_regime = x_units = None
        self.custom_x         = x_var
        self.custom_x_regime  = x_regime
        if y1_var:
            y1_units = self._verify_custom_var(y1_regime, y1_var, 'y1')
        else:
            y1_var = y1_regime = y1_units = None
        self.custom_y1        = y1_var
        self.custom_y1_regime = y1_regime
        self.units.custom1    = {"x": x_units, "y": y1_units}
        if y2_var:
            y2_units = self._verify_custom_var(y2_regime, y2_var, 'y2')
        else:
            y2_var = y2_regime = y2_units = None
        self.custom_y2_regime = y2_regime
        self.custom_y2        = y2_var
        self.units.custom2    = {"x": x_units, "y": y2_units}

        self.custom_graph = self.custom_x and (self.custom_y1 or self.custom_y2)

    def set_org_dir(self, element, v_ambient):
        """ Set the projection plane direction for profile graphs.
        Args:
            element: Element instance
            v_ambient: ambient current vector
        """
        message = None
        match self.projplane:
            case ProjectionPlane.EFFLUENT:
                if abs(np.dot(
                    element.v_velocity,
                    make_vector(x=element.v_velocity[0], y=element.v_velocity[1])
                )) < 1E-8:
                    if magnitude(v_ambient) > 0:
                        if self.show_vert_dis_uam_msg:
                            message = "Vertical discharge, projection plane defined by current"
                        self.org_dir = unit_vector(v_ambient)
                    elif self.show_vert_dis_const_msg:
                        message = "Vertical discharge, projection plane not changed"
                else:
                    self.org_dir = unit_vector(project_vector(element.v_velocity, GRAVITY_VECTOR))
            case ProjectionPlane.CURRENT:
                if magnitude(v_ambient) < 1E-8:
                    if self.show_no_cur_const_msg:
                        message = "No current case"
                else:
                    self.org_dir = unit_vector(v_ambient)
            case ProjectionPlane.SET:
                self.org_dir = unit_vector(make_vector(
                    x=math.cos(self.set_plan_rad),
                    y=math.sin(self.set_plan_rad)
                ))
            case _:
                pass
        return message

    def _convert(self, unit_type, in_units, value, max_length=2, max_depth=1, multiplier=None):
        """ Convert a value or values to desired units. Assumes value in default unit type (=1).
        Args:
            unit_type: the desired unit type value
            in_units: the unit subclass for the value and conversion
            value: the value to confirm from (in base/default units) (can be list)
            max_length: if value is list, maximum expected list length (default=2)
            max_depth: if value is list, maximum expected list depth/nesting (default=1)
            multiplier: optional additional multiplier after conversion
        """
        if isinstance(value, (list, tuple, np.ndarray)):
            max_depth -= 1
            if max_depth < 0:
                raise Exception("Invalid graph output coordinates (exceeding dimensionality).")
            if len(value) < max_length:
                raise Exception("Invalid graph output coordinates (below expected length).")
            if len(value) > max_length:
                value = value[:max_length]
            return tuple(map(
                lambda v: self._convert(unit_type, in_units, v, max_length, max_depth, multiplier),
                value
            ))
        if value is None or np.isnan(value):
            return value
        if multiplier is None:
            return unit_type.convert(value, ufrom=1, uto=in_units)
        else:
            return multiplier*unit_type.convert(value, ufrom=1, uto=in_units)

    def x(self, x, max_depth=0):
        """ Shorthand for calling _convert() on length values in the set length units for plan/profile graphs.
        Args:
            x: value to convert
            max_depth: maximum expected list depth/nesting (default=0)
        """
        return self._convert(units.Length, self.length_units, x, max_depth=max_depth)

    def z(self, z, max_depth=0):
        """ Shorthand for calling _convert() on depth values in the set length units for profile graphs.
        Args:
            z: value to convert
            max_depth: maximum expected list depth/nesting (default=0)
        """
        multiplier = -1 if self.graph_by_depth and self.negative_depth_axis else None
        return self._convert(units.Length, self.depth_units, z, max_depth=max_depth, multiplier=multiplier)

    def xz(self, xz):
        """ Shorthand for calling _convert() on x/z values in the set length units for profile graphs.
        Args:
            xz: list of x and z values
        """
        return tuple((self.x(xz[0]), self.z(xz[1])))

    def add_nulls(self):
        """ Add nulls to graphs (to indicate break and new line). """
        for name in (
            # not cleanest but following o.c., these are series mainly tracked in graph() where this is called from
            "trajectory", "boundary1", "boundary2",
            "path", "out1", "out2",
            "density", "ambdensity",
            "dilution", "cldilution",
            "custom1", "custom2"
        ):
            if name.startswith("custom") and not self.custom_graph:
                continue
            self.series[name].append((np.NaN, np.NaN))

    def graph_vector(self, start, end):
        """ Add vector to pathv series. """
        if len(self.series['pathv']):
            self.series['pathv'].append((np.NaN, np.NaN))
        self.series['pathv'].append(self.x(start, 1))
        self.series['pathv'].append(self.x(end, 1))

    def graph(self, umunit, casecount, net_dilution):
        """ Add next series to graphs.
        Args:
            umunit: instance of UMUnit
            casecount: case number
            net_dilution: net dilution
        """
        if self.stop_planprofile_graphs:
            return

        xy_pos = self.x(umunit.element.v_surface_tdsp, 1)
        if self.graph_by_depth:
            z_pos = self.z(umunit.element.depth)
        else:
            z_pos = self.z(umunit.diff_params.depth - umunit.element.depth)

        if umunit.step == 0:
            # break series lines for new case
            if casecount > 1:
                self.add_nulls()
            # plot ambient density profile
            for amb in umunit.ambient_stack:
                 self.series['ambdensity'].append((
                     seawater_density(
                         ambient_cond=amb,
                         at_equilibrium=umunit.model_params.at_equilibrium,
                         in_sigma=True
                     ),
                     self.z(amb.z)
                ))
            # plot initial trajectories at origin
            self.plot_profile(trajectory=(0, z_pos), convert_xz=False)
            # original was _d(n.R_, _sx(-1, n.A)) -- why difference of inverse? why not just add?
            element_span_1 = self.x(umunit.element.v_surface_tdsp + umunit.element.v_radius, 1)
            element_span_2 = self.x(umunit.element.v_surface_tdsp - umunit.element.v_radius, 1)
            self.plot_plan(
                path=xy_pos,
                out_1=element_span_1,
                out_2=element_span_2,
                convert_x=False
            )
            # initialize these vars
            # if self.graph_by_depth:
            #     self._cl = [0, z_pos]
            #     self._up = [0, z_pos]
            #     self._ark = [0, z_pos]
            # boundary thingy
            # self.bdycount += 1
            # if self.bdycount+1 < self.bdysize:
            #     self.bdy[self.bdycount][0] = missing
            #     self.bdy[self.bdycount][1] = missing
            #     self.bdycount += 1
            #     # this was a separate if statement with exact same condition
            #     # here actually difference, not sure why it differs yet
            #     element_span = umunit.element.v_surface_tdsp - umunit.element.v_radius
            #     element_span = self.x(element_span, 1)
            #     self.bdy[self.bdycount][0] = element_span[0]
            #     self.bdy[self.bdycount][1] = element_span[1]
            # exit function here for step=0
            return

        if umunit.element.cl_dilution < 1:
            self.plot_dilution(np.NaN, np.NaN, cl_series=True, convert_x=False)
        x_travel = self.x(magnitude(project_vector(
            (umunit.v_source - umunit.element.v_surface_tdsp),
            GRAVITY_VECTOR
        )))
        self.plot_dilution(x_travel, umunit.element.cl_dilution, convert_x=True, cl_series=True)
        self.plot_dilution(x_travel, net_dilution, convert_x=True)

        if umunit.model_params == Model.DOS3 and umunit.element.density > 800:
            the_density = umunit.element.density - 1000.0  # in sigma
        else:
            the_density = seawater_density(
                salinity=umunit.element.salinity,
                temperature=umunit.element.temperature,
                at_equilibrium=umunit.model_params.at_equilibrium,
                depth=umunit.element.depth,
                in_sigma=True
            )
        self.plot_density(the_density, z_pos, convert_z=False)

        dumx_r    = self.x(np.dot(umunit.element.v_surface_tdsp, self.org_dir))
        dumx_n    = self.x(np.dot(umunit.N_, self.org_dir))
        dumx_nopp = self.x(np.dot(umunit.Nopp_, self.org_dir))
        # dumx_max  = max(dumx_n, dumx_nopp)
        # if dumx_max > self.gtmax:
        #     self.gtmax = dumx_max
        nz = np.array([-umunit.N_[2], -umunit.Nopp_[2]])
        if self.graph_by_depth:
            # self.plot_profile(trajectory=self._cl, convert_xz=False)  # TODO: added in clwhite?
            # self._cl[0] = dumx_r
            # self._cl[1] = z_pos
            pass
        else:
            nz += umunit.element.depth
        nz = self.z(nz, 1)
        self.plot_profile(
            trajectory=(dumx_r, z_pos),
            boundary_1=(dumx_n, nz[0]),
            boundary_2=(dumx_nopp, nz[1]),
            convert_xz=False
        )

        # self.plot_profile(trajectory=self._up, convert_xz=False)  # TODO: added in clwhite?
        # self._up[0] = dumx_n
        # self._up[1] = -self.z(umunit.N_[2])
        # self.plot_profile(trajectory=self._up, convert_xz=False)
        #
        # self.plot_profile(trajectory=self._ark, convert_xz=False)  # TODO: added in clwhite?
        # self._ark[0] = dumx_nopp
        # self._ark[1] = -self.z(umunit.Nopp_[2])
        # self.plot_profile(trajectory=self._ark, convert_xz=False)

        # continue trajectory in plan view
        if math.sqrt(abs(umunit.element.cl_dilution / umunit.element.dilution)) > 1:  # (oc) Curt Dalton changes
            dumie = 1  # (oc) decay fix
        else:
            dumie = umunit.um3isoplet
        if umunit.element.diameter == 0:
            self.plot_plan(
                path=xy_pos,
                out_1=(np.NaN, np.NaN),
                out_2=(np.NaN, np.NaN),
                convert_x=False
            )
        else:
            v_radius_adj = (dumie*umunit.element.v_radius)
            element_span_1 = self.x(umunit.element.v_surface_tdsp + v_radius_adj, 1)
            element_span_2 = self.x(umunit.element.v_surface_tdsp - v_radius_adj, 1)
            self.plot_plan(
                path=xy_pos,
                out_1=element_span_1,
                out_2=element_span_2,
                convert_x=False
            )
            # self.bdycount += 1
            # if self.bdycount+1 < self.bdysize:
            #     self.bdy[self.bdycount][0] = calc_vector[0]
            #     self.bdy[self.bdycount][1] = calc_vector[1]

        if not self.custom_graph:
            return
        cvars = [np.nan, np.nan, np.nan]
        for i, (regime, varname, seriesname) in enumerate((
            (self.custom_x_regime, self.custom_x, None),
            (self.custom_y1_regime, self.custom_y1, "custom1"),
            (self.custom_y2_regime, self.custom_y2, "custom2")
        )):
            if not regime or not varname:
                continue
            custom_convert = False
            if varname == "decay_rate" and self.ambient_store.decay_rate.units == units.DecayRate.LY_PER_HOUR:
                custom_convert = True
            cvars[i] = self._get_custom_var(regime, varname, umunit, do_convert=(not custom_convert))
            if custom_convert:
                cvars[i] = mancini(
                    umunit.model_params.bacteria_model,
                    umunit.ambient_stack[0].decay_rate,
                    umunit.ambient_stack[0].salinity,
                    umunit.ambient_stack[0].temperature,
                    0,
                    False
                )
                cvars[i] *= math.exp(-umunit.model_params.light_absorb_coeff*umunit.ambient.depth)
            if seriesname is not None:
                self.series[seriesname].append((cvars[0], cvars[i]))

    def graph_mixing_zone(self, umunit, casecount, netdilution):
        """ Graph mixing zone dilution.
        Args:
            umunit: instance of UMUnit
            casecount: case number
            netdilution: net dilution
        """
        # from vmunit.mixgraf
        if umunit.element.concentration == 0:
            return
        self.series["mzdilution"].append((casecount, netdilution))
        if umunit.model_params.model != Model.PDS:
            self.series["mzconcentration"].append((
                casecount,
                units.Concentration.convert(umunit.element.concentration, self.concentration_units)
            ))
        if umunit.element.dilution < 1.0:
            self.plot_dilution(np.nan, np.nan, convert_x=False)
        if self.custom_graph:
            y1_is_dilution = self.custom_y1 == "dilution"
            y2_is_dilution = self.custom_y2 == "dilution"
            if y1_is_dilution or y2_is_dilution:
                custom_x = self.x((umunit.element.x_displacement**2 + umunit.element.y_displacement**2)**0.5)
                if y1_is_dilution:
                    self.series["custom1"].append((
                        custom_x,
                        self._get_custom_var(self.custom_y1_regime, self.custom_y1, umunit, do_convert=True)
                    ))
                if y2_is_dilution:
                    self.series["custom2"].append((
                        custom_x,
                        self._get_custom_var(self.custom_y2_regime, self.custom_y2, umunit, do_convert=True)
                    ))
        self.plot_dilution(
            magnitude(project_vector(umunit.v_source - umunit.element.v_surface_tdsp, GRAVITY_VECTOR)),
            netdilution,
            convert_x=True
        )
        if umunit.model_params.model != Model.PDS:
            self.plot_plan(path=(umunit.element.v_velocity[0], umunit.element.v_velocity[1]))

    def graph_end(self, umunit, casecount, netdilution):  # (oc) intended for plotting custom casecount graph
        """ Graph end of plume / stop condition dilution.
        Args:
            umunit: instance of UMUnit
            casecount: case number
            netdilution: net dilution
        """
        # from vmunit.graph4panel
        self.series["enddilution"].append((casecount, netdilution))
        self.series["concentration"].append((
            casecount,
            units.Concentration.convert(umunit.element.concentration, self.concentration_units)
        ))
        if self.stop_planprofile_graphs:
            return
        # from vmunit.graph39
        if self.custom_graph:
            if self.custom_y1:
                if self.custom_y1 == 'depth' and not self.plot_centerline:
                    arg = umunit.element.depth - magnitude(project_vector(umunit.element.v_radius, I_UNIT_VECTOR))
                    arg = units.Length.convert(arg, ufrom=1, uto=self.length_units)
                else:
                    arg = self._get_custom_var(self.custom_y1_regime, self.custom_y1, umunit, do_convert=True)
                self.series["custom1"].append((casecount, arg))
            if self.custom_y2:
                if self.custom_y2 == 'depth' and not self.plot_centerline:
                    arg = umunit.element.depth - magnitude(project_vector(umunit.element.v_radius, I_UNIT_VECTOR))
                    arg = units.Length.convert(arg, ufrom=1, uto=self.length_units)
                else:
                    arg = self._get_custom_var(self.custom_y2_regime, self.custom_y2, umunit, do_convert=True)
                self.series["custom2"].append((casecount, arg))

    def plot_ff_concentration(self, v_displacement, add_bdy):
        """ Plot far-field concentration? Copied name from o.c. actually plots displacement.
        Args:
            v_displacement: vector displacement
            add_bdy: add boundary buffer
        """
        # from wfar.plotconc() -- some parts moved into BrooksFarField.py
        self.plot_plan(
            path=v_displacement,
            out_1=(v_displacement+add_bdy),
            out_2=(v_displacement-add_bdy)
        )

    def plot_profile(self, trajectory=None, boundary_1=None, boundary_2=None, convert_xz=True):
        """ Plot various profile graphs.
        Args:
            trajectory: trajectory coordinates (as x,z) (or None to skip)
            boundary_1: boundary coordinates (as x,z) (or None to skip)
            boundary_2: second boundary coordinates (as x,z) (or None to skip)
            convert_xz: convert units or false if already converted (default=True)
        """
        if self.stop_planprofile_graphs:
            return
        if trajectory is not None:
            self.series["trajectory"].append(self.xz(trajectory) if convert_xz else trajectory)
        if boundary_1 is not None:
            self.series["boundary1"].append(self.xz(boundary_1) if convert_xz else boundary_1)
        if boundary_2 is not None:
            self.series["boundary2"].append(self.xz(boundary_2) if convert_xz else boundary_2)

    def plot_plan(self, path=None, out_1=None, out_2=None, convert_x=True):
        """ Plot various plan graphs.
        Args:
            path: path coordinates (as x,y) (or None to skip)
            out_1: boundary coordinates (as x,y) (or None to skip)
            out_2: second boundary coordinates (as x,y) (or None to skip)
            convert_x: convert units or false if already converted (default=True)
        """
        if self.stop_planprofile_graphs:
            return
        if path is not None:
            self.series["path"].append(self.x(path, 1) if convert_x else path)
        if out_1 is not None:
            self.series["out1"].append(self.x(out_1, 1) if convert_x else out_1)
        if out_2 is not None:
            self.series["out2"].append(self.x(out_2, 1) if convert_x else out_2)

    def plot_dilution(self, displacement, dilution, convert_x=True, cl_series=False):
        """ Plot dilution profile graph.
        Args:
            displacement: x-displacement
            dilution: dilution value
            convert_x: convert displacement units or false if already converted (default=True)
            cl_series: true if cl-dilution otherwise normal dilution (default=False)
        """
        if self.stop_planprofile_graphs:
            return
        sname = "cldilution" if cl_series else "dilution"
        self.series[sname].append((
            self.x(displacement) if convert_x else displacement,
            dilution
        ))

    def plot_density(self, density, depth, convert_z=True):
        """ Plot density profile graph.
        Args:
            density: density value
            depth: depth value
            convert_z: convert depth units or false if already converted (default=True)
        """
        if self.stop_planprofile_graphs:
            return
        self.series["density"].append((
            density,
            self.z(depth) if convert_z else depth
        ))

    def _get_custom_units(self, regime, varname, umunit=None):
        """ Get units for a custom variable.
        Args:
            regime: regime name
            varname: var name
            umunit: instance of UMUnit or None
        """
        match regime:
            case 'diffuser':
                template = umunit.diff_params if umunit else None
                to_units = self.diffuser_store.get(varname).units
            case 'ambient':
                template = umunit.ambient if umunit else None
                if varname == 'density':
                    # for now, density standardized in kg/m3
                    to_units = units.Density.KILOGRAMS_PER_CUBIC_METER
                else:
                    to_units = self.ambient_store.get(varname).units
            case 'element':
                template = umunit.element if umunit else None
                if varname in ('height', 'diameter'):
                    to_units = self.diffuser_store.diameter.units
                elif varname in ('depth', 'x_displacement', 'y_displacement'):
                    to_units = self.length_units
                elif varname in ('concentration', 'v4o3'):
                    to_units = self.concentration_units
                elif varname in ('vertical_angle', 'horizontal_angle', 'temperature', 'salinity'):
                    to_units = self.diffuser_store.get(varname).units
                elif varname == 'density':
                    # for now, density standardized in kg/m3
                    to_units = units.Density.KILOGRAMS_PER_CUBIC_METER
                elif varname == 'speed':
                    to_units = self.ambient_store.current_speed.units
                elif varname == 'total_time':
                    to_units = self.diffuser_store.start_time.units
                else:
                    # mass/buildup isn't handled, left kilograms
                    # dilution/vcld left as unitless
                    to_units = None
            case _:
                template = to_units = None
        return template, to_units

    def _get_custom_var(self, regime, varname, umunit, do_convert=True):
        """ Get value for custom variable.
        Args:
            regime: regime name
            varname: var name
            umunit: instance of UMUnit or None
            do_convert: convert units (default=True)
        """
        template, to_units = self._get_custom_units(regime, varname, umunit)
        if not template:
            return np.NaN
        val = template.get(varname)
        if do_convert and to_units:
            val = units.convert(
                val,
                varname,
                from_units=1,
                to_units=to_units,
                model_params=umunit.model_params,
                celsius=umunit.element.temperature,
                psu=umunit.element.salinity,
                depth=umunit.element.depth
            )
        return val
