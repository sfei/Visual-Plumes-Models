import math
from .. import units
from ..globals import missing, UserInputError
from .. import helpers
from ..params.ModelParameters import ModelParameters
from .Ambient import Ambient
from .AmbientStore import AmbientStore, Interpolation
from .calculations import seawater_density, mancini
from ..vectors import make_vector, DEGREES_TO_RADIAN, vector_average


class AmbientHandler:

    def __init__(self, model_params, ambient_store):
        assert isinstance(model_params, ModelParameters)
        assert isinstance(ambient_store, AmbientStore)
        self.model_params  = model_params
        self.ambient_store = ambient_store
        self.ambient_stack = None
        self.ts_amb_stacks = {}
        # self.bottom        = None
        self.Qstream       = 0.0
        self._was_filled   = False

    def fill(self, model_params, ambient_stack, ambient_ts_stacks, diff_params=None, orig_ambient=None):
        """
        Set up the handler by providing the ambient stack of ambient conditions by depth layers. Conditions at the depth
        layers will be appropriately processed after copying.
        :param model_params:       The ModelParameters.
        :param ambient_stack:      The ambient stack (as a list of Ambient ordered in increasing depth).
        :param ambient_ts_stacks:  Additional ambient conditions pulled from timeseries (as map of stacks per var type)
        :param diff_params:        Optional DiffuserParameters (just used to pull depth for QStream calc).
        :param orig_ambient:       Optional Ambient to copy initial values.
        """
        # This has been reworked somewhat to separate the UI from model. Assumes records is a list of Ambient instances sent
        # by UI, for a given ambient case, in order of top to bottom.
        self.ambient_stack = [r.copy() for r in ambient_stack]
        # convert z to depths
        last_z = -9999
        for i, amb in enumerate(self.ambient_stack):
            amb.depth = self.get_z(amb.z, model_params.bottom_depth)
            if amb.depth <= last_z:
                raise UserInputError(f"Ambient layer out of order (f{i}:f{amb.z}). Ambient layers must be in order of increasing depth.")
            last_z = amb.depth
        record_count = self._validate_records(self.ambient_stack)

        self.ts_amb_stacks = {}
        if ambient_ts_stacks:
            for vkey in self.ambient_store._input_vars_:
                store = self.ambient_store.get(vkey)
                if not store.from_time_series:
                    continue
                if vkey not in ambient_ts_stacks:
                    raise UserInputError(f"Ambient timeseries data indicated but not supplied for {vkey}")
                self.ts_amb_stacks[vkey] = [r.copy() for r in ambient_ts_stacks[vkey]]
                # convert z to depths
                for amb in self.ts_amb_stacks[vkey]:
                    amb.depth = self.get_z(
                        amb.z,
                        model_params.bottom_depth,
                        z_is_depth=store.z_is_depth,
                        z_units=store.ts_depth_units
                    )
                self._validate_records(self.ts_amb_stacks[vkey])

        # first pass fills/interpolates in original units
        for ambient_cond in self.ambient_stack:
            self._fill_level_1(model_params, ambient_cond)
        # second pass converts to standard units, calcs density, sets bottom/orig conditions
        for ambient_cond in self.ambient_stack:
            self._fill_level_2(model_params, ambient_cond, orig_ambient)
            ambient_cond.density = seawater_density(
                at_equilibrium=self.model_params.at_equilibrium,
                ambient_cond=ambient_cond,
                in_sigma=False
            )
        # self._fill_level_2(model_params, self.bottom)
        # self.bottom.density = seawater_density(
        #     at_equilibrium=self.model_params.at_equilibrium,
        #     ambient_cond=self.bottom,
        #     in_sigma=False
        # )

        # for i, ambient_cond in enumerate(self.ambient_stack):
        #     # keep overwritting bottom conditions from `below` layer returned
        #     self.bottom = self._fill_level(model_params, ambient_cond, orig_ambient)
        #     ambient_cond.density = seawater_density(
        #         at_equilibrium=self.model_params.at_equilibrium,
        #         ambient_cond=ambient_cond,
        #         in_sigma=False
        #     )
        #     # TODO: what does this achieve?
        #     # pre_z = z
        #     # while ini_z == pre_z:
        #     #     if ambient_cond.z is not None and ambient_cond.z != missing:
        #     #         ini_z = units.Length.convert(ambient_cond.z, self.ambient_store.z.units, units.Length.METERS)

        # TODO: Qstream was in loop but just kept getting overwritten until the last..
        # TODO: handler timeseries current speeds
        self.Qstream = model_params.tpb_channel_width
        if record_count > 1:
            self.Qstream *= 0.5*(self.ambient_stack[-1].current_speed + self.ambient_stack[-2].current_speed)
        else:
            self.Qstream *= self.ambient_stack[0].current_speed
        if diff_params and diff_params.depth > self.ambient_stack[-1].depth:
            self.Qstream *= diff_params.depth
        else:
            self.Qstream *= self.ambient_stack[-1].depth
        #ambbufi = i
        self._was_filled = True

    # height adjustment by depth or height option
    def get_z(self, depth_value, bottom_depth=1e9, z_is_depth=None, z_units=None):
        if depth_value is None:
            return None
        if z_units is None:
            z_units = self.ambient_store.z.units
        if z_is_depth is None:
            z_is_depth = self.ambient_store.z.z_is_depth
        z = units.Length.convert(depth_value, ufrom=z_units)
        return z if z_is_depth else bottom_depth - z

    def _validate_records(self, ambient_records):
        record_count = len(ambient_records)
        if not record_count:
            raise UserInputError("No ambient conditions defined")
        if record_count < 2:
            raise UserInputError("Ambient water column must have at least 2 levels")
        # if self.get_z(ambient_records[0].z) != 0:
        #     raise UserInputError("Top row depth should be zero, please correct")
        if ambient_records[0].depth == ambient_records[-1].depth:
            raise UserInputError("Surface and bottom depths are equal")
        # new check to assure depths are properly ordered
        last_depth = 0
        for ambient_row in ambient_records:
            if ambient_row.z is None:
                raise UserInputError("Invalid ambient depth. No value provided.")
            if ambient_row.depth is None:
                raise UserInputError("Invalid ambient depth. Not depth converted.")
            if ambient_row.depth < 0:
                raise UserInputError("Invalid ambient depth. Negative depth.")
            if ambient_row.depth < last_depth:
                raise UserInputError("Depth layers must be supplied in increasing order from surface to bottom")
            last_depth = ambient_row.depth
        return record_count

    def _fill_level_1(self, model_params, ambient_cond):
        """
        Fill ambient condition values at the given ambient layer. Interpolations are handled for values at the
        sandwiching depth layers (for values which are null/None), then the exact value at the target depth linearly
        interpolated between. Though since this is only used by fill, which processes the ambient stack, I think the
        linear interpolation is a moot process since either the above or below layer should be exact. NO UNIT
        CONVERSIONS YET. That's handled in part 2. Must be split to avoid converting before blanks are filled.
        """
        depth = ambient_cond.depth
        bottom_depth = model_params.bottom_depth
        # get the values (interpolated if necessary) for the nearest sandwiching depth layers
        above, below = self._ambient_limits(depth, bottom_depth)
        if below.depth == above.depth:
            raise UserInputError("Same bracketing water levels!")

        # set the ambient values at the diffuser origin (pulled in from ambientvalues)
        # while also checking for any changed values from original (part of fillit)
        # TODO: may have to be split again later
        ambient_cond.depth = depth
        ambient_cond.z = ambient_cond.depth  # should not use unconverted `z` after fill() was called, but just in case
        z_factor = (depth - below.depth) / (above.depth - below.depth)

        for vkey in self.ambient_store._input_vars_:
            ambient_store_at = self.ambient_store.get(vkey)
            if ambient_store_at.from_time_series:
                ts_above, ts_below = self._series_limits(
                    ambient_store_at,
                    self.ts_amb_stacks[vkey],
                    depth,
                    bottom_depth
                )
                above_val = ts_above.get(vkey)
                below_val = ts_below.get(vkey)
                use_z_factor = (depth - ts_below.depth) / (ts_above.depth - ts_below.depth)
            else:
                above_val = above.get(vkey)
                below_val = below.get(vkey)
                use_z_factor = z_factor
            # convert all input values to default units
            # exception for decay rate (not in ly/hr) -- left as is, conversion done at level-specific value pull
            if vkey != "decay_rate" or self.ambient_store.decay_rate.units != units.DecayRate.LY_PER_HOUR:
                ambient_cond.set(vkey, (below_val + (above_val - below_val) * use_z_factor))

    def _fill_level_2(self, model_params, ambient_cond, orig_ambient=None):
        """
        Converts ambient condition values at the given ambient layer into standard units. Exception for decay rate not
        in ly/hr -- left as is, conversion done at level-specific value pull.
        """
        for vkey in self.ambient_store._input_vars_:
            # exception for decay rate (not in ly/hr) -- left as is, conversion done at level-specific value pull
            if vkey != "decay_rate" or self.ambient_store.decay_rate.units != units.DecayRate.LY_PER_HOUR:
                from_units = self.ambient_store.get(vkey).units
                converted = units.convert(
                    ambient_cond.get(vkey),
                    units_or_var_name=vkey,
                    from_units=from_units,
                    model_params=model_params,
                    celsius=ambient_cond.temperature,
                    psu=ambient_cond.salinity,
                    depth=ambient_cond.depth
                )
                if orig_ambient:
                    orig_ambient.set(vkey, converted)
                ambient_cond.set(vkey, converted)

    def _fill_level(self, model_params, ambient_cond, orig_ambient):
        """
        Fill ambient condition values at the given depth. Interpolations are handled for values at the sandwiching depth
        layers (for values which are null/None), then the exact value at the target depth linearly interpolated between.
        Though since this is only used by fill, which processes the ambient stack, I think the linear interpolation is a
        moot process since either the above or below layer should be exact.
        :param umunit:          The UMUnit instance (for pulling various model and diffuser parameters).
        :param ambient_cond:    The ambient conditions for the given depth (copied in handler).
        :param depth:           The depth.
        :param bottom_depth:    The bottom depth/bathymetry.
        :return:                The values for the bottom/below sandwiching depth layer
        """
        depth = ambient_cond.depth
        bottom_depth = model_params.bottom_depth
        # get the values (interpolated if necessary) for the nearest sandwiching depth layers
        above, below = self._ambient_limits(depth, bottom_depth)
        if below.depth == above.depth:
            raise UserInputError("Same bracketing water levels!")

        # set the ambient values at the diffuser origin (pulled in from ambientvalues)
        # while also checking for any changed values from original (part of fillit)
        # TODO: may have to be split again later
        ambient_cond.depth = depth
        ambient_cond.z = ambient_cond.depth  # should not use unconverted `z` after fill() was called, but just in case
        z_factor = (depth - below.depth) / (above.depth - below.depth)
        # changed = False

        def fill_var(vkey):
            nonlocal self, depth, bottom_depth, z_factor, ambient_cond, orig_ambient
            ambient_store_at = self.ambient_store.get(vkey)
            if ambient_store_at.from_time_series:
                ts_above, ts_below = self._series_limits(
                    ambient_store_at,
                    self.ts_amb_stacks[vkey],
                    depth,
                    bottom_depth
                )
                above_val = ts_above.get(vkey)
                below_val = ts_below.get(vkey)
                from_units = ambient_store_at.units
                use_z_factor = (depth - ts_below.depth) / (ts_above.depth - ts_below.depth)
            else:
                above_val = above.get(vkey)
                below_val = below.get(vkey)
                from_units = ambient_store_at.units
                use_z_factor = z_factor
            # convert all input values to default units
            # exception for decay rate (not in ly/hr) -- left as is, conversion done at level-specific value pull
            if vkey != "decay_rate" or self.ambient_store.decay_rate.units != units.DecayRate.LY_PER_HOUR:
                converted = units.convert(
                    (below_val + (above_val - below_val) * use_z_factor),
                    units_or_var_name=vkey,
                    from_units=from_units,
                    model_params=model_params,
                    celsius=ambient_cond.temperature,
                    psu=ambient_cond.salinity,
                    depth=depth
                )
                if orig_ambient:
                    orig_ambient.set(vkey, converted)
                ambient_cond.set(vkey, converted)
            # check if changed versus original (I don't think this is used anywhere..)
            # if not changed:
            #     og_val = units.convert(
            #         og_ambient_cond.get(vkey),
            #         var_name=vkey,
            #         from_units=from_units,
            #         umunit=umunit,
            #         z=z
            #     )
            #     if abs(og_val - converted) > 1e-6:
            #         changed = True

        fill_var('temperature')  # cause salinity conversion requires temperature to be in celsius
        fill_var('salinity')     # cause decay rate also requires salinity to be in PSU
        # do rest of conversions
        for vkey in self.ambient_store._input_vars_:
            if vkey not in ('salinity', 'temperature'):
                fill_var(vkey)

        # TODO: my understanding is that `bottom` keeps getting overwritten as fillit is called, but since it processes
        # only once from surface to bottom, the last filled value becomes the permanent 'bottom' as used in interpolate
        # self.bottom = below
        # TODO: move bottom being set outside of this function, also is changed actually even used anywhere?
        # return changed
        return below

    def _ambient_limits(self, depth, bottom_depth):
        """
        Get the nearest over/under ambient condition layers with valid values that sandwich a given z-value.
        :param depth:        he target depth.
        :param bottom_depth: The bottom depth/bathymetry.
        :return: Pair of Ambient conditions for values under and over
        """
        ambient_records = [
            r.to_dict() if isinstance(r, Ambient) else r
            for r in self.ambient_stack
        ]

        above = Ambient()  # formerly oben
        below = Ambient()  # formerly unten
        # these represent the depth of nearest sandwiching layers with non-null values, whereas above['z'] and
        # below['z'] are simply the nearest sandwiching layers with may be null for value of interest
        above_value_dep = below_value_dep = 0  # formerly obva and unva

        # get ambient layers above and below target depth (or if exact exists, at and above -- or and below if at surf)
        for i, ambient_row in enumerate(self.ambient_stack):
            above.depth = below.depth
            below.depth = ambient_row.depth
            if i and depth <= below.depth:
                break

        # fill in ambient values
        for vkey in self.ambient_store._input_vars_:
            if self.ambient_store.get(vkey).from_time_series:
                continue

            record_iterator = iter(ambient_records)
            ambient_row = next(record_iterator)

            # get first valid value
            if ambient_row.get(vkey) is None:
                while (ambient_row := next(record_iterator, None)) is not None:
                    if ambient_row.get(vkey) is not None:
                        break
            if ambient_row is None:
                # TODO: there is a checkset (see procedure makecheckset in main.pas) that determines which vars are
                # required by model selected. Not implementing for now. UM3 needs all vars so not needed if only
                # considering that
                raise UserInputError(f"Ambient values for col={vkey} does not have a single, valid value")
            above_value = below_value = ambient_row[vkey]
            above_value_dep = below_value_dep = ambient_row['depth']
            if below_value_dep <= depth:
                # get sandwiching (under/over) non-null values and depths of each
                while (ambient_row := next(record_iterator, None)) is not None:
                    value_at = ambient_row.get(vkey)
                    if value_at is not None:
                        above_value = below_value
                        below_value = value_at
                        above_value_dep = below_value_dep
                        below_value_dep = ambient_row['depth']
                        if depth == below_value_dep:
                            # exact depth match
                            above_value = below_value
                            above_value_dep = below_value_dep
                            break
                        if depth < below_value_dep:
                            # break when target depth is between
                            break

            # (oc) Find ambient values at oben[valen] (& unten)
            above_value, below_value = self._interpolate_value(
                depth,
                bottom_depth,
                self.ambient_store.get(vkey),
                above,
                below,
                above_value,
                below_value,
                above_value_dep,
                below_value_dep
            )
            above.set(vkey, above_value)
            below.set(vkey, below_value)

        # pin above/below to surface/bottom if not sandwiching target depth
        if depth < above.depth:
            below.depth = above.depth
            above.depth = 0
        if depth > below.depth:
            above.depth = below.depth
            below.depth = bottom_depth

        return above, below

    def _series_limits(self, ambient_store_at, ambient_ts_stack, depth, bottom_depth):
        """
        Get the nearest over/under ambient condition layers with valid values that sandwich a given z-value.
        :param ambient_store_at: The AmbientStoreSubset for the ambient variable.
        :param ambient_ts_stack: The stack of Ambient conditions for the ambient variable.
        :param depth:            The target depth.
        :param bottom_depth:     The bottom depth/bathymetry.
        :return: Pair of Ambient conditions for values under and over
        """
        # I reused the same logic mostly from ambient limits, not sure why it would have to be any different except
        # using a different ambient and depth stacks
        above = Ambient()  # formerly tsob
        below = Ambient()  # formerly tsunt
        # these represent the depth of nearest sandwiching layers with non-null values, whereas above['z'] and
        # below['z'] are simply the nearest sandwiching layers with may be null for value of interest
        above_value_dep = below_value_dep = 0  # formerly obva and unva

        # get ambient layers above and below target depth (or if exact exists, at and above -- or and below if at surf)
        for i, ambient_row in enumerate(ambient_ts_stack):
            above.depth = below.depth
            below.depth = ambient_row.depth
            if i and depth <= below.depth:
                break

        # fill in ambient values
        above.set(ambient_store_at.var_key, missing)
        below.set(ambient_store_at.var_key, missing)

        record_iterator = iter(ambient_ts_stack)
        ambient_row = next(record_iterator)

        # get first valid value
        if ambient_row.get(ambient_store_at.var_key) is None:
            while (ambient_row := next(record_iterator, None)) is not None:
                if ambient_row.get(ambient_store_at.var_key) is not None:
                    break
        if ambient_row is None:
            raise UserInputError(f"Empty ambient timeseries data ({ambient_store_at.var_key})")
        value_at = ambient_row.get(ambient_store_at.var_key)
        above_value = below_value = value_at
        above_value_dep = below_value_dep = ambient_row.depth

        if below_value_dep <= depth:
            # get sandwiching (under/over) non-null values and depths of each
            while (ambient_row := next(record_iterator, None)) is not None:
                value_at = ambient_row.get(ambient_store_at.var_key)
                if helpers.is_truthy(value_at, zero_is_true=True):
                    assert isinstance(value_at, (int, float))
                    above_value = below_value
                    below_value = value_at
                    above_value_dep = below_value_dep
                    below_value_dep = ambient_row.depth
                    if depth == below_value_dep:
                        # exact depth match
                        above_value = below_value
                        above_value_dep = below_value_dep
                        break
                    if depth < below_value_dep:
                        # break when target depth is between
                        break

        above_value, below_value = self._interpolate_value(
            depth,
            bottom_depth,
            ambient_store_at,
            above,
            below,
            above_value,
            below_value,
            above_value_dep,
            below_value_dep
        )
        above.set(ambient_store_at.var_key, above_value)
        below.set(ambient_store_at.var_key, below_value)

        # pin above/below to surface/bottom if not sandwiching target depth
        if depth < above.depth:
            below.depth = above.depth
            above.depth = 0
        if depth > below.depth:
            above.depth = below.depth
            below.depth = bottom_depth

        return above, below

    def _interpolate_value(self, depth, bottom_depth, ambient_store_at, above, below, above_value, below_value,
                           above_value_dep, below_value_dep):
        # exact level match for values
        if above_value_dep == below_value_dep:  # (oc) constant ambient except linear to zero
            if depth < above_value_dep and ambient_store_at.extrapolation_sfc == Interpolation.LINEAR_TO_ZERO:  # (oc) to zero
                # both matches deeper than target depth, interpolate to surface 0
                above_value *= above.depth / above_value_dep
                below_value *= below.depth / below_value_dep
                if depth < above.depth:
                    below_value = above_value
                    above_value = 0
            elif depth > below_value_dep and ambient_store_at.extrapolation_btm == Interpolation.LINEAR_TO_ZERO:  # (oc) to zero
                # both matches shallower than target depth, interpolate to bottom 0
                rate = above_value / (above_value_dep - bottom_depth)
                below_value = rate * (below.depth - bottom_depth)
                above_value = rate * (above.depth - bottom_depth)
                if depth > below.depth:
                    above_value = below_value
                    below_value = 0
            else:  # (oc) multiple levels
                # use exact match values
                pass
        else:
            ratio = (above_value - below_value) / (above_value_dep - below_value_dep)
            if depth < above_value_dep:  # (oc) extrapolate
                # both matches deeper than target depth, interpolate to surface
                match ambient_store_at.extrapolation_sfc:
                    case Interpolation.CONSTANT:
                        below_value = above_value  # (oc) constant
                    case Interpolation.EXTRAPOLATED:  # (oc) extrapolate
                        dumie = above_value + ratio * (above.depth - above_value_dep)
                        below_value = above_value + ratio * (below.depth - above_value_dep)
                        above_value = dumie
                    case Interpolation.LINEAR_TO_ZERO:
                        dumie = above_value * above.depth / above_value_dep
                        below_value = above_value * below.depth / above_value_dep
                        above_value = dumie
            elif depth > below_value_dep:  # (oc) extrapolate
                # both matches shallower than target depth, interpolate to bottom
                match ambient_store_at.extrapolation_btm:
                    case Interpolation.CONSTANT:
                        above_value = below_value  # (oc) constant
                    case Interpolation.EXTRAPOLATED:  # (oc) extrapolate
                        dumie = above_value + ratio * (above.depth - above_value_dep)
                        below_value = above_value + ratio * (below.depth - above_value_dep)
                        above_value = dumie
                    case Interpolation.LINEAR_TO_ZERO:
                        dumie = below_value * (above.depth - bottom_depth) / (below_value_dep - bottom_depth)
                        below_value *= (below.depth - bottom_depth) / (below_value_dep - bottom_depth)
                        above_value = dumie
            else:  # (oc) interpolate
                # match value to nearest sandwiching depths, linearly interpolate where necessary
                if above.depth != above_value_dep:
                    above_value += (below_value - above_value) * (above.depth - above_value_dep) / (
                            below_value_dep - above_value_dep)
                    above_value_dep = above.depth
                if below.depth != below_value_dep:
                    below_value += ratio * (below.depth - below_value_dep)
        return above_value, below_value

    def ambient_level(self, umunit, ambient_cond, depth, bottom_depth):
        """
        Fills a given ambient condition object for a given depth.
        :param ambient_cond: The Ambient object to fill values for.
        :param depth:        The depth at which to interpolate values for.
        :Param bottom_depth: The bottom/bathymetry depth.
        :return:
        """
        assert self._was_filled
        assert isinstance(ambient_cond, Ambient)
        if depth is None:
            depth = ambient_cond.depth
        ambient_cond.z = ambient_cond.depth  # should not use unconverted `z` after fill() was called, but just in case
        self.interpolate_ambient(
            umunit=umunit,
            ambient_cond=ambient_cond,
            depth=depth,
            bottom_depth=bottom_depth
        )
        ambient_cond.density = seawater_density(
            at_equilibrium=self.model_params.at_equilibrium,
            ambient_cond=ambient_cond,
            in_sigma=False
        )
        v_amb_current = make_vector(
            x=ambient_cond.current_speed * math.cos(DEGREES_TO_RADIAN * ambient_cond.current_dir),
            y=ambient_cond.current_speed * math.sin(DEGREES_TO_RADIAN * ambient_cond.current_dir),
            z=0
        )
        if ambient_cond.decay_rate < 0:
            ambient_cond.decay_rate = 0
        if self.ambient_store.decay_rate.units != units.DecayRate.LY_PER_HOUR:
            ambient_cond.kt = ambient_cond.decay_rate
        else:
            # conversion from ly/hr is done here (as opposite to values in stack)
            # TODO: why ly per hour doesn't consider depth layers, instead based only off surface rate?
            ambient_cond.kt = mancini(
                self.model_params.bacteria_model,
                mancini(
                    self.model_params.bacteria_model,
                    self.ambient_stack[0].decay_rate,
                    self.ambient_stack[0].salinity,
                    self.ambient_stack[0].temperature,
                    0,
                    False
                ),
                ambient_cond.salinity,
                ambient_cond.temperature,
                depth,
                True
            )
            if ambient_cond.kt < 0:
                ambient_cond.kt = 0
        return v_amb_current

    def interpolate_ambient(self, umunit, ambient_cond, depth, bottom_depth):
        """
        Fills a given ambient condition object for a given depth. This subfunction in particular handles interpolating
        values between the under/over layers.
        :param ambient_cond: The Ambient object to fill values for.
        :param depth:        The depth at which to interpolate values for.
        :Param bottom_depth: The bottom/bathymetry depth.
        """
        assert self._was_filled
        assert isinstance(ambient_cond, Ambient)
        over_above = True
        surfacing = False
        i = 0
        if depth <= 0:
            if umunit is not None:
                surfacing = umunit.status['surfacing'] = True
            depth = 0
        elif depth >= self.ambient_stack[0].depth:
            over_above = False
            i = len(self.ambient_stack)-1
            for n, amb in enumerate(self.ambient_stack):
                # TODO: I don't understand need for the second condition exactly
                if depth <= amb.depth or (n > 0 and amb.depth == 0):
                    i = n
                    break
        below        = self.ambient_stack[i if not over_above else 1]
        below_under  = not over_above and depth > below.depth
        below_bottom = below_under and depth >= bottom_depth
        above        = self.ambient_stack[i-1 if not over_above else 1]
        z_factor     = (depth - above.depth)/(below.depth - above.depth) if below.depth != above.depth else 0
        # if temp_only:
        #     # TODO: why doesn't temp only handle different extrapolation types? or timeseries data? is it even used?
        #     if depth <= 0:  # (oc) whitman bug
        #         ambient_cond.temperature = self.ambient_stack[0].temperature  # (oc) 4 Dec 2005 correction of Whitman bug change
        #     elif z_factor:
        #         ambient_cond.temperature = above.temperature + z_factor*(below.temperature - above.temperature)
        #     else:
        #         ambient_cond.temperature = below.temperature
        #     return

        def fill_var(vkey):
            nonlocal self, above, below, over_above, below_under, surfacing, z_factor, depth, umunit, ambient_cond
            do_conversion = False
            use_above     = above
            use_below     = below
            is_over       = over_above
            is_below      = below_under
            is_bottom     = below_bottom
            use_z_factor  = z_factor

            # timeseries may have their own/different depth stacks
            ambient_store_at = self.ambient_store.get(vkey)
            if ambient_store_at.from_time_series:
                do_conversion = ambient_store_at.units if ambient_store_at.units != 1 else 0
                ts_stack = self.ts_amb_stacks[vkey]
                if surfacing or depth < ts_stack[0].depth:
                    is_over   = True
                    is_below  = is_bottom = False
                    use_above = ts_stack[0]
                    use_below = ts_stack[1]
                else:
                    is_over = False
                    for amb in ts_stack:
                        use_above = use_below
                        use_below = amb
                        if depth <= use_below.depth:
                            break
                    is_below = depth > use_below.depth
                    is_bottom = is_below and depth > bottom_depth
                use_z_factor = (
                    (depth - use_above.depth) / (use_below.depth - use_above.depth)
                    if use_below.depth != use_above.depth
                    else 0
                )

            above_val = use_above.get(vkey)
            below_val = use_below.get(vkey)

            if is_below:
                match ambient_store_at.extrapolation_btm:  # (oc) cerucci addition; need z < oben[valen]?
                    case Interpolation.CONSTANT:
                        set_val = below_val
                    case Interpolation.EXTRAPOLATED:
                        set_val = above_val + use_z_factor * (below_val - above_val)
                    case Interpolation.LINEAR_TO_ZERO:
                        set_val = 0 if is_bottom else below_val * (use_below.depth - bottom_depth) / (depth - bottom_depth)
                    case _:
                        set_val = below_val
            elif is_over:
                match ambient_store_at.extrapolation_sfc:
                    case Interpolation.CONSTANT:
                        set_val = above_val
                    case Interpolation.EXTRAPOLATED:
                        set_val = above_val + use_z_factor * (below_val - above_val)
                    case Interpolation.LINEAR_TO_ZERO:
                        set_val = 0 if surfacing else above_val * (use_above.depth - depth) / use_above.depth
                    case _:
                        set_val = above_val
            elif self.model_params.current_vector_averaging and vkey == 'current_dir':
                # (oc) direction averaging
                set_val = vector_average(
                    depth,
                    below_val,
                    above_val,
                    use_below.depth,
                    use_above.depth
                )
            else:
                set_val = above_val + use_z_factor * (below_val - above_val)
            if do_conversion:
                set_val = units.convert(
                    set_val,
                    vkey,
                    from_units=do_conversion,
                    model_params=umunit.model_params,
                    celsius=ambient_cond.temperature,
                    psu=ambient_cond.salinity,
                    depth=depth
                )
            ambient_cond.set(vkey, set_val)

        fill_var('temperature')  # cause salinity conversion requires temperature to be in celsius
        fill_var('salinity')     # cause decay rate also requires salinity to be in PSU
        # do rest of conversions
        for vkey in self.ambient_store._input_vars_:
            if vkey not in ('salinity', 'temperature'):
                fill_var(vkey)
