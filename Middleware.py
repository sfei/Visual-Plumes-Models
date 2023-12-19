import math, traceback
import numpy as np
from .globals import UserInputError
from . import units
from .helpers import num_format
from .vectors import DEGREES_TO_RADIAN, GRAVITY_VECTOR, make_vector, unit_vector, rescale_vector, magnitude
from .UMUnit import UMUnit
from .ambient.Ambient import Ambient
from .ambient.AmbientStore import AmbientStore
from .farfield.BrooksFarField import BrooksFarField, run_far_field
from .farfield.TidalPollutionBuildup import TidalPollutionBuildup, run_tidal_pollution_buildup
from .params.DiffuserParameters import DiffuserParameters
from .params.DiffuserStore import DiffuserStore
from .params.ModelParameters import ModelParameters, Model, FarfieldDiffusivity, BacteriaModel, SimilarityProfile, \
    MaxVerticalReversals
from .Output import OutputUM3
from .timeseries.TimeseriesHandler import TimeseriesHandler


def run(model_params, diffuser_params, diffuser_store, timeseries_handler, ambient_stack, ambient_store, output_handler):
    """ Middleware handler per diffuser-ambient combination. Handles multiple cases for timeseries-based cases. That is,
    for cases created as time-slices from the start to end time in the diffuser settings. Conditions can change based on
    timeseries files provided for ambient and/or diffuser parameters. However the diffuser and ambient settings
    themselves are constant. It does not handle changing combination cases of different diffuser/ambient settings.

    Currently runs UM3 plume model and optionally a far-field model. Far-field may be either Brooks far-field or tidal
    pollution buildup, but not both.

    Roughly transcribed from TMainform.wUM1Click in main.pas

    Args:
        model_params       (ModelParameters):    Model parameters.
        diffuser_params    (DiffuserParameters): Diffuser parameters.
        diffuser_store     (DiffuserStore):      Diffuser parameter meta/info.
        timeseries_handler (TimeseriesHandler):  Timeseries handler (or None) loaded with timeseries params and data. If
                                                 None, just runs a single case.
        ambient_stack      (Ambient[]):          Stack of ambient conditions in order of increasing depth.
        ambient_store      (AmbientStore):       Ambient paramter meta/info.
        output_handler     (OutputUM3):          Output handler (or None) for defining tracked parameters in UM3 model.
                                                 If None, loads default set of output parameters to track.
    Returns:
        A dictionary of outputs. See README.md for full description.
    """
    try:
        return run_unsafe(model_params, diffuser_params, diffuser_store, timeseries_handler, ambient_stack, ambient_store, output_handler)
    except AssertionError as e:
        return {
            "success": False,
            "error": e,
            "errortraceback": traceback.format_exc()
        }
    except UserInputError as e:
        return {
            "success": False,
            "error": e,
            "errortraceback": traceback.format_exc()
        }


def run_unsafe(model_params, diffuser_params, diffuser_store, timeseries_handler, ambient_stack, ambient_store, output_handler):
    assert isinstance(model_params, ModelParameters)
    assert isinstance(diffuser_params, DiffuserParameters)
    assert isinstance(diffuser_store, DiffuserStore)
    assert isinstance(diffuser_params, DiffuserParameters)
    for stack in ambient_stack:
        assert isinstance(stack, Ambient)
    assert isinstance(ambient_store, AmbientStore)
    if output_handler:
        assert isinstance(output_handler, OutputUM3)

    if model_params.model != Model.UM3:
        raise UserInputError("Currently only the UM3 plume model is supported.")
    # TODO: if other supported models, some config and warnings need to be transcribed, makecheckset, modellabel

    if model_params.brooks_far_field and model_params.tidal_pollution_buildup:
        raise UserInputError("Cannot select both Brooks Far-Field  and Tidal Pollution Buildup. Select one or neither.")

    # otherwise tool will complain about no valid far-field current values, when not needed to be supplied
    if not model_params.brooks_far_field and len(ambient_stack):
        has_ff_vel = False
        has_ff_dir = False
        has_ff_coeff = False
        for stack in ambient_stack:
            if not has_ff_vel and stack.ff_velocity is not None:
                has_ff_vel = True
            if not has_ff_dir and stack.ff_dir is not None:
                has_ff_dir = True
            if not has_ff_coeff and stack.ff_diff_coeff is not None:
                has_ff_coeff = True
            if min(has_ff_vel, has_ff_dir, has_ff_coeff):
                break
        if not min(has_ff_vel, has_ff_dir, has_ff_coeff):
            ambient_stack = list(ambient_stack)
            ambient_stack[0] = ambient_stack[0].copy()
            if not has_ff_vel:
                ambient_stack[0].ff_velocity = 0.0
            if not has_ff_dir:
                ambient_stack[0].ff_dir = 0.0
            if not has_ff_coeff:
                ambient_stack[0].ff_diff_coeff = 0.0
            ambient_stack = tuple(ambient_stack)

    # calculate intervals for time-slice cases
    if not timeseries_handler:
        current_time_secs   = 0
        end_time_secs       = 0
        increment_time_secs = 1
    else:
        assert isinstance(timeseries_handler, TimeseriesHandler)
        current_time_secs = units.Time.convert(
            timeseries_handler.start_time,
            timeseries_handler.units.start_time,
            units.Time.SECONDS
        )
        end_time_secs = units.Time.convert(
            timeseries_handler.end_time,
            timeseries_handler.units.end_time,
            units.Time.SECONDS
        )
        increment_time_secs = units.Time.convert(
            timeseries_handler.time_increment,
            timeseries_handler.units.time_increment,
            units.Time.SECONDS
        )
    if end_time_secs < current_time_secs:
        raise UserInputError("End time must be same as or after start time")

    match diffuser_store.isopleth.units:
        case units.Isopleth.CONCENTRATION:
            isopleth_utype = units.Concentration
            isopleth_units = diffuser_store.concentration.units
        case units.Isopleth.SALINITY:
            isopleth_utype = units.Salinity
            isopleth_units = diffuser_store.salinity.units
        case units.Isopleth.TEMPERATURE:
            isopleth_utype = units.Temperature
            isopleth_units = diffuser_store.temperature.units
        case units.Isopleth.SPEED:
            isopleth_utype = units.Speed
            isopleth_units = ambient_store.current_speed.units
        case _:
            raise UserInputError("Invalid isopleth units")

    # model parameters output
    memos_model_params = []
    memos_model_params.append("---------------------------------------------------")
    memos_model_params.append("Model Parameters")
    memos_model_params.append("---------------------------------------------------")
    memos_model_params.append(f"Report effective dilution:    {bool(model_params.report_effective_dillution)}")
    memos_model_params.append(f"Current vector averaging:     {bool(model_params.current_vector_averaging)}")
    memos_model_params.append(f"Allow induced current:        {bool(model_params.allow_induced_current)}")
    memos_model_params.append(f"Bacteria model:               {BacteriaModel.label(model_params.bacteria_model)}")
    memos_model_params.append(f"Equation of state:            {'p,' if not model_params.at_equilibrium else ''}S,T")
    memos_model_params.append(f"Similarity profile:           {SimilarityProfile.label(model_params.similarity_profile)}")
    memos_model_params.append("Coefficients")
    memos_model_params.append(f"   Diffuser port contraction: {model_params.contraction_coeff}")
    memos_model_params.append(f"   Light absorption:          {model_params.light_absorb_coeff}")
    memos_model_params.append(f"   UM3 aspiration:            {model_params.aspiration_coeff}")
    memos_model_params.append("Stop conditions")
    memos_model_params.append(f"   Max. reported dilution:    {model_params.max_dilution}")
    memos_model_params.append(f"   Max. vertical reversals:   {MaxVerticalReversals.label(model_params.max_reversals)}")
    memos_model_params.append(f"   Stop on bottom hit:        {bool(model_params.stop_on_bottom_hit)}")
    memos_model_params.append(f"   Don't stop on surface hit: {bool(model_params.dont_stop_on_surface_hit)}")
    if model_params.brooks_far_field:
        memos_model_params.append("Brooks far-field solution")
        ff_increment = "ALL" if model_params.output_all_ff_increments else f"{model_params.ff_increment} m"
        memos_model_params.append(f"   Far-field increment:       {ff_increment}")
        memos_model_params.append(f"   Estimate far-field bg.:    {bool(model_params.estimate_ff_background)}")
        memos_model_params.append(f"   Far-field diffusivity:     {FarfieldDiffusivity.label(model_params.farfield_diffusivity)}")
    elif model_params.tidal_pollution_buildup:
        memos_model_params.append("Tidal pollution buildup")
        memos_model_params.append(f"   Channel width:             {model_params.tpb_channel_width} m")
        memos_model_params.append(f"   Channel seg. length:       {model_params.tpb_segment_length} m")
        memos_model_params.append(f"   Upstream dir.:             {model_params.tpb_upstream_dir}°")
        memos_model_params.append(f"   Coast bin:                 {model_params.tpb_coast_bin}")
        memos_model_params.append(f"   Coast concentration:       {model_params.tpb_coast_concentration}")
        memos_model_params.append(f"   Mixing zone depth:         {model_params.tpb_mixing_zone_ceil} m")
    else:
        memos_model_params.append("No far-field model selection")
    if model_params.use_shore_vector:
        memos_model_params.append(f"Shore vector {model_params.dist_to_shore} m at {model_params.dir_to_shore}°")
    memos_model_params.append(f"Elevation projection plane at {model_params.elevation_proj_plane}°")
    memos_model_params.append(f"Output each {model_params.write_step_freq} steps")

    if timeseries_handler:
        memos_model_params.append("")
        memos_model_params.append("---------------------------------------------------")
        memos_model_params.append("Timeseries Data")
        memos_model_params.append("---------------------------------------------------")
        memos_model_params.append(f"Start/end time: {timeseries_handler.start_time} {units.Time.label(timeseries_handler.units.start_time)} / {timeseries_handler.end_time} {units.Time.label(timeseries_handler.units.end_time)}")
        memos_model_params.append(f"Time increment: {timeseries_handler.time_increment} {units.Time.label(timeseries_handler.units.time_increment)}")
        i = 0
        diffuser_ts_info = timeseries_handler.get_diffuser_info()
        if len(diffuser_ts_info):
            for vkey, info in diffuser_ts_info.items():
                i += 1
                in_units = units.from_var_name(vkey)
                if in_units:
                    units_label = in_units.label(info['store'].units)
                else:
                    units_label = ""
                memos_model_params.append(f"{i}. Diffuser timeseries ({vkey})")
                memos_model_params.append(f"   units:          {units_label}")
                memos_model_params.append(f"   increment:      {info['store'].ts_increment} hrs")
                memos_model_params.append(f"   cycling period: {info['length']} -> {info['store'].ts_increment*info['length']} hrs")
        ambient_ts_info = timeseries_handler.get_ambient_info()
        if len(ambient_ts_info):
            for vkey, info in ambient_ts_info.items():
                i += 1
                in_units = units.from_var_name(vkey)
                if in_units:
                    units_label = in_units.label(info['store'].units)
                else:
                    units_label = ""
                depth_units = units.Length.label(info['store'].ts_depth_units)
                memos_model_params.append(f"{i}. Ambient timeseries ({vkey})")
                memos_model_params.append(f"   units:          {units_label}")
                memos_model_params.append(f"   depth(s):       {', '.join([str(d) for d in info['depths']])} {depth_units}")
                memos_model_params.append(f"   increment:      {info['store'].ts_increment} hrs")
                memos_model_params.append(f"   cycling period: {info['length']} -> {info['store'].ts_increment*info['length']} hrs")

    # for diffuser table output
    diff_out_format = [
        ("Port diameter",        'diameter',          units.Length),
        ("Vertical angle",       'vertical_angle',    units.Angle),
        ("Horizontal angle",     'horizontal_angle',  units.Angle),
        ("X-offset",             'offset_x',          units.Length),
        ("Y-offset",             'offset_y',          units.Length),
        ("Num. ports",           'num_ports',         units.Unitless)
    ]
    if diffuser_params.num_ports > 1:
        diff_out_format.append(
            ("Port spacing", 'port_spacing', units.Length)
        )
    diff_out_format += [
        ("Mixing zone distance", 'acute_mixing_zone', units.Length),
        ("Isopleth",             'isopleth',          isopleth_utype),
        ("Port depth",           'depth',             units.Length),
        ("Total flow",           'effluent_flow',     units.FlowRate),
        ("Salinity",             'salinity',          units.Salinity),
        ("Temperature",          'temperature',       units.Temperature),
        ("Pollutant",            'concentration',     units.Concentration)
    ]
    diff_out_units = tuple(
        (diffuser_store.get(d[1]).units if d[1] != 'isopleth' else isopleth_units)
        for d in diff_out_format
    )
    headers_diff = []
    outputs_diff = []
    for i, d in enumerate(diff_out_format):
        headers_diff.append({
            'name':        d[1],
            'label':       d[0],
            'units':       d[2],
            'in_units':    diff_out_units[i],
            'units_label': "" if d[2] == units.Unitless else d[2].label(diff_out_units[i])
        })

    # for ambient table output
    amb_out_format = (
        ("Depth",                      'depth',         units.Length,           True),
        ("Current speed",              'current_speed', units.Speed,            True),
        ("Current direction",          'current_dir',   units.Angle,            True),
        ("Salinity",                   'salinity',      units.Salinity,         True),
        ("Temperature",                'temperature',   units.Temperature,      True),
        ("Background concentration",   'bg_conc',       units.Concentration,    True),
        ("Decay rate",                 'decay_rate',    units.DecayRate,        True),
        ("Far-field speed",            'ff_velocity',   units.Speed,            model_params.brooks_far_field),
        ("Far-field direction",        'ff_dir',        units.Angle,            model_params.brooks_far_field),
        ("Dispersion",                 'ff_diff_coeff', units.EddyDiffusivity,  model_params.brooks_far_field),
        ("Density",                    'density',       units.Density,          True)
    )
    amb_out_format = tuple(a for a in amb_out_format if a[-1])
    amb_out_units = tuple(
        (ambient_store.z.units if a[1] == "depth" else (
            units.Density.SIGMA_T if a[1] == "density" else ambient_store.get(a[1]).units
        ))
        for a in amb_out_format
    )
    headers_amb = []
    outputs_amb = []
    for i, a in enumerate(amb_out_format):
        if a[2] == units.Unitless:
            units_label = ""
        else:
            units_label = a[2].label(amb_out_units[i])
        headers_amb.append({
            'name':        a[1],
            'label':       a[0],
            'units':       a[2],
            'in_units':    amb_out_units[i],
            'units_label': units_label
        })

    # timeseries indices by case
    memos_ts_indices = []
    # plume model outputs
    memos_um      = []
    memos_um_post = []
    headers_um    = tuple(hdr for hdr in output_handler.headers()) if output_handler else None
    outputs_um    = []
    # brooks ff model outputs
    ff_was_run    = False
    memos_ff      = []
    headers_ff    = []
    outputs_ff    = []
    # tidal pollution buildup model outputs
    tpb_was_run   = False
    memos_tpb     = []

    graph_handler = None

    # prepare far-field model
    tpb_model = brooks_ff_model = None
    if model_params.tidal_pollution_buildup:
        tpb_model = TidalPollutionBuildup(model_params)
    elif model_params.brooks_far_field:
        if model_params.model in (Model.DKH, Model.NRFIELD):
            model_params.tpb_channel_width = 1e6
        brooks_ff_model = BrooksFarField(model_params)
        model_params.tpb_channel_width = 1e6
    else:
        model_params.tpb_channel_width = 1e6

    # TODO: other necessary processes? setvisiblecomps?

    model_params.casecount = 0

    # loop from start time to end time, each iteration being a case
    casetimes = []
    while current_time_secs <= end_time_secs:
        model_params.casecount += 1
        casetimes.append(current_time_secs)

        diffuser_ts_indices = {}
        ambient_ts_indices = {}
        if timeseries_handler:
            timeseries_handler.set_time(current_time_secs)
            # diffuser settings are single value in dictionary by var name
            diffuser_ts_values = timeseries_handler.get_diffuser()
            # diffuser timeseries can just overwrite the values in diffuser params since it's just a set of values
            for vkey, (value, ts_index) in diffuser_ts_values.items():
                diffuser_params.set(vkey, value)
                diffuser_ts_indices[vkey] = ts_index
            # ambient settings are an ambient stack of Ambient by var name, have to be handled specially because of depth layers
            ambient_ts_stacks, ambient_ts_indices = timeseries_handler.get_ambient()
        else:
            diffuser_ts_values = None
            ambient_ts_stacks = None

        for key in diffuser_store._vars_:
            diffuser_store.get(key).from_time_series = diffuser_ts_values and (key in diffuser_ts_values)
        for key in ambient_store._input_vars_:
            ambient_store.get(key).from_time_series = ambient_ts_stacks and (key in ambient_ts_stacks)

        ts_indices = list(diffuser_ts_indices.values()) + list(ambient_ts_indices.values())
        if len(ts_indices):
            memos_ts_indices.append([f"Timeseries indices: {', '.join([str(i) for i in ts_indices])}"])
        else:
            memos_ts_indices.append([])

        # output handler copy per case
        um_output = OutputUM3(output_handler) if output_handler else None
        # set up model for this time slice
        um_model = UMUnit(
            model_parameters=model_params,
            diffuser_parameters=diffuser_params,
            diffuser_store=diffuser_store,
            ambient_stack=ambient_stack,
            ambient_store=ambient_store,
            ambient_ts_stacks=ambient_ts_stacks,
            output_handler=um_output,
            graph_handler=graph_handler
        )
        # run the model
        um_model.merge()
        memos_um.append(um_output.memos)
        # if some other fatal error occured
        if um_model.status['nogo']:
            raise UserInputError("Correct deficiencies and try again")

        # post run memos
        post_memos = []
        if um_model.status['close_isoplet']:
            post_memos.append("Isopleth closed. Extend simulation by reducing the isopleth value.")
            post_memos.append("")
        # (oc) 2012 Lmz calculation
        v_c1_vg    = unit_vector(np.cross(um_model.v_align, GRAVITY_VECTOR))
        Horradius_ = rescale_vector(0.5*um_model.element.diameter, np.cross(GRAVITY_VECTOR, um_model.element.v_velocity))
        Lmz_       = np.dot(um_model.element.v_surface_tdsp - um_model.v_source, v_c1_vg) * v_c1_vg
        util       = abs(np.dot(Horradius_, v_c1_vg))
        util2      = abs(np.dot(Lmz_, v_c1_vg))
        # precise stuff removed, annoying just number-to-string formatting
        post_memos.append(
            f"Horiz plane projections in effluent direction: radius(m): {num_format(util)};\n" +
            f"CL(m):      {num_format(util2)}\n" +
            f"Lmz(m):     {num_format(util + util2)}\n" +
            f"Rate sec-1: {num_format(um_model.ambient.decay_rate)} s-1  /  {num_format(um_model.ambient.decay_rate*86400)} d-1  /  {num_format(um_model.ambient.kt)} kt"
            # f"Amb Sal:    {um_model.ambient.salinity:.3f}"
        )
        # TODO: what are these outputs? They look like just for debugging purposes, do we need to keep?
        # z_travel = um_model.diff_params.depth - um_model.element.depth
        # forced_entrain_outputs = [
        #     "{:>9.3f}".format(2 * um_model.entrainment_spco2 * z_travel * um_model.ambient.current_speed / um_model.diff_params.effluent_flow * um_model.diff_params.num_ports),
        #     "{:>9.3f}".format(z_travel),
        #     "{:>9.3f}".format(um_model.element.diameter),
        #     "{:>9.3f}".format(um_model.cth)
        # ]
        # post_memos.append(f"\nforced entrain:")
        # post_memos.append(''.join(forced_entrain_outputs))
        # mag_A = magnitude(um_model.orig_element.v_radius)
        # if mag_A == 0:
        #     post_memos.append(
        #         "                      division by zero avoided"
        #     )
        #     post_memos.append(
        #         "      999" +
        #         "{:>9.3f}".format(um_model.ambient.current_speed) +
        #         "{:>9.3f}".format(um_model.element.diameter)
        #     )
        # else:
        #     post_memos.append(
        #         "      999" +
        #         "{:>9.3f}".format(um_model.diff_params.effluent_flow / um_model.diff_params.num_ports / mag_A / mag_A / math.pi) +
        #         "{:>9.3f}".format(um_model.ambient.current_speed) +
        #         "{:>9.3f}".format(um_model.element.diameter)
        #     )
        # if um_model.dil100mx > 1.0:
        #     post_memos.append(
        #         "{:>9.3f}".format(um_model.ambient.depth) +
        #         "{:>9.3f}".format(um_model.dil100mx if um_model.dil100mx < um_model.dil100mp else um_model.dil100mp) +
        #         "  -999.99" +
        #         "{:>9.3f}".format(um_model.diff_params.acute_mixing_zone) +
        #         "{:>9.3f}".format(um_model.element.total_time) +
        #         "{:>9.3f}".format(um_model.travel)
        #     )
        # else:
        #     # arg_vector = rescale_vector(1.0, project_vector(um_model.orig_element.v_velocity, K_UNIT_VECTOR))
        #     # um_model.dis100mx = abs(np.dot(um_model.element.v_surface_tdsp, arg_vector))
        #     post_memos.append(
        #         "{:>9.3f}".format(um_model.ambient.depth) +
        #         "{:>9.3f}".format(um_model.element.dilution) +
        #         "{:>9.3f}".format(units.Concentration.convert(um_model.element.concentration, 1, um_model.diffuser_store.concentration.units)) +
        #         "{:>9.3f}".format(um_model.element.total_time) +
        #         "{:>9.3f}".format(um_model.travel)
        #     )
        memos_um_post.append(post_memos)

        # save diffuser and ambient values
        outputs_diff.append(tuple(
            d[2].convert(
                um_model.diff_params.get(d[1]),
                ufrom=1,
                uto=diff_out_units[i]
            )
            for i, d in enumerate(diff_out_format)
        ))
        outputs_amb.append(tuple(
            tuple(
                a[2].convert(
                    amb.get(a[1]),
                    ufrom=1,
                    uto=amb_out_units[i]
                )
                for i, a in enumerate(amb_out_format)
            )
            for amb in um_model.ambient_stack
        ))

        # default outputs if not supplied
        if not output_handler:
            graph_handler = um_model.graphit
        if not graph_handler:
            graph_handler = um_model.graphit
        # append the outputs from output handler
        outputs_um.append(tuple(outs for outs in um_output.outputs()))

        # run the farfield
        if tpb_model:
            tpb_was_run = True
            model_params.tpb_bincon = run_tidal_pollution_buildup(um_model, tpb_model, increment_time_secs)
        elif brooks_ff_model:
            ff_was_run = True
            ff_output = run_far_field(um_model, brooks_ff_model, timeseries_handler)
            memos_ff.append(ff_output.memos)
            if ff_output.length:
                if not headers_ff:
                    headers_ff = tuple(hdr for hdr in ff_output.headers())
                outputs_ff.append(tuple(outs for outs in ff_output.outputs()))
            else:
                outputs_ff.append([])

        if model_params.casecount == 1:
            verify_rad = model_params.elevation_proj_plane * DEGREES_TO_RADIAN
            graph_handler.graph_vector(
                um_model.v_source,
                um_model.v_source + 0.5*um_model.element.diameter*make_vector(math.cos(verify_rad), math.sin(verify_rad))
            )

        # increment time
        current_time_secs += increment_time_secs

        # TODO: other stop condition for stopbtnhit?

    # TODO: some graphing stuff for other models not handled

    if tpb_model:
        memos_tpb.append(f"Highest buildup of concentration: {num_format(tpb_model.worstback)}")
        if tpb_model.ibmax >= 200:
            memos_tpb.append("Internal array cells used >= 200 perhaps signaling undesirable loss of pollutant upstream. Increase Channel seg. length on the Settings tab")
        elif tpb_model.ibmax < 190 and model_params.casecount > 200:
            memos_tpb.append("Internal array cells used < 190 Consider decreasing Channel seg. length on the Settings tab to optimize the precison of the pollutant buildup method.")

    return {
        "success": True,
        "error": None,
        "errortraceback": None,
        "cases": model_params.casecount,
        "casetime": casetimes,
        "modelparams": {
           "memos":      memos_model_params
        },
        "timeseries": {
            "memos":     memos_ts_indices,
        },
        "diffuser": {
            "headers":   headers_diff,
            "outputs":   outputs_diff
        },
        "ambient": {
            "headers":   headers_amb,
            "outputs":   outputs_amb
        },
        "plume": {
            "headers":   headers_um,
            "outputs":   outputs_um,
            "memos":     memos_um,
            "postmemos": memos_um_post
        },
        "farfield": {
            "was_run":   ff_was_run,
            "headers":   headers_ff,
            "outputs":   outputs_ff,
            "memos":     memos_ff
        },
        "tpb": {
            "was_run":   tpb_was_run,
            "memos":     memos_tpb
        },
        "graphs": {
            series_name: {
                'coords': to_list_nan_to_none(series_coords),
                'units':  graph_handler.units[series_name]
            }
            for series_name, series_coords in graph_handler.series.items()
        }
    }


def to_list_nan_to_none(arr):
    return list(map(_list_nan_to_none, arr))


def _list_nan_to_none(val):
    if isinstance(val, np.ndarray):
        return list(map(_list_nan_to_none, val.tolist()))
    if isinstance(val, (list, tuple)):
        return list(map(_list_nan_to_none, val))
    return None if np.isnan(val) else val
