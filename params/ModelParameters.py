from enum import Enum
from ..globals import UserInputError
from .AbstractParameters import AbstractParameters


# MODEL enumerated type (sorta)
class Model(Enum):
    UM3, DOS3, PDS, DKH, NRFIELD = range(1, 6)
    @staticmethod
    def parse(input_str):
        input_str = input_str.trim().upper() if input_str else ""
        match input_str:
            case "UM3":
                return Model.UM3
            case "DOS3":
                # TODO: not supported
                return Model.DOS3
            case "PDS":
                # TODO: not supported
                return Model.PDS
            case "DKH":
                # TODO: not supported
                return Model.DKH
            case "NRFIELD":
                # TODO: not yet transcribed
                return Model.NRFIELD
            case _:
                return None


# BACTERIA MODEL enumerated type (sorta)
class BacteriaModel(Enum):
    COLIFORM_MANCINI = "coliform_mancini"
    COLIFORM_301H    = "coliform_301h"
    ENTEROCCOUS_301H = "enterococcus_301h"
    @staticmethod
    def label(input):
        match input:
            case BacteriaModel.COLIFORM_MANCINI:
                return "Mancini (1978) coliform model"
            case BacteriaModel.COLIFORM_301H:
                return "301(h) TSD (1994) coliform (for saltwater, Eqn B-68)"
            case BacteriaModel.ENTEROCCOUS_301H:
                return "301(h) TSD (1994) enterococcus (for saltwater, Eqn B-69)"
            case _:
                return ""


class MaxVerticalReversals(Enum):
    INITIAL_TRAP_LEVEL      = 0
    MAX_RISE_OR_FALL        = 1
    SECOND_TRAP_LEVEL       = 2
    SECOND_MAX_RISE_OR_FALL = 3
    @staticmethod
    def label(input):
        match input:
            case MaxVerticalReversals.INITIAL_TRAP_LEVEL:
                return "to initial trap level"
            case MaxVerticalReversals.MAX_RISE_OR_FALL:
                return "to max rise or fall"
            case MaxVerticalReversals.SECOND_TRAP_LEVEL:
                return "to 2nd trap level"
            case MaxVerticalReversals.SECOND_MAX_RISE_OR_FALL:
                return "to 2nd max rise or fall"
            case _:
                return ""


class SimilarityProfile(Enum):
    DEFAULT   = 0
    POWER_3_2 = 1
    GAUSSIAN  = 2
    @staticmethod
    def label(input):
        match input:
            case SimilarityProfile.DEFAULT:
                return "default profile"
            case SimilarityProfile.POWER_3_2:
                return "3/2 power profile"
            case SimilarityProfile.GAUSSIAN:
                return "gaussian profile"
            case _:
                return ""


class FarfieldDiffusivity(Enum):
    CONSTANT  = 0
    POWER_4_3 = 1
    @staticmethod
    def label(input):
        match input:
            case FarfieldDiffusivity.CONSTANT:
                return "constant diffusivity"
            case FarfieldDiffusivity.POWER_4_3:
                return "4/3 power diffusivity"
            case _:
                return ""


class ModelParameters(AbstractParameters):

    def __init__(self):
        super().__init__()

        # model type
        self.model                      = Model.UM3
        self.casecount                  = 1

        ################################################################################################################
        # Model configuration
        ################################################################################################################
        # current vector averaging for ambient interpolation >> configstr.checked[3]
        self.current_vector_averaging   = False
        # allow induced current (multiport) >> >> UM3optionscb.checked[4]
        self.allow_induced_current      = False
        # report effective dillution >> configstr.checked[1]
        self.report_effective_dillution = False
        # was isopleth boudary >> configstr.checked[2] (TODO: unused)
        self.was_isopleth_boundary      = False

        ################################################################################################################
        # Tidal pollution buildup parameters
        ################################################################################################################
        # whether to calc tidal pollution buildup >> configstr.checked[4]
        self.tidal_pollution_buildup    = False
        # (formerly binwid) >> mainform.binwided -- required, assumed meters
        self.tpb_channel_width          = 1e6
        # (formerly binlen) >> mainform.lenbinded -- assumed meters
        self.tpb_segment_length         = 50
        # (formerly bindir) >> mainform.upstrmed -- assumed degrees
        self.tpb_upstream_dir           = 90
        # (formerly icst) >> mainform.coasted
        self.tpb_coast_bin              = 0
        # (formerly ccst) >> mainform.coastced
        self.tpb_coast_concentration    = 0
        # (formerly cig) >> mainform.mixciged (in UI this is labeled as mixing zone depth)
        self.tpb_mixing_zone_ceil       = 0
        # (formerly bincon) >> from tidal pollution buildup parameters
        self.tpb_bincon                 = 0  # calculated from inputs

        ################################################################################################################
        # Brooks far-field solution parameters
        ################################################################################################################
        # whether to run far field model >> configstr.checked[0]
        self.brooks_far_field           = False
        # far field increment >> farfieldinced
        self.ff_increment               = 100
        # estimate farfield background >> UM3optionscb.checked[0]
        self.estimate_ff_background     = False
        # output all farfield time increments >> UM3optionscb.checked[1]
        self.output_all_ff_increments   = False
        # farfield diffusivity option >> from special settings on middle-right as button (in oc checked as v4o3 in exset/outputs)
        # output option enable/disabled '4/3Eddy' option based on this
        self.farfield_diffusivity       = FarfieldDiffusivity.CONSTANT

        ################################################################################################################
        # Additional model inputs
        ################################################################################################################
        # diffuser port contraction coefficient
        self.contraction_coeff          = 1.0
        # (formerly ke) light absorption coefficient >> mainform.absorbed
        self.light_absorb_coeff         = 0.16
        # (formerly asp) UM3 aspiration coefficient >> mainform.aspirationed
        self.aspiration_coeff           = 0.1
        # (formerly bacteriaarg) bacteria model >> mainform.bacteriarg
        self.bacteria_model             = BacteriaModel.COLIFORM_MANCINI
        # equation of state >> mainform.EqnStateRg (option "S, T" is true, options "p, S, T" is false)
        self.at_equilibrium             = True
        # similarity profile >> mainform.simprofilerg (enum 0-2 inclusive)
        self.similarity_profile         = SimilarityProfile.DEFAULT
        # TODO: dynamic DKH spacing (needed?)
        self.dynamic_dkh_spacing        = False

        ################################################################################################################
        # Stop condition settings
        ################################################################################################################
        # maximum number of reversals before stop >> mainform.maxreversalsrg (enum 0-3 inclusive)
        # this must get converted to an integer before model run based on enum type
        self.max_reversals              = MaxVerticalReversals.SECOND_MAX_RISE_OR_FALL
        self.max_reversals_limit        = 5
        self.max_traps_limit            = 5
        # stop on bottom hit >> UM3optionscb.checked[2]
        self.stop_on_bottom_hit         = False
        # do not stop on surface hit >> UM3optionscb.checked[3]
        self.dont_stop_on_surface_hit   = False
        # maximum dilution reported >> mainform.steped
        self.max_dilution               = 10000

        ################################################################################################################
        # Output settings
        ################################################################################################################
        # (formerly ifrq) write each # of steps >> mainform.frqed
        self.write_step_freq            = 100
        # elevation projection plane >> from graphic settings -- assumed degrees
        self.elevation_proj_plane       = 0
        # Shore vector >> from graphic settings -- mainform.ShoreDisEd.text and mainform.ShoreDirEd.text
        self.use_shore_vector           = False    # the checkbox
        self.dist_to_shore              = 7250     # first input -- assumed meters
        self.dir_to_shore               = 58       # second input -- assumed degrees

        ################################################################################################################
        # calculated vars
        ################################################################################################################
        # (formerly c.vbot) bottom depth
        self.bottom_depth = 1e9

        self._vars_ = list(vars(self).keys())

    def copy(self):
        c = ModelParameters()
        for key in self._vars_:
            setattr(c, key, getattr(self, key))
        return c

    def validate(self):
        self._validate_enum(self.model, "model type", Model)
        self._validate_enum(self.bacteria_model, "bacteria model", BacteriaModel)
        self._validate_enum(self.similarity_profile, "similarity profile", SimilarityProfile)
        self.contraction_coeff  = self._validate(self.contraction_coeff, "contraction coefficient")
        self.light_absorb_coeff = self._validate(self.light_absorb_coeff, "light absorption coefficient")
        self.aspiration_coeff   = self._validate(self.aspiration_coeff, "aspiration coefficient")
        if self.brooks_far_field:
            self.ff_increment   = self._validate(self.ff_increment, "far-field increment", allow_zero=False, allow_negative=False)
            self._validate_enum(self.farfield_diffusivity, "far-field diffusivity", FarfieldDiffusivity)
        if self.tidal_pollution_buildup:
            if self.brooks_far_field:
                raise UserInputError("Cannot select both Brooks Far Field and Tidal Pollution Buildup models")
            self.tpb_channel_width       = self._validate(self.tpb_channel_width, "channel width (TPB)", allow_zero=False, allow_negative=False)
            self.tpb_segment_length      = self._validate(self.tpb_segment_length, "segment length (TPB)", allow_zero=False, allow_negative=False)
            self.tpb_upstream_dir        = self._validate(self.tpb_upstream_dir, "upstream direction (TPB)")
            self.tpb_coast_bin           = self._validate(self.tpb_coast_bin, "coast bin (TPB)", allow_zero=False, allow_negative=False, as_integer=True, min=10, max=99)
            self.tpb_coast_concentration = self._validate(self.tpb_coast_concentration, "coast concentration (TPB)", allow_negative=False)
            self.tpb_mixing_zone_ceil    = self._validate(self.tpb_mixing_zone_ceil, "mixing zone depth (TPB)", allow_negative=False)
        self._validate_enum(self.max_reversals, "max vertical reversals", MaxVerticalReversals)
        self.max_dilution                = self._validate(self.max_dilution, "maximum dilution", allow_zero=False, allow_negative=False)
        self.write_step_freq             = self._validate(self.write_step_freq, "write step frequency", allow_zero=False, allow_negative=False)
        if self.use_shore_vector:
            self.dist_to_shore           = self._validate(self.dist_to_shore, "distance to shore", allow_zero=False, allow_negative=False)
            self.dir_to_shore            = self._validate(self.dir_to_shore, "direction to shore")
