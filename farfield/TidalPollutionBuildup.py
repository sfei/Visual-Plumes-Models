import math
import numpy as np
from .. import helpers, UserInputError
from ..vectors import DEGREES_TO_RADIAN
from ..params.ModelParameters import ModelParameters


NUMBER_BINS  = 200
BLOCK_LENGTH = 10


def run_tidal_pollution_buildup(umunit, tpb_model, increment_time_secs):
    if not tpb_model.bin_volume:
        tpb_model.bottom_depth = umunit.model_params.bottom_depth
        tpb_model.bin_volume = tpb_model.segment_length*tpb_model.channel_width*tpb_model.bottom_depth
    message = tpb_model.run(
        casecount=umunit.model_params.casecount,
        effluent_flow=umunit.diff_params.effluent_flow,
        time_increment=increment_time_secs,
        orig_mass=umunit.orig_element.mass,
        orig_mass_pollutant=umunit.ini_mass_pollutant,
        ambient=umunit.ambient
    )
    if message:
        umunit.outputit.memo(message)
    return tpb_model.bincon


class TidalPollutionBuildup:

    def __init__(self, model_params):
        assert isinstance(model_params, ModelParameters)
        self.channel_width       = model_params.tpb_channel_width        # (formerly binwid)
        self.segment_length      = model_params.tpb_segment_length       # (formerly binlen)
        self.upstream_dir        = model_params.tpb_upstream_dir         # (formerly bindir)
        self.coast_bin           = model_params.tpb_coast_bin            # (formerly icst)
        self.coast_concentration = model_params.tpb_coast_concentration  # (formerly ccst)
        self.mixing_zone_ceil    = model_params.tpb_mixing_zone_ceil     # (formerly cig)
        self.bottom_depth        = model_params.bottom_depth             # (formerly vbot)
        self.bin_volume          = None  # bottom depth and volume calculated after plume model sets bottom depth

        if not (self.channel_width > 0):
            raise UserInputError("Invalid TPB channel width (must be >0).")
        if not (self.segment_length > 0):
            raise UserInputError("Invalid TPB segment length (must be >0).")
        if not (10 <= self.coast_bin <= 99):
            raise UserInputError("Invalid TPB coast bins (must be between 10 and 99, inclusive).")
        if not (self.coast_concentration >= 0):
            raise UserInputError("Invalid TPB segment length (must be ≥0).")
        if not (self.mixing_zone_ceil >= 0):
            raise UserInputError("Invalid TPB mixing zone depth (must be ≥0).")

        self.bin        = np.zeros((NUMBER_BINS,))
        self.binvol     = np.zeros((NUMBER_BINS,))
        self.bincon     = 0
        self.ibin       = NUMBER_BINS//2
        self.ibin1      = self.ibin
        self.binx       = self.ibin*self.segment_length
        self.binx0      = 0
        self.wasneg     = False
        self.worstback  = 0
        self.ibmax      = 0
        self.seg_len_10 = BLOCK_LENGTH*self.segment_length

    def run(self, casecount, effluent_flow, time_increment, orig_mass, orig_mass_pollutant, ambient):
        while self.ibin1 >= NUMBER_BINS:
            # shift left by block length, filling in with zeros on right
            self.bin    = helpers.shift(self.bin, -BLOCK_LENGTH, 0)
            self.ibin1 -= BLOCK_LENGTH
            self.ibin  -= BLOCK_LENGTH
            self.binx  -= self.seg_len_10
            self.binx0 -= self.seg_len_10

        addvol = effluent_flow*time_increment
        addition = addvol*orig_mass_pollutant/orig_mass
        from_bin, to_bin = sorted((self.ibin, self.ibin1))
        if from_bin != to_bin:
            to_bin   += 1
            divide    = to_bin - from_bin
            addvol   /= divide
            addition /= divide
        else:
            to_bin += 1
        self.bin[from_bin:to_bin] += addition
        self.binvol[from_bin:to_bin] += addvol

        self.binx0 = self.binx
        self.binx -= (
            ambient.current_speed
            *math.cos(DEGREES_TO_RADIAN*(ambient.current_dir - self.upstream_dir))
            *time_increment
        )
        self.ibin  = self.ibin1
        self.ibin1 = int(self.binx//self.segment_length)
        if self.ibin1 < 0:
            # if nseries == 0:
            #     showmessage("Are time-series files linked?")
            return "Outside of upstream storage limits, increase channel segment length"

        if self.ibin1 > 100 and self.ibin1 > self.ibin:
            end_i = self.ibin1 - (100 - self.coast_bin) + 1
            self.bin[:end_i] = self.coast_concentration*self.bin_volume
            self.binvol[:end_i] = 0

        # decay
        if ambient.decay_rate and time_increment:
            self.bin -= self.bin*ambient.decay_rate*time_increment
        # TODO: used anywhere?
        # binsum = np.sum(self.bin)

        if casecount == 1:
            self.wasneg = self.binx > self.binx0
        bi = self.ibin+1 if self.wasneg else self.ibin-1
        self.bincon = self.bin[bi]/(self.bin_volume + self.binvol[bi])

        if self.worstback < self.bincon:
            self.worstback = self.bincon

        try:
            # get index of last, positive, nonzero value (or first in reversed order)
            ibmax = next(len(self.bin)-i for i, v in enumerate(reversed(self.bin)) if v > 0)
            if ibmax > self.ibmax:
                self.ibmax = ibmax
        except StopIteration:
            pass
