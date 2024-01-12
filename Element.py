import math
import numpy as np
from .vectors import ZERO_VECTOR, GRAVITY_VECTOR, angle, magnitude, rescale_vector, make_vector


class Element:
    """ Represents the control volume element that is simulating the bouyant jet from the diffuser. """

    def __init__(self):
        self.v_displace       = ZERO_VECTOR        # (formerly dR)  movement per timestep (displacement)
        self.v_surface_tdsp   = ZERO_VECTOR        # (formerly R_)  total displacement, but from x/y of source on water surface
        self.v_velocity       = ZERO_VECTOR        # (formerly V)   velocity vector
        self.v_radius         = ZERO_VECTOR        # (formerly A)   radius vector aligned to plane of control volume face
        self.Cc               = ZERO_VECTOR        #                Rc rescaled to 0.05% of diameter ?? maybe shear entrainment related TODO: ever used even?
        self.Rc               = ZERO_VECTOR        #                Something aligned to plane of control volume face
        self.v_radii          = np.empty((37, 3))  # (formerly B)   rotated versions of the radius vector on plane aligned to plane of control volume face
        self.height           = 0.0                # (formerly h)   height
        self.mass             = 0.0                # (formerly m)   mass
        self.d_mass           = 0.0                # (formerly dm)  change in mass
        self.density          = 1000.0             # (formerly rho / c.vden) density

        self.diameter         = 0.0  # (formerly c.vdia) diameter
        self.depth            = 0.0  # (formerly c.vdep) port depth
        self.vertical_angle   = 0.0  # (formerly c.vang) vertical angle (in degrees)
        self.horizontal_angle = 0.0  # (formerly c.vhor) horizontal angle
        self.salinity         = 0.0  # (formerly c.vsal) salinity
        self.temperature      = 0.0  # (formerly c.vtem) temperature
        self.concentration    = 0.0  # (formerly c.vpol) pollutant concentration
        # self.v4o3             = 0.0  # (formerly c.v4o3) something related to pollutant
        self.speed            = 0.0  # (formerly c.vvel) magnitude of velocity
        self.dilution         = 0.0  # (formerly c.vdil) dilution
        self.cl_dilution      = 0.0  # (formerly c.vcld) CL-dilution
        self.total_surf_dsp   = 0.0  # (formerly Rf -- in far-field model) total surface displacement distance
        self.x_displacement   = 0.0  # (formerly c.vx) total displacement from origin in x-axis
        self.y_displacement   = 0.0  # (formerly c.vy) total displacement from origin in y-axis
        self.buildup          = 0.0  # (formerly c.vbup) tidal pollution buildup
        self.total_time       = 0.0  # (formerly c.vtime) total simulated time passed in model run

        self._flag_init_      = False

    def copy(self):
        copy = Element()
        copy.v_displace       = np.copy(self.v_displace)
        copy.v_surface_tdsp   = np.copy(self.v_surface_tdsp)
        copy.v_velocity       = np.copy(self.v_velocity)
        copy.v_radius         = np.copy(self.v_radius)
        copy.Cc               = np.copy(self.Cc)
        copy.Rc               = np.copy(self.Rc)
        copy.v_radii          = np.copy(self.v_radii)
        copy.height           = self.height
        copy.mass             = self.mass
        copy.d_mass           = self.d_mass
        copy.density          = self.density
        copy.diameter         = self.diameter
        copy.depth            = self.depth
        copy.vertical_angle   = self.vertical_angle
        copy.horizontal_angle = self.horizontal_angle
        copy.salinity         = self.salinity
        copy.temperature      = self.temperature
        copy.concentration    = self.concentration
        # copy.v4o3           = self.v4o3
        copy.speed            = self.speed
        copy.dilution         = self.dilution
        copy.cl_dilution      = self.cl_dilution
        copy.total_surf_dsp   = self.total_surf_dsp
        copy.x_displacement   = self.x_displacement
        copy.y_displacement   = self.y_displacement
        copy.buildup          = self.buildup
        copy.total_time       = self.total_time
        return copy

    def get(self, key):
        return getattr(self, key)

    def set(self, key, value):
        return setattr(self, key, value)

    def body_calc(self, last_element, dt, um3iso=1, in_init=False):
        """ Collection of body calcs.
        Args:
            last_element: The previous element conditions.
            dt: The timestep.
            um3iso: I don't know exactly, comes from isoplet calcs. Defaults to 1.
            in_init: True if from initialization.
        """
        radius  = self.diameter*0.5
        d_angle = angle(self.v_velocity, last_element.v_velocity)
        # TODO: not sure what Rc represents exactly, but all these cross products produce a vector aligned to the faces
        # of the control volume
        if in_init:
            if d_angle == 0:
                # cross product of velocity against gravity, so a vector aligned to horizontal plane results
                # why rescale twice? why not just once? why to random large value?
                self.Rc = rescale_vector(
                    1e15,
                    rescale_vector(
                        self.diameter,
                        np.cross(self.v_velocity, GRAVITY_VECTOR)
                    )
                )
            else:
                self.Rc = rescale_vector(
                    magnitude(dt*self.v_velocity) / d_angle,  # no v_displace on init, so next expected displacement
                    np.cross(
                        self.v_velocity,
                        np.cross(self.v_velocity, last_element.v_velocity)
                    )
                )
            self.v_radii[0] = rescale_vector(radius, self.Rc)
            self.v_radius = rescale_vector(radius, np.cross(self.Rc, self.v_velocity))

        else:
            if d_angle != 0:
                # cross product twice just rotates it round the plane -- if it matters, the vector points in directional
                # change of vector from last element to element
                self.Rc = rescale_vector(
                    magnitude(self.v_displace) / d_angle,
                    np.cross(
                        self.v_velocity,
                        np.cross(self.v_velocity, last_element.v_velocity)
                    )
                )
            else:
                self.Rc = rescale_vector(
                    1e15,
                    rescale_vector(self.diameter, np.cross(self.v_velocity, self.v_radius))
                )
            # from part of bdycalcs() but pieced out here
            if angle(self.v_velocity, GRAVITY_VECTOR) in (0, math.pi):
                self.v_radius = make_vector(x=radius)
            else:
                self.v_radius = rescale_vector(radius, np.cross(GRAVITY_VECTOR, self.v_velocity))
            self.v_radii[0] = rescale_vector(radius, np.cross(self.v_radius, self.v_velocity))

        N_ = self.v_surface_tdsp + um3iso*self.v_radii[0]
        Nopp_ = self.v_surface_tdsp - um3iso*self.v_radii[0]
        return N_, Nopp_
