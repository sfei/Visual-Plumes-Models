import math
import numpy as np


#-----------------------------------------------------------------------------------
# VARIABLE REMAPS
#-----------------------------------------------------------------------------------
# zero_     -> ZERO_VECTOR      : vector where (i,j,k) = (0,0,0)
# i_        -> I_UNIT_VECTOR    : vector where (i,j,k) = (1,0,0)
# j_        -> J_UNIT_VECTOR    : vector where (i,j,k) = (0,1,0)
# j_        -> K_UNIT_VECTOR    : vector where (i,j,k) = (0,0,1)
# G         -> GRAVITY_VECTOR   : vector where (i,j,k) = (0,0,9.807)

#-----------------------------------------------------------------------------------
# FUNCTION REMAPS
#-----------------------------------------------------------------------------------
# angle(A, B)       -> angle(A, B)
# _b(i, j, k)       -> make_vector(x, y, z)
# _c(m, A)          -> rescale_vector(new_magnitude, A)
# cosine(A, B)      -> vector_cosine(A, B)
# dot(A, B)         -> NOT RECREATED
#                      this is just the dot product two vectors
#                      instead, with numpy, just do: np.dot(A, B)
# mag(A)            -> magnitude(A)
# _pp(A, B)         -> project_vector(A, onto_B)
# _s(A, B)          -> NOT RECREATED
#                      this is just the vector sum
#                      instead, with numpy, just do: A + B
# _d(A, B)          -> NOT RECREATED
#                      this is just the vector difference
#                      instead, with numpy, just do: A - B
# _sx(f, A)         -> NOT RECREATED
#                      this is just the product of a scalar (f) and a vector (A)
#                      instead, with numpy, just do: f*A
# _x(A, B)          -> NOT RECREATED
#                      this is just the cross product two vectors
#                      instead, with numpy, just do: np.cross(A, B)
# xcomp(A)          -> xcomp(A)
#                      left as is, though can just as well do A[0]
# ycomp(A)          -> ycomp(A)
#                      left as is, though can just as well do A[1]
# zcomp(A)          -> zcomp(A)
#                      left as is, though can just as well do A[2]
# vectoraverage(z, dird, diru, zd, zu)
#                   -> vector_average(z, dird, diru, zd, zu)
# convertvectors(c) -> convert_vectors(params, element)


DEGREES_TO_RADIAN = math.pi/180.0
RADIANS_TO_DEGREE = 1.0/DEGREES_TO_RADIAN


def make_vector(x=0.0, y=0.0, z=0.0):
    # formerly _b()
    # create a vector up to three dimensions
    return np.array([x, y, z])


# formerly zero_
ZERO_VECTOR = make_vector()

# formerly i_, j_, and k_
I_UNIT_VECTOR = make_vector(x=1.0)
J_UNIT_VECTOR = make_vector(y=1.0)
K_UNIT_VECTOR = make_vector(z=1.0)

# formerly G
MAGNITUDE_GRAVITY = 9.807
GRAVITY_VECTOR = make_vector(z=-9.807)


def change_vector(A, x=None, y=None, z=None):
    # change a vector's component(s) to given values
    B = np.copy(A)
    for i, new_value in enumerate((x,y,z)):
        if new_value is not None:
            B[i] = new_value
    return B


def unit_vector(A):
    # convert a vector to a unit vector, that is, rescale components so magnitude=1
    # basically rescale_vector(1, A)
    if is_zero_vector(A):
        return ZERO_VECTOR
    return A / magnitude(A)


def rescale_vector(new_magnitude, A):
    # formerly _c()
    # rescale a vector to a given magnitude, while maintaining angle/proportions
    if is_zero_vector(A):
        return ZERO_VECTOR
    return A*(new_magnitude / magnitude(A))


def magnitude(A):
    # formerly mag()
    # get magnitude of vector
    # equivalent to np.sum(np.square(A))**0.5, but using shorter numpy function
    return np.linalg.norm(A) if not is_zero_vector(A) else 0


def is_zero_vector(A):
    # quick check if zero vector
    return A[0] == A[1] == A[2] == 0


# these functions could be called manually, but will leave in just in case
def xcomp(A):
    return A[0]
def ycomp(A):
    return A[1]
def zcomp(A):
    return A[2]


def angle(A, B):
    # angle between two vectors
    if is_zero_vector(A) or is_zero_vector(B):
        # TODO: technically this is undefined, but (o.c.) returns this to match validation results
        return 90*DEGREES_TO_RADIAN
    # some fanciness to handle vectors on the same line
    return np.arccos(
        np.clip(
            np.dot(unit_vector(A), unit_vector(B)), 
            -1.0, 
            1.0
        )
    )


def vector_cosine(A, B):
    # cosine of angle between two vectors
    # expanded name to be more clear and separate from simple cos() function
    if is_zero_vector(A) or is_zero_vector(B):
        return 0
    return np.dot(A, B) / (magnitude(A) * magnitude(B))


def project_vector(A, onto_B):
    # TODO, this is a little different than typical projection formula -- double check why
    # Traditionally it's dot product A*B divided by magnitude of B
    B = unit_vector(onto_B)
    return A - rescale_vector(np.dot(A, B), B)


def vector_average(z, dird, diru, zd, zu):
    dird_radians = DEGREES_TO_RADIAN*dird
    diru_radians = DEGREES_TO_RADIAN*diru
    vdird = make_vector(math.cos(dird_radians), math.sin(dird_radians), 0)
    vdiru = make_vector(math.cos(diru_radians), math.sin(diru_radians), 0)
    vdirm = np.cross(vdird, vdiru)
    avg   = RADIANS_TO_DEGREE * math.acos(vector_cosine(vdird, vdiru)) * (zd - z)/(zd - zu)
    avg   = dird + (avg * (-1 if vdirm[0] else 1))
    if avg < -180:
        return avg + 360
    if avg > 360:
        return avg - 360
    return avg

