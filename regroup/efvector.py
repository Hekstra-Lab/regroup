from pymol.cgo import *
from pymol import cmd
import gemmi
import numpy as np


def get_orthogonalization_matrix(cell):
    """
    Compute real-space orthogonalization matrix from unit cell parameters.
    """
    a, b, c, alpha, beta, gamma = cell.parameters
    alpha = np.deg2rad(alpha)
    beta = np.deg2rad(beta)
    gamma = np.deg2rad(gamma)

    # Compute unit cell volume
    V = (
        a
        * b
        * c
        * np.sqrt(
            1
            - np.cos(alpha) ** 2
            - np.cos(beta) ** 2
            - np.cos(gamma) ** 2
            + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma)
        )
    )

    # Compute Cartesian orthogonalization matrix (Rupp, Page 746)
    O = np.zeros((3, 3))
    O[0, 0] = a
    O[0, 1] = b * np.cos(gamma)
    O[1, 1] = b * np.sin(gamma)
    O[0, 2] = c * np.cos(beta)
    O[1, 2] = c * (np.cos(alpha) - (np.cos(beta) * np.cos(gamma))) / np.sin(gamma)
    O[2, 2] = V / (a * b * np.sin(gamma))

    return O


def add_efvector(
    obj,
    h,
    k,
    l,
    pos=None,
    invert_polarity=None,
    expansion=None,
    name="efvector",
    red=True,
    color=None,
):
    """
    Add EF vector based on `obj` cell parameters in direction of Miller plane [h, k, l]. 

    Example:
    ```
    pos = [85.2496566772461, 72.81723022460938, 80.97605895996094]
    add_efvector("ls", -1, 1, -2, pos, name="efvector7", expansion = 8)
    ```
    """

    # Get unit vector corresponding to EF direction
    sym = cmd.get_symmetry(obj)
    cell = gemmi.UnitCell(*np.array(sym[:6]).round(3))
    print(cell.parameters)
    O = get_orthogonalization_matrix(cell)

    if isinstance(h, int) and isinstance(k, int) and isinstance(l, int):
        millerplane = -1 * np.array([h, k, l], dtype=int)
        if invert_polarity:
            millerplane *= -1
            

        #Rupp p. 748, eqn. A-48 shows that the facet normal in the Euclidean crystal
        # frame is (np.linalg.inv(O).T@millerplane).
        v = (np.linalg.inv(O).T @ millerplane)
        v = v / np.linalg.norm(v)
        print(v)
    else:
        print("""treating hkl arguments as
        real-space fractional coordinates, 
        not reciprocal-space grid points.""")

        fractional = -1 * np.array([h,k,l])
        if invert_polarity:
            fractional *= -1
        
        #treat as row vector. see Rupp (2018) p. 233, eqn. 5-4
        v = fractional @ O.T
        v = v / np.linalg.norm(v)
        print(v)

    if color is None:
        if red:
            color = [1.0, 0.0, 0.0]
        else:
            color = [0.0, 0.0, 1.0]
    else:
        color = color

    # Position for start of arrow
    if pos is None:
        pos = cmd.get_position()

    # Set up aesthetics of arrow
    w = 0.2  # cylinder width
    l = 2  # cylinder length
    h = 0.66  # cone height
    d = w * 1.618  # cone base diameter

    if expansion:
        w *= expansion
        l *= expansion
        h *= expansion
        d *= expansion

    obj = [
        CYLINDER,
        pos[0],
        pos[1],
        pos[2],
        pos[0] + (v[0] * l),
        pos[1] + (v[1] * l),
        pos[2] + (v[2] * l),
        w,
        *color,
        *color,
        CONE,
        pos[0] + (v[0] * l),
        pos[1] + (v[1] * l),
        pos[2] + (v[2] * l),
        pos[0] + (v[0] * (l + h)),
        pos[1] + (v[1] * (l + h)),
        pos[2] + (v[2] * (l + h)),
        d,
        0.0,
        *color,
        *color,
        1.0,
        1.0,
    ]

    cmd.load_cgo(obj, name)
    return


cmd.extend("add_efvector", add_efvector)
