import gemmi
import numpy as np

##COORDINATE TRANSFORMATIONS


def get_normal_vector(hkl, Astar):
    """
    Get lab-frame normal vector to real-space Miller plane hkl.

    Note
    ----
        The lab-frame normal vector to a real-space Miller plane is collinear with
        the reciprocal dHKL vector. As such, we can use this simpler 
        formula in the reciprocal lattice basis to get the correct 
        orientation in the lab frame. For a graphical explanation of this,
        please look at Rupp, p238.
    """
    return hkl @ Astar.T


def lab_vec_to_crystal(v_lab, Astar):
    """
    Convert a lab-frame vector into crystal fractional coordinates.
    """
    v_lab = np.array(v_lab, dtype=float)
    v_cryst = v_lab @ Astar
    return v_cryst / np.linalg.norm(v_cryst)


def facet_normal_to_crystal_frame(hkl, O):
    """
    Convert a reciprocal-lattice facet normal hkl into a direct-space (fractional coordinates)
    vector parallel to the real-space plane normal.
    """
    hkl = np.array(hkl, dtype=float)
    v = hkl @ np.linalg.inv(O.T @ O)
    return v / np.linalg.norm(v)


## HELPERS 

def fmt_vec(v, ndigits=4):
    """Pretty-print a vector as a rounded tuple."""
    v = np.array(v, dtype=float)
    return tuple(np.round(v, ndigits))

## INNER PRODUCT COMPUTATION IN EUCLIDEAN COORDINATES AND 
## FRACTIONAL COORDINATES

def angle(v1, v2):
    """Compute angle between two vectors in Euclidean space"""
    return np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))



def dot_metric(u, v, G):
    """
    Physical dot product between fractional-coordinate vectors u and v.
    """
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    return float(u @ G @ v)


def norm_metric(v, G):
    """
    Physical norm of a fractional-coordinate vector.
    """
    vv = dot_metric(v, v, G)
    if vv < 0 and abs(vv) < 1e-12:
        vv = 0.0
    if vv < 0:
        raise ValueError(f"Metric norm squared is negative: {vv}")
    return np.sqrt(vv)


def unit_metric(v, G):
    """
    Normalize a fractional-coordinate vector using the physical crystal metric.
    """
    v = np.asarray(v, dtype=float)
    n = norm_metric(v, G)
    if n == 0:
        raise ValueError("Vector has zero physical norm.")
    return v / n


def cosine_metric(u, v, G):
    """
    Cosine of the physical angle between two fractional-coordinate vectors.
    """
    cosang = dot_metric(u, v, G) / (norm_metric(u, G) * norm_metric(v, G))
    return float(np.clip(cosang, -1.0, 1.0))


def angle_metric(u, v, G):
    """
    Physical angle in degrees between two fractional-coordinate vectors.
    """
    return float(np.degrees(np.arccos(cosine_metric(u, v, G))))

def gemmi_to_rot(op, den=None):
    """
    Convert Gemmi symop rotation to a 3x3 float matrix.

    This intentionally matches fsa_broken_symops.py.
    """
    if den is None:
        den = op.DEN if hasattr(op, "DEN") else 24

    return np.array(op.rot, dtype=float).reshape(3, 3) / den
    

def symop_string(op):
    """
    Pretty symop string, with a fallback for older cctbx builds.
    """
    try:
        return op.as_xyz()
    except Exception:
        return str(op)


def symop_key(op):
    """
    Rotation/translation key for comparing parent and subgroup operations.
    """
    return (
        tuple(np.round(op.r().as_double(), 12)),
        tuple(np.round(op.t().as_double(), 12)),
    )


def metric_close(u, v, G, atol=1e-5):
    """
    Compare two fractional-coordinate vectors using the physical metric norm.
    """
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)

    diff_norm = norm_metric(u - v, G)
    scale = max(norm_metric(u, G), norm_metric(v, G), 1.0)

    return diff_norm <= atol * scale

