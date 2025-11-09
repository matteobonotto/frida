import numpy as np

eps = np.finfo(np.double).resolution


def calc_complete_elliptic_integral(m, tol: float = eps):
    """
    Computes the value of the complete elliptic
    integrals of the first and second kinds, evaluated for each
    element of M.  As currently implemented, M is limited to 0 <= M <= 1.
    """
    assert isinstance(m, np.ndarray), "m must be ndarray"
    # if not isinstance(m, np.ndarray):
    if (m == 1).all():
        e = np.ones(m)
        return np.inf + e, e

    a0 = 1
    b0 = np.sqrt(1 - m)
    c0 = np.nan
    s0 = m
    i1 = 0
    mm = np.inf
    while mm > tol:
        a1 = (a0 + b0) / 2
        b1 = np.sqrt(a0 * b0)
        c1 = (a0 - b0) / 2
        i1 = i1 + 1
        w1 = 2**i1 * c1**2
        mm = w1.max()

        # test for stagnation (may happen for TOL < machine precision)
        if np.allclose(c0, c1, atol=tol):
            raise ValueError("Convergence failed.")

        s0 = s0 + w1
        a0 = a1
        b0 = b1
        c0 = c1

    k = np.pi / (2 * a1)
    e = k * (1 - s0 / 2)
    idx = np.argwhere(m == 1)
    if idx.size != 0:
        k[idx] = np.inf
        e[idx] = 1.0

    return k, e


def green_flux_filament(source, target):

    RR, ZZ = target[:, 0], target[:, 1]
    flux = np.zeros(target.shape[0])
    factor = 4e-7 * 2 * np.pi * RR

    idx_axis = np.argwhere(target[:, 0] == 0)
    RR[idx_axis] += 1e-14

    for r, z, current in source:

        kk = np.sqrt(4 * r * RR / ((r + RR) ** 2 + (ZZ - z) ** 2))
        kk[idx_axis] += 1e-14

        kk_square = kk**2
        kk_square[kk_square > 1] = 1

        e, k = calc_complete_elliptic_integral(kk_square)
        flux_ = factor / kk * np.sqrt(r / RR) * ((1 - kk**2 / 2) * e - k)
        flux += flux_ * current

    # impose flux = 0 on target points where r=0
    flux[idx_axis] = 0.0

    return flux


def green_brbz_filament(source, target):

    RR, ZZ = target[:, 0], target[:, 1]
    br = np.zeros(target.shape[0])
    bz = np.zeros(target.shape[0])
    factor = 4e-7 * 2 * np.pi * RR

    idx_axis = np.argwhere(target[:, 0] == 0)
    RR[idx_axis] += 1e-14

    for r, z, current in source:
        dr_squared = (r - RR) ** 2
        dz_squared = (ZZ - z) ** 2

        kk = np.sqrt(4 * r * RR / ((r + RR) ** 2 + dz_squared))
        kk_square = kk**2
        kk_square[kk_square > 1] = 1.0

        e, k = calc_complete_elliptic_integral(kk_square)

        br_ = (
            1e-7
            * kk
            * (ZZ - z)
            / (RR * np.sqrt(r * RR))
            * (-e + (r**2 + RR**2 + dz_squared) / (dr_squared + dz_squared) * k)
        )
        bz_ = (
            1e-7
            * kk
            / np.sqrt(r * RR)
            * (e + (r**2 - RR**2 - dz_squared) / (dr_squared + dz_squared) * k)
        )

        br += br_ * current
        bz += bz_ * current

    # Find points on axis (r=0)
    br[idx_axis] = 0
    if idx_axis.size > 0:
        bz[idx_axis] = 0
        factor = 4e-7 * np.pi
        for r, z, current in source:
            bz[idx_axis] += factor * r / (2 * (r**2 + (ZZ[idx_axis] - z) ** 2) ** 1.5)

    return br, bz
