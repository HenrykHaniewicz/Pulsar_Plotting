# List of utilities (mostly from PRESTO).
# Called 'routines' to differentiates from all other packages where the utilities file is called 'utils', etc.

import numpy as np
import constants

def p_to_f(p, pd, pdd=None):
    """
    p_to_f(p, pd, pdd=None):
       Convert period, period derivative and period second
       derivative to the equivalent frequency counterparts.
       Will also convert from f to p.
    """
    f = 1.0 / p
    fd = -pd / (p * p)
    if (pdd is None):
        return [f, fd]
    else:
        if (pdd == 0.0):
            fdd = 0.0
        else:
            fdd = 2.0 * pd * pd / (p ** 3.0) - pdd / (p * p)
        return [f, fd, fdd]

def pferrs(porf, porferr, pdorfd=None, pdorfderr=None):
    """
    pferrs(porf, porferr, pdorfd=None, pdorfderr=None):
       Calculate the period or frequency errors and
       the pdot or fdot errors from the opposite one.
    """
    if (pdorfd is None):
        return [1.0 / porf, porferr / porf ** 2.0]
    else:
        forperr = porferr / porf ** 2.0
        fdorpderr = np.sqrt((4.0 * pdorfd ** 2.0 * porferr ** 2.0) / porf ** 6.0 + pdorfderr ** 2.0 / porf ** 4.0)
        [forp, fdorpd] = p_to_f(porf, pdorfd)
        return [forp, forperr, fdorpd, fdorpderr]

def rad_to_hms(rad):
    """
    rad_to_hms(rad):
       Convert radians to hours, minutes, and seconds of arc.
    """
    rad = np.fmod(rad, constants.TWOPI)
    if (rad < 0.0): rad = rad + constants.TWOPI
    arc = constants.RADTOHRS * rad
    h = int(arc)
    arc = (arc - h) * 60.0
    m = int(arc)
    s = (arc - m) * 60.0
    return (h, m, s)

def rad_to_dms(rad):
    """
    rad_to_dms(rad):
       Convert radians to degrees, minutes, and seconds of arc.
    """
    if (rad < 0.0):
        sign = -1
    else:
        sign = 1
    arc = constants.RADTODEG * np.fmod(np.fabs(rad), constants.PI)
    d = int(arc)
    arc = (arc - d) * 60.0
    m = int(arc)
    s = (arc - m) * 60.0
    if sign == -1 and d == 0:
        return (sign * d, sign * m, sign * s)
    else:
        return (sign * d, m, s)

def coord_to_string(h_or_d, m, s):
    """
    coord_to_string(h_or_d, m, s):
       Return a formatted string of RA or DEC values as
       'hh:mm:ss.ssss' if RA, or 'dd:mm:ss.ssss' if DEC.
    """
    retstr = ""
    if h_or_d < 0:
        retstr = "-"
    elif abs(h_or_d) == 0:
        if (m < 0.0) or (s < 0.0):
            retstr = "-"
    h_or_d, m, s = abs(h_or_d), abs(m), abs(s)
    if (s >= 9.9995):
        return retstr + f"{h_or_d:.2d}:{m:.2d}:{s:.4f}"
    else:
        return retstr + f"{h_or_d:.2d}:{m:.2d}:0{s:.4f}"

def pdot_from_B(p, B):
    """
    pdot_from_B(p, B):
        Return a pdot (or p, actually) that a pulsar with spin
        period (or pdot) 'p' (in sec) would experience given a
        magnetic field strength 'B' in gauss.
    """
    return (B / 3.2e19) ** 2.0 / p

def pdot_from_edot(p, edot, I = 1.0e45):
    """
    pdot_from_edot(p, edot, I=1.0e45):
        Return the pdot that a pulsar with spin period 'p (in sec)
        would experience given an Edot 'edot' (in ergs/s) and a
        moment of inertia I.
    """
    return (p ** 3.0 * edot) / (4.0 * constants.PI * constants.PI * I)

def pdot_from_age(p, age):
    """
    pdot_from_age(p, age):
        Return the pdot that a pulsar with spin period 'p' (in sec)
        would experience given a characteristic age 'age' (in yrs).
    """
    return p / (2.0 * age * constants.SECPERJULYR)
