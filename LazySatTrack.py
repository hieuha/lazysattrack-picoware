# LazySatTrack (single-file Picoware app) -- offline SGP4 pass tracker.
# Drop THIS ONE FILE into /sd/picoware/apps/ and open it from Library.
# sgp4lite is embedded below; TLEs load from /sd/tles.txt if present,
# else a built-in 2021 sample is used. Near-Earth (LEO) only.
#
#   // weak transmissions from an imagined future -- LazySatTrack
# ===========================================================================

# ============================ embedded sgp4lite ============================
# sgp4lite.py -- compact near-Earth SGP4 (WGS-72), pure Python / MicroPython safe.
# Deep-space (SDP4) NOT included: only for objects with period < ~225 min (LEO).
# Returns TEME position (km) & velocity (km/s). Based on the public SGP4 algorithm
# (Hoots/Vallado, "Revisiting Spacetrack Report #3"), near-Earth branch only.
import math

pi = math.pi
twopi = 2.0 * pi
deg2rad = pi / 180.0
x2o3 = 2.0 / 3.0

# WGS-72 constants (the set TLEs are fit to)
mu = 398600.8
radiusearthkm = 6378.135
xke = 60.0 / math.sqrt(radiusearthkm * radiusearthkm * radiusearthkm / mu)
tumin = 1.0 / xke
j2 = 0.001082616
j3 = -0.00000253881
j4 = -0.00000165597
j3oj2 = j3 / j2


def _expfield(s):
    """Parse TLE exponential field like ' 51432-3' -> 0.51432e-3."""
    s = s.strip()
    if not s or s in ('0', '00000-0', '00000+0'):
        return 0.0
    sign = 1.0
    if s[0] == '-':
        sign = -1.0; s = s[1:]
    elif s[0] == '+':
        s = s[1:]
    ex = 0
    mant = s
    for i in range(len(s) - 1, 0, -1):
        if s[i] in '+-':
            mant = s[:i]; ex = int(s[i:]); break
    return sign * float('0.' + mant) * (10.0 ** ex)


class Sat:
    __slots__ = ('name', 'satnum', 'epochyr', 'epochdays', 'jdepoch', 'jdepochF', 'bstar', 'inclo',
                 'nodeo', 'ecco', 'argpo', 'mo', 'no', 'a', 'alta', 'altp',
                 'error', 'method',
                 'aycof', 'con41', 'cc1', 'cc4', 'cc5', 'd2', 'd3', 'd4',
                 'delmo', 'eta', 'argpdot', 'omgcof', 'sinmao', 't2cof', 't3cof',
                 't4cof', 't5cof', 'x1mth2', 'x7thm1', 'mdot', 'nodedot',
                 'xlcof', 'xmcof', 'nodecf', 'isimp', 'gsto')


def jday(yr, mon, day, hr, minute, sec):
    # Split Julian Date: whole part (an exact X.5 half-integer, so it survives
    # 32-bit float without loss) plus the sub-day fraction in [0,1). Single-
    # precision MicroPython cannot hold an absolute JD (~2.46e6) to sub-day
    # resolution -- a 60 s step is below its ULP -- so time is carried as this
    # pair everywhere and the two are only ever combined via exact integer diffs.
    jdw = (367.0 * yr - int(7 * (yr + int((mon + 9) / 12.0)) * 0.25)
           + int(275 * mon / 9.0) + day + 1721013.5)
    jdf = ((sec / 60.0 + minute) / 60.0 + hr) / 24.0
    return jdw, jdf


def parse_tle(line1, line2):
    s = Sat()
    s.error = 0
    s.satnum = line1[2:7].strip()
    two_yr = int(line1[18:20])
    s.epochyr = two_yr + 2000 if two_yr < 57 else two_yr + 1900
    s.epochdays = float(line1[20:32])
    s.bstar = _expfield(line1[53:61])
    s.inclo = float(line2[8:16]) * deg2rad
    s.nodeo = float(line2[17:25]) * deg2rad
    s.ecco = float('0.' + line2[26:33].strip())
    s.argpo = float(line2[34:42]) * deg2rad
    s.mo = float(line2[43:51]) * deg2rad
    no_revday = float(line2[52:63])
    s.no = no_revday * twopi / 1440.0  # rad/min (kozai)
    # epoch julian date
    yr = s.epochyr
    days = s.epochdays
    mon, day, hr, mi, sc = _days2mdhms(yr, days)
    s.jdepoch, s.jdepochF = jday(yr, mon, day, hr, mi, sc)
    _sgp4init(s)
    return s


def _days2mdhms(year, days):
    lmonth = [31, 29 if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)) else 28,
              31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    dayofyr = int(days)
    i = 0; inttemp = 0
    while dayofyr > inttemp + lmonth[i] and i < 11:
        inttemp += lmonth[i]; i += 1
    mon = i + 1
    day = dayofyr - inttemp
    temp = (days - dayofyr) * 24.0
    hr = int(temp)
    temp = (temp - hr) * 60.0
    mi = int(temp)
    sec = (temp - mi) * 60.0
    return mon, day, hr, mi, sec


def _gstime(jdw, jdf):
    # Greenwich sidereal time from a split JD (whole + fraction). Subtracting the
    # J2000 epoch from the whole part is an exact integer difference in float32.
    tut1 = ((jdw - 2451545.0) + jdf) / 36525.0
    t = (-6.2e-6 * tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1
         + (876600.0 * 3600 + 8640184.812866) * tut1 + 67310.54841)
    t = math.fmod(t * deg2rad / 240.0, twopi)
    if t < 0.0:
        t += twopi
    return t


def _sgp4init(s):
    ecco = s.ecco; epoch = (s.jdepoch - 2433281.5) + s.jdepochF
    inclo = s.inclo
    no = s.no
    # recover original mean motion (no_unkozai)
    cosio = math.cos(inclo); cosio2 = cosio * cosio
    eccsq = ecco * ecco; omeosq = 1.0 - eccsq; rteosq = math.sqrt(omeosq)
    x3thm1 = 3.0 * cosio2 - 1.0
    ak = (xke / no) ** x2o3
    d1 = 0.75 * j2 * x3thm1 / (rteosq * omeosq)
    del_ = d1 / (ak * ak)
    adel = ak * (1.0 - del_ * del_ - del_ * (1.0 / 3.0 + 134.0 * del_ * del_ / 81.0))
    del_ = d1 / (adel * adel)
    no = no / (1.0 + del_)          # no_unkozai
    s.no = no
    ao = (xke / no) ** x2o3
    sinio = math.sin(inclo)
    po = ao * omeosq
    con42 = 1.0 - 5.0 * cosio2
    s.con41 = -con42 - cosio2 - cosio2       # = 3cos^2 - 1 ... check: -(1-5c2)-2c2 = -1+5c2-2c2=3c2-1 OK
    con41 = s.con41
    ainv = 1.0 / ao
    posq = po * po
    rp = ao * (1.0 - ecco)
    s.a = ao
    s.altp = rp - 1.0
    s.alta = ao * (1.0 + ecco) - 1.0
    s.gsto = _gstime(s.jdepoch, s.jdepochF)

    s.isimp = 0
    if rp < (220.0 / radiusearthkm + 1.0):
        s.isimp = 1

    sfour = 78.0 / radiusearthkm + 1.0
    qzms24 = ((120.0 - 78.0) / radiusearthkm) ** 4
    perige = (rp - 1.0) * radiusearthkm
    if perige < 156.0:
        sfour = perige - 78.0
        if perige < 98.0:
            sfour = 20.0
        qzms24 = ((120.0 - sfour) / radiusearthkm) ** 4
        sfour = sfour / radiusearthkm + 1.0

    pinvsq = 1.0 / posq
    tsi = 1.0 / (ao - sfour)
    eta = ao * ecco * tsi
    etasq = eta * eta
    eeta = ecco * eta
    psisq = abs(1.0 - etasq)
    coef = qzms24 * (tsi ** 4)
    coef1 = coef / (psisq ** 3.5)
    cc2 = (coef1 * no * (ao * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq))
           + 0.375 * j2 * tsi / psisq * con41 * (8.0 + 3.0 * etasq * (8.0 + etasq))))
    s.cc1 = s.bstar * cc2
    cc3 = 0.0
    if ecco > 1.0e-4:
        cc3 = -2.0 * coef * tsi * j3oj2 * no * sinio / ecco
    s.x1mth2 = 1.0 - cosio2
    s.cc4 = (2.0 * no * coef1 * ao * omeosq *
             (eta * (2.0 + 0.5 * etasq) + ecco * (0.5 + 2.0 * etasq)
              - j2 * tsi / (ao * psisq) *
              (-3.0 * con41 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta))
               + 0.75 * s.x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) *
               math.cos(2.0 * s.argpo))))
    s.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq)
    cosio4 = cosio2 * cosio2
    temp1 = 1.5 * j2 * pinvsq * no
    temp2 = 0.5 * temp1 * j2 * pinvsq
    temp3 = -0.46875 * j4 * pinvsq * pinvsq * no
    s.mdot = (no + 0.5 * temp1 * rteosq * con41 + 0.0625 * temp2 * rteosq *
              (13.0 - 78.0 * cosio2 + 137.0 * cosio4))
    s.argpdot = (-0.5 * temp1 * con42 + 0.0625 * temp2 *
                 (7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
                 temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4))
    xhdot1 = -temp1 * cosio
    s.nodedot = (xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) +
                 2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio)
    xpidot = s.argpdot + s.nodedot
    s.omgcof = s.bstar * cc3 * math.cos(s.argpo)
    s.xmcof = 0.0
    if ecco > 1.0e-4:
        s.xmcof = -x2o3 * coef * s.bstar / eeta
    s.nodecf = 3.5 * omeosq * xhdot1 * s.cc1
    s.t2cof = 1.5 * s.cc1
    if abs(cosio + 1.0) > 1.5e-12:
        s.xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio)
    else:
        s.xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / 1.5e-12
    s.aycof = -0.5 * j3oj2 * sinio
    s.delmo = (1.0 + eta * math.cos(s.mo)) ** 3
    s.sinmao = math.sin(s.mo)
    s.x7thm1 = 7.0 * cosio2 - 1.0
    s.eta = eta

    # drag terms
    if s.isimp != 1:
        cc1sq = s.cc1 * s.cc1
        s.d2 = 4.0 * ao * tsi * cc1sq
        temp = s.d2 * tsi * s.cc1 / 3.0
        s.d3 = (17.0 * ao + sfour) * temp
        s.d4 = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) * s.cc1
        s.t3cof = s.d2 + 2.0 * cc1sq
        s.t4cof = 0.25 * (3.0 * s.d3 + s.cc1 * (12.0 * s.d2 + 10.0 * cc1sq))
        s.t5cof = (0.2 * (3.0 * s.d4 + 12.0 * s.cc1 * s.d3 + 6.0 * s.d2 * s.d2 +
                   15.0 * cc1sq * (2.0 * s.d2 + cc1sq)))
    else:
        s.d2 = s.d3 = s.d4 = 0.0
        s.t3cof = s.t4cof = s.t5cof = 0.0


def sgp4(s, tsince):
    """Propagate. tsince in minutes from epoch. Returns (r[3] km, v[3] km/s)."""
    s.error = 0
    xmdf = s.mo + s.mdot * tsince
    argpdf = s.argpo + s.argpdot * tsince
    nodedf = s.nodeo + s.nodedot * tsince
    argpm = argpdf
    mm = xmdf
    t2 = tsince * tsince
    nodem = nodedf + s.nodecf * t2
    tempa = 1.0 - s.cc1 * tsince
    tempe = s.bstar * s.cc4 * tsince
    templ = s.t2cof * t2

    if s.isimp != 1:
        delomg = s.omgcof * tsince
        delm = s.xmcof * ((1.0 + s.eta * math.cos(xmdf)) ** 3 - s.delmo)
        temp = delomg + delm
        mm = xmdf + temp
        argpm = argpdf - temp
        t3 = t2 * tsince
        t4 = t3 * tsince
        tempa = tempa - s.d2 * t2 - s.d3 * t3 - s.d4 * t4
        tempe = tempe + s.bstar * s.cc5 * (math.sin(mm) - s.sinmao)
        templ = templ + s.t3cof * t3 + t4 * (s.t4cof + tsince * s.t5cof)

    nm = s.no
    em = s.ecco
    inclm = s.inclo

    am = ((xke / nm) ** x2o3) * tempa * tempa
    nm = xke / (am ** 1.5)
    em = em - tempe
    if em >= 1.0 or em < -0.001:
        s.error = 1
        return None, None
    if em < 1.0e-6:
        em = 1.0e-6

    mm = mm + s.no * templ
    xlm = mm + argpm + nodem
    nodem = math.fmod(nodem, twopi)
    argpm = math.fmod(argpm, twopi)
    xlm = math.fmod(xlm, twopi)
    mm = math.fmod(xlm - argpm - nodem, twopi)

    sinim = math.sin(inclm); cosim = math.cos(inclm)
    ep = em
    xincp = inclm
    argpp = argpm
    nodep = nodem
    mp = mm
    sinip = sinim; cosip = cosim

    # long-period periodics
    axnl = ep * math.cos(argpp)
    temp = 1.0 / (am * (1.0 - ep * ep))
    aynl = ep * math.sin(argpp) + temp * s.aycof
    xl = mp + argpp + nodep + temp * s.xlcof * axnl

    # kepler for E+w
    u = math.fmod(xl - nodep, twopi)
    eo1 = u
    for _ in range(10):
        sineo1 = math.sin(eo1); coseo1 = math.cos(eo1)
        tem5 = 1.0 - coseo1 * axnl - sineo1 * aynl
        tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5
        if abs(tem5) >= 0.95:
            tem5 = 0.95 if tem5 > 0 else -0.95
        eo1 = eo1 + tem5
        if abs(tem5) < 1.0e-12:
            break

    ecose = axnl * coseo1 + aynl * sineo1
    esine = axnl * sineo1 - aynl * coseo1
    el2 = axnl * axnl + aynl * aynl
    pl = am * (1.0 - el2)
    if pl < 0.0:
        s.error = 2
        return None, None
    rl = am * (1.0 - ecose)
    rdotl = math.sqrt(am) * esine / rl
    rvdotl = math.sqrt(pl) / rl
    betal = math.sqrt(1.0 - el2)
    temp = esine / (1.0 + betal)
    sinu = am / rl * (sineo1 - aynl - axnl * temp)
    cosu = am / rl * (coseo1 - axnl + aynl * temp)
    su = math.atan2(sinu, cosu)
    sin2u = (cosu + cosu) * sinu
    cos2u = 1.0 - 2.0 * sinu * sinu
    temp = 1.0 / pl
    temp1 = 0.5 * j2 * temp
    temp2 = temp1 * temp

    # short-period periodics
    cosisq = cosip * cosip
    con41 = 3.0 * cosisq - 1.0
    x1mth2 = 1.0 - cosisq
    x7thm1 = 7.0 * cosisq - 1.0
    mrt = rl * (1.0 - 1.5 * temp2 * betal * con41) + 0.5 * temp1 * x1mth2 * cos2u
    su = su - 0.25 * temp2 * x7thm1 * sin2u
    xnode = nodep + 1.5 * temp2 * cosip * sin2u
    xinc = xincp + 1.5 * temp2 * cosip * sinip * cos2u
    mvt = rdotl - nm * temp1 * x1mth2 * sin2u / xke
    rvdot = rvdotl + nm * temp1 * (x1mth2 * cos2u + 1.5 * con41) / xke

    # orientation vectors
    sinsu = math.sin(su); cossu = math.cos(su)
    snod = math.sin(xnode); cnod = math.cos(xnode)
    sini = math.sin(xinc); cosi = math.cos(xinc)
    xmx = -snod * cosi; xmy = cnod * cosi
    ux = xmx * sinsu + cnod * cossu
    uy = xmy * sinsu + snod * cossu
    uz = sini * sinsu
    vx = xmx * cossu - cnod * sinsu
    vy = xmy * cossu - snod * sinsu
    vz = sini * cossu

    r = [mrt * ux * radiusearthkm, mrt * uy * radiusearthkm, mrt * uz * radiusearthkm]
    vkmpersec = radiusearthkm * xke / 60.0
    v = [(mvt * ux + rvdot * vx) * vkmpersec,
         (mvt * uy + rvdot * vy) * vkmpersec,
         (mvt * uz + rvdot * vz) * vkmpersec]
    return r, v


gstime = _gstime


# ============================ LazySatTrack app ============================
# LazySatTrack -- offline SGP4 satellite pass tracker, ported to Picoware.
# Drop this file + py into /sd/picoware/apps/ ; it appears in Library.
# TLEs live at /sd/tles.txt (3-line sets). Near-Earth (LEO) only.
#
#   // weak transmissions from an imagined future -- LazySatTrack
#
# Controls (PicoCalc keyboard):
#   UP/DOWN  prev/next satellite     LEFT/RIGHT  select pass (updates sky plot)
#   L        find / search picker    C           recompute passes
#   BACK/ESC exit (or leave picker)
# In the picker: type digits=NORAD id / letters=name, UP/DOWN move,
#                ENTER track, BACKSPACE edit, BACK/ESC cancel.
import math
import gc

from picoware.system.vector import Vector
from picoware.system import buttons as B

# ============================ CONFIG (edit me) =============================
LAT = 21.0278        # ground station latitude  (deg, +N)
LON = 105.8342       # ground station longitude (deg, +E)
ALT = 0.01           # altitude (km)
TZ = 7.0             # display timezone offset (hours from UTC)
HORIZON_H = 24       # look-ahead window (hours)
MIN_PEAK = 10.0      # ignore passes below this max elevation (deg)
STEP_S = 60          # pass-search coarse step (s); AOS/LOS bisection-refined after
CACHE_H = 6.0        # reuse a satellite's computed passes for this many hours
                     # before recomputing (press C to force a fresh recompute)
TLE_PATHS = ("/sd/tles.txt", "/sd/picoware/apps/tles.txt", "tles.txt")

# fallback catalog so the app runs even with no /sd/tles.txt (2021 epoch sample)
_DEFAULT_TLES = (
    "ISS (ZARYA)\n"
    "1 25544U 98067A   21122.75616700  .00027980  00000-0  51432-3 0  9994\n"
    "2 25544  51.6442 207.4449 0002769 310.1189 193.6568 15.48993527281553\n"
    "HST\n"
    "1 20580U 90037B   21123.53321839  .00001264  00000-0  63000-4 0  9990\n"
    "2 20580  28.4707 288.8102 0002830 291.9319 244.0674 15.09299865499811\n"
    "NOAA 19\n"
    "1 33591U 09005A   21123.55575694  .00000075  00000-0  67238-4 0  9998\n"
    "2 33591  99.1912 129.3324 0013711 288.9666  71.0177 14.12466021627079\n"
)
# ==========================================================================

pi = math.pi
d2r = pi / 180.0
r2d = 180.0 / pi
WGS_A = 6378.137
WGS_F = 1 / 298.257223563
WGS_E2 = WGS_F * (2 - WGS_F)
RE = radiusearthkm

# Optional ground-station override file on the SD card. If present it wins over the
# LAT/LON/ALT/TZ defaults above, so the station can be moved without recompiling.
_STATION_PATHS = ("/sd/station.txt", "/sd/picoware/apps/station.txt", "station.txt")


def _load_station():
    """Read LAT/LON/ALT/TZ overrides from the first station file found.

    Accepts either 'KEY value' lines (KEY in LAT/LON/ALT/TZ; '=' or ',' also work)
    or one positional line 'lat lon [alt] [tz]'. Blank lines and '#'/';' comments
    are ignored. Returns a dict of the values it managed to parse (may be empty).
    """
    keys = ("LAT", "LON", "ALT", "TZ")
    for p in _STATION_PATHS:
        try:
            f = open(p)
        except Exception:
            continue
        vals = {}
        try:
            for raw in f:
                line = raw.strip()
                if not line or line[0] in "#;":
                    continue
                toks = line.replace(",", " ").replace("=", " ").split()
                if not toks:
                    continue
                k = toks[0].upper()
                if k in keys and len(toks) >= 2:
                    try:
                        vals[k] = float(toks[1])
                    except Exception:
                        pass
                else:                              # positional: lat lon [alt] [tz]
                    for name, t in zip(keys, toks):
                        try:
                            vals[name] = float(t)
                        except Exception:
                            pass
        finally:
            f.close()
        if vals:
            return vals
    return {}


_st = _load_station()
LAT = _st.get("LAT", LAT)
LON = _st.get("LON", LON)
ALT = _st.get("ALT", ALT)
TZ = _st.get("TZ", TZ)

# Ground-station latitude/longitude sin&cos, precomputed once. _look() is called
# thousands of times during a pass search, so hoisting these out of the hot loop
# avoids ~8 trig calls per propagation step.
_SLA = math.sin(LAT * d2r)
_CLA = math.cos(LAT * d2r)
_SLO = math.sin(LON * d2r)
_CLO = math.cos(LON * d2r)

# button->char maps for the picker filter (robust: built from constants)
_DIGITS = {getattr(B, "BUTTON_%d" % i): str(i) for i in range(10)}
_LETTERS = {getattr(B, "BUTTON_%s" % c): c.lower() for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"}

# ---- module-level state (Picoware apps are frame-driven; persist via globals) --
_state = "track"
_sats = []
_si = 0
_passes = []
_arc = None
_sel = 0
_filter = ""
_obs = None
_col = {}
_fs = None            # font size Vector
_last_sec = -1
_need_redraw = True
_msg = None           # transient status line
_JD0 = 0.0            # whole-day reference JD (X.5, exact in float32); times are
                      # carried as small day-offsets relative to this, set per recompute
_cache = {}           # _si -> {"passes","jd0","jf0"}: computed passes per satellite
_clock_ok = False     # True once the RTC holds a real date (>=2020, i.e. NTP-synced)


# ============================ ORBITAL GEOMETRY ============================
def _observer_ecef(lat, lon, alt):
    la = lat * d2r
    lo = lon * d2r
    sN = math.sin(la)
    cN = math.cos(la)
    N = WGS_A / math.sqrt(1 - WGS_E2 * sN * sN)
    return ((N + alt) * cN * math.cos(lo),
            (N + alt) * cN * math.sin(lo),
            (N * (1 - WGS_E2) + alt) * sN)


def _look(sat, f, obs):
    # f is a day-offset relative to the module reference _JD0. Minutes-since-epoch
    # is built from the exact whole-day difference plus the small fractional parts.
    tmin = ((_JD0 - sat.jdepoch) + (f - sat.jdepochF)) * 1440.0
    r, _ = sgp4(sat, tmin)
    if r is None:
        return None
    g = _gstime(_JD0, f)
    c = math.cos(g)
    s = math.sin(g)
    ex = c * r[0] + s * r[1]
    ey = -s * r[0] + c * r[1]
    ez = r[2]
    rx = ex - obs[0]
    ry = ey - obs[1]
    rz = ez - obs[2]
    sla = _SLA
    cla = _CLA
    slo = _SLO
    clo = _CLO
    rS = sla * clo * rx + sla * slo * ry - cla * rz
    rE = -slo * rx + clo * ry
    rZ = cla * clo * rx + cla * slo * ry + sla * rz
    rng = math.sqrt(rS * rS + rE * rE + rZ * rZ)
    el = math.asin(rZ / rng) * r2d
    az = math.atan2(rE, -rS) * r2d
    if az < 0:
        az += 360
    subLon = math.atan2(ey, ex) * r2d
    subLat = math.atan2(ez, math.sqrt(ex * ex + ey * ey)) * r2d
    alt = math.sqrt(ex * ex + ey * ey + ez * ez) - RE
    return {"az": az, "el": el, "rng": rng, "subLat": subLat, "subLon": subLon, "alt": alt}


def _elev(sat, f, obs):
    # Elevation-only path for the pass search hot loop: no dict allocation and
    # skips the azimuth / sub-point / altitude trig that _look() computes.
    # f is a day-offset relative to _JD0 (see _look for the split-JD rationale).
    tmin = ((_JD0 - sat.jdepoch) + (f - sat.jdepochF)) * 1440.0
    r, _ = sgp4(sat, tmin)
    if r is None:
        return -90.0
    g = _gstime(_JD0, f)
    c = math.cos(g)
    s = math.sin(g)
    ex = c * r[0] + s * r[1]
    ey = -s * r[0] + c * r[1]
    rx = ex - obs[0]
    ry = ey - obs[1]
    rz = r[2] - obs[2]
    rS = _SLA * _CLO * rx + _SLA * _SLO * ry - _CLA * rz
    rE = -_SLO * rx + _CLO * ry
    rZ = _CLA * _CLO * rx + _CLA * _SLO * ry + _SLA * rz
    return math.asin(rZ / math.sqrt(rS * rS + rE * rE + rZ * rZ)) * r2d


def _refine(sat, a, b, obs):
    ea = _elev(sat, a, obs)
    for _ in range(16):
        m = (a + b) / 2
        em = _elev(sat, m, obs)
        if (ea < 0) == (em < 0):
            a = m
            ea = em
        else:
            b = m
    return (a + b) / 2


def _find_passes(sat, f0, obs):
    # f0 and every time below are day-offsets relative to _JD0. Stepping a small
    # offset (rather than an absolute ~2.46e6 JD) is what keeps the increment above
    # the float32 ULP, so the loop actually advances instead of stalling forever.
    dstep = STEP_S / 86400.0
    end = f0 + HORIZON_H / 24.0
    out = []
    prev = _elev(sat, f0, obs)
    aos = f0 if prev >= 0 else None
    mx = -90.0
    mt = f0
    f = f0 + dstep
    while f <= end:
        el = _elev(sat, f, obs)
        if prev < 0 and el >= 0:
            aos = _refine(sat, f - dstep, f, obs)
            mx = -90.0
            mt = aos
        if aos is not None and el > mx:
            mx = el
            mt = f
        if prev >= 0 and el < 0:
            los = _refine(sat, f - dstep, f, obs)
            a = mt - dstep
            b = mt + dstep
            for _ in range(18):
                m1 = a + (b - a) / 3
                m2 = b - (b - a) / 3
                if _elev(sat, m1, obs) < _elev(sat, m2, obs):
                    a = m1
                else:
                    b = m2
            pk = (a + b) / 2
            pl = _look(sat, pk, obs)
            al = _look(sat, aos, obs)
            ll = _look(sat, los, obs)
            if pl and pl["el"] >= MIN_PEAK:
                out.append({"aos": aos, "los": los, "peak": pk, "maxel": pl["el"],
                            "aosaz": al["az"], "losaz": ll["az"], "peakaz": pl["az"],
                            "dur": (los - aos) * 86400.0})
            aos = None
            if len(out) >= 8:
                break
        prev = el
        f += dstep
    return out


# ============================ TIME HELPERS ===============================
def _now_jd_parts():
    """Current time as a split JD (whole, fraction). Falls back to 2021-05-03
    when the RTC is unset so the bundled 2021 sample TLEs still give sane passes."""
    try:
        from machine import RTC
        dt = RTC().datetime()   # (Y,M,D,wd,h,m,s,subs)
        return jday(dt[0], dt[1], dt[2], dt[4], dt[5], dt[6])
    except Exception:
        return jday(2021, 5, 3, 0, 0, 0)


def _now_set():
    """Latch the current whole-day JD into _JD0 and return the fractional offset."""
    global _JD0
    jw, jf = _now_jd_parts()
    _JD0 = jw
    return jf


def _now_f():
    """Current time as a day-offset relative to _JD0 (whole-day diff stays exact)."""
    jw, jf = _now_jd_parts()
    return (jw - _JD0) + jf


def _clock_valid():
    """True when the RTC holds a plausible real date (year >= 2020). An unset RTC
    reads back as year 2000/2001, which would make every pass meaningless."""
    try:
        from machine import RTC
        return RTC().datetime()[0] >= 2020
    except Exception:
        return False


def _sync_time(vm):
    """Kick off a one-shot NTP sync (async) if the clock isn't set yet. Needs WiFi
    connected; the RTC is written in UTC (offset 0) to match the Z-time display."""
    try:
        t = vm.time
        if t and not t.is_set and not t.is_fetching and not _clock_valid():
            t.fetch(0)
    except Exception:
        pass


def _jd_to_cal(f):
    # f is a day-offset relative to _JD0. Recombine without adding a small value to
    # the large whole part: fold whole days of f into the (exact integer) day count.
    import math as _m
    fl = _m.floor(f)
    Z = int(_JD0 + 0.5) + int(fl)
    F = f - fl
    if Z < 2299161:
        A = Z
    else:
        al = int((Z - 1867216.25) / 36524.25)
        A = Z + 1 + al - int(al / 4)
    Bb = A + 1524
    C = int((Bb - 122.1) / 365.25)
    D = int(365.25 * C)
    E = int((Bb - D) / 30.6001)
    day = Bb - D - int(30.6001 * E) + F
    mon = E - 1 if E < 14 else E - 13
    yr = C - 4716 if mon > 2 else C - 4715
    d = int(day)
    fr = (day - d) * 24
    h = int(fr)
    fr = (fr - h) * 60
    mi = int(fr)
    s = int((fr - mi) * 60)
    return yr, mon, d, h, mi, s


def _hms(jd, tz):
    _, _, _, h, mi, s = _jd_to_cal(jd + tz / 24.0)
    return "%02d:%02d:%02d" % (h, mi, s)


def _hm(jd, tz):
    _, _, _, h, mi, _ = _jd_to_cal(jd + tz / 24.0)
    return "%02d:%02d" % (h, mi)


def _datestr(jd, tz):
    _, _, d, h, mi, _ = _jd_to_cal(jd + tz / 24.0)
    return "%02d/%02d %02d:%02d" % (d, __mon(jd, tz), h, mi)


def __mon(jd, tz):
    _, mo, _, _, _, _ = _jd_to_cal(jd + tz / 24.0)
    return mo


def _cd(sec):
    if sec < 0:
        return "--:--:--"
    sec = int(sec)
    h = sec // 3600
    sec -= h * 3600
    m = sec // 60
    sec -= m * 60
    return "%02d:%02d:%02d" % (h, m, sec)


_COMP = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
         "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]


def _comp(az):
    return _COMP[int(((az % 360) + 11.25) // 22.5) % 16]


# ============================ TLE LOADING ================================
def _load_tles():
    txt = None
    for p in TLE_PATHS:
        try:
            f = open(p)
            txt = f.read()
            f.close()
            break
        except OSError:
            continue
    if not txt:
        txt = _DEFAULT_TLES     # embedded fallback -> single-file, zero-config
    lines = [l.strip() for l in txt.split("\n") if l.strip()]
    sats = []
    i = 0
    while i < len(lines):
        if i + 2 < len(lines) and lines[i + 1].startswith("1 ") and lines[i + 2].startswith("2 "):
            name = lines[i][:16]
            l1 = lines[i + 1]
            l2 = lines[i + 2]
            i += 3
        elif lines[i].startswith("1 ") and i + 1 < len(lines) and lines[i + 1].startswith("2 "):
            name = "SAT " + lines[i][2:7]
            l1 = lines[i]
            l2 = lines[i + 1]
            i += 2
        else:
            i += 1
            continue
        try:
            s = parse_tle(l1, l2)
            s.name = name
            sats.append(s)
        except Exception:
            pass
    return sats


# ============================ DRAWING ====================================
CX = 78
CY = 124             # polar-plot centre; sits just under the header (was lower, left a top gap)
PR = 70


def _polar_xy(az, el):
    if el < 0:
        el = 0
    r = PR * (1 - el / 90.0)
    return int(CX + r * math.sin(az * d2r)), int(CY - r * math.cos(az * d2r))


def _build_arc(sat, p, obs):
    # Each point carries the polar-plot (x, y) plus the raw elevation (deg), so the
    # same sampling feeds both the sky plot and the elevation-vs-time chart. The 33
    # samples are evenly spaced in time from AOS to LOS, i.e. index i maps linearly
    # to the pass timeline.
    pts = []
    for i in range(33):
        f = p["aos"] + (p["los"] - p["aos"]) * i / 32.0
        l = _look(sat, f, obs)
        if l:
            el = max(0.0, l["el"])
            x, y = _polar_xy(l["az"], el)
            pts.append((x, y, el))
    return pts


def _T(d, x, y, s, c):
    d.text(Vector(int(x), int(y)), s, c)


def _LN(d, x0, y0, x1, y1, c):
    d.line_custom(Vector(int(x0), int(y0)), Vector(int(x1), int(y1)), c)


def _FR(d, x, y, w, h, c):
    d.fill_rectangle(Vector(int(x), int(y)), Vector(int(w), int(h)), c)


def _CI(d, x, y, r, c):
    d.circle(Vector(int(x), int(y)), int(r), c)


def _draw_track(vm):
    global _last_sec
    d = vm.draw
    ch = _fs.y if _fs else 10
    d.clear(Vector(0, 0), d.size, _col["bg"])
    sat = _sats[_si]
    f = _now_f()                       # current time as day-offset from _JD0

    _T(d, 4, 2, "LAZYSATTRACK", _col["amber"])
    if _clock_ok:
        _T(d, 232, 2, _hms(f, 0) + "Z", _col["green"])
    else:
        _T(d, 232, 2, "NTP SYNC..", _col["amber"])
    y0 = 2 + ch + 1
    _LN(d, 0, y0, 319, y0, _col["dim"])

    hid = getattr(sat, "satnum", "?")
    _T(d, 4, y0 + 3, "#" + hid, _col["cyan"])
    _T(d, 4 + (len(hid) + 2) * (_fs.x if _fs else 8), y0 + 3, getattr(sat, "name", "?")[:14], _col["white"])
    age = int((_JD0 - sat.jdepoch) + (f - sat.jdepochF))
    _T(d, 280, y0 + 3, "%dd" % age, _col["red"] if age > 14 else _col["dim"])

    # polar plot (brighter grid so the rings and compass letters read clearly)
    for e in (0, 30, 60):
        _CI(d, CX, CY, PR * (1 - e / 90.0), _col["grid"])
    _LN(d, CX - PR, CY, CX + PR, CY, _col["grid"])
    _LN(d, CX, CY - PR, CX, CY + PR, _col["grid"])
    _T(d, CX - 3, CY - PR - 10, "N", _col["grid"])
    _T(d, CX - 3, CY + PR + 2, "S", _col["grid"])
    _T(d, CX + PR + 3, CY - 4, "E", _col["grid"])
    _T(d, CX - PR - 9, CY - 4, "W", _col["grid"])
    if _arc:
        for i in range(1, len(_arc)):
            _LN(d, _arc[i - 1][0], _arc[i - 1][1], _arc[i][0], _arc[i][1], _col["amber"])
        _FR(d, _arc[0][0] - 1, _arc[0][1] - 1, 3, 3, _col["green"])
        _FR(d, _arc[-1][0] - 1, _arc[-1][1] - 1, 3, 3, _col["red"])
    live = _look(sat, f, _obs)
    if live and live["el"] > 0:
        x, y = _polar_xy(live["az"], live["el"])
        _FR(d, x - 2, y - 2, 5, 5, _col["green"])
        _CI(d, x, y, 5, _col["green"])

    # telemetry column
    tx = 170
    ty = y0 + 3 + ch + 4
    step = ch + 2

    def row(lbl, val, c):
        nonlocal ty
        _T(d, tx, ty, lbl, _col["dim"])
        _T(d, tx + 60, ty, val, c)
        ty += step
    if live:
        vis = live["el"] > 0
        row("STATE", "VISIBLE" if vis else "hidden", _col["green"] if vis else _col["dim"])
        row("AZ", "%5.1f %s" % (live["az"], _comp(live["az"])), _col["white"])
        row("EL", "%5.1f" % live["el"], _col["green"] if vis else _col["white"])
        row("RNG", "%6.0fkm" % live["rng"], _col["white"])
        row("SUB", "%.0f,%.0f" % (live["subLat"], live["subLon"]), _col["white"])
        row("ALT", "%5.0fkm" % live["alt"], _col["white"])
        ty += 2
        nxt = None
        for p in _passes:
            if f < p["aos"]:
                nxt = ("AOS", p["aos"])
                break
            if p["aos"] <= f <= p["los"]:
                nxt = ("LOS", p["los"])
                break
        if nxt:
            _T(d, tx, ty, "NEXT " + nxt[0], _col["cyan"])
            _T(d, tx, ty + step, "T-" + _cd((nxt[1] - f) * 86400), _col["cyan"])

    # selected-pass peak label + elevation-vs-time chart (fills the free lower-right)
    yb = CY + PR + 8
    if _passes and _sel < len(_passes):
        sp = _passes[_sel]
        pky = (y0 + 3 + ch + 4) + 8 * (ch + 2) + 2
        _T(d, tx, pky, "PEAK %2.0fd @%s" % (sp["maxel"], _hm(sp["peak"], TZ)), _col["amber"])
        bx = tx
        bw = 140
        bty = pky + ch + 2                 # chart top
        bby = yb - 4 - ch                  # chart baseline (leave a row for time labels)
        bh = bby - bty
        if _arc and len(_arc) >= 2 and bh > 20:
            for gg in (30, 60):            # elevation gridlines at 30 and 60 deg
                gy = bby - int(gg / 90.0 * bh)
                _LN(d, bx, gy, bx + bw, gy, _col["dim"])
            _LN(d, bx, bby, bx + bw, bby, _col["dim"])   # horizon baseline
            n = len(_arc)
            px = py = None
            for i in range(n):
                el = _arc[i][2]
                x = bx + int(i * bw / (n - 1))
                yv = bby - int((el if el < 90.0 else 90.0) / 90.0 * bh)
                if px is not None:
                    _LN(d, px, py, x, yv, _col["green"])
                px, py = x, yv
            if sp["aos"] <= f <= sp["los"]:              # 'you are here' on the curve
                nx = bx + int((f - sp["aos"]) / (sp["los"] - sp["aos"]) * bw)
                _LN(d, nx, bty, nx, bby, _col["cyan"])
            _T(d, bx, bby + 2, _hm(sp["aos"], TZ), _col["dim"])
            _T(d, bx + bw - 4 * (_fs.x if _fs else 8), bby + 2, _hm(sp["los"], TZ), _col["dim"])

    # pass list
    _LN(d, 0, yb, 319, yb, _col["dim"])
    _T(d, 4, yb + 3, "NEXT PASSES  L=find C=calc", _col["dim"])
    yy = yb + 3 + ch + 2
    if not _passes:
        _T(d, 6, yy, "none >%d deg in %dh" % (int(MIN_PEAK), HORIZON_H), _col["red"])
    else:
        for i, p in enumerate(_passes[:4]):
            c = _col["amber"] if i == _sel else _col["white"]
            if i == _sel:
                _FR(d, 0, yy - 1, 3, ch, _col["amber"])
            _T(d, 6, yy, "%d %s %2.0fd %s>%s" % (
                i + 1, _datestr(p["aos"], TZ), p["maxel"],
                _comp(p["aosaz"]), _comp(p["losaz"])), c)
            yy += ch + 1
    if _msg:
        _T(d, 4, 320 - ch - 1, _msg[:40], _col["dim"])
    d.swap()
    _last_sec = int(f * 86400)


def _draw_pick(vm):
    d = vm.draw
    ch = _fs.y if _fs else 10
    cw = _fs.x if _fs else 8
    d.clear(Vector(0, 0), d.size, _col["bg"])
    bl = _filter.lower()
    fl = [i for i, s in enumerate(_sats)
          if _filter == "" or _filter in getattr(s, "satnum", "")
          or bl in getattr(s, "name", "").lower()]
    _T(d, 4, 2, "TRACK / SELECT", _col["amber"])
    _T(d, 248, 2, "%d/%d" % (len(fl), len(_sats)), _col["dim"])
    _LN(d, 0, 2 + ch + 1, 319, 2 + ch + 1, _col["dim"])
    _T(d, 4, 4 + ch, "find: " + _filter + "_", _col["green"])
    yy = 4 + ch + ch + 4
    if not fl:
        _T(d, 6, yy, "no match", _col["dim"])
    else:
        sel = _sel if _sel in fl else fl[0]
        for i in fl[:18]:
            s = _sats[i]
            arrow = ">" if i == sel else " "
            cur = "*" if i == _si else " "
            nm = getattr(s, "name", "?")[:13]
            inc = getattr(s, "inclo", 0.0) * r2d
            if i == sel:
                _FR(d, 0, yy - 1, 3, ch, _col["amber"])
            _T(d, 4, yy, "%s%s#%-5s %-13s %2.0f" % (
                arrow, cur, getattr(s, "satnum", "?"), nm, inc),
                _col["amber"] if i == sel else _col["white"])
            yy += ch + 1
    yb = 320 - ch - 2
    _LN(d, 0, yb, 319, yb, _col["dim"])
    _T(d, 4, yb + 2, "type=find ENTER=track ^v BKSP ESC", _col["dim"])
    d.swap()


# ============================ RECOMPUTE ==================================
def _recompute(vm):
    global _passes, _arc, _sel
    d = vm.draw
    ch = _fs.y if _fs else 10
    d.clear(Vector(0, 0), d.size, _col["bg"])
    _T(d, 6, 40, "COMPUTING PASSES...", _col["amber"])
    _T(d, 6, 40 + ch + 2, getattr(_sats[_si], "name", "?"), _col["dim"])
    d.swap()
    f0 = _now_set()                    # latch reference day into _JD0, get start offset
    _passes = _find_passes(_sats[_si], f0, _obs)
    _sel = 0
    _arc = _build_arc(_sats[_si], _passes[0], _obs) if _passes else None
    # cache the result so revisiting this satellite skips the (multi-second) search
    _cache[_si] = {"passes": _passes, "jd0": _JD0, "jf0": f0}


def _cache_fresh(ent):
    """A cached pass list is reusable until it ages past CACHE_H or all its passes
    have already ended (whichever comes first)."""
    nw, nf = _now_jd_parts()
    off = (nw - ent["jd0"]) + nf              # current time as offset vs the cache's _JD0
    if (off - ent["jf0"]) * 24.0 > CACHE_H:
        return False
    ps = ent["passes"]
    if ps and off >= ps[-1]["los"]:
        return False
    return True


def _show_sat(vm):
    """Switch the view to the current satellite, reusing cached passes when fresh."""
    global _JD0, _passes, _arc, _sel
    ent = _cache.get(_si)
    if ent and _cache_fresh(ent):
        _JD0 = ent["jd0"]
        _passes = ent["passes"]
        _sel = 0
        _arc = _build_arc(_sats[_si], _passes[0], _obs) if _passes else None
        return
    _recompute(vm)


# ============================ APP ENTRY POINTS ===========================
def start(view_manager):
    global _sats, _obs, _col, _fs, _state, _si, _sel, _filter, _need_redraw, _msg, _clock_ok
    d = view_manager.draw
    _fs = d.font_size

    def C(r, g, b):
        # Picoware Draw (lcd.LCD) primitives take RGB565 color ints -- convert directly.
        return ((r & 0xF8) << 8) | ((g & 0xFC) << 3) | (b >> 3)
    _col = {"bg": C(0, 0, 0), "dim": C(120, 90, 20), "red": C(255, 60, 50),
            "green": C(0, 255, 120), "amber": C(255, 176, 0), "cyan": C(90, 210, 255),
            "white": C(230, 230, 200), "grid": C(175, 170, 130)}

    _obs = _observer_ecef(LAT, LON, ALT)
    _sats = _load_tles()
    _cache.clear()                     # fresh session -> drop any stale cached passes
    _sync_time(view_manager)           # kick off NTP (async) if the clock isn't set
    _clock_ok = _clock_valid()
    _state = "track"
    _si = 0
    _sel = 0
    _filter = ""
    _msg = None
    _need_redraw = True

    if not _sats:
        d.clear(Vector(0, 0), d.size, _col["bg"])
        _T(d, 6, 30, "NO TLE FILE", _col["red"])
        _T(d, 6, 30 + (_fs.y if _fs else 10) + 2, "put /sd/tles.txt", _col["dim"])
        _T(d, 6, 30 + 2 * ((_fs.y if _fs else 10) + 2), "(3-line TLE sets)", _col["dim"])
        d.swap()
        return True

    _recompute(view_manager)
    _draw_track(view_manager)
    return True


def run(view_manager):
    global _state, _si, _sel, _filter, _need_redraw, _last_sec, _clock_ok
    im = view_manager.input_manager
    b = im.button

    if not _sats:
        if b in (B.BUTTON_BACK, B.BUTTON_ESCAPE):
            im.reset()
            view_manager.back()
        return

    handled = False
    if b != -1:
        handled = True
        im.reset()
        if _state == "track":
            if b in (B.BUTTON_BACK, B.BUTTON_ESCAPE):
                view_manager.back()
                return
            elif b == B.BUTTON_DOWN and _passes:
                _sel = (_sel + 1) % len(_passes)
                globals()["_arc"] = _build_arc(_sats[_si], _passes[_sel], _obs)
            elif b == B.BUTTON_UP and _passes:
                _sel = (_sel - 1) % len(_passes)
                globals()["_arc"] = _build_arc(_sats[_si], _passes[_sel], _obs)
            elif b == B.BUTTON_RIGHT:
                _si = (_si + 1) % len(_sats)
                _show_sat(view_manager)
            elif b == B.BUTTON_LEFT:
                _si = (_si - 1) % len(_sats)
                _show_sat(view_manager)
            elif b == B.BUTTON_L:
                _state = "pick"
                _filter = ""
                _sel = _si
            elif b == B.BUTTON_C:
                _sync_time(view_manager)   # retry NTP (e.g. after connecting WiFi)
                _recompute(view_manager)
        else:  # picker
            bl = _filter.lower()
            fl = [i for i, s in enumerate(_sats)
                  if _filter == "" or _filter in getattr(s, "satnum", "")
                  or bl in getattr(s, "name", "").lower()]
            if b in (B.BUTTON_BACK, B.BUTTON_ESCAPE):
                _state = "track"
            elif b in (B.BUTTON_ENTER, B.BUTTON_CENTER):
                if fl:
                    _si = _sel if _sel in fl else fl[0]
                    _show_sat(view_manager)
                    _state = "track"
            elif b == B.BUTTON_BACKSPACE:
                _filter = _filter[:-1]
            elif b == B.BUTTON_UP and fl:
                p = fl.index(_sel) if _sel in fl else 0
                _sel = fl[max(0, p - 1)]
            elif b == B.BUTTON_DOWN and fl:
                p = fl.index(_sel) if _sel in fl else 0
                _sel = fl[min(len(fl) - 1, p + 1)]
            elif b in _DIGITS:
                _filter += _DIGITS[b]
            elif b in _LETTERS:
                _filter += _LETTERS[b]

    if _state == "pick":
        if handled or _need_redraw:
            _draw_pick(view_manager)
            _need_redraw = False
    else:
        if not _clock_ok and _clock_valid():
            # NTP just landed -> recompute everything against the now-correct clock
            _clock_ok = True
            _cache.clear()
            _recompute(view_manager)
            _need_redraw = True
        sec = int(_now_f() * 86400)
        if handled or sec != _last_sec or _need_redraw:
            _draw_track(view_manager)
            _need_redraw = False


def stop(view_manager):
    global _sats, _passes, _arc, _obs, _need_redraw
    _sats = []
    _passes = []
    _arc = None
    _obs = None
    _cache.clear()
    _need_redraw = True
    gc.collect()
