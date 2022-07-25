import numpy as np
import string

# Copyright (C) 2022  Yingbo Li
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# REFERENCE FOR MATHS: A Guide to Coordinate Systems in Great Britain
# V3.6 © OS 2020 (original Copyright Ordnance Survey 2018)
# https://www.ordnancesurvey.co.uk/documents/resources/guide-coordinate-systems-great-britain.pdf

def letter_grid2numbers(letters):
    ''' convert the first two letters of an OS Grid Ref into numbers '''
    alphabet = list(string.ascii_uppercase)
    del(alphabet[8])
    a,b = list(map(lambda letter : alphabet.index(letter) + 1,letters))
    E_500,N_500 = (a % 5) - 1, (25 - a)//5
    if (a%5) == 0:
        # Strictly speaking, valid OS Grid Refs don't require this as none of the 500's are in the right hand column
        E_500 = 5 - 1
    E_100,N_100 =  (b % 5) - 1, (25 - b)//5
    if (b%5) == 0:
        E_100 = 5 - 1

    # ensure coords are relative to false origin
    transform_E_500 = -2
    transform_N_500 = -1
    return (E_500 + transform_E_500)*5 + E_100, (N_500 + transform_N_500)*5 + N_100

def grid_ref2osgb(grid_ref):
    ''' convert grid reference to (E,N) tuple.
    Note: decimals in the grid reference are not currently supported
    '''
    letters = grid_ref[:2]
    try:
        E_pre = str(int(letters[0]))
        N_pre = str(int(letters[1]))
    except ValueError:
        E_pre,N_pre = letter_grid2numbers(letters)

    numbers = grid_ref[2:].strip()
    if len(numbers) == 0:
        raw_E,raw_N = "",""
    else:
        raw_E,raw_N = numbers.split(" ")
        if len(raw_E) > 5 or len(raw_N) > 5:
            raise Exception("Each number in the numerical part of the grid reference must have 5 or less digits, with no decimal digits.")

    E_zeros_append = "0" * (int((5 - len(raw_E))))
    N_zeros_append = "0" * (int((5 - len(raw_N))))

    E = str(E_pre) + raw_E + E_zeros_append
    N = str(N_pre) + raw_N + N_zeros_append
    return int(E),int(N)

# Maximum number of iterations for any steps requiring it (mainly any type of Cartesian -> LatLon inc. EN -> LatLon and XYZ -> LatLon)
max_iters = 10 # this shouldn't need to be set too high

deg2rad = lambda deg : (deg/360) * 2 * np.pi
rad2deg = lambda rad : (rad/(2*np.pi)) * 360
sec2rad = lambda sec : deg2rad(sec * 1/3600)

# Airy 1830 for OSGB36, National Grid
osgb36_data = {"a":6377563.396,
               "b":6356256.909,
               "F_0":0.9996012717}

# WGS84
wgs84_data = {"a":6378137.000,
              "b":6356752.3141}

def osgbEN2osgbLatLon(E,N,H=0):
    '''
    Convert Northing (N), Easting (E) and optional altitude (H) on an OSGB36 ellipsoid to Latitude (osgb_phi) and Longitude (osgb_lam) on an OSGB36 ellipsoid
    '''

    a = osgb36_data["a"]
    b = osgb36_data["b"]
    F_0 = osgb36_data["F_0"]

    e_squared_osgb = (a**2 - b**2)/(a**2)
    E_0 = 400000
    N_0 = -100000

    # positive phi: N, positive lambda: E
    phi_0 = deg2rad(49)
    lambda_0 = deg2rad(-2)

    n = (a-b)/(a+b)
    M = lambda phi : b*F_0*(
        (1 + n + (5/4)*n**2 + (5/4)*n**3) * (phi - phi_0)
        - (3*n + 3*n**2 + (21/8)*n**3) * np.sin(phi - phi_0) * np.cos(phi + phi_0)
        + ((15/8)*n**2 + (15/8)*n**3) * np.sin(2*(phi - phi_0)) * np.cos(2*(phi + phi_0))
        - (35/24)*n**3 * np.sin(3*(phi - phi_0)) * np.cos(3*(phi + phi_0))
    )

    phi_prime = ((N-N_0)/(a*F_0)) + phi_0

    test_val = N - N_0 - M(phi_prime)
    c = 0
    while (abs(test_val) > 1e-5) and (c < max_iters):
        phi_prime = ((test_val)/(a*F_0)) + phi_prime
        test_val = N - N_0 - M(phi_prime)
        c += 1

    nu = a*F_0 * (1 - e_squared_osgb*np.sin(phi_prime)**2)**-0.5
    rho = a*F_0 * (1 - e_squared_osgb) * (1 - e_squared_osgb*np.sin(phi_prime)**2)**-1.5
    eta_squared = (nu/rho) - 1

    sec = lambda theta : 1/np.cos(theta)
    VII = np.tan(phi_prime)/(2*rho*nu)
    VIII = (np.tan(phi_prime)/(24*rho*nu**3)) * (5 + 3 * np.tan(phi_prime)**2 + eta_squared - 9 * (np.tan(phi_prime)**2) * eta_squared)
    IX = (np.tan(phi_prime)/(720*rho*nu**5)) * (61 + 90 * np.tan(phi_prime)**2 + 45 * np.tan(phi_prime)**4)
    X = sec(phi_prime)/nu
    XI = (sec(phi_prime)/(6*nu**3)) * ((nu/rho) + 2 * np.tan(phi_prime)**2)
    XII = (sec(phi_prime)/(120*nu**5)) * (5 + 28 * np.tan(phi_prime)**2 + 24 * np.tan(phi_prime)**4)
    XIIA = (sec(phi_prime)/(5040*nu**7)) * (61 + 662 * np.tan(phi_prime)**2 + 1320 * np.tan(phi_prime)**4 + 720 * np.tan(phi_prime)**6)

    phi_osgb = phi_prime - VII*(E-E_0)**2 + VIII*(E-E_0)**4 - IX*(E-E_0)**6
    lam_osgb = lambda_0 + X*(E-E_0) - XI*(E-E_0)**3 + XII*(E-E_0)**5 - XIIA*(E-E_0)**7
    return phi_osgb,lam_osgb

def osgbLatLon2osgbCartesian(phi_osgb,lam_osgb,H=0):
    '''
    Convert Latitude (phi_osgb) and Longitude (lam_osgb) on an OSGB36 ellipsoid to a vector of Cartesian coordinates [x,y,z] (v_osgb) on an OSGB36 ellipsoid
    '''

    a = osgb36_data["a"]
    b = osgb36_data["b"]
    F_0 = osgb36_data["F_0"]

    e_squared_osgb = (a**2 - b**2)/(a**2)
    nu_osgb = a * (1 - e_squared_osgb*np.sin(phi_osgb)**2)**-0.5

    x_osgb = (nu_osgb + H) * np.cos(phi_osgb) * np.cos(lam_osgb)
    y_osgb = (nu_osgb + H) * np.cos(phi_osgb) * np.sin(lam_osgb)
    z_osgb = ((1-e_squared_osgb)*nu_osgb + H) * np.sin(phi_osgb)

    v_osgb = np.array([x_osgb,y_osgb,z_osgb])
    return v_osgb

def osgbCartesian2wgsCartesian(v_osgb):
    '''
    Convert a vector of Cartesian coordinates [x,y,z] (v_osgb) on an OSGB36 ellipsoid to a vector of Cartesian coordinates [x,y,z] (v_wgs) on a WGS84 ellipsoid
    '''

    # transformations
    t_x = -446.448
    t_y = 125.157
    t_z = -542.060

    # scale
    s = 20.4894 * 1e-6

    # rotations
    r_x = sec2rad(-0.1502)
    r_y = sec2rad(-0.2470)
    r_z = sec2rad(-0.8421)

    helmert_matrix = np.array([[1 + s,-r_z,r_y],
                               [r_z,1 + s,-r_x],
                               [-r_y,r_x,1 + s]])

    v_t = np.array([t_x,t_y,t_z])

    a = wgs84_data["a"]
    b = wgs84_data["b"]

    v_wgs = np.linalg.solve(helmert_matrix,v_osgb - v_t)
    return v_wgs

def wgsCartesian2wgsLatLon(v_wgs):
    '''
    Convert Cartesian coordinates [x,y,z] (v_wgs) on a WGS84 ellipsoid to Latitude (wgs_phi) and Longitude (wgs_lam) on a WGS84 ellipsoid
    '''

    a = wgs84_data["a"]
    b = wgs84_data["b"]

    x_wgs = v_wgs[0]
    y_wgs = v_wgs[1]
    z_wgs = v_wgs[2]

    wgs_e_squared = (a**2 - b**2)/(a**2)

    lam_wgs = np.arctan2(y_wgs,x_wgs)
    p = (x_wgs**2+y_wgs**2)**0.5
    phi_wgs = np.arctan2(z_wgs,(p*(1-wgs_e_squared)))
    wgs_nu = a * (1 - wgs_e_squared*np.sin(phi_wgs)**2)**-0.5
    reached_precision = False
    c = 0
    while (not reached_precision) and (c < max_iters):
        new_phi_wgs = np.arctan2((z_wgs + wgs_e_squared*wgs_nu*np.sin(phi_wgs)),p)
        wgs_nu = a * (1 - wgs_e_squared*np.sin(phi_wgs)**2)**-0.5
        if abs(phi_wgs - new_phi_wgs) < 1e-30:
            reached_precision = True
        phi_wgs = new_phi_wgs
        c += 1

    H = (p/np.cos(phi_wgs)) - wgs_nu

    return phi_wgs,lam_wgs,H

def osgbEN2wgsLatLon(E,N,H=0):
    ''' Convert osgb36 (E,N,H) to wgs84 (lat,lon,H) with lat and lon in degrees (not radians)
            E: Easting
            N: Northing
            H: Altitude AOD
    '''
    phi_osgb,lam_osgb = osgbEN2osgbLatLon(E,N,H)
    v_osgb = osgbLatLon2osgbCartesian(phi_osgb,lam_osgb,H)
    v_wgs = osgbCartesian2wgsCartesian(v_osgb)
    phi_wgs,lam_wgs,H = wgsCartesian2wgsLatLon(v_wgs)
    lat = rad2deg(phi_wgs)
    lon = rad2deg(lam_wgs)
    return lat,lon,H
