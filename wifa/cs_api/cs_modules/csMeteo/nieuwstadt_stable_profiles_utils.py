#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 10:13:07 2024

@author: j26483
"""

import numpy as np


def phim_hogstrom_s(z, dlmo):
    zeta = z * dlmo
    if zeta < 0.5:
        return 1 + 4.8 * zeta
    elif zeta < 10:
        return 7.9 - 4.25 / zeta + (1 / zeta) ** 2
    else:
        return 0.7485 * zeta


def phih_hogstrom_s(z, dlmo):
    zeta = z * dlmo
    return 0.95 + 7.8 * zeta


def phim_u_nieuwstadt(z, zi, dlmo, alpha):
    a = 0.5 * alpha - 1
    b = 0.5 * pow(alpha * (alpha + 2), 0.5)
    coef = (
        (phim_hogstrom_s(z, dlmo) / z)
        * pow(max(1 - z / zi, 0), a)
        * np.cos(b * np.log(max(1 - z / zi, 1e-12)))
    )
    return coef


def phim_v_nieuwstadt(z, zi, dlmo, alpha):
    a = 0.5 * alpha - 1
    b = 0.5 * pow(alpha * (alpha + 2), 0.5)
    coef = (
        (phim_hogstrom_s(z, dlmo) / z)
        * pow(max(1 - z / zi, 0), a)
        * np.sin(b * np.log(max(1 - z / zi, 1e-12)))
    )
    return coef


def phih_nieuwstadt(z, zi, dlmo, alpha):
    a = alpha - 2
    coef = (phih_hogstrom_s(z, dlmo) / z) * pow(max(1 - z / zi, 0), a)
    return coef


def psim_u_nieuwstadt(z, zi, dlmo, z0, alpha):
    psi = 0
    dz0 = z0
    k = 0
    z_id = z0
    z_idp1 = z0 + dz0
    while z_idp1 < z:
        dz = max(dz0 * pow(1.01, k), dz0)
        derivee = phim_u_nieuwstadt(0.5 * (z_id + z_idp1), zi, dlmo, alpha)
        psi = psi + derivee * dz
        k = k + 1
        z_id = z_idp1
        z_idp1 = z_idp1 + dz
    derivee = phim_u_nieuwstadt(0.5 * (z_id + z), zi, dlmo, alpha)
    psi = psi + derivee * max(z - z_id, 0)
    return psi


def psim_v_nieuwstadt(z, zi, dlmo, z0, alpha):
    psi = 0
    dz0 = z0
    k = 0
    z_id = z0
    z_idp1 = z0 + dz0
    while z_idp1 < z:
        dz = max(dz0 * pow(1.01, k), dz0)
        derivee = phim_v_nieuwstadt(0.5 * (z_id + z_idp1), zi, dlmo, alpha)
        psi = psi + derivee * dz
        k = k + 1
        z_id = z_idp1
        z_idp1 = z_idp1 + dz
    derivee = phim_v_nieuwstadt(0.5 * (z_id + z), zi, dlmo, alpha)
    psi = psi + derivee * max(z - z_id, 0)
    return psi


def psih_nieuwstadt(z, zi, dlmo, z0, alpha):
    psi = 0
    dz0 = z0
    k = 0
    z_id = z0
    z_idp1 = z0 + dz0
    while z_idp1 < z:
        dz = max(dz0 * pow(1.01, k), dz0)
        derivee = phih_nieuwstadt(0.5 * (z_id + z_idp1), zi, dlmo, alpha)
        psi = psi + derivee * dz
        k = k + 1
        z_id = z_idp1
        z_idp1 = z_idp1 + dz
    derivee = phih_nieuwstadt(0.5 * (z_id + z), zi, dlmo, alpha)
    psi = psi + derivee * max(z - z_id, 0)
    return psi


def tke_nieuwstadt(z, zi, dlmo, alpha):
    d = 9.7
    c = 6
    c2 = 0.1
    c3 = -0.8
    ctheta = 1.4
    w0 = 1 / 3 - 2 / (3 * c) + 2 / c + 4 * c3 / (3 * c)
    w1 = 2 * c2 / (3 * c) - 2 / c - 4 * c3 / (3 * c)
    zeta = dlmo * z
    rich = phih_hogstrom_s(z, dlmo) / pow(phim_hogstrom_s(z, dlmo), 2) * zeta
    a1 = 0.5 + 1.5 * pow(rich, 2) - pow(rich, 3)
    tke = (
        0.5
        * np.sqrt(
            d
            * (phim_hogstrom_s(z, dlmo) - zeta)
            / (
                phih_hogstrom_s(z, dlmo)
                * (
                    w0
                    + w1 * phim_hogstrom_s(z, dlmo) / (phim_hogstrom_s(z, dlmo) - zeta)
                    - (1 - a1) * zeta / (ctheta * (phim_hogstrom_s(z, dlmo) - zeta))
                )
            )
        )
        * pow(1 - z / zi, alpha / 2 + 1)
    )
    return tke


def epsilon_nieuwstadt(z, zi, dlmo, alpha):
    zeta = dlmo * z
    epsilon = (phim_hogstrom_s(z, dlmo) - zeta) * pow(1 - z / zi, alpha) / z
    return epsilon
