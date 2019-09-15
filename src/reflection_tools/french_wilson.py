# Clipper plugin to UCSF ChimeraX
# Copyright (C) 2016-2019 Tristan Croll, University of Cambridge

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

# Note that this software makes use of modified versions of the Clipper, LibCCP4
# and MMDB libraries, as well as portions of the Intel Math Kernel Library. Each
# of these is redistributed under its own license terms.

# Conversion of intensities to amplitudes, based on
# Read & McCoy (2016) Acty Cryst D72(3): 375-387

_MIN_REFLECTIONS_PER_BIN = 40

def french_wilson_analytical(i_sigi, max_bins = 60):
    '''
    Perform French & Wilson scaling (conversion of intensities to amplitudes)
    using the analytical approach described in
        Read & McCoy (2016), Acta Cryst D72(3): 375-387
    Args:
        i_sigi: either a HKL_data_I_sigI or HKL_data_I_sigI_ano (anomalous data
            will be automatically merged)
        max_bins: maximum number of bins (resolution shells)for calculation of
        mean intensity

    Returns:
        A HKL_data_F_sigF
    '''
    from chimerax.clipper import (
        HKL_data_I_sigI, HKL_data_I_sigI_ano, HKL_data_F_sigF, F_sigF, I_sigI
    )
    from chimerax.clipper.clipper_python import (
        TargetFn_meanInth_I_sigI_float,
        TargetFn_meanInth_I_sigI_ano_float,
        TargetFn_scaleI1I2_I_sigI_float,
        BasisFn_spline, BasisFn_binner,
        ResolutionFn
    )
    from math import sqrt

    hkls = i_sigi.base_hkl_info

    if isinstance(i_sigi, HKL_data_I_sigI_ano):
        i_sigi_merged = HKL_data_I_sigI(hkls)
        ih = i_sigi.first_data
        while not ih.last():
            isigi_data = i_sigi[ih]
            i_sigi_merged[ih.hkl] = I_sigI(isigi_data.i, isigi_data.sigi)
            i_sigi.next_data(ih)
        i_sigi = i_sigi_merged

    hkls = i_sigi.base_hkl_info
    n_reflections = hkls.num_reflections
    if n_reflections < _MIN_REFLECTIONS_PER_BIN * 3:
        err_str = ('Intensity dataset contains only {} reflections, which is '
            'insufficient to provide reliable statistics for French & Wilson '
            'scaling. Are you sure this is a macromolecular dataset?')
        err_str.format(n_reflections)
        raise TypeError(err_str)
    reflections_per_bin = max(n_reflections//max_bins, _MIN_REFLECTIONS_PER_BIN)
    n_bins = n_reflections // reflections_per_bin

    import numpy
    params = [1.0]*n_bins #    numpy.ones(n_bins, numpy.float32)
    #basisfn = BasisFn_spline(hkls, n_bins, 1.0)
    basisfn = BasisFn_binner(i_sigi, n_bins, 1.0)
    target = TargetFn_meanInth_I_sigI_float(i_sigi, 1.0)
    rfn = ResolutionFn(hkls, basisfn, target, params)

    scaled_is = HKL_data_I_sigI(hkls)
    s_isigi_data = I_sigI()
    ih = i_sigi.first
    count = 0
    while not ih.last():
        isigi_data = i_sigi[ih]
        #fsigf_data = f_sigf[ih]
        if isigi_data.missing:
            s_isigi_data.set_null()
            scaled_is[ih.hkl] = s_isigi_data
            ih.next()
            continue
        imean = rfn.f(ih)
        eps = ih.hkl_class.epsilon
        eobs_sq = isigi_data.i / eps / imean
        sig_eobs_sq = isigi_data.sigi / eps / imean
        try:
            if ih.hkl_class.centric:
                e_xpct = _expected_E_FW_cen(eobs_sq, sig_eobs_sq)
                esq_xpct = _expected_Esq_FW_cen(eobs_sq, sig_eobs_sq)
            else:
                e_xpct = _expected_E_FW_acen(eobs_sq, sig_eobs_sq)
                esq_xpct = _expected_Esq_FW_acen(eobs_sq, sig_eobs_sq)
        except ValueError:
            print ('Failed on {}'.format(str(ih.hkl)))
            print('N_bins: {} Imean: {} I: {} sigI: {} eobs_sq: {} sig_eobs_sq: {}'.format(
                n_bins, imean, isigi_data.i, isigi_data.sigi, eobs_sq, sig_eobs_sq
            ))
            raise
        sig_e = sqrt(esq_xpct - e_xpct**2)

        sqrt_eps = sqrt(eps)
        s_isigi_data.i = e_xpct**2 * eps
        s_isigi_data.sigi = sig_e * eps * e_xpct
        scaled_is[ih.hkl] = s_isigi_data

        ih.next()
        count+=1

    scale_fn = TargetFn_scaleI1I2_I_sigI_float(scaled_is, i_sigi)
    f_sigf = HKL_data_F_sigF(hkls)
    fsigf_data = F_sigF()
    rfn = ResolutionFn(hkls, basisfn, scale_fn, params)
    ih = scaled_is.first_data
    while not ih.last():
        s_isigi_data = scaled_is[ih]
        s_isigi_data.scale(sqrt(rfn.f(ih)))
        f = sqrt(s_isigi_data.i)
        fsigf_data.f = f
        fsigf_data.sigf = s_isigi_data.sigi / f
        f_sigf[ih.hkl] = fsigf_data
        scaled_is.next_data(ih)


    return f_sigf







def _expected_E_FW_acen(eobs_sq, sig_eobs_sq):
    ''' First moment of the expected E distribution i.e. <E>'''
    from math import sqrt, pi, exp
    CROSSOVER1 = -12.5
    CROSSOVER2 = 18
    x = (eobs_sq - sig_eobs_sq**2)/sig_eobs_sq
    xsqr = x**2

    if x < CROSSOVER1: # Large negative argument: asymptotic approximation
        e_xpct = sqrt(-pi*sig_eobs_sq/x) * (
            (-916620705 + xsqr * (
                91891800 + xsqr * (
                    -11531520 + xsqr * (
                        1935360 + xsqr * (
                            -491520 + xsqr * 262144
                        )
                    )
                )
            )) /
            (-495452160 + xsqr * (
                55050240 + xsqr * (
                    -7864320 + xsqr * (
                        1572864 + xsqr * (
                            -524288 + xsqr * 524288
                        )
                    )
                )
            ))
        )
    elif x > CROSSOVER2: # Large positive argument: asymptotic approximation
        e_xpct = sqrt(sig_eobs_sq) * (
            (-45045 + 32*xsqr * (
                -315 + 8*xsqr * (
                    -15 - 16*xsqr + 128 * xsqr**2
                )
            )) /
            (32768 * x**7.5)
        )
    else:
        from scipy.special import pbdv, erfc # parabolic cylinder function D, complementary error function
        e_xpct = (sqrt(sig_eobs_sq/2) * exp(-xsqr/4) *
            pbdv(-1.5, -x)[0] / # Scipy parabolic cylinder function returns value *and* derivative
            erfc(-x/sqrt(2))
        )
    return e_xpct


def _expected_Esq_FW_acen(eobs_sq, sig_eobs_sq):
    ''' Second moment of the expected E distribution, i.e. <E**2> '''
    from math import sqrt, pi, exp
    CROSSOVER1 = -8.9
    CROSSOVER2 = 5.7
    esq_xpct_base = eobs_sq - sig_eobs_sq**2
    x = esq_xpct_base/(sqrt(2)*sig_eobs_sq)
    xsqr = x**2

    if x < CROSSOVER1:
        esq_xpct = esq_xpct_base * ((
            -135135 + xsqr * (
                20790 + xsqr * (
                    -3780 + xsqr * (
                        840 + xsqr * (
                            -240 + xsqr * (
                                96 - xsqr * 64
                            )
                        )
                    )
                )
            )
        ) / (
            -135135 + xsqr * (
                20790 + xsqr * (
                    -3780 + xsqr * (
                        840 + xsqr * (
                            -240 + xsqr * (
                                96 + xsqr * (
                                    -64     + xsqr * 128
                                )
                            )
                        )
                    )
                )
            )
        )
        )
    elif x > CROSSOVER2:
        esq_xpct = esq_xpct_base
    else:
        from scipy.special import erfc
        esq_xpct = esq_xpct_base + sqrt(2/pi) * sig_eobs_sq / (exp(xsqr) * erfc(-x))
    return esq_xpct

def _expected_E_FW_cen(eobs_sq, sig_eobs_sq):
    ''' First moment of the expected E distribution i.e. <E>'''
    CROSSOVER1 = -17.5
    CROSSOVER2 = 17.5
    from math import pi, sqrt
    x = sig_eobs_sq/2 - eobs_sq/sig_eobs_sq
    xsqr = x**2
    if x < CROSSOVER1:
        pcd_ratio = ((1024 * sqrt(pi)*(-x)**6.5) /
            (3465 + xsqr * (
                840 + xsqr * (
                    384 + xsqr * 1024
                )
            ))
        )
    elif x > CROSSOVER2:
        pcd_ratio = (
        (3440640 + xsqr *(
            -491520 + xsqr *(
                98304 + xsqr *(
                    -32768 + xsqr * 32768
                )
            )
        )) /
        (675675 + xsqr *(
            -110880 + xsqr *(
                26880 + xsqr * (
                    -12288 + xsqr * 32768
                )
            )
        ))
        ) / sqrt(x)
    else:
        from scipy.special import pbdv
        pcd_ratio = pbdv(-1,x)[0] / pbdv(-0.5,x)[0]
    e_xpct = sqrt(sig_eobs_sq /pi)*pcd_ratio
    return e_xpct

def _expected_Esq_FW_cen(eobs_sq, sig_eobs_sq):
    ''' Second moment of the expected E distribution, i.e. <E**2> '''
    CROSSOVER1 = -17.5
    CROSSOVER2 = 17.5
    x = sig_eobs_sq/2 - eobs_sq/sig_eobs_sq
    xsqr = x**2
    if x < CROSSOVER1:
        pcd_ratio = (
        ( 45045 + xsqr * (
            10080 + xsqr * (
                3840 + xsqr * (
                    4096 - xsqr * 32768
                )
            )
        )) /
        ( x * (
            55440 + xsqr * (
                13440 + xsqr * (
                    6144 + xsqr * 16384
                )
            )
        ))
        )
    elif x > CROSSOVER2:
        pcd_ratio = (
        ( 11486475 + xsqr * (
            -1441440 + xsqr * (
                241920 + xsqr * (
                    -61440 + xsqr * 32768
                )
            )
        )) /
        ( x * (
            675675 + xsqr * (
                -110880 + xsqr * (
                    26880 + xsqr * (
                        -12288 + xsqr * 32768
                    )
                )
            )
        ))
        )
    else:
        from scipy.special import pbdv
        pcd_ratio = pbdv(-1.5,x)[0] / pbdv(-0.5,x)[0]
    esq_xpct = sig_eobs_sq*pcd_ratio/2
    return esq_xpct
