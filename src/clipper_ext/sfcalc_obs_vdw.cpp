// ChimeraX-Clipper
// Copyright (C) 2016-2019 Tristan Croll, University of Cambridge
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
// Note that this software makes use of modified versions of the Clipper, LibCCP4
// and MMDB libraries, as well as portions of the Intel Math Kernel Library. Each
// of these is redistributed under its own license terms.

#include "sfcalc_obs_vdw.h"
#include "edcalc_ext.h"
#include "scaling.h"

using namespace clipper;
using namespace clipper::datatypes;

namespace clipper_cx
{

template<class dtype> class Compute_add_scaled_fphi
{
public:
    Compute_add_scaled_fphi(const dtype& scale) : scale_(scale) {}
    const F_phi<dtype> operator () (const HKL_info::HKL_reference_index& ih,
            const F_phi<dtype>& fphi1, const F_phi<dtype>& fphi2) const
    {
        clipper::datatypes::F_phi<dtype> fphi;
        if (!fphi1.missing() && !fphi2.missing() )
            fphi = F_phi<dtype>(std::complex<dtype>(fphi1) + scale_*std::complex<dtype>(fphi2));
        return fphi;
    }
private:
    dtype scale_;
};


template<class T>
T optimize_k_sol(HKL_data<datatypes::F_phi<T>>& fphi,
        const HKL_data<datatypes::F_phi<T>>& fphi_atom,
        const HKL_data<datatypes::F_phi<T>>& fphi_mask,
        HKL_data<datatypes::F_phi<T>>& fphi_mask_final,
        const HKL_data<datatypes::F_sigF<T>>& fsig,
        const HKL_info& hkls,
        std::vector<ftype>& params, T k, T dk, T& min_r)
{
    std::cout << "Recalculating bulk solvent B-factor and scale..." << std::endl;
    // try some different scale factors
    //std::vector<double> params( nparams, 1.0 );
    //BasisFn_spline basisfn( hkls, nparams, 1.0 );
    BasisFn_aniso_gaussian basisfn;
    TargetFn_scaleFobsFcalc<T> targetfn(fsig, fphi);
    // TargetFn_scaleF1F2<datatypes::F_phi<T>,datatypes::F_sigF<T> > targetfn( fphi, fsig );
    T x1 = k, dx = dk, x;
    ftype y[3] = { 0.0, 0.0, 0.0 };
    for ( int i = 0; i < 8; i++ ) {
      // take 3 samples of function
      for ( int d = -1; d <= 1; d++ ) if ( y[d+1] == 0.0 ) {
        x = x1 + T(d)*dx;
        x = (x < 0 ? 0.0 : x);
        fphi.compute(fphi_atom, fphi_mask, Compute_add_scaled_fphi<T>(x));
        // for ( ih = fphi.first(); !ih.last(); ih.next() )
        //   fphi[ih] = std::complex<T>(fphi_atom[ih]) +
        //          x * std::complex<T>(fphi_mask[ih]);
        // auto params = guess_initial_aniso_gaussian_params(fsig, fphi);
        ResolutionFn_nonlinear rfn( hkls, basisfn, targetfn, params, 20.0 );
        auto r = quick_r(fphi, fsig, rfn);
        // std::cerr << "Bulk solvent scale: " << x << "R: " << r << std::endl;
        y[d+1] = r;
        if (r < min_r)
        {
            params = rfn.params();
            min_r = r;
        }
        //std::cout << d << "\t" << x << "\t" << r << "\n";
      }
      // find minimum of current 3 samples
      if      ( y[0] < y[1] && y[0] < y[2] ) { y[1] = y[0]; x1 -= dx; }
      else if ( y[2] < y[1] && y[2] < y[0] ) { y[1] = y[2]; x1 += dx; }
      x1 = (x1 < 0 ? 0.0 : x1 );
      // reduce step and search again
      y[0] = y[2] = 0.0;
      dx /= 2.0;
    }
    return x1;

}

template <class T>
T optimize_b_sol(HKL_data<datatypes::F_phi<T>>& fphi,
        const HKL_data<datatypes::F_phi<T>>& fphi_atom,
        const HKL_data<datatypes::F_phi<T>>& fphi_mask,
        HKL_data<datatypes::F_phi<T>>& fphi_mask_final,
        const HKL_data<datatypes::F_sigF<T>>& fsig,
        const HKL_info& hkls,
        std::vector<ftype>& params, T k_sol, T ua1, T dua, T& min_r)
{
    T ua;
    // BasisFn_spline basisfn( hkls, nparams, 1.0 );
    BasisFn_aniso_gaussian basisfn;
    // TargetFn_scaleF1F2<datatypes::F_phi<T>,datatypes::F_sigF<T> > targetfn( fphi, fsig );
    TargetFn_scaleFobsFcalc<T> targetfn(fsig, fphi);
    // HKL_data<F_phi<T> > fphi_mask_final (hkls, cell);
    ftype y[3] = { 0.0, 0.0, 0.0 };
    for (int i=0; i<8; ++i) {
        for (int d=-1; d<=1; ++d ) if (y[d+1] == 0.0 ) {
            HKL_data<data32::F_phi>::HKL_reference_index ih;
            ua = ua1+T(d)*dua;
            fphi_mask_final.compute(fphi_mask, datatypes::Compute_scale_u_iso<datatypes::F_phi<T> >(1.0, -ua));
            fphi.compute(fphi_atom, fphi_mask_final, Compute_add_scaled_fphi<T>(k_sol));
            // for ( HKL_data<data32::F_phi>::HKL_reference_index ih = fphi.first();
            //       !ih.last(); ih.next() )
            //     fphi[ih] = std::complex<T>(fphi_atom[ih]) +
            //         k_sol * std::complex<T>(fphi_mask_final[ih]);
            // auto params = guess_initial_aniso_gaussian_params(fsig, fphi);
            ResolutionFn_nonlinear rfn( hkls, basisfn, targetfn, params, 20.0 );
            auto r = quick_r(fphi, fsig, rfn);
            //std::cerr << "B_add: " << Util::u2b(ua) << " R: " << r << std::endl;
            y[d+1] = r;
            if (r < min_r)
            {
                params = rfn.params();
                min_r = r;
            }
            //std::cout << d << "\t" << x << "\t" << r << "\n";
        }
        // find minimum of current 3 samples
        if      ( y[0] < y[1] && y[0] < y[2] ) { y[1] = y[0]; ua1 -= dua; }
        else if ( y[2] < y[1] && y[2] < y[0] ) { y[1] = y[2]; ua1 += dua; }
        //else params = working_params;
        // reduce step and search again
        y[0] = y[2] = 0.0;
        dua /= 2.0;

    }
    return ua1;

}

template<class T>
bool SFcalc_obs_bulk_vdw<T>::operator() ( HKL_data<datatypes::F_phi<T> >& fphi,
    const HKL_data<datatypes::F_sigF<T> >& fsig, const Atom_list& atoms)
{
  // std::cout << "Starting bulk solvent calculation..." << std::endl << std::flush;
  // set U value constants
  //double u_atom = Util::b2u( 0 ); // 20.0 );
  double u_mask = Util::b2u( 50.0 );

  // increase the U values
  // Atom_list atomu = atoms;
  // U_aniso_orth uadd( u_atom ), u;
  // for ( int i = 0; i < atomu.size(); i++ ) if ( !atomu[i].is_null() ) {
  //   u = atomu[i].u_aniso_orth();
  //   if ( u.is_null() ) u = U_aniso_orth( atomu[i].u_iso() );
  //   atomu[i].set_u_aniso_orth( u + uadd );
  // }

  // now make the map for ed calcs
  const HKL_info&   hkls = fsig.base_hkl_info();
  const Spacegroup& spgr = hkls.spacegroup();
  const Cell&       cell = fsig.base_cell();
  HKL_data<datatypes::F_phi<T> >
    fphi_atom( hkls, cell ), fphi_mask( hkls, cell );
  const Grid_sampling grid( spgr, cell, hkls.resolution() );
  Xmap<float> xmap( spgr, cell, grid );

  // do ed calc from atomu
  // EDcalc_aniso<ftype32> edcalc;
  // auto start = std::chrono::steady_clock  ::now();
  EDcalc_aniso_thread<ftype32> edcalc(nthreads);
  // std::cout << "Number locked before: " << xmap.count_locked() << std::endl;
  edcalc( xmap, atoms);
  // std::cout << "Number locked after: " << xmap.count_locked() << std::endl;
  // if (!xmap.all_unlocked())
    // std::cerr << "ERROR: not all grid points are unlocked!" << std::endl;
  // auto end = std::chrono::steady_clock  ::now();

  // std::chrono::duration<double> elapsed = end-start;
  // std::cout << "EDcalc with " << nthreads << " threads took " << elapsed.count() << " seconds." << std::endl;
  // start = std::chrono::steady_clock  ::now();
  xmap.fft_to( fphi_atom, nthreads );
  // end = std::chrono::steady_clock  ::now();
  // elapsed = end-start;
  // std::cout << "Single x-to-h FFT took " << elapsed.count() << " seconds." << std::endl;

  // fphi_atom.compute( fphi_atom, datatypes::Compute_scale_u_iso<datatypes::F_phi<T> >( 1.0, u_atom ) );

  // do density calc from mask

  // start = std::chrono::steady_clock  ::now();
  auto emcalc = EDcalc_mask_vdw<ftype32>();
  emcalc.set_num_threads(nthreads);
  emcalc( xmap, atoms );
  for ( Xmap<ftype32>::Map_reference_index ix = xmap.first();
        !ix.last(); ix.next() )
    xmap[ix] = 1.0 - xmap[ix];
  xmap.fft_to( fphi_mask, nthreads );

  // end = std::chrono::steady_clock  ::now();
  // elapsed = end-start;
  // std::cout << "Single-threaded bulk solvent mask + FFT took " << elapsed.count() << " seconds." << std::endl;


  HKL_data<F_phi<T> > fphi_mask_final (hkls, cell);

  // set (0,0,0) terms to null
  fphi_mask.compute( fphi_mask, datatypes::Compute_scale_u_iso<datatypes::F_phi<T> >( 1.0, -u_mask ) );
  // std::cerr << "First HKL: " << fphi_atom.first().hkl().format() << std::endl;
  // fphi_atom[fphi_atom.first()].set_null();
  // fphi_mask[fphi_mask.first()].set_null();
  auto c000 = HKL_data<data32::F_phi>::HKL_reference_coord(hkls, HKL(0,0,0)).index();
  fphi_atom[c000].set_null();
  fphi_mask[c000].set_null();

  if (bulk_solvent_optimization_needed_)
  {
      T x1 = 0.35;
      T min_r;
      fphi.compute(fphi_atom, fphi_mask, Compute_add_scaled_fphi<T>(0.5));
      *params_ = guess_initial_aniso_gaussian_params(fsig, fphi, min_r);
      // std::cerr << "Initial anisotropic scaling params: ";
      // for (auto p: *params_)
      //   std::cerr << p << ",";
      // std::cerr << std::endl;
      x1 = optimize_k_sol<T>(fphi, fphi_atom, fphi_mask, fphi_mask_final, fsig, hkls, *params_, x1, x1, min_r);
      auto ua1 = optimize_b_sol<T>(fphi, fphi_atom, fphi_mask, fphi_mask_final, fsig, hkls, *params_, x1, Util::b2u(0), Util::b2u(50), min_r);

      // adopt final scale and B-factor
      bulk_u = ua1;
      bulkscl = x1;
      std::cout << "Final solvent B-factor: " << Util::u2b(u_mask + ua1) << " scale: " << x1 << std::endl;
      fphi_mask_final.compute(fphi_mask, datatypes::Compute_scale_u_iso<datatypes::F_phi<T> >(1.0, -ua1));
      bulk_solvent_optimization_needed_ = false;
  } else {
      // Just use the stored values
      fphi_mask_final.compute(fphi_mask, datatypes::Compute_scale_u_iso<datatypes::F_phi<T>>(1.0, -bulk_u));
  }
  for ( HKL_data<data32::F_phi>::HKL_reference_index ih = fphi.first();
        !ih.last(); ih.next() )
    fphi[ih] = std::complex<T>(fphi_atom[ih]) +
          bulkscl * std::complex<T>(fphi_mask_final[ih]);
  // end = std::chrono::steady_clock  ::now();
  // elapsed = end-start;
  // std::cout << "Complete single-threaded bulk solvent calc took " << elapsed.count() << " seconds." << std::endl;

  // store stats
  ftype64 w, s0 = 0.0, s1 = 0.0;
  for ( Xmap<ftype32>::Map_reference_index ix = xmap.first();
        !ix.last(); ix.next() ) {
    w = 1.0/ftype64( xmap.multiplicity( ix ) );
    s0 += w;
    s1 += w*xmap[ix];
  }
  bulkfrc = s1/s0;
  // bulkscl = x1;

  return true;
}

// compile templates

template class CLIPPER_CX_IMEX SFcalc_obs_bulk_vdw<ftype32>;

template class CLIPPER_CX_IMEX SFcalc_obs_bulk_vdw<ftype64>;




} // namespace clipper_cx
