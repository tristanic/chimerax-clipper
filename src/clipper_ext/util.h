/**
 * @Author: Tristan Croll <tic20>
 * @Date:   25-Nov-2020
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 25-Nov-2020
 * @License: Lesser GNU Public License version 3.0 (see LICENSE.md)
 * @Copyright: 2016-2019 Tristan Croll
 */
#pragma once
#include <algorithm>
#include <clipper/clipper.h>
#include <clipper/core/resol_basisfn.h>

using namespace clipper;
using namespace clipper::datatypes;

namespace clipper_cx
{


template<template<class> class dtype, typename T>
HKL_info select_random_reflections(const HKL_data<dtype<T>>& data, size_t num_reflections)
{
    std::vector<HKL> valid_hkls;
    for (auto ih=data.first_data(); !ih.last(); ih=data.next_data(ih))
    {
        valid_hkls.push_back(ih.hkl());
    }
    std::random_shuffle(valid_hkls.begin(), valid_hkls.end());
    if (valid_hkls.size() > num_reflections)
    {
        valid_hkls.erase(valid_hkls.begin()+num_reflections, valid_hkls.end());
    }
    auto base_hkl = data.base_hkl_info();
    HKL_info selected_hkls = HKL_info(base_hkl.spacegroup(), base_hkl.cell(), base_hkl.resolution());
    selected_hkls.add_hkl_list(valid_hkls);
    return selected_hkls;
}

template<template<class> class dtype, typename T>
HKL_info select_random_reflections_in_bins(const HKL_data<dtype<T>>& data, size_t reflections_per_bin, size_t num_bins)
{
    if (data.hkl_info().num_reflections() <= reflections_per_bin*num_bins) return data.hkl_info();
    std::map< int, std::vector<HKL> > bins;
    for (size_t i=0; i<num_bins; ++i)
    {
        bins[i] = std::vector<HKL>();
    }
    Resolution_ordinal s_ord;
    s_ord.init(data.hkl_info(), 1.0);
    for (auto ih = data.first_data(); !ih.last(); ih=data.next_data(ih))
    {
        int bin = Util::bound( 0, Util::intf( ftype(num_bins) * s_ord.ordinal( ih.invresolsq() ) ), (int)num_bins-1 );
        bins[bin].push_back(ih.hkl());
    }
    std::vector<HKL> selected_hkls_;
    for (auto& pair: bins)
    {
        auto& hkls = pair.second;
        auto count = std::min(reflections_per_bin, hkls.size());
        std::random_shuffle(hkls.begin(), hkls.end());
        for (size_t i=0; i<count; ++i)
        {
            selected_hkls_.push_back(hkls[i]);
        }
    }
    auto base_hkl = data.base_hkl_info();
    HKL_info selected_hkls = HKL_info(base_hkl.spacegroup(), base_hkl.cell(), base_hkl.resolution());
    selected_hkls.add_hkl_list(selected_hkls_);
    return selected_hkls;
    


}


template<class DType>
HKL_data<DType> reflection_subset(const HKL_data<DType>& data_in,
    const HKL_info& selected_hkls)
{
    HKL_data<DType> data_out = HKL_data<DType>(selected_hkls, selected_hkls.cell());
    for (auto ih = selected_hkls.first(); !ih.last(); ih.next() )
    {
        data_out[ih] = data_in[ih.hkl()];
    }
    return data_out;
}


}
