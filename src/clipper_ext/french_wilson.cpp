/**
 * @Author: Tristan Croll <tic20>
 * @Date:   30-May-2019
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 30-May-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2017-2018 Tristan Croll
 */
#include <vector>

namespace clipper_cx { namespace french_wilson {

// static const double acentric_zj[71] = {
//      0.226, 0.230, 0.235, 0.240, 0.246, 0.251, 0.257, 0.263, 0.270,  0.276,
//      0.283, 0.290, 0.298, 0.306, 0.314, 0.323, 0.332, 0.341, 0.351,  0.362,
//      0.373, 0.385, 0.397, 0.410, 0.424, 0.439, 0.454, 0.470, 0.487,  0.505,
//      0.525, 0.545, 0.567, 0.590, 0.615, 0.641, 0.668, 0.698, 0.729,  0.762,
//      0.798, 0.835, 0.875, 0.917, 0.962, 1.009, 1.059, 1.112, 1.167,  1.226,
//      1.287, 1.352, 1.419, 1.490, 1.563, 1.639, 1.717, 1.798, 1.882,  1.967,
//      2.055, 2.145, 2.236, 2.329, 2.422, 2.518, 2.614, 2.710, 2.808,  2.906,
//      3.004
//  };
//
// static const double acentric_sig_zj [71] = {
//     0.217, 0.221, 0.226, 0.230, 0.235, 0.240, 0.245, 0.250, 0.255,  0.261,
//     0.267, 0.273, 0.279, 0.286, 0.292, 0.299, 0.307, 0.314, 0.322,  0.330,
//     0.339, 0.348, 0.357, 0.367, 0.377, 0.387, 0.398, 0.409, 0.421,  0.433,
//     0.446, 0.459, 0.473, 0.488, 0.503, 0.518, 0.535, 0.551, 0.568,  0.586,
//     0.604, 0.622, 0.641, 0.660, 0.679, 0.698, 0.718, 0.737, 0.757,  0.776,
//     0.795, 0.813, 0.831, 0.848, 0.865, 0.881, 0.895, 0.909, 0.921,  0.933,
//     0.943, 0.953, 0.961, 0.968, 0.974, 0.980, 0.984, 0.988, 0.991,  0.994,
//     0.996
// };

static const double acentric_zf[71] = {
    0.423, 0.428, 0.432, 0.437, 0.442, 0.447, 0.453, 0.458, 0.464,  0.469,
    0.475, 0.482, 0.488, 0.495, 0.502, 0.509, 0.516, 0.524, 0.532,  0.540,
    0.549, 0.557, 0.567, 0.576, 0.586, 0.597, 0.608, 0.619, 0.631,  0.643,
    0.656, 0.670, 0.684, 0.699, 0.714, 0.730, 0.747, 0.765, 0.783,  0.802,
    0.822, 0.843, 0.865, 0.887, 0.911, 0.935, 0.960, 0.987, 1.014,  1.042,
    1.070, 1.100, 1.130, 1.161, 1.192, 1.224, 1.257, 1.289, 1.322,  1.355,
    1.388, 1.421, 1.454, 1.487, 1.519, 1.551, 1.583, 1.615, 1.646,  1.676,
    1.706
};

static const double acentric_sig_zf[71] = {
    0.216, 0.218, 0.220, 0.222, 0.224, 0.226, 0.229, 0.231, 0.234,  0.236,
    0.239, 0.241, 0.244, 0.247, 0.250, 0.253, 0.256, 0.259, 0.262,  0.266,
    0.269, 0.272, 0.276, 0.279, 0.283, 0.287, 0.291, 0.295, 0.298,  0.302,
    0.307, 0.311, 0.315, 0.319, 0.324, 0.328, 0.332, 0.337, 0.341,  0.345,
    0.349, 0.353, 0.357, 0.360, 0.364, 0.367, 0.369, 0.372, 0.374,  0.375,
    0.376, 0.377, 0.377, 0.377, 0.376, 0.374, 0.372, 0.369, 0.366,  0.362,
    0.358, 0.353, 0.348, 0.343, 0.338, 0.332, 0.327, 0.321, 0.315,  0.310,
    0.304
};

// static const double centric_zj[71] = {
//     0.114, 0.116, 0.119, 0.122, 0.124, 0.127, 0.130, 0.134, 0.137,  0.141,
//     0.145, 0.148, 0.153, 0.157, 0.162, 0.166, 0.172, 0.177, 0.183,  0.189,
//     0.195, 0.202, 0.209, 0.217, 0.225, 0.234, 0.243, 0.253, 0.263,  0.275,
//     0.287, 0.300, 0.314, 0.329, 0.345, 0.363, 0.382, 0.402, 0.425,  0.449,
//     0.475, 0.503, 0.534, 0.567, 0.603, 0.642, 0.684, 0.730, 0.779,  0.833,
//     0.890, 0.952, 1.018, 1.089, 1.164, 1.244, 1.327, 1.416, 1.508,  1.603,
//     1.703, 1.805, 1.909, 2.015, 2.123, 2.233, 2.343, 2.453, 2.564,  2.674,
//     2.784, 2.894, 3.003, 3.112, 3.220, 3.328, 3.435, 3.541, 3.647,  3.753,
//     3.962
// };
//
// static const double centric_sig_zj[71] = {
//     0.158, 0.161, 0.165, 0.168, 0.172, 0.176, 0.179, 0.184, 0.188,  0.192,
//     0.197, 0.202, 0.207, 0.212, 0.218, 0.224, 0.230, 0.236, 0.243,  0.250,
//     0.257, 0.265, 0.273, 0.282, 0.291, 0.300, 0.310, 0.321, 0.332,  0.343,
//     0.355, 0.368, 0.382, 0.397, 0.412, 0.428, 0.445, 0.463, 0.481,  0.501,
//     0.521, 0.543, 0.565, 0.589, 0.613, 0.638, 0.664, 0.691, 0.718,  0.745,
//     0.773, 0.801, 0.828, 0.855, 0.881, 0.906, 0.929, 0.951, 0.971,  0.989,
//     1.004, 1.018, 1.029, 1.038, 1.044, 1.049, 1.052, 1.054, 1.054,  1.053,
//     1.051, 1.049, 1.047, 1.044, 1.041, 1.039, 1.036, 1.034, 1.031,  1.029,
//     1.028
// };

static const double centric_zf[81] = {
    0.269, 0.272, 0.276, 0.279, 0.282, 0.286, 0.289, 0.293, 0.297,  0.301,
    0.305, 0.309, 0.314, 0.318, 0.323, 0.328, 0.333, 0.339, 0.344,  0.350,
    0.356, 0.363, 0.370, 0.377, 0.384, 0.392, 0.400, 0.409, 0.418,  0.427,
    0.438, 0.448, 0.460, 0.471, 0.484, 0.498, 0.512, 0.527, 0.543,  0.560,
    0.578, 0.597, 0.618, 0.639, 0.662, 0.687, 0.713, 0.740, 0.769,  0.800,
    0.832, 0.866, 0.901, 0.938, 0.976, 1.016, 1.057, 1.098, 1.140,  1.183,
    1.227, 1.270, 1.313, 1.356, 1.398, 1.439, 1.480, 1.519, 1.558,  1.595,
    1.632, 1.667, 1.701, 1.735, 1.767, 1.799, 1.829, 1.859, 1.889,  1.917,
    1.945
};

static const double centric_sig_zf[81] = {
    0.203, 0.205, 0.207, 0.209, 0.211, 0.214, 0.216, 0.219, 0.222,  0.224,
    0.227, 0.230, 0.233, 0.236, 0.239, 0.243, 0.246, 0.250, 0.253,  0.257,
    0.261, 0.265, 0.269, 0.273, 0.278, 0.283, 0.288, 0.293, 0.298,  0.303,
    0.309, 0.314, 0.320, 0.327, 0.333, 0.340, 0.346, 0.353, 0.361,  0.368,
    0.375, 0.383, 0.390, 0.398, 0.405, 0.413, 0.420, 0.427, 0.433,  0.440,
    0.445, 0.450, 0.454, 0.457, 0.459, 0.460, 0.460, 0.458, 0.455,  0.451,
    0.445, 0.438, 0.431, 0.422, 0.412, 0.402, 0.392, 0.381, 0.370,  0.360,
    0.349, 0.339, 0.330, 0.321, 0.312, 0.304, 0.297, 0.290, 0.284,  0.278,
    0.272
};

}}
