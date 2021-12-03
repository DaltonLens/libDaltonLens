/* libDaltonLens - public domain library - http://daltonlens.org
                                  no warranty implied; use at your own risk

    Author: Nicolas Burrus <nicolas@burrus.name>

    LICENSE information at the end of the file.
*/

/*
    For more information about the math of these two algorithms
    see https://daltonlens.org/understanding-cvd-simulation/

    The notebook used to dump the precomputed matrices can be found
    in the DaltonLens-Python project:
    https://github.com/DaltonLens/DaltonLens-Python/blob/master/research/for-desktop/precomputed_matrices.ipynb
*/

#include "libDaltonLens.h"

#include <math.h>

/*
    From https://en.wikipedia.org/wiki/SRGB

    This is the standard used by most images and displays. This is the textbook
    implementation that uses powf. 

    Optimizations are beyond the scope of this library, but it is possible to
    make it faster by using a lookup-table (only 256 values), but then it's less
    SIMD-friendly. It's also unclear if it'll be much faster on recent CPUs as
    they are often memory-bound.

    Another option is to use a fast pow implementation, e.g. with a Chebychev
    approximation:
    https://stackoverflow.com/questions/6475373/optimizations-for-pow-with-const-non-integer-exponent
*/
static inline float linearRGB_from_sRGB(unsigned char v)
{
    float fv = v / 255.f;
    if (fv < 0.04045f) return fv / 12.92f;
    return powf((fv + 0.055f) / 1.055f, 2.4f);
}

static inline unsigned char sRGB_from_linearRGB(float v)
{
    if (v <= 0.f) return 0;
    if (v >= 1.f) return 255;
    if (v < 0.0031308f) return 0.5f + (v * 12.92 * 255.f);
    return 0.f + 255.f * (powf(v, 1.f / 2.4f) * 1.055f - 0.055f);
}

/*
    Brettel 1997 precomputed parameters.

    LMS model
    =========

    These values use the sRGB standard to go from linearRGB to CIE XYZ, and the
    Smith & Pokorny 1975 model for CIE XYZ to LMS. This is the LMS model used by
    Viénot, Brettel and Mollon, but upgraded to use the sRGB standard used by
    modern monitors.

    Projection Planes
    =================

    These were computed using RGB white as the neutral element, not the
    equal-energy E. This option is commonly chosen by the Brettel
    implementations (including Vischeck), and it increases the range of colors
    that project within the sRGB gamut and avoids clipping too much.

    DaltonLens-Python Code to regenerate
    ====================================

    simulator = simulate.Simulator_Brettel1997(convert.LMSModel_sRGB_SmithPokorny75())
    simulator.dumpPrecomputedValues = True
    simulator.simulate_cvd(np.zeros((1,1,3), dtype=np.uint8), simulate.Deficiency.PROTAN, severity=1.0)
    simulator.simulate_cvd(np.zeros((1,1,3), dtype=np.uint8), simulate.Deficiency.DEUTAN, severity=1.0)
    simulator.simulate_cvd(np.zeros((1,1,3), dtype=np.uint8), simulate.Deficiency.TRITAN, severity=1.0)

    Alternative code to get the same output as Vischeck (as implemented in GIMP display filters):
    simulator = simulate.Simulator_Vischeck()
    simulator.dumpPrecomputedValues = True
    simulator.simulate_cvd(np.zeros((1,1,3), dtype=np.uint8), simulate.Deficiency.PROTAN, severity=1.0)
    simulator.simulate_cvd(np.zeros((1,1,3), dtype=np.uint8), simulate.Deficiency.DEUTAN, severity=1.0)
    simulator.simulate_cvd(np.zeros((1,1,3), dtype=np.uint8), simulate.Deficiency.TRITAN, severity=1.0)


    The implementation got simplified to minimize compute (with the
    exact same output), refer to 
    https://daltonlens.org/cvd-simulation-svg-filters/ for more details.
*/

/*
    These parameters combine the full projection pipeline for each half-plane so
    we don't need an explicit transform to the LMS space.

    To check on which plane the LMS point should project we also don't need to
    actually compute the LMS coordinates, and can do it directly in the RGB space.
*/
struct DLBrettel1997Params
{
    // Transformation using plane 1 == rgbFromLms . projection1 . lmsFromRgb
    float rgbCvdFromRgb_1[9];
    
    // Full transformation using plane 2 == rgbFromLms . projection2 . lmsFromRgb
    float rgbCvdFromRgb_2[9];

    // Normal of the separation plane to pick the right transform, already in the RGB space.
    // == normalInLms . lmsFromRgb
    float separationPlaneNormalInRgb[3];
};

static struct DLBrettel1997Params brettel_protan_params = {
    {
        0.14980, 1.19548, -0.34528,
        0.10764, 0.84864, 0.04372,
        0.00384, -0.00540, 1.00156,
    },
    {
        0.14570, 1.16172, -0.30742,
        0.10816, 0.85291, 0.03892,
        0.00386, -0.00524, 1.00139,
    },
    { 0.00048, 0.00393, -0.00441 }
};
static struct DLBrettel1997Params brettel_deutan_params = {
    {
        0.36477, 0.86381, -0.22858,
        0.26294, 0.64245, 0.09462,
        -0.02006, 0.02728, 0.99278,
    },
    {
        0.37298, 0.88166, -0.25464,
        0.25954, 0.63506, 0.10540,
        -0.01980, 0.02784, 0.99196,
    },
    { -0.00281, -0.00611, 0.00892 }
};
static struct DLBrettel1997Params brettel_tritan_params = {
    {
        1.01277, 0.13548, -0.14826,
        -0.01243, 0.86812, 0.14431,
        0.07589, 0.80500, 0.11911,
    },
    {
        0.93678, 0.18979, -0.12657,
        0.06154, 0.81526, 0.12320,
        -0.37562, 1.12767, 0.24796,
    },
    { 0.03901, -0.02788, -0.01113 }
};

void dl_simulate_cvd_brettel1997 (enum DLDeficiency deficiency, float severity, unsigned char* srgba_image, size_t width, size_t height, size_t bytesPerRow)
{
    // Compute a default bytesPerRow if it wasn't specified.
    if (bytesPerRow == 0)
    {
        bytesPerRow = width * 4;
    }

    const struct DLBrettel1997Params* params = NULL;
    switch (deficiency)
    {
        case DLDeficiency_Protan: params = &brettel_protan_params; break;
        case DLDeficiency_Deutan: params = &brettel_deutan_params; break;
        case DLDeficiency_Tritan: params = &brettel_tritan_params; break;
    }

    for (size_t row = 0; row < height; ++row)
    {
        unsigned char* rowPtr = srgba_image + bytesPerRow*row;
        for (size_t col = 0; col < width*4; col += 4)
        {
            // rgb = linearRGB_from_sRGB(srgb)
            // alpha is untouched.
            const float rgb[3] = {
                linearRGB_from_sRGB(rowPtr[col + 0]),
                linearRGB_from_sRGB(rowPtr[col + 1]),
                linearRGB_from_sRGB(rowPtr[col + 2])
            };
            
            // Check on which plane we should project by comparing wih the separation plane normal.
            const float* n = params->separationPlaneNormalInRgb;
            const float dotWithSepPlane = rgb[0]*n[0] + rgb[1]*n[1] + rgb[2]*n[2];
            const float* rgbCvdFromRgb = (dotWithSepPlane >= 0 ? params->rgbCvdFromRgb_1 : params->rgbCvdFromRgb_2);

            float rgb_cvd[3] = {
                rgbCvdFromRgb[0]*rgb[0] + rgbCvdFromRgb[1]*rgb[1] + rgbCvdFromRgb[2]*rgb[2],
                rgbCvdFromRgb[3]*rgb[0] + rgbCvdFromRgb[4]*rgb[1] + rgbCvdFromRgb[5]*rgb[2],
                rgbCvdFromRgb[6]*rgb[0] + rgbCvdFromRgb[7]*rgb[1] + rgbCvdFromRgb[8]*rgb[2]
            };

            // Apply the severity factor as a linear interpolation.
            // It's the same to do it in the RGB space or in the LMS
            // space since it's a linear transform.
            rgb_cvd[0] = rgb_cvd[0]*severity + rgb[0]*(1.f-severity);
            rgb_cvd[1] = rgb_cvd[1]*severity + rgb[1]*(1.f-severity);
            rgb_cvd[2] = rgb_cvd[2]*severity + rgb[2]*(1.f-severity);

            // Encode as sRGB and write the result.
            rowPtr[col + 0] = sRGB_from_linearRGB(rgb_cvd[0]);
            rowPtr[col + 1] = sRGB_from_linearRGB(rgb_cvd[1]);
            rowPtr[col + 2] = sRGB_from_linearRGB(rgb_cvd[2]);
        }
    }
}

/*
    Viénot 1999 precomputed parameters.

    This follows the paper exactly, but using the modern sRGB standard to decode
    the input RGB values.

    Since there is only one projection plane, the entire pipeline can be reduced
    to a single 3x3 matrix multiplication in the linearRGB space.

    DaltonLens-Python Code to regenerate
    ====================================

    simulator =
    simulate.Simulator_Vienot1999(convert.LMSModel_sRGB_SmithPokorny75())
    simulator.dumpPrecomputedValues = True
    simulator.simulate_cvd(np.zeros((1,1,3), dtype=np.uint8), simulate.Deficiency.PROTAN, severity=1.0)
    simulator.simulate_cvd(np.zeros((1,1,3), dtype=np.uint8), simulate.Deficiency.DEUTAN, severity=1.0)
    simulator.simulate_cvd(np.zeros((1,1,3), dtype=np.uint8), simulate.Deficiency.TRITAN, severity=1.0)
*/

static float dl_vienot_protan_rgbCvd_from_rgb[] = {
    0.11238, 0.88762, 0.00000,
    0.11238, 0.88762, -0.00000,
    0.00401, -0.00401, 1.00000
};

static float dl_vienot_deutan_rgbCvd_from_rgb[] = {
    0.29275, 0.70725, 0.00000,
    0.29275, 0.70725, -0.00000,
    -0.02234, 0.02234, 1.00000
};

// WARNING: Viénot 1999 is not accurate for tritanopia. Use Brettel 1997 instead.
static float dl_vienot_tritan_rgbCvd_from_rgb[] = {
    1.00000, 0.14461, -0.14461,
    0.00000, 0.85924, 0.14076,
    -0.00000, 0.85924, 0.14076
};

void dl_simulate_cvd_vienot1999 (enum DLDeficiency deficiency, float severity, unsigned char* srgba_image, size_t width, size_t height, size_t bytesPerRow)
{
    // Compute a default bytesPerRow if it wasn't specified.
    if (bytesPerRow == 0)
    {
        bytesPerRow = width * 4;
    }

    const float* rgbCvd_from_rgb = NULL;
    switch (deficiency)
    {
        case DLDeficiency_Protan: rgbCvd_from_rgb = dl_vienot_protan_rgbCvd_from_rgb; break;
        case DLDeficiency_Deutan: rgbCvd_from_rgb = dl_vienot_deutan_rgbCvd_from_rgb; break;
        case DLDeficiency_Tritan: rgbCvd_from_rgb = dl_vienot_tritan_rgbCvd_from_rgb; break;
    }

    for (size_t row = 0; row < height; ++row)
    {
        unsigned char* rowPtr = srgba_image + bytesPerRow*row;
        for (size_t col = 0; col < width*4; col += 4)
        {
            // rgb = linearRGB_from_sRGB(srgb)
            // alpha is untouched.
            float rgb[3] = {
                linearRGB_from_sRGB(rowPtr[col + 0]),
                linearRGB_from_sRGB(rowPtr[col + 1]),
                linearRGB_from_sRGB(rowPtr[col + 2])
            };
            
            // rgb_cvd = rgbCvd_from_rgb * rgb
            float rgb_cvd[3] = {
                rgbCvd_from_rgb[0]*rgb[0] + rgbCvd_from_rgb[1]*rgb[1] + rgbCvd_from_rgb[2]*rgb[2],
                rgbCvd_from_rgb[3]*rgb[0] + rgbCvd_from_rgb[4]*rgb[1] + rgbCvd_from_rgb[5]*rgb[2],
                rgbCvd_from_rgb[6]*rgb[0] + rgbCvd_from_rgb[7]*rgb[1] + rgbCvd_from_rgb[8]*rgb[2]
            };

            // Implement the severity factor as a linear interpolation.
            if (severity < 0.999f)
            {
                rgb_cvd[0] = severity*rgb_cvd[0] + (1.f - severity)*rgb[0];
                rgb_cvd[1] = severity*rgb_cvd[1] + (1.f - severity)*rgb[1];
                rgb_cvd[2] = severity*rgb_cvd[2] + (1.f - severity)*rgb[2];
            }

            // Write the result, encoded to sRGB
            rowPtr[col + 0] = sRGB_from_linearRGB(rgb_cvd[0]);
            rowPtr[col + 1] = sRGB_from_linearRGB(rgb_cvd[1]);
            rowPtr[col + 2] = sRGB_from_linearRGB(rgb_cvd[2]);
        }
    }
}

void dl_simulate_cvd (enum DLDeficiency deficiency, float severity, unsigned char* srgba_image, size_t width, size_t height, size_t bytesPerRow)
{
    // Viénot 1999 is not accurate for tritanopia, so use Brettel in that case.
    // Otherwise use Viénot 1999 because it's a bit faster.
    if (deficiency == DLDeficiency_Tritan)
    {
        dl_simulate_cvd_brettel1997(deficiency, severity, srgba_image, width, height, bytesPerRow);
    }
    else
    {
        dl_simulate_cvd_vienot1999(deficiency, severity, srgba_image, width, height, bytesPerRow);
    }
}

/*
    LICENSE

    This is free and unencumbered software released into the public domain.

    Anyone is free to copy, modify, publish, use, compile, sell, or
    distribute this software, either in source code form or as a compiled
    binary, for any purpose, commercial or non-commercial, and by any
    means.

    In jurisdictions that recognize copyright laws, the author or authors
    of this software dedicate any and all copyright interest in the
    software to the public domain. We make this dedication for the benefit
    of the public at large and to the detriment of our heirs and
    successors. We intend this dedication to be an overt act of
    relinquishment in perpetuity of all present and future rights to this
    software under copyright law.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
    IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
    OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
    ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
    OTHER DEALINGS IN THE SOFTWARE.

    For more information, please refer to <https://unlicense.org>
*/
