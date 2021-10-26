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

*/

static const float LMS_from_linearRGB[] = {
    0.17886, 0.43997, 0.03597,
    0.03380, 0.27515, 0.03621,
    0.00031, 0.00192, 0.01528
};

static const float linearRGB_from_LMS[] = {
    8.00533, -12.88195, 11.68065,
    -0.97821, 5.26945, -10.18300,
    -0.04017, -0.39885, 66.48079
};

/*
    The projections are always parallel to one of the LMS axes. So the
    projection matrix only modifies one coordinate, depending on the deficiency
    type. So here we only store the matrix row that actually does something,
    along with the coordinate it applies to.

    Also to check on which plane the LMS point projects, we store the normal of
    the separation plane and use a dot product to get the side instead of
    comparing the ratios like the original paper. It's equivalent.
*/
struct DLBrettel1997Params
{
    int lmsElementToProject;
    float projectionOnPlane1[3];
    float projectionOnPlane2[3];
    float separationPlaneNormal[3];
};

static struct DLBrettel1997Params brettel_protan_params = {
    0, // only this LMS coordinate is affected for protan
    { 0.00000, 2.18394, -5.65554 }, // Projection to plane 1
    { 0.00000, 2.16614, -5.30455 }, // Projection to plane 2
    { 0.00000, 0.01751, -0.34516 }  // Normal of the separation plane to pick the projection plane.
};

static struct DLBrettel1997Params brettel_deutan_params = {
    1, // only this LMS coordinate is affected for deutan
    { 0.46165, 0.00000, 2.44885 }, // Projection to plane 1
    { 0.45789, 0.00000, 2.58960 }, // Projection to plane 2
    { -0.01751, 0.00000, 0.65480 }  // Normal of the separation plane to pick the projection plane.
};

static struct DLBrettel1997Params brettel_tritan_params = {
    2, // only this LMS coordinate is affected for tritan
    { -0.00213, 0.05477, 0.00000 }, // Projection to plane 1
    { -0.06195, 0.16826, 0.00000 }, // Projection to plane 2
    { 0.34516, -0.65480, 0.00000 }  // Normal of the separation plane to pick the projection plane.
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
            float rgb[3] = {
                linearRGB_from_sRGB(rowPtr[col + 0]),
                linearRGB_from_sRGB(rowPtr[col + 1]),
                linearRGB_from_sRGB(rowPtr[col + 2])
            };
            
            // lms = dl_LMS_from_linearRGB * rgb
            float lms[3] = {
                LMS_from_linearRGB[0]*rgb[0] + LMS_from_linearRGB[1]*rgb[1] + LMS_from_linearRGB[2]*rgb[2],
                LMS_from_linearRGB[3]*rgb[0] + LMS_from_linearRGB[4]*rgb[1] + LMS_from_linearRGB[5]*rgb[2],
                LMS_from_linearRGB[6]*rgb[0] + LMS_from_linearRGB[7]*rgb[1] + LMS_from_linearRGB[8]*rgb[2]
            };

            // Check on which plane we should project by comparing wih the separation plane normal.
            float dotWithSepPlane = lms[0]*params->separationPlaneNormal[0] + lms[1]*params->separationPlaneNormal[1] + lms[2]*params->separationPlaneNormal[2];
            const float* projectionOnPlane = (dotWithSepPlane >= 0 ? params->projectionOnPlane1 : params->projectionOnPlane2);

            // Project on the plane. Only one coordinate changes (the axis
            // corresponding to the missing cone cells), so no need to perform a
            // full 3x3 multiplication.
            // The severity factor is implemented as a linear interpolation with the original value.
            float projected_element = projectionOnPlane[0]*lms[0] + projectionOnPlane[1]*lms[1] + projectionOnPlane[2]*lms[2];
            lms[params->lmsElementToProject] = (projected_element*severity) + (lms[params->lmsElementToProject]*(1.f-severity));

            // Go back to linear RGB
            // rgb_cvd = dl_linearRGB_from_LMS * lms
            float rgb_cvd[3] = {
                linearRGB_from_LMS[0]*lms[0] + linearRGB_from_LMS[1]*lms[1] + linearRGB_from_LMS[2]*lms[2],
                linearRGB_from_LMS[3]*lms[0] + linearRGB_from_LMS[4]*lms[1] + linearRGB_from_LMS[5]*lms[2],
                linearRGB_from_LMS[6]*lms[0] + linearRGB_from_LMS[7]*lms[1] + linearRGB_from_LMS[8]*lms[2]
            };

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

float dl_vienot_protan_rgbCvd_from_rgb[] = {
    0.10889, 0.89111, -0.00000,
    0.10889, 0.89111, 0.00000,
    0.00447, -0.00447, 1.00000
};

float dl_vienot_deutan_rgbCvd_from_rgb[] = {
    0.29031, 0.70969, -0.00000,
    0.29031, 0.70969, -0.00000,
    -0.02197, 0.02197, 1.00000
};

// WARNING: Viénot 1999 is not accurate for tritanopia. Use Brettel 1997 instead.
float dl_vienot_tritan_rgbCvd_from_rgb[] = {
    1.00000, 0.15236, -0.15236,
    0.00000, 0.86717, 0.13283,
    -0.00000, 0.86717, 0.13283
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
