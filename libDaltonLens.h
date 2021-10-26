/* libDaltonLens - public domain library - http://daltonlens.org
                                  no warranty implied; use at your own risk

    Author: Nicolas Burrus <nicolas@burrus.name>

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

#include <stddef.h>

enum DLDeficiency
{
    DLDeficiency_Protan,
    DLDeficiency_Deutan,
    DLDeficiency_Tritan
};

/*
    Automatically picks the best CVD simulation algorithm (Brettel 1997 for
    tritanopia, Viénot 1999 for protanopia and deuteranopia).

    For a comparison of the available algorithms see:
    https://daltonlens.org/opensource-cvd-simulation/

    For more information about the math of the chosen algorithms see
    https://daltonlens.org/understanding-cvd-simulation/

    'srgba_image' is expected to be in the RGBA32 format, 8 bits per channel,
    encoded as sRGB.

    'severity' should be between 0.0 and 1.0 and is implemented via a linear
    interpolation with the original image.

    'bytesPerRow' can be 0, in which case it'll be automatically computed as
    width*4.
*/
void dl_simulate_cvd (enum DLDeficiency deficiency, float severity, unsigned char *srgba_image, size_t width, size_t height, size_t bytesPerRow);

/*
    Implements the algorithm proposed in 1997 by 
    Brettel, H. and Viénot, F. and Mollon, J. D.

    'Computerized simulation of color appearance for dichromats'
    Journal of the Optical Society of America. A, Optics, Image Science, and Vision    

    It works well for all kinds of dichromacies, but is a bit more expensive
    than Viénot 1999.

    This version is adapted to modern sRGB monitors.
*/
void dl_simulate_cvd_brettel1997 (enum DLDeficiency deficiency, float severity, unsigned char* srgba_image, size_t width, size_t height, size_t bytesPerRow);

/*
    Implements the algorithm proposed in 1999 by
    Viénot, F. and Brettel, H. and Mollon, J. D.

    Digital video colourmaps for checking the legibility of displays by dichromats
	Color Research & Application

    It works well for protanopia and deuteranopia, but NOT for tritanopia.

    This version is adapted to modern sRGB monitors.
*/
void dl_simulate_cvd_vienot1999 (enum DLDeficiency deficiency, float severity, unsigned char* srgba_image, size_t width, size_t height, size_t bytesPerRow);
