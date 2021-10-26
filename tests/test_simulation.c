#include <libDaltonLens.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define SOKOL_TIME_IMPL
#include "sokol_time.h"

#include <stdio.h>

/*

    The "ground truth" images come from DaltonLens-Python, which is itself
    tested against external references like Vischeck.

*/

int compare_images (const char* gt_name, unsigned char* inputImage, int width, int height)
{
    char fullGroundTruthPath[1024];
    snprintf(fullGroundTruthPath, 1024, "%s%s", TEST_IMAGES_DIR, gt_name);

    int gt_width = 0, gt_height = 0, gt_channels = 0;
    unsigned char* gtImage = stbi_load (fullGroundTruthPath, &gt_width, &gt_height, &gt_channels, 4 /* rgba */);

    if (gt_width != width || gt_height != height)
    {
        fprintf (stderr, "FAIL: (%s) size does not match\n", gt_name);
        return 1;
    }

    for (int r = 0; r < height; ++r)
    {
        const unsigned char* inputRowPtr = inputImage + width*r;
        const unsigned char* gtRowPtr = gtImage + width*r;
        for (int c = 0; c < width*4; ++c)
        {
            int diff = abs(inputRowPtr[c] - gtRowPtr[c]);
            if (diff > 1)
            {
                fprintf (stderr, "FAIL: (%s) pixel differs at (%d,%d)[%d] diff=%d\n", gt_name, c/4, r, c%4, diff);
                char outputPath[1024];
                snprintf(outputPath, 1024, "output_%s", gt_name);
                stbi_write_png (outputPath, width, height, 4, inputImage, 0 /* default stride */);
                return 2;
            }
        }
    }

    fprintf (stderr, "GOOD: (%s)\n", gt_name);
    return 0;
}

struct Context
{
    unsigned char* inputBuffer;
    unsigned char* tmpImageBuffer;
    int width;
    int height;
};

int test_Vienot1999 (struct Context* context)
{
    int cumulatedComparison = 0;

    int w = context->width;
    int h = context->height;    

    {
        memcpy(context->tmpImageBuffer, context->inputBuffer, w * h * 4);
        uint64_t timeStart = stm_now();
        dl_simulate_cvd_vienot1999(DLDeficiency_Protan, 1.0, context->tmpImageBuffer, w, h, 0 /* default */);
        fprintf (stderr, "TIMING dl_simulate_cvd_vienot1999 = %.1f ms\n", stm_ms(stm_since(timeStart)));
        cumulatedComparison += compare_images("vienot1999_protan_1.0.png", context->tmpImageBuffer, w, h);
    }

    {
        memcpy(context->tmpImageBuffer, context->inputBuffer, w * h * 4);
        dl_simulate_cvd_vienot1999(DLDeficiency_Deutan, 1.0, context->tmpImageBuffer, w, h, 0 /* default */);
        cumulatedComparison += compare_images("vienot1999_deutan_1.0.png", context->tmpImageBuffer, w, h);
    }

    {
        memcpy(context->tmpImageBuffer, context->inputBuffer, w * h * 4);
        dl_simulate_cvd_vienot1999(DLDeficiency_Tritan, 1.0, context->tmpImageBuffer, w, h, 0 /* default */);
        cumulatedComparison += compare_images("vienot1999_tritan_1.0.png", context->tmpImageBuffer, w, h);
    }

    {
        memcpy(context->tmpImageBuffer, context->inputBuffer, w * h * 4);
        dl_simulate_cvd_vienot1999(DLDeficiency_Protan, 0.55, context->tmpImageBuffer, w, h, 0 /* default */);
        cumulatedComparison += compare_images("vienot1999_protan_0.55.png", context->tmpImageBuffer, w, h);
    }

    {
        memcpy(context->tmpImageBuffer, context->inputBuffer, w * h * 4);
        dl_simulate_cvd_vienot1999(DLDeficiency_Deutan, 0.55, context->tmpImageBuffer, w, h, 0 /* default */);
        cumulatedComparison += compare_images("vienot1999_deutan_0.55.png", context->tmpImageBuffer, w, h);
    }

    {
        memcpy(context->tmpImageBuffer, context->inputBuffer, w * h * 4);
        dl_simulate_cvd_vienot1999(DLDeficiency_Tritan, 0.55, context->tmpImageBuffer, w, h, 0 /* default */);
        cumulatedComparison += compare_images("vienot1999_tritan_0.55.png", context->tmpImageBuffer, w, h);
    }

    return cumulatedComparison;
}

int test_Brettel1997 (struct Context* context)
{
    int cumulatedComparison = 0;

    int w = context->width;
    int h = context->height;    

    {
        memcpy(context->tmpImageBuffer, context->inputBuffer, w * h * 4);
        uint64_t timeStart = stm_now();
        dl_simulate_cvd_brettel1997(DLDeficiency_Protan, 1.0, context->tmpImageBuffer, w, h, 0 /* default */);
        fprintf (stderr, "TIMING dl_simulate_cvd_brettel1997 = %.1f ms\n", stm_ms(stm_since(timeStart)));
        cumulatedComparison += compare_images("brettel1997_protan_wn_1.0.png", context->tmpImageBuffer, w, h);
    }

    {
        memcpy(context->tmpImageBuffer, context->inputBuffer, w * h * 4);
        dl_simulate_cvd_brettel1997(DLDeficiency_Deutan, 1.0, context->tmpImageBuffer, w, h, 0 /* default */);
        cumulatedComparison += compare_images("brettel1997_deutan_wn_1.0.png", context->tmpImageBuffer, w, h);
    }

    {
        memcpy(context->tmpImageBuffer, context->inputBuffer, w * h * 4);
        dl_simulate_cvd_brettel1997(DLDeficiency_Tritan, 1.0, context->tmpImageBuffer, w, h, 0 /* default */);
        cumulatedComparison += compare_images("brettel1997_tritan_wn_1.0.png", context->tmpImageBuffer, w, h);
    }

    {
        memcpy(context->tmpImageBuffer, context->inputBuffer, w * h * 4);
        dl_simulate_cvd_brettel1997(DLDeficiency_Protan, 0.55, context->tmpImageBuffer, w, h, 0 /* default */);
        cumulatedComparison += compare_images("brettel1997_protan_wn_0.55.png", context->tmpImageBuffer, w, h);
    }

    {
        memcpy(context->tmpImageBuffer, context->inputBuffer, w * h * 4);
        dl_simulate_cvd_brettel1997(DLDeficiency_Deutan, 0.55, context->tmpImageBuffer, w, h, 0 /* default */);
        cumulatedComparison += compare_images("brettel1997_deutan_wn_0.55.png", context->tmpImageBuffer, w, h);
    }

    {
        memcpy(context->tmpImageBuffer, context->inputBuffer, w * h * 4);
        dl_simulate_cvd_brettel1997(DLDeficiency_Tritan, 0.55, context->tmpImageBuffer, w, h, 0 /* default */);
        cumulatedComparison += compare_images("brettel1997_tritan_wn_0.55.png", context->tmpImageBuffer, w, h);
    }

    return cumulatedComparison;
}

int test_automaticDispatch (struct Context* context)
{
    int cumulatedComparison = 0;

    int w = context->width;
    int h = context->height;    

    {
        memcpy(context->tmpImageBuffer, context->inputBuffer, w * h * 4);
        dl_simulate_cvd(DLDeficiency_Protan, 1.0, context->tmpImageBuffer, w, h, 0 /* default */);
        cumulatedComparison += compare_images("vienot1999_protan_1.0.png", context->tmpImageBuffer, w, h);
    }

    {
        memcpy(context->tmpImageBuffer, context->inputBuffer, w * h * 4);
        dl_simulate_cvd(DLDeficiency_Deutan, 1.0, context->tmpImageBuffer, w, h, 0 /* default */);
        cumulatedComparison += compare_images("vienot1999_deutan_1.0.png", context->tmpImageBuffer, w, h);
    }

    {
        memcpy(context->tmpImageBuffer, context->inputBuffer, w * h * 4);
        dl_simulate_cvd(DLDeficiency_Tritan, 1.0, context->tmpImageBuffer, w, h, 0 /* default */);
        cumulatedComparison += compare_images("brettel1997_tritan_wn_1.0.png", context->tmpImageBuffer, w, h);
    }

    return cumulatedComparison;
}

int main ()
{
    stm_setup();

    char inputImagePath[1024];
    snprintf (inputImagePath, 1024, "%s%s", TEST_IMAGES_DIR, "input.png");

    struct Context context;
    int channels = 0;
    context.inputBuffer = stbi_load (inputImagePath, &context.width, &context.height, &channels, 4 /* rgba */);

    if (context.inputBuffer == NULL)
    {
        fprintf (stderr, "Could not read image %s", inputImagePath);
        return 2;
    }

    context.tmpImageBuffer = malloc(context.width * context.height * 4);

    int numFailed = 0;

    fprintf (stderr, ">> Testing Vienot 1999\n");
    if (test_Vienot1999 (&context) != 0)
    {
        fprintf (stderr, "TEST FAILED: Vienot 1999\n");
        ++numFailed;
    }

    fprintf (stderr, ">> Testing Brettel 1997\n");
    if (test_Brettel1997 (&context) != 0)
    {
        fprintf (stderr, "TEST FAILED: Brettel 1997\n");
        ++numFailed;
    }

    fprintf (stderr, ">> Testing Automatic Dispatch\n");
    if (test_automaticDispatch (&context) != 0)
    {
        fprintf (stderr, "TEST FAILED: automatic dispatch\n");
        ++numFailed;
    }

    // int success = stbi_write_png(argv[2], width, height, 4, imageBuffer, 0 /* default stride */);
    // if (!success)
    // {
    //     fprintf(stderr, "Could not write image to %s", argv[2]);
    //     return 4;
    // }

    stbi_image_free (context.inputBuffer);
    free (context.tmpImageBuffer);

    return numFailed;
}
