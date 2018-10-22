#include "camera_pipeline.hpp"
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <stdlib.h>     /* malloc, calloc, realloc, free */
#include <math.h>

// mutates inputs.
float** CameraPipeline::GaussianChannelBlur(float** channel, int height, int width) const {
    float** expandedChannel = (float**) malloc(sizeof(float*) * (height + 2));
    for (int i = 0; i < height + 2; i++) {
        expandedChannel[i] = (float*) malloc(sizeof(float) * (width+2));
    }
    float** blurredChannel = (float**) malloc(sizeof(float*) * height);
    for (int i = 0; i < height + 2; i++) {
        blurredChannel[i] = (float*) malloc(sizeof(float) * width);
    }
    std::cout << "Allocated" << std::endl;
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            expandedChannel[row + 1][col + 1] = channel[row][col];
        }
    }
    std::cout << "Expanded allocation" << std::endl;
    // mirror the edges
    for (int z = 0; z < height + 2; z++) {
        expandedChannel[z][0] = expandedChannel[z][1];
        expandedChannel[z][width+1] = expandedChannel[z][width];
        expandedChannel[0][z] = expandedChannel[1][z];
        expandedChannel[height+1][z] = expandedChannel[height][z];
    }
    std::cout << "Mirrored edges" << std::endl;
    float weights[] = {1.f/64, 3.f/64, 3.f/64, 1.f/64,
        3.f/64, 9.f/64, 9.f/64, 3.f/64,
        3.f/64, 9.f/64, 9.f/64, 3.f/64,
        1.f/64, 3.f/64, 3.f/64, 1.f/64
    };
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            float tmp = 0.f;
            for (int jj=0; jj < 3; jj++) {
                for (int ii=0; ii<3; ii++) {
                    tmp += expandedChannel[j + jj][i + ii] * weights[jj * 3 + ii];
                }
            }
            blurredChannel[j][i] = tmp;
        }
    }
    for(int i = 0; i < height+2; ++i) {
        free(expandedChannel[i]);
    }
    free(expandedChannel);
    return blurredChannel;
}

std::unique_ptr<Image<RgbPixel>> CameraPipeline::GaussianBlur(Image<RgbPixel>* originalImage) const {
    int height = originalImage->height();
    int width = originalImage->width();
    // copy into float array
    float** new_img_r = new float*[height];
    float** new_img_g = new float*[height];
    float** new_img_b = new float*[height];
    for (int i = 0; i < height; i++) {
        new_img_r[i] = new float[width];
        new_img_g[i] = new float[width];
        new_img_b[i] = new float[width];
    }
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            auto& pixel = (*originalImage)(row, col);
            new_img_r[row][col] = pixel.r;
            new_img_g[row][col] = pixel.g;
            new_img_b[row][col] = pixel.b;
        }
    }
    CameraPipeline::GaussianChannelBlur(new_img_r, height, width);
    CameraPipeline::GaussianChannelBlur(new_img_g, height, width);
    CameraPipeline::GaussianChannelBlur(new_img_b, height, width);
    std::unique_ptr<Image<RgbPixel>> blurredImage(new Image<RgbPixel>(width, height));
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            auto& pixel = (*blurredImage)(row, col);
            pixel.r = new_img_r[row][col];
            pixel.g = new_img_g[row][col];
            pixel.b = new_img_b[row][col];
        }
    }
    for(int i = 0; i < height+2; ++i) {
        delete [] new_img_r[i];
        delete [] new_img_g[i];
        delete [] new_img_b[i];
    }
    delete [] new_img_r;
    delete [] new_img_g;
    delete [] new_img_b;
    return blurredImage;
}

std::unique_ptr<Image<RgbPixel>> CameraPipeline::DownSample(std::unique_ptr<Image<RgbPixel>> &originalImage) const {
    std::unique_ptr<Image<RgbPixel>> blurredImage = CameraPipeline::GaussianBlur(originalImage.get());
    std::unique_ptr<Image<RgbPixel>> downImg(new Image<RgbPixel>(originalImage->width()/2, originalImage->height()/2));
    int newRowIdx = 0;
    int newColIdx = 0;
    for (int row = 0; row < originalImage->height(); row++) {
        if (row % 2 == 0) {
            for (int col = 0; col < originalImage->width(); col++) {
                if (col % 2 == 0) {
                    auto& newPixel = (*downImg)(newRowIdx, newColIdx);
                    auto& originalPixel = (*blurredImage)(row, col);
                    newPixel.r = originalPixel.r;
                    newPixel.g = originalPixel.g;
                    newPixel.b = originalPixel.b;
                    newColIdx++;
                }
            }
            newRowIdx++;
        }
    }
    return downImg;
}

float** CameraPipeline::UpSampleChannel(float** originalImage, int height, int width) const {
    std::cout << "Sizes:" << std::endl;
    std::cout << height << std::endl;
    std::cout << width << std::endl;
    float** upImg = (float**) malloc(sizeof(float*) * 2 * height);
    for (int z = 0; z < height*2; z++) {
        upImg[z] = (float*) malloc(sizeof(float) * 2 * width);
    }
    std::cout << "upsamplealloc" << std::endl;
    for (int j = 0; j < 2*height; j++) {
        for (int i = 0; i < 2*width; i++) {
            int row = j/2;
            int col = i/2;
            if (row == height-1) {
                std::cout << "Issued correction" << std::endl;
                row = height - 2;
            }
            if (col == width-1) {
                std::cout << "Issued correction" << std::endl;
                col = width - 2;
            }
            float w1 = (i%2) ? .75f : .25f;
            float w2 = (j%2) ? .75f : .25f;
//            std::cout << "Results:" << std::endl;
//            std::cout << upImg[j][i] << std::endl;
//            std::cout << originalImage[row][col] << std::endl;
//            std::cout << originalImage[row][col+1] << std::endl;
//            std::cout << originalImage[row+1][col] << std::endl;
//            std::cout << originalImage[row+1][col+1] << std::endl;
            upImg[j][i] = w1 * w2 * originalImage[row][col] + (1.0-w1) * w2 * originalImage[row][col+1] + w1 * (1 - w2) * originalImage[row+1][col] + (1.0-w1) * (1.0-w2) * originalImage[row+1][col+1];
        }
    }
    return upImg;
}

std::unique_ptr<Image<RgbPixel>> CameraPipeline::UpSample(std::unique_ptr<Image<RgbPixel>> &originalImage) const {
    std::unique_ptr<Image<RgbPixel>> upImg(new Image<RgbPixel>(2 * originalImage->height(), 2 * originalImage->width()));
    for (int j = 0; j < 2*originalImage->height(); j++) {
        for (int i = 0; i < 2* originalImage->width(); i++) {
            int row = j/2;
            int col = i/2;
            float w1 = (i%2) ? .75f : .25f;
            float w2 = (j%2) ? .75f : .25f;
            (*upImg)(j, i).r = w1 * w2 * (*originalImage)(row,col).r + (1.0-w1) * w2 * (*originalImage)(row,col+1).r + w1 * (1 - w2) * (*originalImage)(row+1,col).r + (1.0-w1) * (1.0-w2) * (*originalImage)(row+1,col+1).r;
            (*upImg)(j, i).g = w1 * w2 * (*originalImage)(row,col).g + (1.0-w1) * w2 * (*originalImage)(row,col+1).g + w1 * (1 - w2) * (*originalImage)(row+1,col).g + (1.0-w1) * (1.0-w2) * (*originalImage)(row+1,col+1).g;
            (*upImg)(j, i).b = w1 * w2 * (*originalImage)(row,col).b + (1.0-w1) * w2 * (*originalImage)(row,col+1).b + w1 * (1 - w2) * (*originalImage)(row+1,col).b + (1.0-w1) * (1.0-w2) * (*originalImage)(row+1,col+1).b;
        }
    }
    return upImg;
}

void CameraPipeline::Denoise(std::unique_ptr<Image<RgbPixel>> &originalImage) const {
    int height = originalImage->height();
    int width = originalImage->width();
    float weights[] = {1.f/16, 1.f/8, 1.f/16, 1.f/8, 1.f/4, 1.f/8, 1.f/16, 1.f/8, 1.f/16};
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            auto& pixel = (*originalImage)(row, col);
            float bilateral_val_r = 0.f;
            float norm_r = 0.f;
            float bilateral_val_g = 0.f;
            float norm_g = 0.f;
            float bilateral_val_b = 0.f;
            float norm_b = 0.f;
            for (int i = -1; i <= 1; i++) {
                for (int j = -1; j <= 1; j++) {
                    if ((row + i >= 0 && row + i < height) && (col + j >= 0 && col + j < width)) {
                        auto& other = (*originalImage)(row + i, col + j);
                        float f_val_r = 1 - std::abs(other.r - pixel.r);
                        float f_val_g = 1 - std::abs(other.g - pixel.g);
                        float f_val_b = 1 - std::abs(other.b - pixel.b);
                        bilateral_val_r += f_val_r * weights[(i + 1) * 3 + j] * other.r;
                        norm_r += f_val_r * weights[(i + 1) * 3 + j];
                        bilateral_val_g += f_val_g * weights[(i + 1) * 3 + j] * other.g;
                        norm_g += f_val_g * weights[(i + 1) * 3 + j];
                        bilateral_val_b += f_val_b * weights[(i + 1) * 3 + j] * other.b;
                        norm_b += f_val_b * weights[(i + 1) * 3 + j];
                    }
                }
            }
            pixel.r = bilateral_val_r / norm_r;
            pixel.g = bilateral_val_g / norm_g;
            pixel.b = bilateral_val_b / norm_b;
        }
    }
}

std::unique_ptr<Image<RgbPixel>> CameraPipeline::subImage(Image<RgbPixel>* im1, Image<RgbPixel>* im2) const {
    int height = im1->height();
    int width = im1->width();
    std::unique_ptr<Image<RgbPixel>> subImg(new Image<RgbPixel>(width, height));
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            auto& pixel1 = (*im1)(row, col);
            auto& pixel2 = (*im2)(row, col);
            (*subImg)(row, col) = pixel1 - pixel2;
        }
    }
    return subImg;
}

float** CameraPipeline::subImageChannel(float** im1, float** im2, int height, int width) const {
    float** subImg = (float**) malloc(sizeof(float*) * height);
    for (int i = 0; i < height; i++) {
        subImg[i] = (float*) malloc(sizeof(float) * width);
    }
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            float pixel1 = im1[row][col];
            float pixel2 = im2[row][col];
            subImg[row][col] = pixel1 - pixel2;
        }
    }
    return subImg;
}


std::unique_ptr<Image<RgbPixel>> CameraPipeline::exposureFusion(std::unique_ptr<Image<RgbPixel>> &originalImage, int levels) const {
    int height = originalImage->height();
    int width = originalImage->width();
    // convert to YUV space
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            auto& pixel = (*originalImage)(i, j);
            Float3Pixel newPixel = Float3Pixel::RgbToYuv(pixel);
            pixel.y = newPixel.y;
            pixel.u = newPixel.u;
            pixel.v = newPixel.v;
        }
    }
    std::cout << "Convert to YUV" << std::endl;
    // generate virtual exposure
    std::unique_ptr<Image<RgbPixel>> exposedImage(new Image<RgbPixel>(width, height));
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            auto& pixel = (*originalImage)(i, j);
            auto& exposedPixel = (*exposedImage)(i, j);
            exposedPixel = pixel * 1.2;
        }
    }
    std::cout << "Virtually exposed" << std::endl;
    // make a gaussian pyramid
    float** darkIm = (float**) malloc(sizeof(float*) * height);
    float** brightIm = (float**) malloc(sizeof(float*) * height);
    for (int i = 0; i < height; i++) {
        darkIm[i] = (float*) malloc(sizeof(float) * width);
        brightIm[i] = (float*) malloc(sizeof(float) * width);
    }
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            darkIm[row][col] = (*originalImage)(row, col).y;
            brightIm[row][col] = (*exposedImage)(row, col).y;
        }
    }
    std::cout << "Generated image copies." << std::endl;
    // generate the blending mask for the dark image.
    float** blend_mask_dark = (float**) malloc(sizeof(float*) * height);
    float** blend_mask_bright = (float**) malloc(sizeof(float*) * height);
    for (int i = 0; i < height; i++) {
        blend_mask_dark[i] = (float*) malloc(sizeof(float) * width);
        blend_mask_bright[i] = (float*) malloc(sizeof(float) * width);
    }
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            float pixel = darkIm[row][col];
            float brightPixel = brightIm[row][col];
            float weight = exp(-pow(pixel - 0.5, 2)/ 0.08);
            float brightWeight = exp(-pow(brightPixel - 0.5, 2)/ 0.08);
            blend_mask_dark[row][col] = weight;
            blend_mask_bright[row][col] = brightWeight;
        }
    }
    std::cout << "Finished masking" << std::endl;
    float** gaussianPyramidLevelsDark[levels];
    float** gaussianPyramidLevelsBright[levels];
    float** currentDark = darkIm;
    float** currentBright = brightIm;
    for (int i = 0; i < levels; i++) {
        std::cout << "Preblur" << std::endl;
        float** blurredDark = CameraPipeline::GaussianChannelBlur(currentDark, height, width);
        float** blurredBright = CameraPipeline::GaussianChannelBlur(currentBright, height, width);
        std::cout << "blurred" << std::endl;
        float** downSampledDark = (float**) malloc(sizeof(float*) * (height/2));
        float** downSampledBright = (float**) malloc(sizeof(float*) * (height/2));
        for (int z = 0; z < height/2; z++) {
            downSampledDark[z] = (float*) malloc(sizeof(float) * (width/2));
            downSampledBright[z] = (float*) malloc(sizeof(float) * (width/2));
        }
        std::cout << "downsample alloc" << std::endl;
        int rowIdx = 0;
        int colIdx = 0;
        for (int row = 0; row < height; row++) {
            if (row % 2 == 0) {
                for (int col = 0; col < width; col++) {
                    if (col % 2 == 0) {
                        downSampledDark[rowIdx][colIdx] = blurredDark[row][col];
                        downSampledBright[rowIdx][colIdx] = blurredBright[row][col];
                        colIdx++;
                    }
                }
                rowIdx++;
            }
        }
        std::cout << "pre-free" << std::endl;
        std::cout << "downsampled" << std::endl;
        float** nextLevelDark = CameraPipeline::UpSampleChannel(downSampledDark, height/2, width/2);
        std::cout << "after dark" << std::endl;
        float** nextLevelBright = CameraPipeline::UpSampleChannel(downSampledBright, height/2, width/2);
        std::cout << "upsampled" << std::endl;
        for (int z = 0; z < height/2; z++) {
            free(downSampledDark[z]);
            free(downSampledBright[z]);
        }
        free(downSampledDark);
        free(downSampledBright);
        for (int z = 0; z < height; z++) {
            free(blurredDark[z]);
            free(blurredBright[z]);
        }
        free(blurredDark);
        free(blurredBright);
        std::cout << i << std::endl;
        gaussianPyramidLevelsDark[i] = nextLevelDark;
        gaussianPyramidLevelsBright[i] = nextLevelDark;
        currentDark = nextLevelDark;
        currentBright = nextLevelBright;
        std::cout << "going to next iteration" << std::endl;
    }
    std::cout << "Finished Gaussian pyramid" << std::endl;
    // make a laplacian pyramid
    float** laplacianPyramidLevelsDark[levels];
    float** laplacianPyramidLevelsBright[levels];
    currentDark = darkIm;
    currentBright = brightIm;
    std::cout << "Starting Laplacian pyramid" << std::endl;
    for (int j = 0; j < levels; j++) {
        std::cout << "Before Sub" << std::endl;
        float** nextDark = CameraPipeline::subImageChannel(currentDark, gaussianPyramidLevelsDark[j], height, width);
        float** nextBright = CameraPipeline::subImageChannel(currentBright, gaussianPyramidLevelsBright[j], height, width);
        std::cout << "After Sub" << std::endl;
        laplacianPyramidLevelsDark[j] = nextDark;
        laplacianPyramidLevelsBright[j] = nextBright;
        currentBright = nextBright;
        currentDark = nextDark;
    }
    // for debugging, what does our mask even look like?
    std::unique_ptr<Image<YuvPixel>> mask(new Image<YuvPixel>(width, height));
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            auto& pixel = (*mask)(row, col);
            pixel.y = blend_mask_dark[row][col];
        }
    }
    // cleanup
    for (int i = 0; i < levels; i++) {
        float** curImgDark = gaussianPyramidLevelsDark[i];
        float** curImgBright = gaussianPyramidLevelsBright[i];
        float** curImgLplBright = laplacianPyramidLevelsBright[i];
        float** curImgLplDark = laplacianPyramidLevelsDark[i];
        for (int j = 0; j < height; j++) {
            free(curImgDark[j]);
            free(curImgBright[j]);
            free(curImgLplBright[j]);
            free(curImgLplDark[j]);
        }
        free(curImgDark);
        free(curImgBright);
        free(curImgLplBright);
        free(curImgLplDark);
    }
    for (int i = 0; i < height; i++) {
        free(blend_mask_dark[i]);
        free(blend_mask_bright[i]);
        free(darkIm[i]);
        free(brightIm[i]);
    }
    free(blend_mask_dark);
    free(blend_mask_bright);
    free(darkIm);
    free(brightIm);
    exposedImage.reset();
    return mask;
    // perform the blending.
}

std::unique_ptr<Image<RgbPixel>> CameraPipeline::ProcessShot() const {
    
    // BEGIN: CS348K STUDENTS MODIFY THIS CODE
    const bool BILINEAR_DEMOSAIC = true;
    // put the lens cap on if you'd like to measure a "dark frame"
    sensor_->SetLensCap(true);
    int width = sensor_->GetSensorWidth();
    int height = sensor_->GetSensorHeight();
    auto raw_data = sensor_->GetSensorData(0, 0, width, height);
    std::set<int> stuckPixels;
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            // find pixels that are on in a way more than noise.
            if (raw_data->data(row, col) > 0.5) {
                stuckPixels.insert(width * row + col);
            }
        }
    }
    sensor_->SetLensCap(false);
    
    // grab RAW pixel data from sensor
    width = sensor_->GetSensorWidth();
    height = sensor_->GetSensorHeight();
    raw_data = sensor_->GetSensorData(0, 0, width, height);
    std::set<int> :: iterator itr;
    for (itr = stuckPixels.begin(); itr != stuckPixels.end(); ++itr)
    {
        int linearized_index = *itr;
        int row = linearized_index / width;
        int col = linearized_index - (row * width);
        if ((row % 2 == 0 && col % 2 == 0) || (row % 2 == 1 and col % 2 == 1)) {
            // pattern is:
            // G, B/R, G
            // R/B, G, R/B
            // G, B/R, G
            raw_data->data(row, col) = 0;
            int g_count = 0;
            if (row - 1 > 0) {
                if (col - 1 >= 0) {
                    raw_data->data(row, col) += raw_data->data(row-1, col-1);;
                    g_count += 1;
                }
                if (col + 1 < width) {
                    raw_data->data(row, col) += raw_data->data(row-1, col+1);
                    g_count += 1;
                }
            }
            if (row + 1 < height) {
                if (col - 1 >= 0) {
                    raw_data->data(row, col) += raw_data->data(row+1, col-1);;
                    g_count += 1;
                }
                if (col + 1 < width) {
                    raw_data->data(row, col) += raw_data->data(row+1, col+1);
                    g_count += 1;
                }
            }
            raw_data->data(row, col) = (raw_data->data(row, col) / g_count);
        } else if ((row % 2 == 0 && col % 2 == 1) || (row % 2 == 1 && col % 2 == 0)) {
            // pattern is:
            // B/R, G, B/R
            // G, R/B, G
            // B/R, G, B/R
            raw_data->data(row, col) = 0;
            int r_count = 0;
            if (row - 2 >= 0) {
                raw_data->data(row, col) += raw_data->data(row-2, col);
                r_count += 1;
            }
            if (row + 2 < height) {
                raw_data->data(row, col) += raw_data->data(row+2, col);
                r_count += 1;
            }
            if (col - 2 >= 0) {
                raw_data->data(row, col) += raw_data->data(row, col-2);
                r_count += 1;
            }
            if (col + 2 < width) {
                raw_data->data(row, col) += raw_data->data(row, col+2);
                r_count += 1;
            }
            raw_data->data(row, col) = (raw_data->data(row, col) / r_count);
        }
    }
    // For denoising, could use bilateral filter.
    std::unique_ptr<Image<RgbPixel>> newImage(new Image<RgbPixel>(width, height));
    // In this function you should implement your full RAW image processing pipeline.
    //   (1) Demosaicing
    //   (2) Address sensing defects such as bad pixels and image noise.
    //   (3) Apply local tone mapping based on the local laplacian filter or exposure fusion.
    //   (4) gamma correction
    // allocate 3-channel RGB output buffer to hold the results after processing
    std::unique_ptr<Image<RgbPixel>> image(new Image<RgbPixel>(width, height));
    if (BILINEAR_DEMOSAIC) {
        for (int row = 0; row < height; row++) {
            for (int col = 0; col < width; col++) {
                int linearized_index = row * width + col;
                auto &pixel = (*image)(row, col);
                const auto pixel_val = raw_data->data(row, col);
                if (row % 2 == 0 && col % 2 == 0) {
                    // pattern is:
                    // G, B, G
                    // R, G, R
                    // G, B, G
                    pixel.g = pixel_val;
                    int b_count = 0;
                    int r_count = 0;
                    if (row - 1 >= 0) {
                        pixel.b += raw_data->data(row-1, col);
                        b_count += 1;
                    }
                    if (row + 1 < height) {
                        pixel.b += raw_data->data(row+1, col);
                        b_count += 1;
                    }
                    if (col - 1 >= 0) {
                        pixel.r += raw_data->data(row, col-1);
                        r_count += 1;
                    }
                    if (col + 1 < width) {
                        pixel.r += raw_data->data(row, col+1);
                        r_count += 1;
                    }
                    pixel.b = (pixel.b / b_count);
                    pixel.r = (pixel.r / r_count);
                } else if (row % 2 == 0 && col % 2 == 1) {
                    // pattern is:
                    // B, G, B
                    // G, R, G
                    // B, G, B
                    pixel.r = pixel_val;
                    int b_count = 0;
                    int g_count = 0;
                    if (row - 1 >= 0) {
                        if (col - 1 >= 0) {
                            pixel.b += raw_data->data(row-1, col-1);
                            b_count += 1;
                        }
                        if (col + 1 < width) {
                            pixel.b += raw_data->data(row-1, col+1);
                            b_count += 1;
                        }
                        pixel.g += raw_data->data(row-1, col);
                        g_count += 1;
                    }
                    if (col - 1 >= 0) {
                        pixel.g += raw_data->data(row, col-1);
                        g_count += 1;
                    }
                    if (col + 1 < width) {
                        pixel.g += raw_data->data(row, col+1);
                        g_count += 1;
                    }
                    if (row + 1 < height) {
                        if (col - 1 >= 0) {
                            pixel.b += raw_data->data(row+1, col-1);
                            b_count += 1;
                        }
                        if (col + 1 < width) {
                            pixel.b += raw_data->data(row+1, col+1);
                            b_count += 1;
                        }
                        pixel.g += raw_data->data(row+1, col);
                        g_count += 1;
                    }
                    pixel.b = (pixel.b / b_count);
                    pixel.g = (pixel.g / g_count);
                } else if (row % 2 == 1 && col % 2 == 0) {
                    // pattern is:
                    // R, G, R
                    // G, B, G
                    // R, G, R
                    pixel.b = pixel_val;
                    int r_count = 0;
                    int g_count = 0;
                    if (row - 1 >= 0) {
                        if (col - 1 >= 0) {
                            pixel.r += raw_data->data(row-1, col-1);
                            r_count += 1;
                        }
                        if (col + 1 < width) {
                            pixel.r += raw_data->data(row-1, col+1);
                            r_count += 1;
                        }
                        pixel.g += raw_data->data(row-1, col);
                        g_count += 1;
                    }
                    if (col - 1 >= 0) {
                        pixel.g += raw_data->data(row, col-1);
                        g_count += 1;
                    }
                    if (col + 1 < width) {
                        pixel.g += raw_data->data(row, col+1);
                        g_count += 1;
                    }
                    if (row + 1 < height) {
                        if (col - 1 >= 0) {
                            pixel.r += raw_data->data(row+1, col-1);
                            r_count += 1;
                        }
                        if (col + 1 < width) {
                            pixel.r += raw_data->data(row+1, col+1);
                            r_count += 1;
                        }
                        pixel.g += raw_data->data(row+1, col);
                        g_count += 1;
                    }
                    pixel.r = (pixel.r / r_count);
                    pixel.g = (pixel.g / g_count);
                } else {
                    // this is when row % 2 == 1 and col % 2 == 1
                    // pattern is:
                    // G, R, G
                    // B, G, B
                    // G, R, G
                    pixel.g = pixel_val;
                    int b_count = 0;
                    int r_count = 0;
                    if (row - 1 >= 0) {
                        pixel.r += raw_data->data(row-1, col);
                        r_count += 1;
                    }
                    if (row + 1 < height) {
                        pixel.r += raw_data->data(row+1, col);
                        r_count += 1;
                    }
                    if (col - 1 >= 0) {
                        pixel.b += raw_data->data(row, col-1);
                        b_count += 1;
                    }
                    if (col + 1 < width) {
                        pixel.b += raw_data->data(row, col+1);
                        b_count += 1;
                    }
                    pixel.b = (pixel.b / b_count);
                    pixel.r = (pixel.r / r_count);
                }
            }
        }
    }
//    // convert to YUV space
//    for (int i = 0; i < height; i++) {
//        for (int j = 0; j < width; j++) {
//            auto& pixel = (*image)(i, j);
//            Float3Pixel newPixel = Float3Pixel::RgbToYuv(pixel);
//            pixel.y = newPixel.y;
//            pixel.u = newPixel.u;
//            pixel.v = newPixel.v;
//        }
//    }
//    // low-pass filter on image
//    std::unique_ptr<Image<RgbPixel>> blurredImage = CameraPipeline::GaussianBlur(image.get());
//    // take the Cb and Cr blurred levels and put those in.
//    for (int i = 0; i < height; i++) {
//        for (int j = 0; j < width; j++) {
//            auto& blurPixel = (*blurredImage)(i, j);
//            auto& pixel = (*image)(i, j);
//            pixel.y = pixel.y;
//            pixel.u = blurPixel.u;
//            pixel.v = blurPixel.v;
//        }
//    }
//    // convert the image back to RGB space.
//    for (int i = 0; i < height; i++) {
//        for (int j = 0; j < width; j++) {
//            auto& pixel = (*image)(i, j);
//            Float3Pixel newPixel = Float3Pixel::YuvToRgb(pixel);
//            pixel.r = newPixel.r;
//            pixel.g = newPixel.g;
//            pixel.b = newPixel.b;
//        }
//    }
    CameraPipeline::Denoise(image);
//    // apply some gain because brighter images tend to look nicer.
//    for (int i = 0; i < height; i++) {
//        for (int j = 0; j < width; j++) {
//            auto& pixel = (*image)(i, j);
//            pixel.r = std::min(pixel.r + 0.15 * pixel.r, 1.0);
//            pixel.g = std::min(pixel.g + 0.15 * pixel.g, 1.0);
//            pixel.b = std::min(pixel.b + 0.15 * pixel.b, 1.0);
//        }
//    }
    
    // gamma correct to use more bit space.
    image->GammaCorrect(0.41);
    std::cout << "Calling exposure fusion" << std::endl;
    image = CameraPipeline::exposureFusion(image, 4);
    // map to 255 for 8-bit output.
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            auto& pixel = (*image)(row, col);
            // pixel data from the sensor is normalized to the 0-1 range, so
            // scale by 255 for the final image output.  Output image pixels
            // should be in the 0-255 range.
            pixel.r = pixel.r * 255.f;
            pixel.g = pixel.g * 255.f;
            pixel.b = pixel.b * 255.f;
        }
    }
    // return processed image output
    return image;
    
    // END: CS348K STUDENTS MODIFY THIS CODE
}
