#include "camera_pipeline.hpp"
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <limits>

std::unique_ptr<Image<FloatPixel>> CameraPipeline::GaussianBlurChannel(Image<FloatPixel>* originalImage) const {
    int height = originalImage->height();
    int width = originalImage->width();
    std::unique_ptr<Image<FloatPixel>> blurredImage(new Image<FloatPixel>(width, height));
    // copy into float array
    float** new_img_r = new float*[height + 2];
    for (int i = 0; i < height + 2; i++) {
        new_img_r[i] = new float[width+2];
    }
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            auto& pixel = (*originalImage)(row, col);
            new_img_r[row + 1][col + 1] = pixel.i;
        }
    }
    // mirror the edges
    for (int z = 0; z < height + 2; z++) {
        new_img_r[z][0] = new_img_r[z][1];
        new_img_r[z][width+1] = new_img_r[z][width];
        new_img_r[0][z] = new_img_r[1][z];
        new_img_r[height+1][z] = new_img_r[height][z];
    }
    float weights[] = {1.f/64, 3.f/64, 3.f/64, 1.f/64,
        3.f/64, 9.f/64, 9.f/64, 3.f/64,
        3.f/64, 9.f/64, 9.f/64, 3.f/64,
        1.f/64, 3.f/64, 3.f/64, 1.f/64
    };
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            float tmpr = 0.f;
            for (int jj=0; jj < 3; jj++) {
                for (int ii=0; ii<3; ii++) {
                    tmpr += new_img_r[j + jj][i + ii] * weights[jj * 3 + ii];
                }
            }
            (*blurredImage)(j, i).i = tmpr;
        }
    }
    for(int i = 0; i < height+2; ++i) {
        delete [] new_img_r[i];
    }
    delete [] new_img_r;
    return blurredImage;
}

std::unique_ptr<Image<RgbPixel>> CameraPipeline::GaussianBlur(Image<RgbPixel>* originalImage) const {
    int height = originalImage->height();
    int width = originalImage->width();
    std::unique_ptr<Image<RgbPixel>> blurredImage(new Image<RgbPixel>(width, height));
    // copy into float array
    float** new_img_r = new float*[height + 2];
    float** new_img_g = new float*[height + 2];
    float** new_img_b = new float*[height + 2];
    for (int i = 0; i < height + 2; i++) {
        new_img_r[i] = new float[width+2];
        new_img_g[i] = new float[width+2];
        new_img_b[i] = new float[width+2];
    }
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            auto& pixel = (*originalImage)(row, col);
            new_img_r[row + 1][col + 1] = pixel.r;
            new_img_g[row + 1][col + 1] = pixel.g;
            new_img_b[row + 1][col + 1] = pixel.b;
        }
    }
    // mirror the edges
    for (int z = 0; z < height + 2; z++) {
        new_img_r[z][0] = new_img_r[z][1];
        new_img_r[z][width+1] = new_img_r[z][width];
        new_img_r[0][z] = new_img_r[1][z];
        new_img_r[height+1][z] = new_img_r[height][z];
        new_img_g[z][0] = new_img_g[z][1];
        new_img_g[z][width+1] = new_img_g[z][width];
        new_img_g[0][z] = new_img_g[1][z];
        new_img_g[height+1][z] = new_img_g[height][z];
        new_img_b[z][0] = new_img_b[z][1];
        new_img_b[z][width+1] = new_img_b[z][width];
        new_img_b[0][z] = new_img_b[1][z];
        new_img_b[height+1][z] = new_img_b[height][z];
    }
    float weights[] = {1.f/64, 3.f/64, 3.f/64, 1.f/64,
        3.f/64, 9.f/64, 9.f/64, 3.f/64,
        3.f/64, 9.f/64, 9.f/64, 3.f/64,
        1.f/64, 3.f/64, 3.f/64, 1.f/64
    };
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            float tmpr = 0.f;
            float tmpg = 0.f;
            float tmpb = 0.f;
            for (int jj=0; jj < 3; jj++) {
                for (int ii=0; ii<3; ii++) {
                    tmpr += new_img_r[j + jj][i + ii] * weights[jj * 3 + ii];
                    tmpg += new_img_g[j + jj][i + ii] * weights[jj * 3 + ii];
                    tmpb += new_img_b[j + jj][i + ii] * weights[jj * 3 + ii];
                }
            }
            (*blurredImage)(j, i).r = tmpr;
            (*blurredImage)(j, i).g = tmpg;
            (*blurredImage)(j, i).b = tmpb;
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

std::unique_ptr<Image<FloatPixel>> CameraPipeline::DownSampleChannel(Image<FloatPixel> *originalImage) const {
    std::unique_ptr<Image<FloatPixel>> blurredImage = CameraPipeline::GaussianBlurChannel(originalImage);
    std::unique_ptr<Image<FloatPixel>> downImg(new Image<FloatPixel>(originalImage->width()/2, originalImage->height()/2));
    int newRowIdx = 0;
    int newColIdx = 0;
    for (int row = 0; row < blurredImage->height(); row++) {
        if (row % 2 == 0) {
            newColIdx = 0;
            for (int col = 0; col < blurredImage->width(); col++) {
                if (col % 2 == 0) {
//                    std::cout << newRowIdx << ", " << newColIdx << std::endl;
                    auto& newPixel = (*downImg)(newRowIdx, newColIdx);
                    auto& originalPixel = (*blurredImage)(row, col);
                    newPixel.i = (*blurredImage)(row, col).i;
                    newColIdx++;
                }
            }
            newRowIdx++;
        }
    }
    return downImg;
}

std::unique_ptr<Image<RgbPixel>> CameraPipeline::DownSample(Image<RgbPixel> *originalImage) const {
    std::unique_ptr<Image<RgbPixel>> blurredImage = CameraPipeline::GaussianBlur(originalImage);
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

std::unique_ptr<Image<FloatPixel>> CameraPipeline::UpSampleChannel(Image<FloatPixel>* originalImage) const {
    std::unique_ptr<Image<FloatPixel>> upImg(new Image<FloatPixel>(2 * originalImage->width(), 2 * originalImage->height()));
    for (int j = 0; j < 2*originalImage->height(); j++) {
        for (int i = 0; i < 2* originalImage->width(); i++) {
            int row = j/2;
            int col = i/2;
            float w1 = (i%2) ? .75f : .25f;
            float w2 = (j%2) ? .75f : .25f;
            (*upImg)(j, i).i = w1 * w2 * (*originalImage)(row,col).i + (1.0-w1) * w2 * (*originalImage)(row,col+1).i + w1 * (1 - w2) * (*originalImage)(row+1,col).i + (1.0-w1) * (1.0-w2) * (*originalImage)(row+1,col+1).i;
        }
    }
    return upImg;
}

std::unique_ptr<Image<RgbPixel>> CameraPipeline::UpSample(Image<RgbPixel>* originalImage) const {
    std::unique_ptr<Image<RgbPixel>> upImg(new Image<RgbPixel>(2 * originalImage->width(), 2 * originalImage->height()));
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

std::unique_ptr<Image<FloatPixel>> CameraPipeline::subImageChannel(Image<FloatPixel>* im1, Image<FloatPixel>* im2) const {
    std::unique_ptr<Image<FloatPixel>> subImg(new Image<FloatPixel>(im1->width(),im1->height()));
    for (int row = 0; row < im1->height(); row++) {
        for (int col = 0; col < im1->width(); col++) {
            float val = (*im1)(row, col).i - (*im2)(row, col).i;
            // clip subtraction.
            if (val < 0) {
                val = 0;
            } else if (val > 1) {
                val = 1;
            }
            (*subImg)(row, col).i = val;
        }
    }
    return subImg;
}

std::unique_ptr<Image<RgbPixel>> CameraPipeline::exposureFusion(Image<RgbPixel> *originalImage, int nLevels) const {
    int height = originalImage->height();
    int width = originalImage->width();
    std::unique_ptr<Image<FloatPixel>> darkImage(new Image<FloatPixel>(width, height));
    // convert to grayscale space
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            auto& pixel = (*originalImage)(i, j);
            auto& darkPixel = (*darkImage)(i, j);
            Float3Pixel newPixel = Float3Pixel::RgbToYuv(pixel);
            darkPixel.i = newPixel.y;
        }
    }
    // generate virtual exposure
    std::unique_ptr<Image<FloatPixel>> exposedImage(new Image<FloatPixel>(width, height));
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            auto& pixel = (*originalImage)(i, j);
            auto& exposedPixel = (*exposedImage)(i, j);
            exposedPixel.i = std::min(pixel.y * 1.4, 1.0);
        }
    }
    std::unique_ptr<Image<FloatPixel>> blend_mask_dark(new Image<FloatPixel>(width, height));
    std::unique_ptr<Image<FloatPixel>> blend_mask_bright(new Image<FloatPixel>(width, height));
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            auto& pixel = (*originalImage)(row, col).y;
            auto& brightPixel = (*exposedImage)(row, col).i;
            float weight = exp(-pow(pixel - 0.5, 2)/ 0.08);
            float brightWeight = exp(-pow(brightPixel - 0.5, 2)/ 0.08);
            (*blend_mask_dark)(row, col).i = weight;
            (*blend_mask_bright)(row, col).i = brightWeight;
        }
    }
    // generate G1 for dark, bright, mask1, mask2
    std::cout << "Downsampling starts now." << std::endl;
    // unrolled 4-level Laplacian Pyramid blending.
    std::unique_ptr<Image<FloatPixel>> downDarkG1 = CameraPipeline::DownSampleChannel(darkImage.get());
    std::unique_ptr<Image<FloatPixel>> downBrightG1 = CameraPipeline::DownSampleChannel(exposedImage.get());
    std::unique_ptr<Image<FloatPixel>> downDarkMaskG1 = CameraPipeline::DownSampleChannel(blend_mask_dark.get());
    std::unique_ptr<Image<FloatPixel>> downBrightMaskG1 = CameraPipeline::DownSampleChannel(blend_mask_bright.get());
    std::unique_ptr<Image<FloatPixel>> DarkG1 = CameraPipeline::UpSampleChannel(downDarkG1.get());
    std::unique_ptr<Image<FloatPixel>> BrightG1 = CameraPipeline::UpSampleChannel(downBrightG1.get());
    std::unique_ptr<Image<FloatPixel>> DarkMaskG1 = CameraPipeline::UpSampleChannel(downDarkMaskG1.get());
    std::unique_ptr<Image<FloatPixel>> BrightMaskG1 = CameraPipeline::UpSampleChannel(downBrightMaskG1.get());
    downDarkG1.reset();
    downBrightG1.reset();
    downDarkMaskG1.reset();
    downBrightMaskG1.reset();
    std::unique_ptr<Image<FloatPixel>> downDarkG2 = CameraPipeline::DownSampleChannel(DarkG1.get());
    std::unique_ptr<Image<FloatPixel>> downBrightG2 = CameraPipeline::DownSampleChannel(BrightG1.get());
    std::unique_ptr<Image<FloatPixel>> downDarkMaskG2 = CameraPipeline::DownSampleChannel(DarkMaskG1.get());
    std::unique_ptr<Image<FloatPixel>> downBrightMaskG2 = CameraPipeline::DownSampleChannel(BrightMaskG1.get());
    std::unique_ptr<Image<FloatPixel>> DarkG2 = CameraPipeline::UpSampleChannel(downDarkG2.get());
    std::unique_ptr<Image<FloatPixel>> BrightG2 = CameraPipeline::UpSampleChannel(downBrightG2.get());
    std::unique_ptr<Image<FloatPixel>> DarkMaskG2 = CameraPipeline::UpSampleChannel(downDarkMaskG2.get());
    std::unique_ptr<Image<FloatPixel>> BrightMaskG2 = CameraPipeline::UpSampleChannel(downBrightMaskG2.get());
    downDarkG2.reset();
    downBrightG2.reset();
    downDarkMaskG2.reset();
    downBrightMaskG2.reset();
    std::unique_ptr<Image<FloatPixel>> DarkL1 = CameraPipeline::subImageChannel(darkImage.get(), DarkG1.get());
    std::unique_ptr<Image<FloatPixel>> BrightL1 = CameraPipeline::subImageChannel(exposedImage.get(), BrightG1.get());
    // generate L1 for images-> G0 - G1
    std::unique_ptr<Image<FloatPixel>> blendedImage(new Image<FloatPixel>(width, height));
    // merge the two.
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            // scale weights to add up to 1.
            float factor = (*blend_mask_dark)(row, col).i/(*blend_mask_bright)(row, col).i;
            float scaledDarkWeight = factor / (factor + 1);
            float scaledBrightWeight = 1 / (factor + 1);
            (*blendedImage)(row, col).i += scaledDarkWeight * (*DarkL1)(row, col).i + scaledBrightWeight * (*BrightL1)(row, col).i;
        }
    }
    DarkL1.reset();
    BrightL1.reset();
    std::unique_ptr<Image<FloatPixel>> DarkL2 = CameraPipeline::subImageChannel(DarkG1.get(), DarkG2.get());
    std::unique_ptr<Image<FloatPixel>> BrightL2 = CameraPipeline::subImageChannel(BrightG1.get(), BrightG2.get());
    // merge the two.
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            // scale weights to add up to 1.
            float factor = (*DarkMaskG1)(row, col).i/(*BrightMaskG1)(row, col).i;
            float scaledDarkWeight = factor / (factor + 1);
            float scaledBrightWeight = 1 / (factor + 1);
            (*blendedImage)(row, col).i += scaledDarkWeight * (*DarkL2)(row, col).i + scaledBrightWeight * (*BrightL2)(row, col).i;
        }
    }
    std::unique_ptr<Image<FloatPixel>> downDarkG3 = CameraPipeline::DownSampleChannel(DarkG2.get());
    std::unique_ptr<Image<FloatPixel>> downBrightG3 = CameraPipeline::DownSampleChannel(BrightG2.get());
    std::unique_ptr<Image<FloatPixel>> downDarkMaskG3 = CameraPipeline::DownSampleChannel(DarkMaskG2.get());
    std::unique_ptr<Image<FloatPixel>> downBrightMaskG3 = CameraPipeline::DownSampleChannel(BrightMaskG2.get());
    std::unique_ptr<Image<FloatPixel>> DarkG3 = CameraPipeline::UpSampleChannel(downDarkG3.get());
    std::unique_ptr<Image<FloatPixel>> BrightG3 = CameraPipeline::UpSampleChannel(downBrightG3.get());
    std::unique_ptr<Image<FloatPixel>> DarkMaskG3 = CameraPipeline::UpSampleChannel(downDarkMaskG3.get());
    std::unique_ptr<Image<FloatPixel>> BrightMaskG3 = CameraPipeline::UpSampleChannel(downBrightMaskG3.get());
    downDarkG3.reset();
    downBrightG3.reset();
    downDarkMaskG3.reset();
    downBrightMaskG3.reset();
    std::unique_ptr<Image<FloatPixel>> DarkL3 = CameraPipeline::subImageChannel(DarkG2.get(), DarkG3.get());
    std::unique_ptr<Image<FloatPixel>> BrightL3 = CameraPipeline::subImageChannel(BrightG2.get(), BrightG3.get());
    // merge the two.
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            // scale weights to add up to 1.
            float factor = (*DarkMaskG2)(row, col).i/(*BrightMaskG2)(row, col).i;
            float scaledDarkWeight = factor / (factor + 1);
            float scaledBrightWeight = 1 / (factor + 1);
            (*blendedImage)(row, col).i += scaledDarkWeight * (*DarkL3)(row, col).i + scaledBrightWeight * (*BrightL3)(row, col).i;
        }
    }
    std::unique_ptr<Image<FloatPixel>> downDarkG4 = CameraPipeline::DownSampleChannel(DarkG3.get());
    std::unique_ptr<Image<FloatPixel>> downBrightG4 = CameraPipeline::DownSampleChannel(BrightG3.get());
    std::unique_ptr<Image<FloatPixel>> downDarkMaskG4 = CameraPipeline::DownSampleChannel(DarkMaskG3.get());
    std::unique_ptr<Image<FloatPixel>> downBrightMaskG4 = CameraPipeline::DownSampleChannel(BrightMaskG3.get());
    std::unique_ptr<Image<FloatPixel>> DarkG4 = CameraPipeline::UpSampleChannel(downDarkG4.get());
    std::unique_ptr<Image<FloatPixel>> BrightG4 = CameraPipeline::UpSampleChannel(downBrightG4.get());
    std::unique_ptr<Image<FloatPixel>> DarkMaskG4 = CameraPipeline::UpSampleChannel(downDarkMaskG4.get());
    std::unique_ptr<Image<FloatPixel>> BrightMaskG4 = CameraPipeline::UpSampleChannel(downBrightMaskG4.get());
    downDarkG3.reset();
    downBrightG3.reset();
    downDarkMaskG3.reset();
    downBrightMaskG3.reset();
    std::unique_ptr<Image<FloatPixel>> DarkL4 = CameraPipeline::subImageChannel(DarkG3.get(), DarkG4.get());
    std::unique_ptr<Image<FloatPixel>> BrightL4 = CameraPipeline::subImageChannel(BrightG3.get(), BrightG4.get());
    // merge the two.
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            // scale weights to add up to 1.
            float factor = (*DarkMaskG3)(row, col).i/(*BrightMaskG3)(row, col).i;
            float scaledDarkWeight = factor / (factor + 1);
            float scaledBrightWeight = 1 / (factor + 1);
            (*blendedImage)(row, col).i += scaledDarkWeight * (*DarkL4)(row, col).i + scaledBrightWeight * (*BrightL4)(row, col).i;
        }
    }
//    // finish up and add back the G4 blend.
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            // scale weights to add up to 1.
            float factor = (*DarkMaskG4)(row, col).i/(*BrightMaskG4)(row, col).i;
            float scaledDarkWeight = factor / (factor + 1);
            float scaledBrightWeight = 1 / (factor + 1);
            (*blendedImage)(row, col).i += scaledDarkWeight * (*DarkG4)(row, col).i + scaledBrightWeight * (*BrightG4)(row, col).i;
        }
    }
    // remove the levels of the Gaussian pyramid
    DarkMaskG1.reset();
    BrightMaskG1.reset();
    DarkMaskG2.reset();
    BrightMaskG2.reset();
    DarkMaskG3.reset();
    BrightMaskG3.reset();
    DarkMaskG4.reset();
    BrightMaskG4.reset();
    DarkG1.reset();
    BrightG1.reset();
    DarkG2.reset();
    BrightG2.reset();
    DarkG3.reset();
    BrightG3.reset();
    DarkG4.reset();
    BrightG4.reset();
    DarkL1.reset();
    BrightL1.reset();
    DarkL2.reset();
    BrightL2.reset();
    DarkL3.reset();
    BrightL3.reset();
    DarkL4.reset();
    BrightL4.reset();
    std::cout << "Resetting is done" << std::endl;
    // clip - we could do a rescale if we wanted to later.
//    for (int row = 0; row < height; row++) {
//        for (int col = 0; col < width; col++) {
//            if ((*blendedImage)(row, col).i > 1) {
//                (*blendedImage)(row, col).i = 1;
//            } else if ((*blendedImage)(row, col).i < 0) {
//                (*blendedImage)(row, col).i = 0;
//            }
//        }
//    }
    std::cout << "Clipping is done" << std::endl;
    std::unique_ptr<Image<YuvPixel>> blendedResult(new Image<YuvPixel>(width, height));
//    // convert to RGB space
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            auto &originalPixel = (*originalImage)(i, j);
            Float3Pixel oldUV = Float3Pixel::RgbToYuv(originalPixel);
            auto& pixel = (*blendedImage)(i, j);
            auto& blendedPixel = (*blendedResult)(i, j);
            blendedPixel.y = pixel.i;
            blendedPixel.u = oldUV.u;
            blendedPixel.v = oldUV.v;
            Float3Pixel newPixel = Float3Pixel::YuvToRgb(blendedPixel);
            blendedPixel.r = newPixel.r;
            blendedPixel.g = newPixel.g;
            blendedPixel.b = newPixel.b;
        }
    }
    // test masks
//    for (int i = 0; i < height; i++) {
//        for (int j = 0; j < width; j++) {
//            auto& pixel = (*blend_mask_bright)(i, j);
//            auto& blendedPixel = (*blendedResult)(i, j);
//            blendedPixel.y = pixel.i;
//            blendedPixel.u = 0;
//            blendedPixel.v = 0;
//            Float3Pixel newPixel = Float3Pixel::YuvToRgb(blendedPixel);
//            blendedPixel.r = newPixel.r;
//            blendedPixel.g = newPixel.g;
//            blendedPixel.b = newPixel.b;
//        }
//    }
    std::cout << "Conversion is done" << std::endl;
    return blendedResult;
}

//std::vector<Image<RgbPixel>> CameraPipeline::alignImages(std::vector<std::unique_ptr<Image<RgbPixel>>> images) {
//
//    return
//}

std::unique_ptr<Image<RgbPixel>> CameraPipeline::ProcessShot() const {
    // BEGIN: CS348K STUDENTS MODIFY THIS CODE
    const bool BILINEAR_DEMOSAIC = true;
    // put the lens cap on if you'd like to measure a "dark frame"
    sensor_->SetLensCap(true);
    int width = sensor_->GetSensorWidth();
    int height = sensor_->GetSensorHeight();
    std::cout << "Sensor Width" << width << std::endl;
    std::cout << "Sensor Height" << height << std::endl;
    auto raw_data_list = sensor_->GetBurstSensorData(0, 0, width, height);
    std::vector<std::unique_ptr<Image<FloatPixel>>> grayDownImages;
    // alignment code.
    for (int i = 0; i < raw_data_list.size(); i++) {
        int newRowIdx = 0;
        int newColIdx = 0;
        auto& raw_data = raw_data_list[i];
        std::unique_ptr<Image<FloatPixel>> newImage(new Image<FloatPixel>(width/2, height/2));
        for (int row = 0; row < height; row += 2) {
            newColIdx = 0;
            for (int col = 0; col < width; col += 2) {
//                std::cout << newRowIdx << ", " << newColIdx << std::endl;
                (*newImage)(newRowIdx, newColIdx).i = (raw_data->data(row, col) + raw_data->data(row, col+1) + raw_data->data(row+1, col) + raw_data->data(row+1, col+1)) / 4.0;
                newColIdx++;
            }
            newRowIdx++;
        }
        grayDownImages.emplace_back(std::move(newImage));
    }
    std::vector<std::vector<std::unique_ptr<Image<FloatPixel>>>> gaussianPyramids;
    for (int i = 0; i < grayDownImages.size(); i++) {
        int nLevels = 4;
        std::vector<std::unique_ptr<Image<FloatPixel>>> pyramid;
        pyramid.emplace_back(std::move(grayDownImages[i]));
        std::unique_ptr<Image<FloatPixel>> downG1 = CameraPipeline::DownSampleChannel(pyramid[0].get());
        std::unique_ptr<Image<FloatPixel>> downG2 = CameraPipeline::DownSampleChannel(downG1.get());
        std::unique_ptr<Image<FloatPixel>> downG3 = CameraPipeline::DownSampleChannel(downG2.get());
        std::unique_ptr<Image<FloatPixel>> downG4 = CameraPipeline::DownSampleChannel(downG3.get());
        pyramid.emplace_back(std::move(downG1));
        pyramid.emplace_back(std::move(downG2));
        pyramid.emplace_back(std::move(downG3));
        pyramid.emplace_back(std::move(downG4));
        gaussianPyramids.emplace_back(std::move(pyramid));
    }
    int tile_size = 16;
    int stride = 8;
    // alignment step
    std::vector<std::vector<int>> final_offsets;
    for (int level = 4; level >= 0; level--) {
        // reference image is the first image's image at this level
        std::cout << "image level" << level << std::endl;
        std::cout << gaussianPyramids[0][level]->height() << std::endl;
        std::cout << gaussianPyramids[0][level]->width() << std::endl;
        // tile size 16 - stride
        std::vector<std::vector<int>> offsets;
        int window = 10; // tiles
        for (int i = 0; i < gaussianPyramids[0][level]->height(); i+=tile_size) {
            for (int j = 0; j < gaussianPyramids[0][level]->width(); j+=tile_size) {
                // we're now iterating over blocks in the original image
                // copy into a buffer representing the current tile.
                std::unique_ptr<Image<FloatPixel>> currentBlock(new Image<FloatPixel>(tile_size, tile_size));
                int rowIdx = 0;
                int colIdx = 0;
                for (int row = i; row < i+16; row++) {
                    colIdx = 0;
                    for (int col = j; col < j+16; col++) {
                        (*currentBlock)(rowIdx, colIdx).i = (*gaussianPyramids[0][level])(row, col).i;
                        colIdx++;
                    }
                    rowIdx++;
                }
                float bestDiff = std::numeric_limits<float>::max();
                int bestVerticalOffset = -1;
                int bestHorizontalOffset = -1;
                // explore the window around the tile for alternative matches.
                for (int k = i - (window * stride); k < i + (window * stride); k += stride) {
                    if ((k < 0) || (k + tile_size > gaussianPyramids[1][level]->height())) {
                        continue;
                    }
                    for (int m = j - (window * stride); m < j + (window * stride); m += stride) {
                        if ((m < 0) || (m + tile_size > gaussianPyramids[1][level]->height())) {
                            continue;
                        }
                        // subtract the tile from the tile in the original image.
                        // generate a candidate tile object.
                        std::unique_ptr<Image<FloatPixel>> otherBlock(new Image<FloatPixel>(tile_size, tile_size));
                        int rowIdx = 0;
                        int colIdx = 0;
                        for (int h = k; h < k+16; h++) {
                            colIdx = 0;
                            for (int n = m; n < m+16; n++) {
                                (*otherBlock)(rowIdx, colIdx).i = (*gaussianPyramids[1][level])(h, n).i;
                                colIdx++;
                            }
                            rowIdx++;
                        }
                        // tile obj is populated.
                        std::unique_ptr<Image<FloatPixel>> difference = CameraPipeline::subImageChannel(currentBlock.get(), otherBlock.get());
                        // L2 difference
                        float diff = 0.0;
                        for (int r = 0; r < tile_size; r++) {
                            for (int c = 0; c < tile_size; c++) {
                                diff += (*difference)(r, c).i * (*difference)(r, c).i;
                            }
                        }
                        difference.reset();
                        otherBlock.reset();
                        if (diff < bestDiff) {
                            bestDiff = diff;
                            bestVerticalOffset = i - k;
                            bestHorizontalOffset = j - m;
                        }
                    }
                }
                if (bestVerticalOffset != -1) {
//                    std::cout << "bestDiff:" << bestDiff << std::endl;
                    std::vector<int> offset = {bestVerticalOffset, bestHorizontalOffset, i, i+tile_size, j, j+tile_size};
                    if (level == 0) {
                        final_offsets.push_back(offset);
                    }
                    offsets.push_back(offset);
                }
                currentBlock.reset();
            }
        }
        std::cout << "Offset size:" << offsets.size() << std::endl;
        // we now have offsets. We want to apply these as starting points for the next level.
    }
    // merging step
    std::cout << "Begin merge" << std::endl;
    std::unique_ptr<CameraSensorData<T>> mergedData(new CameraSensorData<T>(width, height));
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            mergedData->data(row, col) = raw_data_list[0]->data(row, col);
        }
    }
    std::cout << "Copied into a new camsensordata" << std::endl;
    raw_data_list = sensor_->GetBurstSensorData(0, 0, width, height);
    for (int i = 0; i < raw_data_list[0]->height(); i+=tile_size) {
        for (int j = 0; j < raw_data_list[0]->width(); j+=tile_size) {
            std::unique_ptr<Image<FloatPixel>> currentBlock(new Image<FloatPixel>(tile_size, tile_size));
            int rowIdx = 0;
            int colIdx = 0;
            for (int row = i; row < i+16; row++) {
                colIdx = 0;
                for (int col = j; col < j+16; col++) {
                    (*currentBlock)(rowIdx, colIdx).i = raw_data_list[0]->data(row, col);
                    colIdx++;
                }
                rowIdx++;
            }
//            std::cout << i << "," << j << std::endl;
            // for this tile, what's a relevant offset?
            int scaledOffsetHeight = i/2;
            int scaledOffsetWidth = j/2;
            for (std::vector<int> offset : final_offsets) {
                if ((scaledOffsetHeight > offset[2]) && (scaledOffsetHeight < offset[3]) && (scaledOffsetWidth > offset[4]) && (scaledOffsetWidth < offset[5])) {
//                    std::cout << "successful offset pass" << std::endl;
                    int verticalOffset = offset[0] * 2;
                    int horizontalOffset = offset[1] * 2;
                    if ((i + verticalOffset > 0) && (i + verticalOffset < raw_data_list[0]->height()) && (j + verticalOffset > 0) && (j + verticalOffset > 0) && (j + verticalOffset < raw_data_list[0]->height())) {
//                        std::cout << "satisfied bounds" << std::endl;
                        std::unique_ptr<Image<FloatPixel>> otherBlock(new Image<FloatPixel>(tile_size, tile_size));
                        int rowIdx = 0;
                        int colIdx = 0;
                        for (int h = i; h < i+16; h++) {
                            colIdx = 0;
                            for (int n = j; n < j+16; n++) {
                                (*otherBlock)(rowIdx, colIdx).i = raw_data_list[1]->data(h, n);
                                colIdx++;
                            }
                            rowIdx++;
                        }
                        std::unique_ptr<Image<FloatPixel>> difference = CameraPipeline::subImageChannel(currentBlock.get(), otherBlock.get());
                        // square differences
                        float diff = 0.0;
                        for (int r = 0; r < tile_size; r++) {
                            for (int c = 0; c < tile_size; c++) {
                                diff += (*difference)(r, c).i * (*difference)(r, c).i;
                            }
                        }
//                        std::cout << "difference computed" << std::endl;
                        difference.reset();
                        float weight = 0;
                        if (diff > 0.04) {
                            std::cout << "best offset is shit" << std::endl;
                            weight = 0;
                        } else if (diff < 0.01) {
                            weight = 1;
                        } else {
                            // linear ramp between these two.
                            weight = 1 + ((diff - 0.01) / (-0.03));
                        }
//                        std::cout << weight << std::endl;
//                        std::cout << "weight computed" << std::endl;
                        rowIdx = 0;
                        for (int h = i; h < i + 16; h++) {
                            colIdx = 0;
                            for (int n = j; n < j + 16; n++) {
                                mergedData->data(h, n) += weight * (*otherBlock)(rowIdx, colIdx).i;
                                colIdx++;
                            }
                            rowIdx++;
                        }
//                        std::cout << "merged in" << std::endl;
                        otherBlock.reset();
                    }
                    break;
                }
            }
            //done with a block
            currentBlock.reset();
        }
    }
    // done with all blocks
    auto& raw_data = mergedData;
//    auto raw_data = sensor_->GetSensorData(0, 0, width, height);
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
    // convert to YUV space
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            auto& pixel = (*image)(i, j);
            Float3Pixel newPixel = Float3Pixel::RgbToYuv(pixel);
            pixel.y = newPixel.y;
            pixel.u = newPixel.u;
            pixel.v = newPixel.v;
        }
    }
    // low-pass filter on image
    std::unique_ptr<Image<RgbPixel>> blurredImage = CameraPipeline::GaussianBlur(image.get());
    // take the Cb and Cr blurred levels and put those in.
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            auto& blurPixel = (*blurredImage)(i, j);
            auto& pixel = (*image)(i, j);
            pixel.y = pixel.y;
            pixel.u = blurPixel.u;
            pixel.v = blurPixel.v;
        }
    }
    // convert the image back to RGB space.
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            auto& pixel = (*image)(i, j);
            Float3Pixel newPixel = Float3Pixel::YuvToRgb(pixel);
            pixel.r = newPixel.r;
            pixel.g = newPixel.g;
            pixel.b = newPixel.b;
        }
    }
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
    std::unique_ptr<Image<RgbPixel>> finalImage = exposureFusion(image.get(), 2);
    // map to 255 for 8-bit output.
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            auto& pixel = (*finalImage)(row, col);
            // pixel data from the sensor is normalized to the 0-1 range, so
            // scale by 255 for the final image output.  Output image pixels
            // should be in the 0-255 range.
            pixel.r = pixel.r * 255.f;
            pixel.g = pixel.g * 255.f;
            pixel.b = pixel.b * 255.f;
        }
    }
    // return processed image output
    return finalImage;
    
    // END: CS348K STUDENTS MODIFY THIS CODE
}
