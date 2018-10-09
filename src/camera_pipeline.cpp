#include "camera_pipeline.hpp"
#include <iostream>
#include <algorithm>
#include <cstdlib>

std::unique_ptr<Image<RgbPixel>> CameraPipeline::GaussianBlur(std::unique_ptr<Image<RgbPixel>> &originalImage) const {
    int height = originalImage->height();
    int width = originalImage->width();
    std::cout << "Obtained the height and width" << std::endl;
    std::unique_ptr<Image<RgbPixel>> blurredImage(new Image<RgbPixel>(width, height));
    std::cout << "Initialize blur img" << std::endl;
    // copy into float array
    float** new_img_r = new float*[height + 2];
    float** new_img_g = new float*[height + 2];
    float** new_img_b = new float*[height + 2];
    for (int i = 0; i < height + 2; i++) {
        new_img_r[i] = new float[width+2];
        new_img_g[i] = new float[width+2];
        new_img_b[i] = new float[width+2];
    }
    std::cout << "Initialized the arrays" << std::endl;
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            auto& pixel = (*originalImage)(row, col);
            new_img_r[row + 1][col + 1] = pixel.r;
            new_img_g[row + 1][col + 1] = pixel.g;
            new_img_b[row + 1][col + 1] = pixel.b;
        }
    }
    //    std::cout << "Performed copy." << std::endl;
    // mirror the edges
    for (int z = 0; z < height + 2; z++) {
        new_img_r[z][0] = new_img_r[z][1];
        new_img_r[z][width+1] = new_img_r[z][width];
        new_img_r[0][z] = new_img_r[1][z];
        new_img_r[height+1][z] = new_img_r[width][z];
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

std::unique_ptr<Image<RgbPixel>> CameraPipeline::DownSample(std::unique_ptr<Image<RgbPixel>> &originalImage) const {
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
    std::unique_ptr<Image<RgbPixel>> blurredImage = CameraPipeline::GaussianBlur(image);
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
    // apply some gain because brighter images tend to look nicer.
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            auto& pixel = (*image)(i, j);
            pixel.r = std::min(pixel.r + 0.15 * pixel.r, 1.0);
            pixel.g = std::min(pixel.g + 0.15 * pixel.g, 1.0);
            pixel.b = std::min(pixel.b + 0.15 * pixel.b, 1.0);
        }
    }
    // gamma correct to use more bit space.
    image->GammaCorrect(0.41);
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
