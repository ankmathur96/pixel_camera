#pragma once

#include <algorithm>
#include <limits>
#include <memory>
#include "camera_pipeline_interface.hpp"
#include "image.hpp"
#include "pixel.hpp"

class CameraPipeline : public CameraPipelineInterface {
 public:
    
  explicit CameraPipeline(CameraSensor* sensor)
    : CameraPipelineInterface(sensor) {}
    
 private:
  using T = typename CameraSensor::T;
  using CameraPipelineInterface::sensor_;

  std::unique_ptr<Image<RgbPixel>> ProcessShot() const override;
  // BEGIN: CS348K STUDENTS MODIFY THIS CODE
  //
  // You can add any necessary private member variables or functions.
  //
  float** GaussianChannelBlur(float** channel, int height, int width) const;
  std::unique_ptr<Image<RgbPixel>> GaussianBlur(Image<RgbPixel>* originalImage) const;
  std::unique_ptr<Image<RgbPixel>> DownSample(std::unique_ptr<Image<RgbPixel>> &originalImage) const;
  float** UpSampleChannel(float** originalImage, int height, int width) const;
  std::unique_ptr<Image<RgbPixel>> UpSample(std::unique_ptr<Image<RgbPixel>> &originalImage) const;
  std::unique_ptr<Image<RgbPixel>> subImage(Image<RgbPixel>* im1, Image<RgbPixel>* im2) const;
  float** subImageChannel(float** im1, float** im2, int height, int width) const;
  std::unique_ptr<Image<RgbPixel>> exposureFusion(std::unique_ptr<Image<RgbPixel>> &originalImage, int levels) const;
  void Denoise(std::unique_ptr<Image<RgbPixel>> &originalImage) const;
  // END: CS348K STUDENTS MODIFY THIS CODE  
};
