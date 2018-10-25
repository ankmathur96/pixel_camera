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
  std::unique_ptr<Image<FloatPixel>> GaussianBlurChannel(Image<FloatPixel>* originalImage) const;
  std::unique_ptr<Image<RgbPixel>> GaussianBlur(Image<RgbPixel>* originalImage) const;
  std::unique_ptr<Image<FloatPixel>> DownSampleChannel(Image<FloatPixel>* originalImage) const;
  std::unique_ptr<Image<RgbPixel>> DownSample(Image<RgbPixel>* originalImage) const;
  std::unique_ptr<Image<FloatPixel>> UpSampleChannel(Image<FloatPixel>* originalImage) const;
  std::unique_ptr<Image<RgbPixel>> UpSample(Image<RgbPixel>* originalImage) const;
  void Denoise(std::unique_ptr<Image<RgbPixel>> &originalImage) const;
  std::unique_ptr<Image<FloatPixel>> subImageChannel(Image<FloatPixel>* im1, Image<FloatPixel>* im2) const;
  std::unique_ptr<Image<RgbPixel>> exposureFusion(Image<RgbPixel> *originalImage, int nLevels) const;
  // END: CS348K STUDENTS MODIFY THIS CODE  
};
