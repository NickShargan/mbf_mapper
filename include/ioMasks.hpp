#pragma once
#include <string>
#include <utility>

#include "itkImage.h"

using Mask2DType = itk::Image<unsigned char, 2>;
using Mask3DType = itk::Image<unsigned char, 3>;

std::pair<Mask2DType::Pointer, Mask2DType::Pointer>
ReadTMaskTiff(const std::string &tiffPath);
