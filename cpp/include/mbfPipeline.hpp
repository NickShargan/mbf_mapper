#pragma once
#include <vector>
#include <string>
#include <Eigen/Dense>
#include "itkImage.h"

using PerfImageType = itk::Image<unsigned short, 3>;
using Mask2DType = itk::Image<unsigned char, 2>;

Eigen::MatrixXd ComputeMbfMap(
    PerfImageType::Pointer perf,
    const std::vector<double>& times,
    Mask2DType::Pointer myoMask,
    const std::vector<double>& aifCurve
);

std::vector<double> EstimateAIF(
    PerfImageType::Pointer perf,
    Mask2DType::Pointer lvMask
);
