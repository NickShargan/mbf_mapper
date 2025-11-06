#include "itkImage.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkPNGImageIO.h"
#include "itkNiftiImageIO.h"
#include <Eigen/Dense>

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using MbfImageType = itk::Image<float, 2>;
using OutputImageType = itk::Image<unsigned char, 2>;

void SaveMbfMapAsPng(const Eigen::MatrixXd& mbfMap, const std::string& path) {
    const int H = mbfMap.rows();
    const int W = mbfMap.cols();

    // Create ITK image
    MbfImageType::Pointer image = MbfImageType::New();
    MbfImageType::RegionType region;
    MbfImageType::IndexType start;
    start.Fill(0);
    MbfImageType::SizeType size;
    size[0] = W;
    size[1] = H;
    region.SetSize(size);
    region.SetIndex(start);
    image->SetRegions(region);
    image->Allocate();

    // Fill pixels from Eigen matrix
    MbfImageType::IndexType idx;
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            idx[0] = x;
            idx[1] = y;
            image->SetPixel(idx, static_cast<float>(mbfMap(y, x)));
        }
    }

    // Rescale to 0–255 for PNG
    using RescaleFilterType = itk::RescaleIntensityImageFilter<MbfImageType, OutputImageType>;
    auto rescale = RescaleFilterType::New();
    rescale->SetInput(image);
    rescale->SetOutputMinimum(0);
    rescale->SetOutputMaximum(255);

    // Write to PNG
    using WriterType = itk::ImageFileWriter<OutputImageType>;
    auto writer = WriterType::New();
    writer->SetFileName(path);
    writer->SetInput(rescale->GetOutput());
    writer->Update();

    std::cout << "MBF map saved to: " << path << std::endl;
}

void SaveMbfMapAsNifti(const Eigen::MatrixXd& mbfMap, const std::string& path)
{
    const int H = mbfMap.rows();
    const int W = mbfMap.cols();

    // Create ITK image
    MbfImageType::Pointer image = MbfImageType::New();
    MbfImageType::RegionType region;
    MbfImageType::IndexType start;
    start.Fill(0);
    MbfImageType::SizeType size;
    size[0] = W;
    size[1] = H;
    region.SetSize(size);
    region.SetIndex(start);
    image->SetRegions(region);
    image->Allocate();

    // Copy Eigen matrix data into ITK image
    MbfImageType::IndexType idx;
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            idx[0] = x;
            idx[1] = y;
            image->SetPixel(idx, static_cast<float>(mbfMap(y, x)));
        }
    }

    // todo: add spacing

    // Write .nii.gz file
    using WriterType = itk::ImageFileWriter<MbfImageType>;
    auto writer = WriterType::New();
    writer->SetFileName(path);
    writer->SetInput(image);
    writer->SetImageIO(itk::NiftiImageIO::New());
    writer->Update();

    std::cout << "MBF NIfTI saved to: " << path << std::endl;
}

void WriteVectorToCSV(const std::vector<double>& vals,
                      const std::vector<double>& timesSec,
                      const std::string& outPath,
                      const std::string& valColName = "value"
                      )
{
    std::ofstream file(outPath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + outPath);
    }

    // header
    file << "Time_s," << valColName << "\n";

    for (std::size_t i = 0; i < vals.size(); ++i) {
        file << timesSec[i] << "," << vals[i] << "\n";
    }
}