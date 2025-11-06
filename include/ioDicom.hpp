#pragma once
#include <string>
#include <vector>
#include "itkImage.h"

using PerfImageType = itk::Image<unsigned short, 3>;

std::vector<std::string> GetDicomFileNames(const std::string& dicomDir);

PerfImageType::Pointer ReadPerfusionSeries(const std::vector<std::string>& fileNames);

std::vector<double> ReadDicomTimestamps(const std::vector<std::string>& fileNames);

double ParseDicomTimeToSeconds(const std::string& tm);
