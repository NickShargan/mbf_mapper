#include "ioDicom.hpp"
#include "itkGDCMSeriesFileNames.h"
#include "itkGDCMImageIO.h"
#include "itkImageSeriesReader.h"
#include "itkMetaDataObject.h"
#include <stdexcept>

double ParseDicomTimeToSeconds(const std::string& tm) {
    if (tm.size() < 6) return 0.0;
    int hh = std::stoi(tm.substr(0, 2));
    int mm = std::stoi(tm.substr(2, 2));
    int ss = std::stoi(tm.substr(4, 2));
    double sec = hh * 3600 + mm * 60 + ss;

    if (tm.size() > 6 && tm[6] == '.') {
        std::string frac = tm.substr(7);
        double f = std::stod("0." + frac);
        sec += f;
    }
    return sec;
}

std::vector<std::string>
GetDicomFileNames(const std::string& dicomDir)
{
    auto nameGenerator = itk::GDCMSeriesFileNames::New();
    nameGenerator->SetDirectory(dicomDir);
    const auto seriesUIDs = nameGenerator->GetSeriesUIDs();
    if (seriesUIDs.empty()) {
        throw std::runtime_error("No DICOM series found");
    }
    return nameGenerator->GetFileNames(seriesUIDs.front());
}

PerfImageType::Pointer
ReadPerfusionSeries(const std::vector<std::string>& fileNames)
{
    auto gdcmIO = itk::GDCMImageIO::New();
    auto reader = itk::ImageSeriesReader<PerfImageType>::New();
    reader->SetImageIO(gdcmIO);
    reader->SetFileNames(fileNames);
    reader->Update();
    return reader->GetOutput();
}

std::vector<double>
ReadDicomTimestamps(const std::vector<std::string>& fileNames)
{
    auto gdcmIO = itk::GDCMImageIO::New();
    std::vector<double> ts;
    ts.reserve(fileNames.size());

    for (const auto& f : fileNames) {
        gdcmIO->SetFileName(f);
        gdcmIO->ReadImageInformation();
        std::string timeStr;
        if (itk::ExposeMetaData<std::string>(gdcmIO->GetMetaDataDictionary(),
                                             "0008|0032", timeStr)) {
            ts.push_back(ParseDicomTimeToSeconds(timeStr));
        } else {
            ts.push_back(0.0);
        }
    }

    if (!ts.empty()) {
        const double t0 = ts.front();
        for (auto& t : ts) t -= t0;
    }

    return ts;
}
