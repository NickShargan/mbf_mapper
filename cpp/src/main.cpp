#include <iostream>

#include "itkImage.h"

#include "utils.hpp"
#include "ioDicom.hpp"
#include "ioMasks.hpp"
#include "mbfPipeline.hpp"


int main(int argc, char* argv[])
{
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <dicom_dir> <mask.tiff> <output.nii>\n";
        return EXIT_FAILURE;
    }

    const std::string dicomDir = argv[1];
    const std::string tiffPath = argv[2];
    const std::string outNii   = argv[3];

    try {
        auto fileNames = GetDicomFileNames(dicomDir);
        auto perf      = ReadPerfusionSeries(fileNames);
        auto times     = ReadDicomTimestamps(fileNames);
        auto [maskLV, maskMyo] = ReadTMaskTiff(tiffPath);

        auto aifCurve = EstimateAIF(perf, maskLV);
        WriteVectorToCSV(aifCurve, times, "aif.csv", "AIF");

        auto mbfMap = ComputeMbfMap(perf, times, maskMyo, aifCurve);
        
        SaveMbfMapAsNifti(mbfMap, outNii);
    }
    catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

// std::vector<double> ExtractTimeCurveAtPixel(
//     PerfImageType::Pointer perf,
//     unsigned int x,
//     unsigned int y)
// {
//     auto region = perf->GetLargestPossibleRegion();
//     auto size   = region.GetSize();

//     if (x >= size[0] || y >= size[1])
//         throw std::out_of_range("Pixel coordinates out of image bounds");

//     const unsigned int T = size[2];
//     std::vector<double> curve(T, 0.0);

//     for (unsigned int t = 0; t < T; ++t)
//     {
//         PerfImageType::IndexType idx;
//         idx[0] = x;
//         idx[1] = y;
//         idx[2] = t;
//         auto val = (perf->GetPixel(idx));
//         curve[t] = val;
//     }

//     return curve;
// }
