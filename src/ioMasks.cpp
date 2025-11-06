#include "ioMasks.hpp"

#include "itkImageFileReader.h"
#include "itkExtractImageFilter.h"

#include <stdexcept>

std::pair<Mask2DType::Pointer, Mask2DType::Pointer>
ReadTMaskTiff(const std::string &tiffPath)
{
    auto reader = itk::ImageFileReader<Mask3DType>::New();
    reader->SetFileName(tiffPath);
    reader->Update();

    auto vol = reader->GetOutput();
    auto region = vol->GetLargestPossibleRegion();
    auto size = region.GetSize();

    if (size[2] < 2)
        throw std::runtime_error("Expected at least 2 pages in TIFF file.");

    auto extractSlice = [&](unsigned int zIndex) -> Mask2DType::Pointer {
        using ExtractFilterType = itk::ExtractImageFilter<Mask3DType, Mask2DType>;
        auto extract = ExtractFilterType::New();

        Mask3DType::RegionType sliceRegion = region;
        Mask3DType::SizeType sliceSize = size;
        sliceSize[2] = 0;
        Mask3DType::IndexType sliceStart = region.GetIndex();
        sliceStart[2] = zIndex;

        sliceRegion.SetSize(sliceSize);
        sliceRegion.SetIndex(sliceStart);

        extract->SetExtractionRegion(sliceRegion);
        extract->SetInput(vol);
        extract->SetDirectionCollapseToSubmatrix();
        extract->Update();
        return extract->GetOutput();
    };

    return {extractSlice(0), extractSlice(1)};
}
