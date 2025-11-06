#include <iostream>

#include "mbfPipeline.hpp"
#include "MbfDeconvolver.hpp"
#include "utils.hpp"

Eigen::MatrixXd ComputeMbfMap(
    PerfImageType::Pointer perf,
    const std::vector<double>& times,
    Mask2DType::Pointer myoMask,
    const std::vector<double>& aifCurve
)
{
    auto region = perf->GetLargestPossibleRegion();
    auto size   = region.GetSize();
    const unsigned int W = size[0];
    const unsigned int H = size[1];

    Eigen::MatrixXd mbfMap = Eigen::MatrixXd::Zero(H, W);

    MbfDeconvolver deconv(aifCurve, 0.02);
    // todo: mode rho to config.hpp
    const double rho = 1.05;
    const double dt  = times.size() > 1 ? (times[1] - times[0]) : 1.0;

    auto isMasked = [&](int x, int y) {
        Mask2DType::IndexType idx = { {x, y} };
        return myoMask->GetPixel(idx) != 0;
    };

    bool isTestSaved = false;

    for (int y = 0; y < static_cast<int>(H); ++y) {
        for (int x = 0; x < static_cast<int>(W); ++x) {
            if (!isMasked(x, y)) continue;

            // extract time curve
            std::vector<double> curve(times.size());
            for (std::size_t t = 0; t < times.size(); ++t) {
                PerfImageType::IndexType idx;
                idx[0] = x; idx[1] = y; idx[2] = t;
                curve[t] = perf->GetPixel(idx);
            }

            double h0 = deconv.computeMBF(curve);
            double mbf = (h0 / dt) * 60.0 / rho;
            mbfMap(y, x) = mbf;

            // std::cout << mbf << "\n";

            // save test
            if (!isTestSaved) {
                auto h = deconv.deconvolve(curve);

                std::vector<double> h_vec(h.data(), h.data() + h.size());

                WriteVectorToCSV(curve, times, "myo_pix.csv", "MYO_pix");
                WriteVectorToCSV(h_vec, times, "h_t.csv", "h_t_svd");

                isTestSaved = true;
            }
        }
    }

    return mbfMap;
}

std::vector<double>
EstimateAIF(PerfImageType::Pointer perf, Mask2DType::Pointer lvMask)
{
    // perf: size[0]=X, size[1]=Y, size[2]=T
    auto perfRegion = perf->GetLargestPossibleRegion();
    auto perfSize   = perfRegion.GetSize();

    auto maskRegion = lvMask->GetLargestPossibleRegion();
    auto maskSize   = maskRegion.GetSize();

    // basic size check on x,y
    if (perfSize[0] != maskSize[0] || perfSize[1] != maskSize[1]) {
        throw std::runtime_error("EstimateAIF: perfusion XY size and mask XY size do not match.");
    }

    const unsigned int X = perfSize[0];
    const unsigned int Y = perfSize[1];
    const unsigned int T = perfSize[2];

    std::cout << "numSlices: " << T << std::endl;

    std::vector<double> aif(T, 0.0);

    for (unsigned int t = 0; t < T; ++t)
    {
        double sum = 0.0;
        std::size_t count = 0;

        for (unsigned int y = 0; y < Y; ++y)
        {
            for (unsigned int x = 0; x < X; ++x)
            {
                Mask2DType::IndexType mIdx;
                mIdx[0] = x;
                mIdx[1] = y;
                auto m = lvMask->GetPixel(mIdx);
                if (m == 0)
                    continue;

                PerfImageType::IndexType pIdx;
                pIdx[0] = x;
                pIdx[1] = y;
                pIdx[2] = t;

                auto val = perf->GetPixel(pIdx);
                sum += static_cast<double>(val);
                ++count;
            }
        }

        if (count > 0) 
            aif[t] = sum / static_cast<double>(count);
        else
            aif[t] = 0.0;
    }

    return aif;
}
