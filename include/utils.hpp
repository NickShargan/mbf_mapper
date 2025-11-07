#pragma once
#include <string>
#include <Eigen/Dense>

void SaveMbfMapAsPng(const Eigen::MatrixXd& mbfMap, const std::string& path);
void SaveMbfMapAsNifti(const Eigen::MatrixXd& mbfMap, const std::string& path);

void WriteVectorToCSV(const std::vector<double>& vals,
                      const std::vector<double>& timesSec,
                      const std::string& outPath,
                      const std::string& valColName);

void zeroLeveling(std::vector<double>& signal);
