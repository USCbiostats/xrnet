#ifndef COORD_DESC_TYPES_H
#define COORD_DESC_TYPES_H

#include <RcppEigen.h>

typedef Eigen::Map<const Eigen::MatrixXd> MapMat;
typedef Eigen::MappedSparseMatrix<double> MapSpMat;
typedef Eigen::Map<const Eigen::VectorXd> MapVec;

#endif // COORD_DESC_TYPES_H
