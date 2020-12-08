#include "biharmonic_distance.h"
#include "Eigen/SparseCholesky"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"
#include "igl/eigs.h"
#include <iostream>
void biharmonic_distance(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    Eigen::MatrixXd &D){
    D.resize(V.rows(), V.rows());
    Eigen::MatrixXd dblA(F.rows(), 1);
    // construct cot matrix
    Eigen::SparseMatrix<double> L, A, A_inv, new_A;
    igl::cotmatrix(V, F, L);
    L *= 1.0;
    // construct A the mass matrix
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, A);
    // using voronoi somehow give incorrect result but the paper says use voronoi
    // very confusing...
    A_inv.resize(A.rows(), A.rows());
    std::vector<Eigen::Triplet<double>> trip;
    for(int i = 0; i < A.rows(); i++){
        trip.push_back(Eigen::Triplet<double>(i, i, 1.0/A.diagonal()[i]));
    }
    // exact distance
    A_inv.setFromTriplets(trip.begin(), trip.end());
    new_A = L.transpose() * A_inv * L;
    new_A.row(0) *= 0;
    new_A.col(0) *= 0;
    new_A.coeffRef(0, 0) += 1.0;
    Eigen::MatrixXd J(V.rows(), V.rows());
    J.setOnes();
    J *= 1.0/double(V.rows());
    J = Eigen::MatrixXd::Identity(V.rows(), V.rows()) - J;
    J.row(0).setZero();
    Eigen::MatrixXd G;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(new_A);
    G = solver.solve(J);
    G -= (1.0 / V.rows()) * G.colwise().sum().replicate(V.rows(), 1);
    Eigen::VectorXd V1 = G.diagonal();
    Eigen::MatrixXd m = V1 * Eigen::VectorXd::Ones(V.rows()).transpose();
    D = sqrt((m + m.transpose() - 2*G).array());
    // approximate distance
    // solve generalized eigen value problem
    

}
