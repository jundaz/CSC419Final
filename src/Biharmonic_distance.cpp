#include "biharmonic_distance.h"
#include "Eigen/SparseCholesky"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"
#include "../spectra/include/Spectra/MatOp/SparseGenMatProd.h"
#include "../spectra/include/Spectra/SymGEigsSolver.h"
#include "../spectra/include/Spectra/MatOp/SparseCholesky.h"
#include <iostream>
void biharmonic_distance(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const int approach,
    const int k,
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
    A_inv.setFromTriplets(trip.begin(), trip.end());
    // exact distance
    if(approach == 0){
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
    }
    else{
        // approximate distance
        // solve generalized eigen value problem
        Eigen::VectorXd x_minus_y, x_minus_y_sqr;
        D.setZero();
        double dist;
        int this_k = k;
        if(k > V.rows()){
            this_k = V.rows();
        }
        Spectra::SparseGenMatProd<double> op(L);
        Spectra::SparseCholesky<double> Bop(A);
        Spectra::SymGEigsSolver<double, Spectra::LARGEST_ALGE, Spectra::SparseGenMatProd<double>, Spectra::SparseCholesky<double>, Spectra::GEIGS_CHOLESKY>geigs(&op, &Bop, this_k-1, this_k * 20);
        geigs.init();
        int nconv = geigs.compute();
        if(geigs.info() == Spectra::SUCCESSFUL){
            // remove first eigenvalue and its eigenvector because its too small
            Eigen::MatrixXd selected_vec = geigs.eigenvectors().rightCols(k-2);
            Eigen::VectorXd selected_val = geigs.eigenvalues().tail(k-2);
            Eigen::VectorXd eig_val_sqr = selected_val.array().pow(2);
            for(int i = 0; i < selected_vec.rows() - 1; i++){
                for(int j = i + 1; j < selected_vec.rows(); j++){
                    x_minus_y = selected_vec.row(i) - selected_vec.row(j);
                    x_minus_y_sqr = x_minus_y.array().pow(2);
                    dist = x_minus_y_sqr.cwiseQuotient(eig_val_sqr).sum();
                    D(i, j) = D(j, i) = sqrt(dist);
                }
            }
        }
        
    }
    
}
