#include "biharmonic_distance.h"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/read_triangle_mesh.h>
#include <Eigen/Core>
#include <iostream>
#include <cstdlib>

// adapted from framework of A5 and tutorial 206 geodesic distance demo.
int main(int argc, char *argv[]) {
    // load V and F
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::opengl::glfw::Viewer viewer;
    igl::read_triangle_mesh(argc > 1 ? argv[1] : "../data/bunny.off", V, F);

    // construct pair distances just for once
    Eigen::MatrixXd D;
    // input 0 for exact dist, 1 for approx dist
    biharmonic_distance(V, F, atoi(argv[2]), atoi(argv[3]), D);

    // bind mouse down call back
    viewer.callback_mouse_down =
        [&](igl::opengl::glfw::Viewer& viewer, int, int)->bool
    {
        // corresponds user's mouse click point to vertex by estimation
        int fid;
        Eigen::Vector3f bc;
        double x_position = viewer.current_mouse_x;
        double y_position = viewer.core().viewport(3) - viewer.current_mouse_y;
        if (igl::unproject_onto_mesh(
            Eigen::Vector2f(x_position, y_position),
            viewer.core().view,
            viewer.core().proj,
            viewer.core().viewport,
            V,
            F,
            fid,
            bc
        )) {
            int max;
            bc.maxCoeff(&max);
            int vid = F(fid, max);

            // retrieve vertex vid's pair distances and use them to color the viewer
            viewer.data().set_data(D.row(vid));
            return true;
        }
        return false;
    };

    viewer.data().set_mesh(V,F);
    viewer.data().show_lines = false;
    return viewer.launch();
}
