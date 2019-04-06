#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/grad.h>
#include <igl/jet.h>
#include <igl/readDMAT.h>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/per_face_normals.h>
#include <igl/segments_to_arrows.h>

#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>


#include <iostream>

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd V_uv;

//********************** STRINGIFY
//http://www.decompile.com/cpp/faq/file_and_line_error_string.htm
#ifndef STRINGIFY
#define STRINGIFY(x) TOSTRING(x)
#define TOSTRING(x) #x
#endif
//*****************************************************************************

int main(int argc, char *argv[])
{
    
    
    bool add_lines = false;
    bool colorize_3d_faces = true;
    bool parametrize = false;

    std::string filepath = STRINGIFY(INPUT_DIR)"/camelhead.off";
    //std::string filepath = STRINGIFY(INPUT_DIR)"/sphere1.obj";

    // Load a mesh in OBJ format    
    //igl::readOBJ(filepath, V, F);
    
    // Load a mesh in OFF format    
    igl::readOFF(filepath, V, F);

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);

    //compute face normal
    Eigen::MatrixXd N_faces;
    igl::per_face_normals(V, F, N_faces);

    if (parametrize) {
        // Find the open boundary
        Eigen::VectorXi bnd;
        igl::boundary_loop(F, bnd);
        // Map the boundary to a circle, preserving edge proportions
        Eigen::MatrixXd bnd_uv;
        igl::map_vertices_to_circle(V, bnd, bnd_uv);
        // Harmonic parametrization for the internal vertices
        igl::harmonic(V, F, bnd, bnd_uv, 1, V_uv);
        // Scale UV to make the texture more clear
        //V_uv *= 5;

        viewer.data().set_mesh(V_uv, F);
        viewer.data().compute_normals();
    }

    if (add_lines) {
        // Average edge length divided by max face normal (for scaling)    
        const double max_size = igl::avg_edge_length(V, F) / N_faces.maxCoeff();
        // Draw a black segment in direction of gradient at face barycenters
        Eigen::MatrixXd BC;
        igl::barycenter(V, F, BC);
        const Eigen::RowVector3d black(0, 0, 0);
        viewer.data().add_edges(BC, BC + max_size*N_faces, black);
    }

    if (colorize_3d_faces) {
        //scale the normal and draw the normal component (x,y,z) as r,g,b for 
        //each face 
        double max_coeff(N_faces.maxCoeff()), min_coeff(N_faces.minCoeff());
        double maxmin_diff = max_coeff - min_coeff;
        for (int i = 0; i < N_faces.rows(); i++) {
            //scale face from [-1,1] to [0,2]
            for (int j = 0; j < 3; j++) {
                N_faces(i, j) -= min_coeff;//from [-1,1] to [0,2]
                N_faces(i, j) = N_faces(i, j) / maxmin_diff; //from [0,2]to[0,1]
                // N_faces(i, j) *= 255; //from [0,1] to [0,255]
            }
        }
        // Compute pseudocolor for original function                               
        Eigen::MatrixXd C;
        viewer.data().set_colors(N_faces);
    }



    // Hide wireframe
    viewer.data().show_lines = false;

    viewer.launch();
}
