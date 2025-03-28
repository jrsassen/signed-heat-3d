#include "geometrycentral/pointcloud/point_position_normal_geometry.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "signed_heat_grid_solver_nps.h"

#include "args/args.hxx"

#include <chrono>

#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;
std::chrono::time_point<high_resolution_clock> t1, t2;
std::chrono::duration<double, std::milli> ms_fp;

using namespace geometrycentral;
using namespace geometrycentral::surface;

int main(int argc, char** argv) {

    // Configure the argument parser
    args::ArgumentParser parser("Solve for generalized signed distance (3D domains).");
    args::HelpFlag help(parser, "help", "Display this help menu", {"help"});
    args::Positional<std::string> meshFilename(parser, "mesh", "A mesh or point cloud file.");
    args::ValueFlag<double> hCoef(parser, "h", "Controls the tet/grid spacing proportional to $2^{-h}$.",
                                  {"h", "hCoef"});

    args::Group group(parser);
    args::Flag grid(group, "grid", "Solve on a background grid (vs. tet mesh).", {"g", "grid"});
    args::Flag fast(group, "fast", "Solve using a less accurate, but significantly faster, method of integration.",
                    {"f", "fast"});
    args::Flag verbose(group, "verbose", "Verbose output", {"V", "verbose"});
    args::Flag headless(group, "headless", "Don't use the GUI.", {"l", "headless"});

    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help&) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    float TCOEF = 1.0;
    float HCOEF = 0.0;
    std::unique_ptr<SignedHeatGridSolver> gridSolver;
    SignedHeat3DOptions SHM_OPTIONS;
    int CONSTRAINT_MODE = static_cast<int>(LevelSetConstraint::ZeroSet);

    // Program variables
    enum MeshMode { Tet = 0, Grid };
    enum InputMode { Mesh = 0, Points };
    int MESH_MODE = MeshMode::Tet;
    int INPUT_MODE = InputMode::Mesh;
    std::string MESHNAME = "input mesh";
    std::string OUTPUT_DIR = "../export";
    std::string OUTPUT_FILENAME;
    int LAST_SOLVER_MODE;
    bool VERBOSE = true;
    bool HEADLESS;
    bool CONTOURED = false;

    std::unique_ptr<SurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geometry;

    Vector<double> PHI;

    if (!meshFilename) {
        std::cerr << "Please specify a mesh file as argument." << std::endl;
        return EXIT_FAILURE;
    }
    if (hCoef) {
        HCOEF = args::get(hCoef);
    }
    HCOEF = 2.;

    // Load mesh
    std::string meshFilepath = args::get(meshFilename);
    MESH_MODE = grid ? MeshMode::Grid : MeshMode::Tet;
    OUTPUT_FILENAME = OUTPUT_DIR + "/GSD.obj";
    HEADLESS = headless;
    SHM_OPTIONS.exportData = true; // always true if in headless mode
    SHM_OPTIONS.meshname = "bunny";
    VERBOSE = true;


    std::tie(mesh, geometry) = readSurfaceMesh(meshFilepath);
    INPUT_MODE = InputMode::Mesh;

    gridSolver = std::unique_ptr<SignedHeatGridSolver>(new SignedHeatGridSolver());
    gridSolver->VERBOSE = true;

    SHM_OPTIONS.levelSetConstraint = static_cast<LevelSetConstraint>(CONSTRAINT_MODE);
    SHM_OPTIONS.tCoef = TCOEF;
    SHM_OPTIONS.hCoef = HCOEF;

    auto t1 = high_resolution_clock::now();
    PHI = gridSolver->computeDistance(*geometry, SHM_OPTIONS);
    auto t2 = high_resolution_clock::now();
    ms_fp = t2 - t1;
    if (VERBOSE) std::cerr << "Solve time (s): " << ms_fp.count() / 1000. << std::endl;




    LAST_SOLVER_MODE = MESH_MODE;
    SHM_OPTIONS.rebuild = false;

    return EXIT_SUCCESS;
}