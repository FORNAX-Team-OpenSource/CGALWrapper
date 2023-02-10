////includes for skeletonization
//#include <CGAL/Simple_cartesian.h>
//#include <CGAL/Polyhedron_3.h>
//#include <CGAL/extract_mean_curvature_flow_skeleton.h>
//#include <CGAL/boost/graph/split_graph_into_polylines.h>
//#include <fstream>
////includes for triangulation feature
//#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
//
//typedef CGAL::Simple_cartesian<double>                        Kernel;
//typedef Kernel::Point_3                                       Point;
//typedef CGAL::Polyhedron_3<Kernel>                            Polyhedron;
//typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
//typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron> Skeletonization;
//typedef Skeletonization::Skeleton                             Skeleton;
//typedef Skeleton::vertex_descriptor                           Skeleton_vertex;
//typedef Skeleton::edge_descriptor                             Skeleton_edge;
//
//std::vector<std::vector<Point>> SkelAllPts; //equivalent to List<List<Point3D>>
//std::vector<Point> SkelSegment; // temporary storage for a segment of the skeleton
//const char* comms1 = "IF EXIST *cline*.json DEL *cline*.json /F /Q";
//struct Display_polylines{//groups the skeleton relevant segments together in the Vertex vector of the skeleton
//  const Skeleton& skeleton;
//  //std::ofstream& out;
//  Display_polylines(const Skeleton& skeleton)//, std::ofstream& out)
//      : skeleton(skeleton) {}//, out(out){}
//  void start_new_polyline(){}
//  void add_node(Skeleton_vertex v){
//      SkelSegment.push_back(skeleton[v].point);
//  }
//  void end_polyline()
//  {
//          SkelAllPts.push_back(SkelSegment);
//          SkelSegment.clear();
//  }
//};
//
//
//// This example extracts a medially centered skeleton from a given mesh.
//int main51(int argc, char* argv[])
//{
//    std::ifstream input(argv[1]);
//    //delete existing cline json files
//    system(comms1);
//    //get base filename of input model
//    std::string modelname="";
//    int i = strlen(argv[1]);
//    char* tempname = &argv[1][i - 1];
//    while (argv[1][i-1] != '\\' && argv[1][i-1]) {
//        i--; tempname--;
//    }
//    ++tempname; i = strlen(tempname);
//    modelname = tempname;
//    modelname.resize(i-4);
//    //process input polyhedron
//    Polyhedron tmesh;
//    input >> tmesh;
//    if (!CGAL::is_triangle_mesh(tmesh))
//    {
//        std::cout << "Input geometry is not triangulated." << std::endl;
//        if (CGAL::Polygon_mesh_processing::triangulate_faces(tmesh))
//            std::cout << "...triangulated the input...\n";
//        else{
//            std::cout << "...could not triangulate...\n";
//            return EXIT_FAILURE;
//        }
//    }
//    Skeleton skeleton;
//    CGAL::extract_mean_curvature_flow_skeleton(tmesh, skeleton);
//    std::cout << "Number of vertices of the skeleton: " << boost::num_vertices(skeleton) << "\n";
//    std::cout << "Number of edges of the skeleton: " << boost::num_edges(skeleton) << "\n";
//    //breaks down the skeleton to segments (List<Point3D> in BIMRL)
//    Display_polylines display(skeleton);//,output);
//    CGAL::split_graph_into_polylines(skeleton, display);
////output json files List<Point3D>
//    for (std::vector<Point> seg : SkelAllPts) {
//        static int j = 0; j++;
//        //char jsonx[100] = "";
//        std::ofstream output(modelname+".cline"+std::to_string(j)+".json");
//        output << "[";
//        for (int i = 0; i < seg.size();i++) {
//            output << "\n\t{\"X\":"<<seg[i].x() << ",\"Y\":"<< seg[i].y() << ",\"Z\":" << seg[i].z() <<"}";
//            if(i<seg.size()-1)
//                output << ',';
//        }
//        output << "\n]";
//        output.close();
//    }
//    return EXIT_SUCCESS;
//}
//
