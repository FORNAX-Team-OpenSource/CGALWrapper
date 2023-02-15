#include <iostream>
#include<fstream>
#include<vector>
#include<string>
#include<algorithm>

#include "Polyhedrons.h"
#include "Polylines.h"
#include "Polyline.h"
#include "Point.h"

// Remesh
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Timer.h>

#include<CGAL/Simple_cartesian.h>
#include<CGAL/Polyhedron_incremental_builder_3.h>
#include<CGAL/Polyhedron_3.h>
#include<CGAL/IO/Polyhedron_iostream.h>

// Skeletonization
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <tuple>
#include <map>

#include <CGAL/facets_in_complex_3_to_triangle_mesh.h>
#include "Polylines.h"

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>


#include <CGAL/Surface_mesh.h>

#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/mesh_segmentation.h>

#include <boost/iterator/function_output_iterator.hpp>

#include <CGAL/Polygon_mesh_processing/detect_features.h>

#include <CGAL/Simple_cartesian.h>

using namespace std;


// Segmentation
typedef CGAL::Exact_predicates_inexact_constructions_kernel KernelSM;
typedef CGAL::Surface_mesh<KernelSM::Point_3> SM;
typedef boost::graph_traits<SM>::face_descriptor face_descriptor;

typedef CGAL::Simple_cartesian<double>     Kernel;
typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
typedef Polyhedron::HalfedgeDS             HalfedgeDS;

typedef Polyhedron::Facet_iterator                   Facet_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;
typedef Polyhedron::Vertex_iterator        Vertex_iterator;


// Remesh ---------------
// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Mesh_polyhedron_3<K>::type MeshPolyhedron;
typedef MeshPolyhedron::HalfedgeDS      MeshHalfedgeDS;

typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, CGAL::Sequential_tag>::type Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr, Mesh_domain::Corner_index, Mesh_domain::Curve_index> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
// Remesh ---------------

//Skeletonization ------------
typedef Kernel::Point_3                                       Point;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron> Skeletonization;
typedef Skeletonization::Skeleton                             Skeleton;
typedef Skeleton::vertex_descriptor                           Skeleton_vertex;
typedef Skeleton::edge_descriptor                             Skeleton_edge;


typedef boost::property_map<SM, CGAL::edge_is_feature_t>::type EIFMap;


typedef boost::graph_traits<SM>::halfedge_descriptor        halfedge_descriptor;
typedef boost::graph_traits<SM>::edge_descriptor            edge_descriptor;


namespace PMP = CGAL::Polygon_mesh_processing;

//Skeletonization ------------

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

// A modifier creating a triangle with the incremental builder.
template<class HDS>
class polyhedron_builder : public CGAL::Modifier_base<HDS> {
public:
    std::vector<double>& coords;
    std::vector<int>& tris;
    polyhedron_builder(std::vector<double>& _coords, std::vector<int>& _tris) : coords(_coords), tris(_tris) {}
    void operator()(HDS& hds) {
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;

        CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
        B.begin_surface(coords.size() / 3, tris.size() / 3);

        for (int i = 0; i < (int)coords.size(); i += 3) {
            B.add_vertex(Point(coords[i + 0], coords[i + 1], coords[i + 2]));
        }

        for (int i = 0; i < (int)tris.size(); i += 3) {
            B.begin_facet();
            B.add_vertex_to_facet(tris[i + 0]);
            B.add_vertex_to_facet(tris[i + 1]);
            B.add_vertex_to_facet(tris[i + 2]);
            B.end_facet();
        }
        B.end_surface();
    }
};

MeshPolyhedron CreateCGALPolyhedron(std::vector<double> coords, std::vector<int>    tris) {
    Polyhedron P;
    polyhedron_builder<MeshHalfedgeDS> builder(coords, tris);
    MeshPolyhedron mP;

    

    mP.delegate(builder);
    //P.delegate(builder);

    //Polyhedron result;

    //CGAL::copy_face_graph(segment_mesh, result);

    //std::ostringstream oss;
    //oss << "CreateCGALPolyhedron.off";
    //std::ofstream os(oss.str().data());
    //os << mP;

    return mP;
}

struct halfedge2edge
{
    halfedge2edge(const SM& m, std::vector<edge_descriptor>& edges)
        : m_mesh(m), m_edges(edges)
    {}
    void operator()(const halfedge_descriptor& h) const
    {
        m_edges.push_back(edge(h, m_mesh));
    }
    const SM& m_mesh;
    std::vector<edge_descriptor>& m_edges;
};

Polyhedron IsotropicRemesh(SM polyhedron, double target_edge_length) {
    try {
        cout << " ISOTROPIC REMESH START" << endl;
        //std::vector<edge_descriptor> border;
        //PMP::border_halfedges(faces(polyhedron), polyhedron, boost::make_function_output_iterator(halfedge2edge(polyhedron, border)));
        //PMP::split_long_edges(border, target_edge_length, polyhedron);

        unsigned int nb_iter = 3;
        //PMP::isotropic_remeshing(
        //    faces(polyhedron),
        //    target_edge_length,
        //    polyhedron,
        //    PMP::parameters::number_of_iterations(nb_iter)
        //    .protect_constraints(true));

        EIFMap eif = get(CGAL::edge_is_feature, polyhedron);
        PMP::detect_sharp_edges(polyhedron, 60, eif);

        PMP::isotropic_remeshing(faces(polyhedron), target_edge_length, polyhedron, CGAL::parameters::edge_is_constrained_map(eif));

        Polyhedron poly;
        CGAL::copy_face_graph(polyhedron, poly);
        cout << " ISOTROPIC REMESH END" << endl;
        return poly;
    }
    catch (const std::exception& e) {

        PMP::isotropic_remeshing(
        faces(polyhedron),
        target_edge_length,
        polyhedron,
        PMP::parameters::number_of_iterations(3)
        .protect_constraints(true));

        Polyhedron poly;
        CGAL::copy_face_graph(polyhedron, poly);
        cout << " ISOTROPIC REMESH END" << endl;
        return poly;
    }
}

Polyhedron RemeshPolyhedron(MeshPolyhedron meshPolyhedron, double edgeSize = 0.025, double facetAngle = 25,
    double facetSize = 0.05, double facetDistance = 0.005, double cellRadiusEdgeRatio = 3, double cellSize = 0.05) {
    try
    {
        //CGAL::Timer t;
        //t.start();
        // Create domain
        Mesh_domain domain(meshPolyhedron);

        // Get sharp features
        domain.detect_features();
       /* double facetSize = edgeSize * 2;
        double facetDistance = edgeSize / 5;*/
        // Mesh criteria
        //Mesh_criteria criteria(edge_size = edgeSize,
        //    facet_angle = 25, facet_size = facetSize, facet_distance = facetDistance,
        //    cell_radius_edge_ratio = 3, cell_size = facetSize);

        Mesh_criteria criteria(edge_size = edgeSize,
            facet_angle = facetAngle, facet_size = facetSize, facet_distance = facetDistance,
            cell_radius_edge_ratio = cellRadiusEdgeRatio, cell_size = cellSize);

        //Mesh_criteria criteria(
        //facet_angle = 30, facet_size = 6, facet_distance = 4,
        //cell_radius_edge_ratio = 3, cell_size = 8);

        //Mesh_criteria criteria(edge_size = 0.5,
        //    facet_angle = 25, facet_size = 0.55, facet_distance = 0.055,
        //    cell_radius_edge_ratio = 3, cell_size = 0.05);

        // Mesh generation
        C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
            no_perturb(), no_exude());
        Polyhedron outputPoly;
        CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, outputPoly);
        //std::cerr << t.time() << " sec. Remeshing ";

        //std::ofstream off_file("new_result_sample.off");
        //c3t3.output_boundary_to_off(off_file);
        //cout << "Remesh done" << endl;
        return outputPoly;
    }
    catch (const std::exception& e) {
        std::cout << e.what();
        Polyhedron outputPoly;
        return outputPoly;
    }
}

std::vector<SM> Segmentation(MeshPolyhedron poly) {
        cout << "SEGMENTATION" << endl;
        SM mesh;
        CGAL::copy_face_graph(poly, mesh);
        typedef SM::Property_map<face_descriptor, double> Facet_double_map;
        Facet_double_map sdf_property_map;
        
        sdf_property_map = mesh.add_property_map<face_descriptor, double>("f:sdf").first;

        CGAL::sdf_values(mesh, sdf_property_map);

        // create a property-map for segment-ids
        typedef SM::Property_map<face_descriptor, std::size_t> Facet_int_map;
        Facet_int_map segment_property_map = mesh.add_property_map<face_descriptor, std::size_t>("f:sid").first;;

        // segment the mesh using default parameters for number of levels, and smoothing lambda
        // Any other scalar values can be used instead of using SDF values computed using the CGAL function
        std::size_t number_of_segments = CGAL::segmentation_from_sdf_values(mesh, sdf_property_map, segment_property_map);

        typedef CGAL::Face_filtered_graph<SM> Filtered_graph;
        Filtered_graph segment_mesh(mesh);
        std::vector<SM> polyhedronResult;
        for (std::size_t id = 0; id < number_of_segments; ++id)
        {
            segment_mesh.set_selected_faces(id, segment_property_map);
            SM result;

            CGAL::copy_face_graph(segment_mesh, result);
            polyhedronResult.push_back(result);
        }
        cout << polyhedronResult.size() << endl;
        if (polyhedronResult.size() == 0) {
            polyhedronResult.push_back(mesh);
        }
        cout << "SEGMENTATION DONE"  << endl;
        return polyhedronResult;
}


std::vector<std::vector<Point>> SkelAllPts;
std::vector<Point> SkelSegment;
struct Display_polylines {
    const Skeleton& skeleton;
    Display_polylines(const Skeleton& skeleton)
        : skeleton(skeleton) {}
    void start_new_polyline() {}
    void add_node(Skeleton_vertex v) {
        SkelSegment.push_back(skeleton[v].point);
    }
    void end_polyline()
    {
        SkelAllPts.push_back(SkelSegment);
        SkelSegment.clear();
    }
};

JSONModels::Polylines SkeletonizePolyhedron(Polyhedron tmesh) {
    try {
        JSONModels::Polylines result;
        if (!CGAL::is_triangle_mesh(tmesh))
        {
            std::cout << "Input geometry is not triangulated." << std::endl;
            if (CGAL::Polygon_mesh_processing::triangulate_faces(tmesh))
                std::cout << "...triangulated the input...\n";
            else {
                std::cout << "...could not triangulate...\n";
                return result;
            }
        }

        Skeleton skeleton;
        CGAL::extract_mean_curvature_flow_skeleton(tmesh, skeleton);
        Display_polylines display(skeleton);
        CGAL::split_graph_into_polylines(skeleton, display);
        for (std::vector<Point> seg : SkelAllPts) {
            JSONModels::Polyline pLine;
            for (int i = 0; i < seg.size(); i++) {
                JSONModels::Point p;
                p.X(seg[i].x());
                p.Y(seg[i].y());
                p.Z(seg[i].z());
                pLine.PointList().push_back(p);
            }
            result.PolylineList().push_back(pLine);
        }
        return result;
    }
    catch (const std::exception& e) {
        std::cout << e.what();
        JSONModels::Polylines result;
        return result;
    }
}

void GeneratePolyhedron(tuple <vector<int>, vector<double>> verticesAndIndices, std::string outputName = "output.json") {
    vector<int> indices = get<0>(verticesAndIndices);
    vector<double> vertices = get<1>(verticesAndIndices);
    vector<JSONModels::Point> points;
    cout << outputName << " OUTPUT NAME" << endl;

    for (size_t i = 0; i < vertices.size(); i+=3)
    {
        JSONModels::Point p;
        p.X(vertices[i]);
        p.Y(vertices[i + 1]);
        p.Z(vertices[i + 2]);
        points.push_back(p);
    }

    vector<JSONModels::Face> faces;
    JSONModels::Polyhedron polyhedron;
    for (size_t i = 0; i < indices.size(); i+=3)
    {
        vector<JSONModels::Polygon> polygons;
        JSONModels::Face face;
            JSONModels::Polygon polygon;
            vector<JSONModels::Point> outputPoints;
            outputPoints.push_back(points[indices[i]]);
            outputPoints.push_back(points[indices[i + 1]]);
            outputPoints.push_back(points[indices[i + 2]]);
            polygon.PointsList(outputPoints);
            polygons.push_back(polygon);
        face.PolygonList(polygons);
        faces.push_back(face);
    }
    polyhedron.FaceList(faces);
    polyhedron.SerializeToFile(outputName);
}

tuple<vector<int>, vector<double>> GetIndicesAndVerticesFromCGALPoly(Polyhedron P) {

    vector<int> indices;
    std::vector<double> overallVector;
    for (Vertex_iterator v = P.vertices_begin(); v != P.vertices_end(); ++v) {
        overallVector.push_back(v->point().x());
        overallVector.push_back(v->point().y());
        overallVector.push_back(v->point().z());
    }

    for (Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
        Halfedge_facet_circulator j = i->facet_begin();
        // Facets in polyhedral surfaces are at least triangles.
        CGAL_assertion(CGAL::circulator_size(j) >= 3);
        //std::cout << CGAL::circulator_size(j) << ' ';
        do {
            indices.push_back(std::distance(P.vertices_begin(), j->vertex()));
        } while (++j != i->facet_begin());
    }
    tuple<vector<int>, vector<double>> verticesIndicies(indices, overallVector);
    return verticesIndicies;
}

tuple <vector<int>, vector<double>> GetVerticesAndIndices(JSONModels::Polyhedron poly) {
	std::map<vector<double>, int> dictionary;
	int index = 0;
	vector<int> indices;
	std::vector<double> overallVector;

	for (size_t i = 0; i < poly.FaceList().size(); i++)
	{
		vector<JSONModels::Polygon> Face = poly.FaceList()[i].PolygonList();
		for (size_t j = 0; j < Face.size(); j++)
		{
			vector<JSONModels::Point> Points = Face[j].PointsList();
			for (size_t k = 0; k < Points.size(); k++)
			{
				std::vector<double> v;
				v.push_back(Points[k].X());
				v.push_back(Points[k].Y());
				v.push_back(Points[k].Z());
				if (dictionary.find(v) == dictionary.end()) {
					dictionary.insert({ v, index });
					overallVector.push_back(Points[k].X());
					overallVector.push_back(Points[k].Y());
					overallVector.push_back(Points[k].Z());
					indices.push_back(index);
					index++;
				}
				else {
					int val = dictionary.at(v);
					indices.push_back(val);
				}
			}
		}
	}
	tuple<vector<int>, vector<double>> verticesIndicies(indices, overallVector);
	return verticesIndicies;
}

int main(int argc, char* argv[]) {

    if (argc < 5) {
        cout << "ERROR: Usage: -GetCenterline -path [input path] -output [output path] -generateSegments [true|false]" << endl;
        return 0;
    }
    char* geometryFunction = argv[0];
    char* inputPath = NULL;
    char* outputPath = NULL;

    bool isSite = false;
    bool test = false;

    for (size_t i = 1; i < argc; i++)
    {
        cout << argv[i] << endl;
        if (strcmp(argv[i], "-path") == 0)
        {
            inputPath = argv[i + 1];
        }
        else if (strcmp(argv[i], "-output") == 0) {
            outputPath = argv[i + 1];
        }
        else if (strcmp(argv[i], "-site") == 0) {
            std::string action(argv[i + 1]);
            if (action == "True" || action == "true") {
                isSite = true;
            }
            else {
                isSite = false;
            }
        }
        else if (strcmp(argv[i], "-test") == 0) {
            std::string action(argv[i + 1]);
            if (action == "True" || action == "true") {
                test = true;
            }
            else {
                test = false;
            }
        }
    }

    if (strcmp(argv[1], "-GetCenterline") == 0) {
        if (((inputPath == NULL)) || ((outputPath == NULL))) {
            cout << "ERROR: -path \"input json path\" and -output \"output path\" should be present as an argument.";
            return 0;
        }
        JSONModels::Polyhedron polyhedrons;
        polyhedrons.DeserializeFromFile(inputPath);
        CGAL::Timer t;
        t.start();
        cout << "Start getting vertices and indices" << endl;
        tuple<vector<int>, vector<double>> vAndI = GetVerticesAndIndices(polyhedrons);
        cout << "Start constructing cgal polyhedron" << endl;
        MeshPolyhedron p = CreateCGALPolyhedron(get<1>(vAndI), get<0>(vAndI));
        cout << "Start remeshing and skeletonization" << endl;
        if (isSite) {
            cout << "site here" << endl;
            int i = 0;
            JSONModels::Polylines result;
            std::string strName = std::string(outputPath);
            std::string inputPathstr = std::string(inputPath);
            std::string s = std::to_string(i);
            strName += "\\segment\\" + s + ".json";
            SM meshSM;
            CGAL::copy_face_graph(p, meshSM);
            Polyhedron seg;
            CGAL::copy_face_graph(meshSM, seg);
            GeneratePolyhedron(GetIndicesAndVerticesFromCGALPoly(seg), strName);

            double volume = CGAL::Polygon_mesh_processing::volume(meshSM);
            double sizeFace = 0;
            sizeFace = volume / 100;
            cout << volume << " Volume" << endl;
            cout << sizeFace << " SIZE OF FACE" << endl;
            if (sizeFace < 0.001) {
                sizeFace = .009;
            }
            Polyhedron remeshed = IsotropicRemesh(meshSM, sizeFace);


            JSONModels::Polylines centerLine = SkeletonizePolyhedron(remeshed);
            result.PolylineList().insert(result.PolylineList().end(), centerLine.PolylineList().begin(), centerLine.PolylineList().end());


            std::string lineStrName = std::string(outputPath);
            lineStrName += "\\centerline\\" + s + ".json";

            if (test) {
                std::string off = std::string(outputPath);
                off += "\\segment\\remeshed_" + s + ".off";
                std::cout << off << endl;
                std::ostringstream oss;
                oss << off;
                std::ofstream os(oss.str().data());
                os << remeshed;

                std::string off1 = std::string(outputPath);
                off1 += "\\segment\\original_" + s + ".off";
                std::cout << off1 << endl;
                std::ostringstream oss1;
                oss1 << off1;
                std::ofstream os1(oss1.str().data());
                os1 << meshSM;
            }

            result.SerializeToFile(lineStrName);
            result.PolylineList().clear();
        }
        else {
            cout << "not site" << endl;
            vector<SM> segmentResult = Segmentation(p);
            for (size_t i = 0; i < segmentResult.size(); i++)
            {
                cout << "NASA LOOP" << endl;
                JSONModels::Polylines result;
                std::string strName = std::string(outputPath);
                std::string inputPathstr = std::string(inputPath);
                std::string s = std::to_string(i);
                strName += "\\segment\\" + s + ".json";
              
                SM meshSM;
                CGAL::copy_face_graph(segmentResult[i], meshSM);
                Polyhedron seg;
                CGAL::copy_face_graph(segmentResult[i], seg);
                GeneratePolyhedron(GetIndicesAndVerticesFromCGALPoly(seg), strName);

                double volume = CGAL::Polygon_mesh_processing::volume(meshSM);
                double sizeFace = 0;
                if (volume < .01) {
                    sizeFace = .00001;
                }
                else if (volume > .01 && volume < .1) {
                    sizeFace = .01;
                }
                else if (volume >= 1 && volume <= 5) {
                    sizeFace = .05;
                }
                else if (volume >= 5 && volume <= 10) {
                    sizeFace = .025;
                }
                else if (volume >= 10 && volume <= 100) {
                    sizeFace = .5;
                }
                else if (volume >= 100) {
                    sizeFace = 1;
                }
                sizeFace = volume / 2;
                if (sizeFace < 0.001) {
                    sizeFace = .009;
                }
                cout << volume << " Volume" << endl;
                cout << sizeFace << " SIZE OF FACE" << endl;
                if (test) {

                    std::string off1 = std::string(outputPath);
                    off1 += "\\segment\\original_" + s + ".off";
                    std::cout << off1 << endl;
                    std::ostringstream oss1;
                    oss1 << off1;
                    std::ofstream os1(oss1.str().data());
                    os1 << meshSM;
                }


                Polyhedron remeshed = IsotropicRemesh(meshSM, sizeFace);

                if (test) {
                    std::string off = std::string(outputPath);
                    off += "\\segment\\remeshed_" + s + ".off";
                    std::cout << off << endl;
                    std::ostringstream oss;
                    oss << off;
                    std::ofstream os(oss.str().data());
                    os << remeshed;
                }

                JSONModels::Polylines centerLine = SkeletonizePolyhedron(remeshed);
                result.PolylineList().insert(result.PolylineList().end(), centerLine.PolylineList().begin(), centerLine.PolylineList().end());

                std::string lineStrName = std::string(outputPath);
                lineStrName += "\\centerline\\" + s + ".json";

               

                result.SerializeToFile(lineStrName);
                result.PolylineList().clear();

            }
        }
        cout << "Finished remeshing and skeletonization" << endl;

        cout << t.time() << " sec. All process done" << std::endl;
        return 1;
    }


    return 0;

}