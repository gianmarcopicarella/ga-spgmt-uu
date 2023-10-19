#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Point_2.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Line_2.h>
#include <CGAL/Direction_3.h>
#include <CGAL/Direction_2.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Sphere_3.h>
//#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Iso_cuboid_3.h>

// Only for debug
//#include <CGAL/Simple_cartesian.h>
//#include <CGAL/Exact_rational.h>

namespace SPGMT
{
    /*typedef double ExactReal;
    typedef CGAL::Simple_cartesian<ExactReal> Kernel;
    typedef Kernel::FT 	FT;
    typedef Kernel::Point_3 Point3;
    typedef CGAL::Point_set_3<Point3> PointSet3;
    typedef Kernel::Plane_3 Plane;
    typedef Kernel::Line_3 Line3;
    typedef Kernel::Direction_3 Dir3;
    typedef Kernel::Vector_3 Vec3;
    typedef Kernel::Line_2 Line2;
    typedef Kernel::Point_2 Point2;
    typedef Kernel::Vector_2 Vec2;*/

    //typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt Kernel;
    typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
    typedef Kernel::FT 	FT;
    typedef Kernel::Point_3 Point3;
    typedef CGAL::Point_set_3<Point3> PointSet3;
    typedef Kernel::Plane_3 Plane;
    typedef Kernel::Line_3 Line3;
    typedef Kernel::Direction_3 Dir3;
    typedef Kernel::Vector_3 Vec3;
    typedef Kernel::Line_2 Line2;
    typedef Kernel::Point_2 Point2;
    typedef Kernel::Vector_2 Vec2;
    typedef Kernel::Direction_2 Dir2;
    typedef CGAL::Segment_3<Kernel> Segment3;
    typedef Kernel::Ray_2 Ray2;
    typedef Kernel::Sphere_3 Sphere3;
    //typedef Kernel::Iso_rectangle_2 Rec2;
    typedef Kernel::Iso_cuboid_3 Cube;

    enum class EdgeType
    {
        LINE,
        HALF_EDGE_EF,
        HALF_EDGE_SF,
        SEGMENT,
        SEGMENT_TRIANGLE
    };

    template<typename T>
    struct Edge
    {
        T myStart, myEnd;
        EdgeType myType{ EdgeType::SEGMENT };
        int myLowestLeftPlane{ -1 };
    };

    using LowerEnvelope3d = std::variant<std::monostate, size_t, std::vector<Edge<Point3>>>;

    enum class ExecutionPolicy
    {
        SEQUENTIAL,
        PARALLEL
    };
}