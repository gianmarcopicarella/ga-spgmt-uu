#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Point_2.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Line_2.h>
#include <CGAL/Direction_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Vector_2.h>

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
}