#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_3.h>

namespace SPGMT
{
    typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
    typedef Kernel::FT FT;
    typedef Kernel::Point_3 Point3;
    typedef CGAL::Point_set_3<Point3> PointSet3;
}