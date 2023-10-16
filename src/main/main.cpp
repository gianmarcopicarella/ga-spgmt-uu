//
//#include "../common/Types.h"
//
//using namespace SPGMT;
//
//int main()
//{
//	/*Plane p{ 10.f, 1.f, 2.03f, -4.f };
//	Plane o{ p.point(), -p.orthogonal_direction() };
//
//	std::cout << p << std::endl;
//	std::cout << o << std::endl;*/
//
//	return 0;
//}

#include "../common/Types.h"

#include <CGAL/boost/graph/Dual.h>
#include <CGAL/boost/graph/helpers.h>

#include <hpx/init.hpp>
#include <hpx/hpx.hpp>

#include <iostream>
#include <chrono>

using namespace SPGMT;

//--hpx:threads = N

int hpx_func()
{
	using std::chrono::high_resolution_clock;
	using std::chrono::duration_cast;
	using std::chrono::duration;
	using std::chrono::milliseconds;

	std::vector<int> values(1000000000);
	std::iota(values.begin(), values.end(), 0);

	auto t1 = high_resolution_clock::now();

	hpx::parallel::sort(
		hpx::execution::par,
		values.begin(), values.end(),
		[](int a, int b)
		{
			return a > b;
		}
	);

	auto t2 = high_resolution_clock::now();

	auto ms_int = duration_cast<milliseconds>(t2 - t1);

	std::cout << "Elapsed time: " << ms_int.count() << std::endl;

	return hpx::local::finalize();
}

int main(int argc, char* argv[]) {

	hpx::local::init(hpx_func, argc, argv);



	//std::cout << "abc" << std::endl;

	/*typedef CGAL::Dual<Plane> DualPlane;

	DualPlane dp{ Plane{10.f, 1.f, 2.03f, -4.f} };

	std::cout << CGAL::num_vertices(dp) << std::endl;*/

	return 0;

	/*unsigned int n = std::thread::hardware_concurrency();
	//std::cout << n << " concurrent threads are supported.\n";*/
}

