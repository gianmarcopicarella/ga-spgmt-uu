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
#include "../common/DebugUtils.h"
//
//#include <CGAL/boost/graph/Dual.h>
//#include <CGAL/boost/graph/helpers.h>

//#include <hpx/init.hpp>
//#include <hpx/hpx.hpp>

//#include <iostream>
//#include <chrono>

using namespace SPGMT;

//--hpx:threads = N

//void testSort(std::vector<Point2>& someOutItems)
//{
//	const auto threadsCount = hpx::resource::get_num_threads();
//	const auto splitCount = threadsCount % 2 == 0 ? threadsCount : (threadsCount + 1);
//	const auto blockSize = someOutItems.size() / splitCount;
//	const auto parLimit = blockSize * splitCount;
//
//	hpx::for_loop_strided(hpx::execution::par, 0, parLimit, blockSize, [&](size_t i) {
//		std::sort(someOutItems.begin() + i, someOutItems.begin() + i + blockSize, CGAL::Less<Point2, Point2>());
//		});
//	std::sort(someOutItems.begin() + parLimit, someOutItems.end(), CGAL::Less<Point2, Point2>());
//
//	auto k = blockSize;
//	while (2 * k <= parLimit)
//	{
//		hpx::for_loop_strided(hpx::execution::par, 0, parLimit, k*2, [&](size_t i) {
//			std::inplace_merge(someOutItems.begin() + i, someOutItems.begin() + i + k, someOutItems.begin() + i + 2 * k, CGAL::Less<Point2, Point2>());
//			});
//		k *= 2;
//	}
//
//	std::inplace_merge(someOutItems.begin(), someOutItems.begin() + parLimit, someOutItems.end(), CGAL::Less<Point2, Point2>());
//}
//
//
//int hpx_func()
//{
//	//using std::chrono::high_resolution_clock;
//	//using std::chrono::duration_cast;
//	//using std::chrono::duration;
//	//using std::chrono::milliseconds;
//	//
//	//std::cout << hpx::resource::get_num_threads() << std::endl;
//
//	//std::vector<Point3> v = SPGMT::Debug::Uniform3DCubeSampling(300.f, 10000);
//	//
//	//std::atomic<size_t> t{ 0 };
//
//	//hpx::for_loop_strided(hpx::execution::seq, 0, v.size(), 100, [](size_t val) {  hpx::util::format_to(std::cout, "{1}\n", val);  });
//
//	///*hpx::execution::experimental::static_chunk_size scs;
//	//hpx::for_loop(
//	//	hpx::execution::par.with(scs),
//	//	0, v.size(),
//	//	[](size_t val) { std::cout << val << std::endl;  });*/
//
//	//auto t1 = high_resolution_clock::now();
//
//	//hpx::for_loop(
//	//	hpx::execution::par,
//	//	0, v.size(), [&](size_t idx) { v[idx] = Point3{ CGAL::Random().get_double(),CGAL::Random().get_double() ,CGAL::Random().get_double() }; }
//	//);
//
//	//auto t2 = high_resolution_clock::now();
//
//	//auto ms_int = duration_cast<milliseconds>(t2 - t1);
//
//	//std::cout << "Elapsed time: " << ms_int.count() << std::endl;
//
//	std::vector<int> test(1034430);
//	CGAL::Random r;
//	for (int i = 0; i < 100; ++i)
//	{
//		test[i] = r.get_int(-200, 200);
//	}
//
//
//	
//
//
//
//
//	return hpx::local::finalize();
//}

int main(int argc, char* argv[]) {

	/*hpx::local::init(hpx_func, argc, argv);*/



	//std::cout << "abc" << std::endl;

	/*typedef CGAL::Dual<Plane> DualPlane;

	DualPlane dp{ Plane{10.f, 1.f, 2.03f, -4.f} };

	std::cout << CGAL::num_vertices(dp) << std::endl;*/

	return 0;

	/*unsigned int n = std::thread::hardware_concurrency();
	//std::cout << n << " concurrent threads are supported.\n";*/
}

