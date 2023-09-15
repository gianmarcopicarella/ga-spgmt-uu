#include <CGAL/CGAL_Ipelet_base.h>
#include "../common/Types.h"

//Function names of the ipelet
const std::string labels[] = {"2D Face Index", "Help" };
//Help message associated to the first function
const std::string hmsg[] = {
  "Find the 2D face index containing a query point"
};
class FaceIndexLookup
  : public CGAL::Ipelet_base<SPGMT::Kernel,2>{
public:
  //declare an ipelet called CGAL Delaunay, with 2 functions (including help message).
  FaceIndexLookup()
    :CGAL::Ipelet_base<SPGMT::Kernel,2>("CGAL Face Index Lookup", labels, hmsg){}
  void protected_run(int);
};
//function called when using the ipelet.
void FaceIndexLookup::protected_run(int fn)
{
  switch (fn){
    case 1:
      show_help(); //print an help message
      return;
    default:
    break;
  };
}

CGAL_IPELET(FaceIndexLookup)