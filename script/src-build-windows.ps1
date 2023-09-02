$SRC_DIR = "../src"
$BUILD_DIR = "../build"
$CGAL_CC_PATH = "../dep/vcpkg/scripts/buildsystems/vcpkg.cmake"
$VCPP_VERSION = '"Visual Studio 17 2022"' # change if needed

Invoke-Expression ("cmake -S " + $SRC_DIR + " -B " + $BUILD_DIR + " -G " + $VCPP_VERSION + " --toolchain " + $CGAL_CC_PATH)