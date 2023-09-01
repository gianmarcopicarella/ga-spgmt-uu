$VCPKG_BASE = "..\dep\vcpkg\"
$VCPKG_EXE = Join-Path $VCPKG_BASE \vcpkg.exe
$CGAL_TRIPLET = "x64-windows" # or x86-windows, based on the build target

Invoke-Expression (Join-Path $VCPKG_BASE \bootstrap-vcpkg.bat)
Invoke-Expression ($VCPKG_EXE + " install yasm-tool:" + $CGAL_TRIPLET)
Invoke-Expression ($VCPKG_EXE + " install cgal --triplet " + $CGAL_TRIPLET)