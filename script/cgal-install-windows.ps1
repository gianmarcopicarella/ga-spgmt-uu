$VCPKG_BASE = "..\dep\vcpkg\"
$VCPKG_EXE = Join-Path $VCPKG_BASE \vcpkg.exe
$CGAL_TRIPLET = "x64-windows" # or x86-windows, based on the build target

Invoke-Expression (Join-Path $VCPKG_BASE \bootstrap-vcpkg.bat)
# because of a bug with gmp in vcpkg for windows the x86 version is needed
Invoke-Expression ($VCPKG_EXE + " install yasm-tool:x86-windows") 
Invoke-Expression ($VCPKG_EXE + " install cgal --triplet " + $CGAL_TRIPLET)