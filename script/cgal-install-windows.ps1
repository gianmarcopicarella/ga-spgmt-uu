$VCPKG_BASE = "..\dep\vcpkg\"
$VCPKG_EXE = Join-Path $VCPKG_BASE \vcpkg.exe
# or x86-windows, based on the build target
$env:VCPKG_DEFAULT_TRIPLET = "x64-windows"

Invoke-Expression (Join-Path $VCPKG_BASE \bootstrap-vcpkg.bat)
# because of a bug with gmp in vcpkg for windows the x86 version is needed
Invoke-Expression ($VCPKG_EXE + " install yasm-tool:x86-windows") 
Invoke-Expression ($VCPKG_EXE + " install cgal")
# install visualization libraries (Corrade and Magnum)
# Invoke-Expression ($VCPKG_EXE + " install --head corrade magnum")