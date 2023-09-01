include_directories(BEFORE ../../include/uBidomain)

# Utilities
add_library(utils INTERFACE)
target_link_libraries(utils INTERFACE Eigen3::Eigen)

add_library(fft fft.cpp FSeries.cpp getV1.cpp)
target_include_directories(fft PUBLIC ${FFTW3_INCLUDE_DIRS})
target_link_directories(fft PUBLIC ${FFTW3_LIBRARY_DIRS})
target_link_libraries(fft PUBLIC Eigen3::Eigen ${FFTW3_LIBRARIES})

add_library(exact_sol exact_sol.cpp)
target_link_libraries(exact_sol PUBLIC Eigen3::Eigen)


# MultiDomains
add_library(side Side.cpp)
target_link_libraries(side PUBLIC Eigen3::Eigen)

add_library(boundary Boundary.cpp)
target_link_libraries(boundary PUBLIC side fft)

add_library(domain Domain.cpp)
target_link_libraries(domain PUBLIC boundary ${GMSH_LIB})
target_include_directories(domain PRIVATE ${GMSH_INC} ../../include/lean_vtk)

add_library(multi_domain MultiDomain.cpp)
target_link_libraries(multi_domain PUBLIC domain)

add_library(cellbycell_model CellByCellModel.cpp)
target_link_libraries(cellbycell_model PUBLIC multi_domain)


# Ionic models
add_library(ionic_model IonicModel.cpp)
target_link_libraries(ionic_model PUBLIC Eigen3::Eigen)

add_library(myokit_ionic_model MyoKitIonicModel.cpp)
target_link_libraries(myokit_ionic_model PUBLIC ionic_model)
