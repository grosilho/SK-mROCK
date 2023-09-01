
include_directories(BEFORE ../include)

add_library(MainHeader INTERFACE)
target_link_libraries(MainHeader INTERFACE Eigen3::Eigen)
#target_include_directories(MainHeader INTERFACE)

add_library(AutodiffLib INTERFACE)
target_link_libraries(AutodiffLib INTERFACE Eigen3::Eigen autodiff::autodiff)
#target_include_directories(AutodiffLib INTERFACE autodiff::autodiff)

add_subdirectory(lean_vtk)
add_subdirectory(uBidomain)

add_library(chebyshev_methods ChebyshevMethods.cpp)
target_include_directories(chebyshev_methods PUBLIC MainHeader)
target_link_libraries(chebyshev_methods PUBLIC MainHeader)

add_library(legendre_methods LegendreMethods.cpp)
target_include_directories(legendre_methods PUBLIC MainHeader)
target_link_libraries(legendre_methods PUBLIC MainHeader)

# ODEs
add_library(ode Ode.cpp)
target_link_libraries(ode PUBLIC MainHeader)
target_include_directories(ode PUBLIC MainHeader)

add_library(ode_problems OdeProblems.cpp OdeProblemsAutodiffRhs.cpp)
target_include_directories(ode_problems PUBLIC ../include/lean_vtk ../include/uBidomain AutodiffLib)
target_link_libraries(ode_problems PUBLIC ode lean_vtk myokit_ionic_model AutodiffLib)

add_library(ode_cellbycell_model ODECellByCellModel.cpp)
target_include_directories(ode_cellbycell_model PUBLIC ../include/uBidomain)
target_link_libraries(ode_cellbycell_model PUBLIC cellbycell_model ode myokit_ionic_model)

add_library(multirate_ode_problems MultirateOdeProblems.cpp)
target_link_libraries(multirate_ode_problems PUBLIC ode_problems)

# DSDEs
add_library(stochastic_integrals StochasticIntegrals.cpp)
target_include_directories(stochastic_integrals PUBLIC MainHeader)
target_link_libraries(stochastic_integrals PUBLIC MainHeader)

add_library(dsde DSde.cpp)
target_link_libraries(dsde PUBLIC ode stochastic_integrals)

add_library(dsde_problems DSdeProblems.cpp)
target_link_libraries(dsde_problems PUBLIC dsde ode_problems)

# Others
add_library(parameters Parameters.cpp)
target_include_directories(parameters PUBLIC ../include/lean_vtk ../include/uBidomain ../include/getpot)
target_link_libraries(parameters PUBLIC 
                            stabilized_dsde_rungekutta_integrator
                            dsde_problems
                            multirate_ode_rungekutta_integrator
                            multirate_ode_problems
                            ode_cellbycell_model
                            classical_dsde_rungekutta_integrator
                            )


# ODE Integrators
add_library(time_integrator TimeIntegrator.cpp)
target_link_libraries(time_integrator PUBLIC ode parameters)

add_library(error_controller ErrorController.cpp)
target_link_libraries(error_controller PUBLIC parameters)

add_library(ode_rungekutta_integrator OdeRungeKuttaIntegrator.cpp)
target_link_libraries(ode_rungekutta_integrator PUBLIC time_integrator error_controller)

add_library(multirate_ode_rungekutta_integrator MultirateOdeRungeKuttaIntegrator.cpp)
target_link_libraries(multirate_ode_rungekutta_integrator PUBLIC ode_rungekutta_integrator)

add_library(stabilized_ode_rungekutta_integrator StabilizedOdeRungeKuttaIntegrators.cpp)
target_link_libraries(stabilized_ode_rungekutta_integrator PUBLIC ode_rungekutta_integrator multirate_ode_rungekutta_integrator chebyshev_methods legendre_methods)

add_library(classical_ode_rungekutta_integrator ClassicalOdeRungeKuttaIntegrators.cpp)
target_link_libraries(classical_ode_rungekutta_integrator PUBLIC ode_rungekutta_integrator)

# DSDE Integrators
add_library(dsde_rungekutta_integrator DSdeRungeKuttaIntegrator.cpp)
target_link_libraries(dsde_rungekutta_integrator PUBLIC ode_rungekutta_integrator dsde chebyshev_methods)

add_library(classical_dsde_rungekutta_integrator ClassicalDSdeRungeKuttaIntegrators.cpp)
target_link_libraries(classical_dsde_rungekutta_integrator PUBLIC dsde_rungekutta_integrator classical_ode_rungekutta_integrator)

add_library(stabilized_dsde_rungekutta_integrator StabilizedDSdeRungeKuttaIntegrators.cpp)
target_link_libraries(stabilized_dsde_rungekutta_integrator PUBLIC  dsde_rungekutta_integrator stabilized_ode_rungekutta_integrator)


# MonteCarlo
add_library(monte_carlo MonteCarlo.cpp)
target_link_libraries(monte_carlo PUBLIC parameters dsde_rungekutta_integrator)

# Parareal
add_library(parareal Parareal.cpp)
target_link_libraries(parareal PUBLIC parameters ode ode_rungekutta_integrator)
