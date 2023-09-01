
include_directories(BEFORE ../include)

add_library(MainHeader INTERFACE)
target_link_libraries(MainHeader INTERFACE Eigen3::Eigen)
#target_include_directories(MainHeader INTERFACE)

add_library(chebyshev_methods ChebyshevMethods.cpp)
target_include_directories(chebyshev_methods PUBLIC MainHeader)
target_link_libraries(chebyshev_methods PUBLIC MainHeader)

# ODEs
add_library(ode Ode.cpp)
target_link_libraries(ode PUBLIC MainHeader)
target_include_directories(ode PUBLIC MainHeader)

add_library(ode_problems OdeProblems.cpp)
target_include_directories(ode_problems PUBLIC MainHeader)
target_link_libraries(ode_problems PUBLIC ode)

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
target_include_directories(parameters PUBLIC ../include/getpot)
target_link_libraries(parameters PUBLIC 
                            stabilized_dsde_rungekutta_integrator
                            dsde_problems
                            multirate_ode_rungekutta_integrator
                            classical_dsde_rungekutta_integrator
                            )


# ODE Integrators
add_library(time_integrator TimeIntegrator.cpp)
target_link_libraries(time_integrator PUBLIC ode parameters)

add_library(ode_rungekutta_integrator OdeRungeKuttaIntegrator.cpp)
target_link_libraries(ode_rungekutta_integrator PUBLIC time_integrator)

add_library(multirate_ode_rungekutta_integrator MultirateOdeRungeKuttaIntegrator.cpp)
target_link_libraries(multirate_ode_rungekutta_integrator PUBLIC ode_rungekutta_integrator)

add_library(stabilized_ode_rungekutta_integrator StabilizedOdeRungeKuttaIntegrators.cpp)
target_link_libraries(stabilized_ode_rungekutta_integrator PUBLIC ode_rungekutta_integrator multirate_ode_rungekutta_integrator chebyshev_methods)

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

