file(GLOB ODE_SRC CONFIGURE_DEPENDS src/*.cpp)
file(GLOB ODE_HEADERS CONFIGURE_DEPENDS include/*.hpp)

# Create ODE solver static library
add_library(ode-solvers STATIC ${ODE_SRC})

# Link the custom linear algebra library
target_link_libraries(ode-solvers PUBLIC linalg)

# Expose ODE headers
target_include_directories(ode-solvers PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
