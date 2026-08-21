# Gather all source and header files in this directory
file(GLOB LINALG_SRC CONFIGURE_DEPENDS src/*.cpp)
file(GLOB LINALG_HEADERS CONFIGURE_DEPENDS include/*.hpp)

# Create a static library
add_library(lin-alg STATIC ${LINALG_SRC})

# Expose the include directory so other projects can find your headers
target_include_directories(lin-alg PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
