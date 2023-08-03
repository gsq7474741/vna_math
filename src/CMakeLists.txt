#add_library(libexample SHARED "")
add_library(vna_math library.cpp)

#message(${CMAKE_SOURCE_DIR})
#include_directories(${CMAKE_SOURCE_DIR}/external/Fast4ier)

aux_source_directory(./fft FFT_SRCS)
target_sources(
        vna_math
        PRIVATE
        library.cpp
        ${CMAKE_CURRENT_LIST_DIR}/library.h # optional
        ${FFT_SRCS}
)

target_include_directories(
        vna_math
        PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}
        ./fft
)

target_link_libraries(vna_math m)

# respective macro definitions
message("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
if(CMAKE_BUILD_TYPE MATCHES Debug)
    target_compile_definitions(${PROJECT_NAME} PUBLIC _DEBUG)
elseif(CMAKE_BUILD_TYPE MATCHES Release)
else()
endif()