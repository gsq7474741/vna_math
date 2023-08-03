cmake_minimum_required(VERSION 3.25)
project(vna_math)

set(CMAKE_CXX_STANDARD 17)



#add_subdirectory(external/googletest)
#target_include_directories(Test PUBLIC
#        ${CMAKE_SOURCE_DIR}/external/googletest/googletest/include
#        )
#target_link_libraries(Test gtest) # 不想自己写main函数的话可以连接gtest_main
## 这里最重要
#if(WIN32)
#    # see https://github.com/google/googletest/issues/4067
#    if (CMAKE_BUILD_TYPE STREQUAL Debug)
#        set_property(TARGET gtest PROPERTY MSVC_RUNTIME_LIBRARY "/MD")
#        set_property(TARGET gtest_main PROPERTY MSVC_RUNTIME_LIBRARY "/MD")
#    else()
#        set_property(TARGET gtest PROPERTY MSVC_RUNTIME_LIBRARY "/MT")
#        set_property(TARGET gtest_main PROPERTY MSVC_RUNTIME_LIBRARY "/MT")
#    endif()
#    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
#endif()
#
#add_test(NAME TEST COMMAND Test)
#
#add_executable(vna_math_test test/test.cpp)
##find_package(GTest REQUIRED)



# enable c++ 11
#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED True)

# take respectives flags for debug & release process
set(CMAKE_CXX_FLAGS_DEBUG_INIT "-std=c++17 -g -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS_RELEASE_INIT "-std=c++17 -g -O2")

# default build type : Debug
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

add_subdirectory(src)

# CTest related
enable_testing()
add_subdirectory(tests)

#export https_proxy=http://172.21.80.1:4780
#export http_proxy=http://172.21.80.1:4780
#export all_proxy=socks5://172.21.80.1:4780