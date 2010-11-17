PROJECT(vpano)
#SET(CMAKE_BUILD_TYPE debug)
ADD_EXECUTABLE(vpano v3dpano.cpp main.cpp geopano.cpp bestfit.cpp)

include_directories(/usr/local/include/jplv)
link_libraries(/usr/local/lib/ljplv.a)

include_directories(/usr/local/include/opencv)
#link_libraries(/home/jaechoon/Desktop/OpenCV-2.1.0/release/lib/libcv.so)
#link_libraries(/home/jaechoon/Desktop/OpenCV-2.1.0/release/lib/libcvaux.so)
FIND_PACKAGE( OpenCV REQUIRED )

find_package(GLUT)
find_package(OpenGL)
include_directories(/usr/include/GL)

set(CORELIBS ${GLUT_LIBRARY} ${OPENGL_LIBRARY} ${OpenCV_LIBS} )
target_link_libraries(vpano ${CORELIBS})
