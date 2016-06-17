# Find pthread library
FIND_PACKAGE(Threads REQUIRED)

LINK_LIBRARIES(gtest gtest_main ${CMAKE_THREAD_LIBS_INIT} ReaDDy)

LIST(APPEND READDY_TEST_SOURCES TestMain.cpp)
LIST(APPEND READDY_TEST_SOURCES TestPlugins.cpp)
LIST(APPEND READDY_TEST_SOURCES TestSimulation.cpp)
LIST(APPEND READDY_TEST_SOURCES TestObservables.cpp)
LIST(APPEND READDY_TEST_SOURCES TestKernelContext.cpp)
ADD_EXECUTABLE(runUnitTests ${READDY_TEST_SOURCES} ${GOOGLETEST_INCLUDE})

INCLUDE_DIRECTORIES(BEFORE ${ReaDDy_INCLUDE_DIRS})

TARGET_LINK_LIBRARIES(runUnitTests ReaDDy ${gtest})

IF (READDY_INSTALL_UNIT_TEST_EXECUTABLE)
    INSTALL(TARGETS runUnitTests
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib)
ENDIF (READDY_INSTALL_UNIT_TEST_EXECUTABLE)