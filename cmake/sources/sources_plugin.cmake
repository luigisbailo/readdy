SET(SOURCES_DIR "${READDY_GLOBAL_DIR}/readdy2/main/plugin")

# includes
SET(PLUGIN_INCLUDE_DIRS "${MODEL_INCLUDE_DIRS};${Boost.DLL_INCLUDE_DIRS}" CACHE INTERNAL "Plugin include dirs" FORCE)

# libraries
SET(READDY_PLUGIN_LIBRARIES "${READDY_MODEL_LIBRARIES};${CMAKE_DL_LIBS}")

# sources
LIST(APPEND READDY_PLUGIN_SOURCES "${SOURCES_DIR}/KernelProvider.cpp")
LIST(APPEND READDY_PLUGIN_SOURCES "${SOURCES_DIR}/KernelPluginDecorator.cpp")

# all sources
LIST(APPEND READDY_ALL_SOURCES ${READDY_PLUGIN_SOURCES})