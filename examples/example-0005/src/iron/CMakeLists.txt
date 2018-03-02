add_executable(${ExampleName} ${ExampleName}.F90 ${IRON_MESH_READER_ROOT}/src/meshReader.F90)

target_link_libraries(${ExampleName} PRIVATE opencmisslibs)

install(TARGETS ${ExampleName} DESTINATION ${PROJECT_SOURCE_DIR}/bin)

