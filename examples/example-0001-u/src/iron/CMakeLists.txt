add_executable(${ExampleName} ${ExampleName}.F90 ~/software/opencmiss/dependencies/meshReader-master/src/meshReader.F90)

target_link_libraries(${ExampleName} opencmisslibs)
install(TARGETS ${ExampleName} DESTINATION ${PROJECT_SOURCE_DIR}/bin)
configure_file(visualise.cmgui visualise.cmgui)


