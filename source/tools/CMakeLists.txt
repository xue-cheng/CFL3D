if(BUILD_TOOLS)
  set(TOOLS_SRC
    grid_perturb_cmplx.F 
    grid_perturb.F
    Get_FD.F
    moovmaker.F
    plot3dg_to_cgns.F
    XINTOUT_to_ovrlp.F
    p3d_to_INGRID.F
    INGRID_to_p3d.F
    cfl3d_to_pegbc.F
    everyother_xyz.F
    initialize_field.F
    p3d_to_cfl3drst.F
    v6_restart_mod.F
    v6inpdoubhalf.F
    cgns_to_cfl3dinput.F
    cgns_readhist.F
    v6_ronnie_mod.F90
    cfl3d_to_nmf.F
    cfl3dinp_to_FVBND.F
    nmf_to_cfl3dinput.F
    gridswitchijk.F
    v6inpswitchijk.F
  )
  set(SHARED_SRC
    ${CFLSRC_D}/umalloc_r.F
    ${CFLSRC_D}/umalloc_c.F
    ${CFLSRC_S}/readkey.F
    ${CFLSRC_S}/global0.F
    ${CFLSRC_S}/parser.F
    ${CFLSRC_S}/outbuf.F
    ${CFLSRC_S}/cgnstools.F 
    ${CFLSRC_S}/my_flush.F
    ${CFLSRC_D}/termn8.F
  )
  add_library(TOOLS_SHARED STATIC ${SHARED_SRC})
  if (WITH_CGNS)
    target_link_libraries(TOOLS_SHARED PUBLIC cgns)
  endif (WITH_CGNS)
  # add dummy target to group the tools
  add_custom_target(tools ALL)
  foreach (tool_src ${TOOLS_SRC})
    get_filename_component(tool ${tool_src} NAME_WE)
    add_executable(${tool} ${tool_src})
    target_link_libraries(${tool} PUBLIC TOOLS_SHARED)
    add_dependencies(tools ${tool})
  endforeach (tool_src)
endif(BUILD_TOOLS)