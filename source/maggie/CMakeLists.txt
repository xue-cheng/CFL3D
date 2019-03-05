
if(BUILD_MAGGIE)
  set(MAG_SRCS cputim.F maggie.F wkstn.F)
  add_executable(maggie_seq ${MAG_SRCS})
  target_include_directories(maggie PUBLIC "${CFL3D_SOURCE_DIR}/header")
endif(BUILD_MAGGIE)
