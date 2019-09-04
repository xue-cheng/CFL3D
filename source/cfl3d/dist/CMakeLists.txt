
set(DIST_SRC
	bc_blkint.F     findmin_new.F   plot3d.F        rrest.F
	bc_patch.F      forceout.F      plot3t.F        rrestg.F
	calyplus.F      pointers.F      setup.F         writ_buf.F
	mgblk.F         qinter.F        prntcp.F        newalpha.F
	cputim.F        patcher.F       qout.F          termn8.F
	dynptch.F       plot3c.F        resp.F          usrint.F
	wrest.F         wrestg.F        pre_bc.F        bc_embed.F
	updateg.F       compg2n.F       resetg.F        bc_period.F
	yplusout.F      sizer.F         cfl3d.F         trnsfr_vals.F
	updatedg.F      ae_corr.F       mgbl.F          setslave.F
	umalloc.F       reass.F         qoutavg.F       plot3davg.F
  qout_coarse.F   qout_2d.F       plot3d_2d.F     plot3d_coarse.F
)

set(CFL3D_SRC ${CFLSRC_S}/ccomplex.F development.F main.F)

# REAL version of CFL3D
if(BUILD_REAL)

  if(BUILD_SEQ)
    add_library(dist ${DIST_SRC})
    target_link_libraries(dist PUBLIC common)
    if(WITH_CGNS)
      target_link_libraries(dist PUBLIC cgns)
    endif(WITH_CGNS)
    add_executable(cfl3d_seq ${CFL3D_SRC})
    target_link_libraries(cfl3d_seq PUBLIC dist)
  endif(BUILD_SEQ)

  if(BUILD_MPI)
    add_library(dist_mpi ${DIST_SRC})
    target_link_libraries(dist_mpi PUBLIC common)
    if(WITH_CGNS)
      target_link_libraries(dist_mpi PUBLIC cgns)
    endif(WITH_CGNS)
    target_link_libraries(dist_mpi PUBLIC MPI::MPI_Fortran)
    add_executable(cfl3d_mpi ${CFL3D_SRC})
    target_link_libraries(cfl3d_mpi PUBLIC dist_mpi)
  endif(BUILD_MPI)
endif(BUILD_REAL)

if(BUILD_CMPLX)
  if(BUILD_SEQ)
    add_library(dist_cmplx ${DIST_SRC})
    target_link_libraries(dist_cmplx PUBLIC common_cmplx)
    if(WITH_CGNS)
      target_link_libraries(dist PUBLIC cgns)
    endif(WITH_CGNS)
    add_executable(cfl3d_seq ${CFL3D_SRC})
    target_link_libraries(cfl3d_seq PUBLIC dist)
  endif(BUILD_SEQ)

  if(BUILD_MPI)
    add_library(dist_cmplx_mpi ${DIST_SRC})
    target_link_libraries(dist_cmplx_mpi PUBLIC common_cmplx)
    if(WITH_CGNS)
      target_link_libraries(dist_cmplx_mpi PUBLIC cgns)
    endif(WITH_CGNS)
    target_link_libraries(dist_cmplx_mpi PUBLIC MPI::MPI_Fortran)
    add_executable(cfl3dcmplx_mpi ${CFL3D_SRC})
    target_link_libraries(cfl3dcmplx_mpi PUBLIC dist_cmplx_mpi)
  endif(BUILD_MPI)

endif(BUILD_CMPLX)