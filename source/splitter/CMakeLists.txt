if(BUILD_SPLITTER)
  
  if(BUILD_REAL)
    set(SPLIT_SRC
      cfl3d_to_tlns3d.F  mapblk.F      movseg.F     spltsg.F
      getdhdr.F          mapin.F       outgrd.F     termn8.F
      ipadd.F            mapout.F      renmbr.F     tlns3d_to_cfl3d.F
      ipfix.F            mapsplt.F     shortinp.F   splitter.F
      iprvsg.F           mkintr.F      sizer.F      
      mkseg.F            splt2.F       ingrd.F
      main.F
      ${CFLSRC_D}/umalloc.F 
      ${CFLSRC_S}/parser.F
      ${CFLSRC_S}/readkey.F 
      ${CFLSRC_S}/outbuf.F
      ${CFLSRC_S}/my_flush.F)
    add_executable(splitter ${SPLIT_SRC})
  endif(BUILD_REAL)
  if(BUILD_CMPLX)
    set(SPLIT_CMPLX_SRC
      cfl3d_to_tlns3d.F  mapblk.F      movseg.F        spltsg.F
      getdhdr.F          mapin.F       outgrd_cmplx.F  termn8.F
      ipadd.F            mapout.F      renmbr.F        tlns3d_to_cfl3d.F
      ipfix.F            mapsplt.F     shortinp.F      splitter_cmplx.F 
      iprvsg.F           mkintr.F      sizer.F         
      mkseg.F            splt2.F       ingrd_cmplx.F 
      main_cmplx.F
      ${CFLSRC_D}/umalloc.F
      ${CFLSRC_D}/umalloc_c.F
      ${CFLSRC_S}/parser.F
      ${CFLSRC_S}/readkey.F
      ${CFLSRC_S}/outbuf.F
      ${CFLSRC_S}/my_flush.F)
    add_executable(splittercmplx ${SPLIT_CMPLX_SRC})
    target_compile_definitions(splittercmplx PUBLIC CMPLX)
  endif(BUILD_CMPLX)
endif(BUILD_SPLITTER)
