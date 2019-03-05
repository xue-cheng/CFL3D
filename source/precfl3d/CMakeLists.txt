if(BUILD_PRECFL3D)
  set(PRECFL_FSRCS
    get_cmd_args.F
    main.F
    usrint.F
    ${CFLSRC_D}/sizer.F
    ${CFLSRC_D}/pointers.F
    ${CFLSRC_D}/compg2n.F
    ${CFLSRC_D}/umalloc.F
    ${CFLSRC_D}/termn8.F
    ${CFLSRC_S}/outbuf.F
    ${CFLSRC_S}/pre_patch.F
    ${CFLSRC_S}/pre_blockbc.F
    ${CFLSRC_S}/pre_period.F
    ${CFLSRC_S}/getdhdr.F
    ${CFLSRC_S}/pre_embed.F
    ${CFLSRC_S}/global.F
    ${CFLSRC_S}/global2.F
    ${CFLSRC_S}/rpatch0.F
    ${CFLSRC_S}/getibk0.F
    ${CFLSRC_S}/cntsurf.F
    ${CFLSRC_S}/lead.F
    ${CFLSRC_S}/global0.F
    ${CFLSRC_S}/readkey.F
    ${CFLSRC_S}/parser.F
    ${CFLSRC_S}/ccomplex.F
    ${CFLSRC_S}/cgnstools.F
    ${CFLSRC_S}/setseg.F
    ${CFLSRC_S}/my_flush.F
    ${CFLSRC_S}/module_kwstm.F90
    ${CFLSRC_S}/module_profileout.F90
    )
  add_executable(precfl3d ${PRECFL_FSRCS})
endif(BUILD_PRECFL3D)


