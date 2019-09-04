c  ---------------------------------------------------------------------------
c  CFL3D is a structured-grid, cell-centered, upwind-biased, Reynolds-averaged
c  Navier-Stokes (RANS) code. It can be run in parallel on multiple grid zones
c  with point-matched, patched, overset, or embedded connectivities. Both
c  multigrid and mesh sequencing are available in time-accurate or
c  steady-state modes.
c
c  Copyright 2001 United States Government as represented by the Administrator
c  of the National Aeronautics and Space Administration. All Rights Reserved.
c 
c  The CFL3D platform is licensed under the Apache License, Version 2.0 
c  (the "License"); you may not use this file except in compliance with the 
c  License. You may obtain a copy of the License at 
c  http://www.apache.org/licenses/LICENSE-2.0. 
c 
c  Unless required by applicable law or agreed to in writing, software 
c  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
c  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
c  License for the specific language governing permissions and limitations 
c  under the License.
c  ---------------------------------------------------------------------------
c
      subroutine bc1000(jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,bcj,bck,
     .                  bci,xtbj,xtbk,xtbi,atbj,atbk,atbi,ista,iend,
     .                  jsta,jend,ksta,kend,nface,tursav,tj0,tk0,
     .                  ti0,vist3d,vj0,vk0,vi0,iuns,nou,bou,nbuf,
     .                  ibufdim,nummem,x,y,z,nbl,myid)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Set freestream boundary conditions 
c***********************************************************************
c
      use module_gust
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension q(jdim,kdim,idim,5), qi0(jdim,kdim,5,4),
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4)
      dimension sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),
     .          sj(jdim,kdim,idim-1,5)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension xtbj(kdim,idim-1,3,2),xtbk(jdim,idim-1,3,2),
     .          xtbi(jdim,kdim,3,2),atbj(kdim,idim-1,3,2),
     .          atbk(jdim,idim-1,3,2),atbi(jdim,kdim,3,2)
      dimension tursav(jdim,kdim,idim,nummem),tj0(kdim,idim-1,nummem,4),
     .          tk0(jdim,idim-1,nummem,4),ti0(jdim,kdim,nummem,4),
     .          vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4),vist3d(jdim,kdim,idim)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension qgv(5)
c
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /gust/ igust
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
c            * * * * * * * * * * * * * * * * * * * * * *
c            * standard boundary condition bctype=1000 *
c            * * * * * * * * * * * * * * * * * * * * * *
c
c******************************************************************************
c      j=1 boundary            freestream                           bctype 1000
c******************************************************************************
c
      if (nface.eq.3) then
      do 100 l=1,5
      do 100 i=ista,iend1
      do 100 k=ksta,kend1
      qj0(k,i,l,1) = qiv(l)
      qj0(k,i,l,2) = qj0(k,i,l,1)
      bcj(k,i,1)   = 0.0
  100 continue
      if (igust .eq. 1) then
        do 112 i=ista,iend1
        do 112 k=ksta,kend1
          xc = (x(1,k,i)+x(1,k+1,i)+x(1,k+1,i+1)+x(1,k,i+1))/4
          yc = (y(1,k,i)+y(1,k+1,i)+y(1,k+1,i+1)+y(1,k,i+1))/4
          zc = (z(1,k,i)+z(1,k+1,i)+z(1,k+1,i+1)+z(1,k,i+1))/4
          call gust_vel(xc,yc,zc,time,qiv,ug,vg,wg,iflag)
          if (iflag.ne.0)then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' stopping in ',
     .        'bc_1000 block=',nbl,' j=1: ', get_err_msg(iflag)
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
          end if
          qj0(k,i,2,1) = qj0(k,i,2,1)+ug
          qj0(k,i,2,2) = qj0(k,i,2,2)+ug
          qj0(k,i,3,1) = qj0(k,i,3,1)+vg
          qj0(k,i,3,2) = qj0(k,i,3,2)+vg
          qj0(k,i,4,1) = qj0(k,i,4,1)+wg
          qj0(k,i,4,2) = qj0(k,i,4,2)+wg
  112   continue
      end if ! igust .eq. 1
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 191 i=ista,iend1
        do 191 k=ksta,kend1
          vj0(k,i,1,1) = vist3d(1,k,i)
          vj0(k,i,1,2) = 0.
  191   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 101 i=ista,iend1
        do 101 k=ksta,kend1
          tj0(k,i,l,1) = tur10(l)
          tj0(k,i,l,2) = tur10(l)
  101   continue
        enddo
      end if
      end if
      end if
c
c******************************************************************************
c      j=jdim boundary             freestream                       bctype 1000
c******************************************************************************
c
      if (nface.eq.4) then
      do 200 l=1,5
      do 200 i=ista,iend1
      do 200 k=ksta,kend1
      qj0(k,i,l,3) = qiv(l)
      qj0(k,i,l,4) = qj0(k,i,l,3)
      bcj(k,i,2)   = 0.0
  200 continue
      if (igust .eq. 1) then
      do 212 i=ista,iend1
      do 212 k=ksta,kend1
        xc = (x(jdim,k,i)+x(jdim,k+1,i)+x(jdim,k+1,i+1)+x(jdim,k,i+1))/4
        yc = (y(jdim,k,i)+y(jdim,k+1,i)+y(jdim,k+1,i+1)+y(jdim,k,i+1))/4
        zc = (z(jdim,k,i)+z(jdim,k+1,i)+z(jdim,k+1,i+1)+z(jdim,k,i+1))/4
        call gust_vel(xc,yc,zc,time,qiv,ug,vg,wg,iflag)
        if (iflag.ne.0)then
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),*)' stopping in ',
     .      'bc_1000 block=',nbl,' jdim: ', get_err_msg(iflag)
          call termn8(myid,-1,ibufdim,nbuf,bou,nou)
        end if
        qj0(k,i,2,3) = qj0(k,i,2,3)+ug
        qj0(k,i,2,4) = qj0(k,i,2,4)+ug
        qj0(k,i,3,3) = qj0(k,i,3,3)+vg
        qj0(k,i,3,4) = qj0(k,i,3,4)+vg
        qj0(k,i,4,3) = qj0(k,i,4,3)+wg
        qj0(k,i,4,4) = qj0(k,i,4,4)+wg
  212 continue
      end if ! igust .eq. 1
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 291 i=ista,iend1
        do 291 k=ksta,kend1
          vj0(k,i,1,3) = vist3d(jdim-1,k,i)
          vj0(k,i,1,4) = 0.0
  291   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 201 i=ista,iend1
        do 201 k=ksta,kend1
          tj0(k,i,l,3) = tur10(l)
          tj0(k,i,l,4) = tur10(l)
  201   continue
        enddo
      end if
      end if
      end if
c
c******************************************************************************
c      k=1 boundary            freestream                           bctype 1000
c******************************************************************************
c
      if (nface.eq.5) then
      do 300 l=1,5
      do 300 i=ista,iend1
      do 300 j=jsta,jend1
      qk0(j,i,l,1) = qiv(l)
      qk0(j,i,l,2) = qk0(j,i,l,1)
      bck(j,i,1)   = 0.0
  300 continue
      if (igust .eq. 1) then
      do 312 i=ista,iend1
      do 312 j=jsta,jend1
        xc = (x(j,1,i)+x(j+1,1,i)+x(j+1,1,i+1)+x(j,1,i+1))/4
        yc = (y(j,1,i)+y(j+1,1,i)+y(j+1,1,i+1)+y(j,1,i+1))/4
        zc = (z(j,1,i)+z(j+1,1,i)+z(j+1,1,i+1)+z(j,1,i+1))/4
        call gust_vel(xc,yc,zc,time,qiv,ug,vg,wg,iflag)
        if (iflag.ne.0)then
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),*)' stopping in ',
     .      'bc_1000 block=',nbl,' k=1: ', get_err_msg(iflag)
          call termn8(myid,-1,ibufdim,nbuf,bou,nou)
        end if
        qk0(j,i,2,1) = qk0(j,i,2,1)+ug
        qk0(j,i,2,2) = qk0(j,i,2,2)+ug
        qk0(j,i,3,1) = qk0(j,i,3,1)+vg
        qk0(j,i,3,2) = qk0(j,i,3,2)+vg
        qk0(j,i,4,1) = qk0(j,i,4,1)+wg
        qk0(j,i,4,2) = qk0(j,i,4,2)+wg
  312 continue
      end if ! igust .eq. 1
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 391 i=ista,iend1
        do 391 j=jsta,jend1
          vk0(j,i,1,1) = vist3d(j,1,i)
          vk0(j,i,1,2) = 0.0
  391   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 301 i=ista,iend1
        do 301 j=jsta,jend1
          tk0(j,i,l,1) = tur10(l)
          tk0(j,i,l,2) = tur10(l)
  301   continue
        enddo
      end if
      end if
      end if
c
c******************************************************************************
c      k=kdim boundary             freestream                       bctype 1000
c******************************************************************************
c
      if (nface.eq.6) then
      do 400 l=1,5
      do 400 i=ista,iend1
      do 400 j=jsta,jend1
      qk0(j,i,l,3) = qiv(l)
      qk0(j,i,l,4) = qk0(j,i,l,3)
      bck(j,i,2)   = 0.0
  400 continue
      if (igust .eq. 1) then
      do 412 i=ista,iend1
      do 412 j=jsta,jend1
        xc = (x(j,kdim,i)+x(j+1,kdim,i)+x(j+1,kdim,i+1)+x(j,kdim,i+1))/4
        yc = (y(j,kdim,i)+y(j+1,kdim,i)+y(j+1,kdim,i+1)+y(j,kdim,i+1))/4
        zc = (z(j,kdim,i)+z(j+1,kdim,i)+z(j+1,kdim,i+1)+z(j,kdim,i+1))/4
        call gust_vel(xc,yc,zc,time,qiv,ug,vg,wg,iflag)
        if (iflag.ne.0)then
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),*)' stopping in ',
     .      'bc_1000 block=',nbl,' kdim: ', get_err_msg(iflag)
          call termn8(myid,-1,ibufdim,nbuf,bou,nou)
        end if
        qk0(j,i,2,3) = qk0(j,i,2,3)+ug
        qk0(j,i,2,4) = qk0(j,i,2,4)+ug
        qk0(j,i,3,3) = qk0(j,i,3,3)+vg
        qk0(j,i,3,4) = qk0(j,i,3,4)+vg
        qk0(j,i,4,3) = qk0(j,i,4,3)+wg
        qk0(j,i,4,4) = qk0(j,i,4,4)+wg
  412 continue
      end if ! igust .eq. 1
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 491 i=ista,iend1
        do 491 j=jsta,jend1
          vk0(j,i,1,3) = vist3d(j,kdim-1,i)
          vk0(j,i,1,4) = 0.0
  491   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 401 i=ista,iend1
        do 401 j=jsta,jend1
          tk0(j,i,l,3) = tur10(l)
          tk0(j,i,l,4) = tur10(l)
  401   continue
        enddo
      end if
      end if
      end if
c
c******************************************************************************
c      i=1 boundary            freestream                           bctype 1000
c******************************************************************************
c
      if (nface.eq.1) then
      do 500 l=1,5
      do 500 k=ksta,kend1
      do 500 j=jsta,jend1
      qi0(j,k,l,1) = qiv(l)
      qi0(j,k,l,2) = qi0(j,k,l,1)
      bci(j,k,1)   = 0.0
  500 continue
      if (igust .eq. 1) then
        do 512 k=ksta,kend1
        do 512 j=jsta,jend1
        xc = (x(j,k,1)+x(j+1,k,1)+x(j+1,k+1,1)+x(j,k+1,1))/4
        yc = (y(j,k,1)+y(j+1,k,1)+y(j+1,k+1,1)+y(j,k+1,1))/4
        zc = (z(j,k,1)+z(j+1,k,1)+z(j+1,k+1,1)+z(j,k+1,1))/4
        call gust_vel(xc,yc,zc,time,qiv,ug,vg,wg,iflag)
        if (iflag.ne.0)then
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),*)' stopping in ',
     .      'bc_1000 block=',nbl,' i=1: ', get_err_msg(iflag)
          call termn8(myid,-1,ibufdim,nbuf,bou,nou)
        end if
        qi0(j,k,2,1) = qi0(j,k,2,1)+ug
        qi0(j,k,2,2) = qi0(j,k,2,2)+ug
        qi0(j,k,3,1) = qi0(j,k,3,1)+vg
        qi0(j,k,3,2) = qi0(j,k,3,2)+vg
        qi0(j,k,4,1) = qi0(j,k,4,1)+wg
        qi0(j,k,4,2) = qi0(j,k,4,2)+wg
  512 continue
      end if ! igust .eq. 1
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 591 k=ksta,kend1
        do 591 j=jsta,jend1
          vi0(j,k,1,1) = vist3d(j,k,1)
          vi0(j,k,1,2) = 0.0
  591   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 501 k=ksta,kend1
        do 501 j=jsta,jend1
          ti0(j,k,l,1) = tur10(l)
          ti0(j,k,l,2) = tur10(l)
  501   continue
        enddo
      end if
      end if
      end if
c
c******************************************************************************
c      i=idim boundary             freestream                       bctype 1000
c******************************************************************************
c
      if (nface.eq.2) then
      do 600 l=1,5
      do 600 k=ksta,kend1
      do 600 j=jsta,jend1
      qi0(j,k,l,3) = qiv(l)
      qi0(j,k,l,4) = qi0(j,k,l,3)
      bci(j,k,2)   = 0.0
  600 continue
      if (igust .eq. 1) then
        do 612 k=ksta,kend1
        do 612 j=jsta,jend1
        xc = (x(j,k,idim)+x(j+1,k,idim)+x(j+1,k+1,idim)+x(j,k+1,idim))/4
        yc = (y(j,k,idim)+y(j+1,k,idim)+y(j+1,k+1,idim)+y(j,k+1,idim))/4
        zc = (z(j,k,idim)+z(j+1,k,idim)+z(j+1,k+1,idim)+z(j,k+1,idim))/4
        call gust_vel(xc,yc,zc,time,qiv,ug,vg,wg,iflag)
        if (iflag.ne.0)then
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),*)' stopping in ',
     .      'bc_1000 block=',nbl,' idim: ', get_err_msg(iflag)
          call termn8(myid,-1,ibufdim,nbuf,bou,nou)
        end if
        qi0(j,k,2,3) = qi0(j,k,2,3)+ug
        qi0(j,k,2,4) = qi0(j,k,2,4)+ug
        qi0(j,k,3,3) = qi0(j,k,3,3)+vg
        qi0(j,k,3,4) = qi0(j,k,3,4)+vg
        qi0(j,k,4,3) = qi0(j,k,4,3)+wg
        qi0(j,k,4,4) = qi0(j,k,4,4)+wg
  612 continue
      end if ! igust .eq. 1
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 691 k=ksta,kend1
        do 691 j=jsta,jend1
          vi0(j,k,1,3) = vist3d(j,k,idim-1)
          vi0(j,k,1,4) = 0.0
  691   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 601 k=ksta,kend1
        do 601 j=jsta,jend1
          ti0(j,k,l,3) = tur10(l)
          ti0(j,k,l,4) = tur10(l)
  601   continue
        enddo
      end if
      end if
      end if
c
      return
      end
