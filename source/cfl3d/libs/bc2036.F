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
      subroutine bc2036(jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,bcj,bck,
     .                  bci,xtbj,xtbk,xtbi,atbj,atbk,atbi,ista,iend,
     .                  jsta,jend,ksta,kend,nface,tursav,tj0,tk0,
     .                  ti0,smin,vist3d,vj0,vk0,vi0,mdim,ndim,bcdata,
     .                  filname,iuns,irelv,snj0,snk0,sni0,ntime,
     .                  snjm,snkm,snim,nou,bou,nbuf,ibufdim,myid,
     .                  nummem)
c
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Set pulsed wall jet BC (viscous wall)
c               with 9 pieces of auxiliary data:
c
c     This is a mod to bc2016, by XUE Cheng
c
c           1) Twtype:
c              Twtype > 0 (fixed wall temperature Tw/Tinf = Twtype)
c              Twtype = 0 (adiabatic wall)
c              Twtype < 0 (fixed wall temperature = stagnation temp)
c
c           2) cq0  = base component of surface jet
c           3) cqu = unsteady component of surface jet
c              Set surface jet boundary conditions via mass-flow 
c              coefficient:
c                                             (rho * u)_jet
c              mass-flow coefficient cqt  = ------------------ 
c                                           (rho * u)_infinity
c           4) sjetx = direction number for surface jet (x-direction)
c           5) sjety = direction number for surface jet (y-direction)
c           6) sjetz = direction number for surface jet (z-direction)
c           7) t_cyc = period or start time, see d_cyc
c           8) t_s99 = time of the pulsed jet achieving its maximum value
c           9) d_cyc :
c              d_cyc > 0 pulsed jet
c              d_cyc = 0 one-shot jet
c              d_cyc < 0 one-minus-cosine
c             a) pulsed jet (d_cyc>0):
c                d_cyc = duty cycle of the pulse between (0, 1)
c                t_cyc = period of unsteady jet
c                t_s99 = time of the pulsed jet achieving its maximum value
c                  t_s99 = 0 => step function
c                  t_s99 > 0 => tanh function
c                  t_s99 < 0 => linear function
c             b) one-shot jet (d_cyc = 0):
c                t < t_cyc       => cqt = cq0
c                t > t_cyc+t_s99 => cqt = cq0+cqu
c                t_cyc = jet start time
c                t_s99 = time of the pulsed jet achieving its maximum value
c                  t_s99 = 0 => step function
c                  t_s99 > 0 => tanh function
c                  t_s99 < 0 => linear function
c             c) one-minus-cosine jet (d_cyc < 0):
c                d_cyc = -1 => continuous
c                d_cyc = -2 => single wave
c                t_s99 = half period
c                t_cyc = jet start time
c                if t < t_cyc:
c                  cqt = cq0
c                else:
c                  cqt = cq0 + 0.5(1-cos)*cqu
c              NOTES:
c              1.  input direction numbers do not have to be normalized
c              2.  Set (sjetx,sjety,sjetz) = (0.0,0.0,0.0) for a
c                  wall-normal surface jet.  In this case,
c                  cq > 0 ==> injection
c                  cq < 0 ==> suction
c
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
      character*80 filname
c
      dimension nou(nbuf)
      dimension q(jdim,kdim,idim,5), qi0(jdim,kdim,5,4),
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),
     .          sj(jdim,kdim,idim-1,5)
      dimension xtbj(kdim,idim-1,3,2),xtbk(jdim,idim-1,3,2),
     .          xtbi(jdim,kdim,3,2),atbj(kdim,idim-1,3,2),
     .          atbk(jdim,idim-1,3,2),atbi(jdim,kdim,3,2)
      dimension bcdata(mdim,ndim,2,12)
      dimension tursav(jdim,kdim,idim,nummem),tj0(kdim,idim-1,nummem,4),
     .          tk0(jdim,idim-1,nummem,4),ti0(jdim,kdim,nummem,4),
     .          smin(jdim-1,kdim-1,idim-1),
     .          vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4),vist3d(jdim,kdim,idim)
      dimension snj0(jdim-1,kdim-1,idim-1),snk0(jdim-1,kdim-1,idim-1),
     .          sni0(jdim-1,kdim-1,idim-1)
      dimension snjm(jdim-1,kdim-1,idim-1),snkm(jdim-1,kdim-1,idim-1),
     .          snim(jdim-1,kdim-1,idim-1)
      dimension a0(7),a1(5),a2(3)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fluid2/ pr,prt,cbar
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /maxiv/ ivmx
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
      common /wallfun/ iwf(3)
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /walljet/ turwj(7)
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)

c   Check inputs

      if (ivmx.ge.2.and.(ivmx.lt.5.or.ivmx.gt.7)) then 
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),*)' stopping in bc2036...',
     .    ' invalid turb. model, ivmx=', ivmx, 
     .    ' must be on of 5, 6 or 7'
        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      endif

      if ((nface.eq.1.or.nface.eq.2).and.(iwf(1).ne.0).or.
     .    (nface.eq.3.or.nface.eq.4).and.(iwf(2).ne.0).or.
     .    (nface.eq.5.or.nface.eq.6).and.(iwf(3).ne.0)) then 
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),*)' stopping in bc2036...',
     .    ' wall function not implemented'
        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      endif

c   *********
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
c     this bc makes use of only one plane of data    
c
      ip    = 1
c
c            * * * * * * * * * * * * * * * * * * * * * *
c            * standard boundary condition bctype=2036 *
c            * * * * * * * * * * * * * * * * * * * * * *
c
c******************************************************************************
c      j=1 boundary        viscous wall with T & cq specified       bctype 2036
c******************************************************************************
      if (nface.eq.3) then
c
      do 400 i=ista,iend1
      ii = i-ista+1
c
      do 300 k=ksta,kend1
      kk = k-ksta+1
c     read data
      cq0           = bcdata(kk,ii,ip,2)
      cqa           = bcdata(kk,ii,ip,3)
      sjetx         = bcdata(kk,ii,ip,4)
      sjety         = bcdata(kk,ii,ip,5)
      sjetz         = bcdata(kk,ii,ip,6)
      t_cyc         = bcdata(kk,ii,ip,7)
      t_s99         = bcdata(kk,ii,ip,8)
      d_cyc         = bcdata(kk,ii,ip,9)

      call bc2036_check(t_cyc,t_s99,d_cyc,nou,bou,ibufdim,nbuf)
      
c
c     surface velocities
c
      uub = 0.
      vvb = 0.
      wwb = 0.
c
c     for dynamic mesh, set velocity at wall to grid velocity at wall
c     if irelv > 0; otherwise, set to zero
c
      if (iuns.gt.0 .and. irelv.gt.0) then
        uub = xtbj(k,i,1,1)
        vvb = xtbj(k,i,2,1)
        wwb = xtbj(k,i,3,1)
      end if

c
c     compute total (steady and unsteady) surface jet
c
      cqt = cq2036(time, cq0, cqa, t_cyc, t_s99, d_cyc)

c     compute direction cosines for surface jet
c
      sjet   = sqrt(sjetx*sjetx+sjety*sjety+sjetz*sjetz)
      snx    = sj(1,k,i,1)
      sny    = sj(1,k,i,2)
      snz    = sj(1,k,i,3)
      if (sjet.ne.0.0) then
        sjetx = sjetx/sjet
        sjety = sjety/sjet
        sjetz = sjetz/sjet
      else
        sjetx = snx 
        sjety = sny 
        sjetz = snz 
      endif

c     wall normal component
c     > 0 means inflow
      cqn = cqt*(sjetx*snx+sjety*sny+sjetz*snz)
c
      pb  = q(1,k,i,5)
      c2  = gamma*q(1,k,i,5)/q(1,k,i,1)
      if (real(cqn).lt.1e-12) then ! without jet
        dpb = q(2,k,i,5)-q(1,k,i,5)
        pb  = pb - dpb/2.0
        if (real(pb).le.0.0) pb = q(1,k,i,5)
        if (real(bcdata(kk,ii,ip,1)) .gt. 0.) then
          c2  = bcdata(kk,ii,ip,1)
        else if (real(bcdata(kk,ii,ip,1)) .lt. 0.) then
          c2=1.e0+gm1*0.5e0*xmach*xmach
        else
          xm2 = (q(1,k,i,2)-uub)**2+
     +          (q(1,k,i,3)-vvb)**2+
     +          (q(1,k,i,4)-wwb)**2
          xm2 = xm2/c2
          c2  = c2*(1.+0.5*gm1*xm2)
        end if
      endif
c     add directional components of surface jet
c
      uub = uub + sjetx*xmach*cqt*c2/(gamma*pb)
      vvb = vvb + sjety*xmach*cqt*c2/(gamma*pb)
      wwb = wwb + sjetz*xmach*cqt*c2/(gamma*pb)
c
      qj0(k,i,1,1) = gamma*pb/c2
      qj0(k,i,2,1) = uub
      qj0(k,i,3,1) = vvb
      qj0(k,i,4,1) = wwb
      qj0(k,i,5,1) = pb
      bcj(k,i,1)   = 1.0
c
c     f23 = 0.0  -  2-point extrapolation
c           1.0  -  3-point extrapolation
c
      f23 = 0.0
c
      j2 = min(2,jdim1)
      if (j2.eq.1) f23 = 0.0
c
      z1  =   2.0 +1.5*f23
      z2  =       -0.5*f23   
      z3  = -(2.0 +    f23)
c
      qj0(k,i,1,2) = z1*q(1,k,i,1) + z2*q(j2,k,i,1) + z3*qj0(k,i,1,1)
      qj0(k,i,2,2) = z1*q(1,k,i,2) + z2*q(j2,k,i,2) + z3*qj0(k,i,2,1)
      qj0(k,i,3,2) = z1*q(1,k,i,3) + z2*q(j2,k,i,3) + z3*qj0(k,i,3,1)
      qj0(k,i,4,2) = z1*q(1,k,i,4) + z2*q(j2,k,i,4) + z3*qj0(k,i,4,1)
      qj0(k,i,5,2) = z1*q(1,k,i,5) + z2*q(j2,k,i,5) + z3*qj0(k,i,5,1)
c 
c     Turbulence
c     only need to do advanced model turbulence B.C.s on finest grid
c
      if (ivmx.ge.2.and.level.ge.lglobal.and.ntime.ne.0) then
        if (abs(cqt).lt.1e-12) then
          vj0(k,i,1,1) = 0.
          vj0(k,i,1,2) = 0.
          if (ivmx.eq.5) then 
            tj0(k,i,1,1) = -tursav(1,k,i,1)
            tj0(k,i,2,1) = -tursav(1,k,i,2)
            tj0(k,i,1,2) = 2.*tj0(k,i,1,1)-tursav(1,k,i,1)
            tj0(k,i,2,2) = 2.*tj0(k,i,2,1)-tursav(1,k,i,2)
          elseif(ivmx.eq.6.or.ivmx.eq.7) then
            c2b=cbar/tinf
            c2bp=c2b+1.0
            re=reue/xmach
            beta1=.075
            j=1
            tt=gamma*qj0(k,i,5,1)/qj0(k,i,1,1)
            fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
            dist=ccabs(smin(j,k,i))
            tj0(k,i,1,1) = 
     =        2.*(60.*fnu/(re**2*q(j,k,i,1)*beta1*dist**2))-
     +        tursav(j,k,i,1)
            tj0(k,i,2,1) = -tursav(j,k,i,2)
            tj0(k,i,1,2) = 2.*tj0(k,i,1,1)-tursav(1,k,i,1)
            tj0(k,i,2,2) = 2.*tj0(k,i,2,1)-tursav(1,k,i,2)
          endif
        else ! use free stream val or user specified via turwj
          vj0(k,i,1,1) = vist3d(1,k,i)
          vj0(k,i,1,0) = 0.0
          do l=1,nummem
            tj0(k,i,l,1) = ccmax(tur10(l), turwj(l))
            tj0(k,i,l,2) = tj0(k,i,l,1)
          enddo
        endif
      endif
  300 continue
  400 continue
c
      end if
c
c******************************************************************************
c      j=jdim boundary     viscous wall with T & cq specified       bctype 2036
c******************************************************************************
      if (nface.eq.4) then
c
      do 800 i=ista,iend1
      ii = i-ista+1
c
      do 700 k=ksta,kend1
      kk = k-ksta+1
c     read data
      cq0           = bcdata(kk,ii,ip,2)
      cqa           = bcdata(kk,ii,ip,3)
      sjetx         = bcdata(kk,ii,ip,4)
      sjety         = bcdata(kk,ii,ip,5)
      sjetz         = bcdata(kk,ii,ip,6)
      t_cyc         = bcdata(kk,ii,ip,7)
      t_s99         = bcdata(kk,ii,ip,8)
      d_cyc         = bcdata(kk,ii,ip,9)
      
      call bc2036_check(t_cyc,t_s99,d_cyc,nou,bou,ibufdim,nbuf)

c
c     surface velocities
c
      uub = 0.
      vvb = 0.
      wwb = 0.
c
c     for dynamic mesh, set velocity at wall to grid velocity at wall
c     if irelv > 0; otherwise, set to zero
c
      if (iuns.gt.0 .and. irelv.gt.0) then
        uub = xtbj(k,i,1,2)
        vvb = xtbj(k,i,2,2)
        wwb = xtbj(k,i,3,2)
      end if

c
c     compute total (steady and unsteady) surface jet
c
      cqt = cq2036(time, cq0, cqa, t_cyc, t_s99, d_cyc)

c     compute direction cosines for surface jet
c
      sjet = sqrt(sjetx*sjetx+sjety*sjety+sjetz*sjetz)
      snx  = -sj(jdim,k,i,1)
      sny  = -sj(jdim,k,i,2)
      snz  = -sj(jdim,k,i,3)
      if (sjet.ne.0.0) then
        sjetx = sjetx/sjet
        sjety = sjety/sjet
        sjetz = sjetz/sjet
      else
        sjetx = snx
        sjety = sny
        sjetz = snz
      endif
c     wall normal component
c     > 0 means inflow
      cqn = cqt*(sjetx*snx+sjety*sny+sjetz*snz)
c
      pb  = q(jdim1,k,i,5)
      c2  = gamma*q(jdim1,k,i,5)/q(jdim1,k,i,1)
      if (real(cqn).lt.1e-12) then ! without jet
        dpb = q(jdim1,k,i,5)-q(jdim1-1,k,i,5)
        pb  = pb + dpb/2.0
        if (real(pb).le.0.0)  pb = q(jdim1,k,i,5)
        if (real(bcdata(kk,ii,ip,1)) .gt. 0.) then
          c2  = bcdata(kk,ii,ip,1)
        else if (real(bcdata(kk,ii,ip,1)) .lt. 0.) then
          c2=1.e0+gm1*0.5e0*xmach*xmach
        else
          xm2 = (q(jdim1,k,i,2)-uub)**2+
     +          (q(jdim1,k,i,3)-vvb)**2+
     +          (q(jdim1,k,i,4)-wwb)**2
          xm2 = xm2/c2
          c2  = c2*(1.+0.5*gm1*xm2)
        end if
      endif

c     add directional components of surface jet
c
      uub = uub + sjetx*xmach*cqt*c2/(gamma*pb)
      vvb = vvb + sjety*xmach*cqt*c2/(gamma*pb)
      wwb = wwb + sjetz*xmach*cqt*c2/(gamma*pb)
c
      qj0(k,i,1,3) = gamma*pb/c2
      qj0(k,i,2,3) = uub
      qj0(k,i,3,3) = vvb
      qj0(k,i,4,3) = wwb
      qj0(k,i,5,3) = pb
      bcj(k,i,2)   = 1.0
c
c     f23 = 0.0  -  2-point extrapolation
c           1.0  -  3-point extrapolation
c
      f23 = 0.0
c
      j2 = max(1,jdim-2)
      if (j2.eq.1) f23 = 0.0
c
      z1  =  -2.0 -1.5*f23
      z2  =       +0.5*f23   
      z3  = +(2.0 +    f23)
c
      qj0(k,i,1,4) = z1*q(jdim1,k,i,1)+z2*q(j2,k,i,1)+z3*qj0(k,i,1,3)
      qj0(k,i,2,4) = z1*q(jdim1,k,i,2)+z2*q(j2,k,i,2)+z3*qj0(k,i,2,3)
      qj0(k,i,3,4) = z1*q(jdim1,k,i,3)+z2*q(j2,k,i,3)+z3*qj0(k,i,3,3)
      qj0(k,i,4,4) = z1*q(jdim1,k,i,4)+z2*q(j2,k,i,4)+z3*qj0(k,i,4,3)
      qj0(k,i,5,4) = z1*q(jdim1,k,i,5)+z2*q(j2,k,i,5)+z3*qj0(k,i,5,3)
c 
c     Turbulence
c     only need to do advanced model turbulence B.C.s on finest grid
c
      if (ivmx.ge.2.and.level.ge.lglobal.and.ntime.ne.0) then 
        if (abs(cqt).lt.1e-12) then
          vj0(k,i,1,3) = 0.
          vj0(k,i,1,4) = 0.
          if (ivmx.eq.5) then 
            tj0(k,i,1,3) = -tursav(jdim-1,k,i,1)
            tj0(k,i,2,3) = -tursav(jdim-1,k,i,2)
            tj0(k,i,1,4) = 2.*tj0(k,i,1,3)-tursav(jdim-1,k,i,1)
            tj0(k,i,2,4) = 2.*tj0(k,i,2,3)-tursav(jdim-1,k,i,2)
          elseif(ivmx.eq.6.or.ivmx.eq.7) then
            c2b=cbar/tinf
            c2bp=c2b+1.0
            re=reue/xmach
            beta1=.075
            j=jdim-1
            tt=gamma*qj0(k,i,5,3)/qj0(k,i,1,3)
            fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
            dist=ccabs(smin(j,k,i))
            tj0(k,i,1,3) = 
     =        2.*(60.*fnu/(re**2*q(j,k,i,1)*beta1*dist**2))-
     +        tursav(j,k,i,1)
            tj0(k,i,2,3) = -tursav(j,k,i,2)
            tj0(k,i,1,4) = 2.*tj0(k,i,1,3)-tursav(jdim-1,k,i,1)
            tj0(k,i,2,4) = 2.*tj0(k,i,2,3)-tursav(jdim-1,k,i,2)
          endif
        else
          vj0(k,i,1,3) = vist3d(jdim-1,k,i)
          vj0(k,i,1,4) = 0.
          do l=1,nummem
            tj0(k,i,l,3) = ccmax(tur10(l), turwj(l))
            tj0(k,i,l,4) = tj0(k,i,l,3)
          enddo
        endif
      endif
  700 continue
  800 continue

c
      end if
c
c******************************************************************************
c      k=1 boundary        viscous wall with T & cq specified       bctype 2036
c******************************************************************************
      if (nface.eq.5) then
c
      do 1200 i=ista,iend1
      ii = i-ista+1
c
      do 1100 j=jsta,jend1
      jj = j-jsta+1
c     read data
      cq0           = bcdata(jj,ii,ip,2)
      cqa           = bcdata(jj,ii,ip,3)
      sjetx         = bcdata(jj,ii,ip,4)
      sjety         = bcdata(jj,ii,ip,5)
      sjetz         = bcdata(jj,ii,ip,6)
      t_cyc         = bcdata(jj,ii,ip,7)
      t_s99         = bcdata(jj,ii,ip,8)
      d_cyc         = bcdata(jj,ii,ip,9)
      
      call bc2036_check(t_cyc,t_s99,d_cyc,nou,bou,ibufdim,nbuf)

c
c     surface velocities
c
      uub = 0.
      vvb = 0.
      wwb = 0.
c
c     for dynamic mesh, set velocity at wall to grid velocity at wall
c     if irelv > 0; otherwise, set to zero
c
      if (iuns.gt.0 .and. irelv.gt.0) then
      uub = xtbk(j,i,1,1)
      vvb = xtbk(j,i,2,1)
      wwb = xtbk(j,i,3,1)
      end if

c
c     compute total (steady and unsteady) surface jet
c
      cqt = cq2036(time, cq0, cqa, t_cyc, t_s99, d_cyc)

c     compute direction cosines for surface jet
c
      sjet   = sqrt(sjetx*sjetx+sjety*sjety+sjetz*sjetz)
      snx    = sk(j,1,i,1)
      sny    = sk(j,1,i,2)
      snz    = sk(j,1,i,3)
      if (sjet.ne.0.0) then
        sjetx = sjetx/sjet
        sjety = sjety/sjet
        sjetz = sjetz/sjet
      else
        sjetx = snx
        sjety = sny
        sjetz = snz
      endif

c     wall normal component
c     > 0 means inflow
      cqn = cqt*(sjetx*snx+sjety*sny+sjetz*snz)

      pb            = q(j,1,i,5)
      c2            = gamma*q(j,1,i,5)/q(j,1,i,1)
      if (real(cqn).lt.1e-12) then ! without jet
        dpb           = q(j,2,i,5)-q(j,1,i,5)
        pb            = pb - dpb/2.0
        if (real(pb).le.0.0)  pb = q(j,1,i,5)
        if (real(bcdata(jj,ii,ip,1)) .gt. 0.) then
          c2  = bcdata(jj,ii,ip,1)
        else if (real(bcdata(jj,ii,ip,1)) .lt. 0.) then
          c2=1.e0+gm1*0.5e0*xmach*xmach
        else
          if (iuns.gt.0 .and. irelv.gt.0) then
          xm2 = (q(j,1,i,2)-xtbk(j,i,1,1))**2+
     +          (q(j,1,i,3)-xtbk(j,i,2,1))**2+
     +          (q(j,1,i,4)-xtbk(j,i,3,1))**2
          else
          xm2 = q(j,1,i,2)**2+q(j,1,i,3)**2+q(j,1,i,4)**2
          end if
          xm2 = xm2/c2
          c2  = c2*(1.+0.5*gm1*xm2)
        end if
      endif
c     add directional components of surface jet
c
      uub = uub + sjetx*xmach*cqt*c2/(gamma*pb)
      vvb = vvb + sjety*xmach*cqt*c2/(gamma*pb)
      wwb = wwb + sjetz*xmach*cqt*c2/(gamma*pb)
c
      qk0(j,i,1,1) = gamma*pb/c2
      qk0(j,i,2,1) = uub
      qk0(j,i,3,1) = vvb
      qk0(j,i,4,1) = wwb
      qk0(j,i,5,1) = pb
      bck(j,i,1)   = 1.0
c
c     f23 = 0.0  -  2-point extrapolation
c           1.0  -  3-point extrapolation
c
      f23 = 0.0
c
      k2 = min(2,kdim1)
      if (k2.eq.1) f23 = 0.0
c
      z1  =   2.0 +1.5*f23
      z2  =       -0.5*f23   
      z3  = -(2.0 +    f23)
c
      qk0(j,i,1,2) = z1*q(j,1,i,1) + z2*q(j,k2,i,1) + z3*qk0(j,i,1,1)
      qk0(j,i,2,2) = z1*q(j,1,i,2) + z2*q(j,k2,i,2) + z3*qk0(j,i,2,1)
      qk0(j,i,3,2) = z1*q(j,1,i,3) + z2*q(j,k2,i,3) + z3*qk0(j,i,3,1)
      qk0(j,i,4,2) = z1*q(j,1,i,4) + z2*q(j,k2,i,4) + z3*qk0(j,i,4,1)
      qk0(j,i,5,2) = z1*q(j,1,i,5) + z2*q(j,k2,i,5) + z3*qk0(j,i,5,1)
c 
c     Turbulence
c     only need to do advanced model turbulence B.C.s on finest grid
c
      if (livmx.ge.2.and.level.ge.lglobal.and.ntime.ne.0) then 
        if (abs(cqt).lt.1e-12) then
          vk0(j,i,1,1) = 0.
          vk0(j,i,1,2) = 0.
          if (ivmx.eq.5) then 
            tk0(j,i,1,1) = -tursav(j,1,i,1)
            tk0(j,i,2,1) = -tursav(j,1,i,2)
            tk0(j,i,1,2) = 2.*tk0(j,i,1,1)-tursav(j,1,i,1)
            tk0(j,i,2,2) = 2.*tk0(j,i,2,1)-tursav(j,1,i,2)
          elseif(ivmx.eq.6.or.ivmx.eq.7) then
            c2b=cbar/tinf
            c2bp=c2b+1.0
            re=reue/xmach
            beta1=.075
            k=1
            tt=gamma*qk0(j,i,5,1)/qk0(j,i,1,1)
            fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
            dist=ccabs(smin(j,k,i))
            tk0(j,i,1,1) = 
     =        2.*(60.*fnu/(re**2*q(j,k,i,1)*beta1*dist**2))-
     +        tursav(j,k,i,1)
            tk0(j,i,2,1) = -tursav(j,k,i,2)
            tk0(j,i,1,2) = 2.*tk0(j,i,1,1)-tursav(j,1,i,1)
            tk0(j,i,2,2) = 2.*tk0(j,i,2,1)-tursav(j,1,i,2)
          endif
        else
          vk0(j,i,1,1) = vist3d(j,1,i)
          vk0(j,i,1,2) = 0.
          do l=1,nummem
            tk0(k,i,l,1) = ccmax(tur10(l), turwj(l))
            tk0(k,i,l,2) = tk0(k,i,l,1)
          enddo
        endif
      endif
 1100 continue
 1200 continue
c
      end if
c
c******************************************************************************
c      k=kdim boundary     viscous wall with T & cq specified       bctype 2036
c******************************************************************************
      if (nface.eq.6) then
c 
      do 1600 i=ista,iend1
      ii = i-ista+1
c
      do 1500 j=jsta,jend1
      jj = j-jsta+1
c     read data
      cq0           = bcdata(jj,ii,ip,2)
      cqa           = bcdata(jj,ii,ip,3)
      sjetx         = bcdata(jj,ii,ip,4)
      sjety         = bcdata(jj,ii,ip,5)
      sjetz         = bcdata(jj,ii,ip,6)
      t_cyc         = bcdata(jj,ii,ip,7)
      t_s99         = bcdata(jj,ii,ip,8)
      d_cyc         = bcdata(jj,ii,ip,9)
      
      call bc2036_check(t_cyc,t_s99,d_cyc,nou,bou,ibufdim,nbuf)

c
c     surface velocities
c
      uub = 0.
      vvb = 0.
      wwb = 0.
c
c     for dynamic mesh, set velocity at wall to grid velocity at wall
c     if irelv > 0; otherwise, set to zero
c
      if (iuns.gt.0 .and. irelv.gt.0) then
      uub = xtbk(j,i,1,2)
      vvb = xtbk(j,i,2,2)
      wwb = xtbk(j,i,3,2)
      end if

c
c     compute total (steady and unsteady) surface jet
c
      cqt = cq2036(time, cq0, cqa, t_cyc, t_s99, d_cyc)

c     compute direction cosines for surface jet
c
      sjet   = sqrt(sjetx*sjetx+sjety*sjety+sjetz*sjetz)
      snx    = -sk(j,kdim,i,1)
      sny    = -sk(j,kdim,i,2)
      snz    = -sk(j,kdim,i,3)
      if (sjet.ne.0.0) then
        sjetx = sjetx/sjet
        sjety = sjety/sjet
        sjetz = sjetz/sjet
      else
        sjetx = snx
        sjety = sny
        sjetz = snz
      endif
      
c     wall normal component
c     > 0 means inflow
      cqn = cqt*(sjetx*snx+sjety*sny+sjetz*snz)
      pb = q(j,kdim1,i,5)
      c2 = gamma*q(j,kdim1,i,5)/q(j,kdim1,i,1)
      
      if (real(cqn).lt.1e-12) then ! without jet
        dpb           = q(j,kdim1,i,5)-q(j,kdim1-1,i,5)
        pb            = pb + dpb/2.0
        if (real(pb).le.0.0)  pb = q(j,kdim1,i,5)
        if (real(bcdata(jj,ii,ip,1)) .gt. 0.) then
          c2  = bcdata(jj,ii,ip,1)
        else if (real(bcdata(jj,ii,ip,1)) .lt. 0.) then
          c2=1.e0+gm1*0.5e0*xmach*xmach
        else
          if (iuns.gt.0 .and. irelv.gt.0) then
          xm2 = (q(j,kdim1,i,2)-xtbk(j,i,1,2))**2+
     +          (q(j,kdim1,i,3)-xtbk(j,i,2,2))**2+
     +          (q(j,kdim1,i,4)-xtbk(j,i,3,2))**2
          else
          xm2 = q(j,kdim1,i,2)**2+q(j,kdim1,i,3)**2+q(j,kdim1,i,4)**2
          end if
          xm2 = xm2/c2
          c2  = c2*(1.+0.5*gm1*xm2)
        end if
      endif
c     add directional components of surface jet
c
      uub = uub + sjetx*xmach*cqt*c2/(gamma*pb)
      vvb = vvb + sjety*xmach*cqt*c2/(gamma*pb)
      wwb = wwb + sjetz*xmach*cqt*c2/(gamma*pb)
c
      qk0(j,i,1,3) = gamma*pb/c2
      qk0(j,i,2,3) = uub
      qk0(j,i,3,3) = vvb
      qk0(j,i,4,3) = wwb
      qk0(j,i,5,3) = pb
      bck(j,i,2)   = 1.0
c
c     f23 = 0.0  -  2-point extrapolation
c           1.0  -  3-point extrapolation
c
      f23 = 0.0
c
      k2 = max(1,kdim-2)
      if (k2.eq.1) f23 = 0.0
c
      z1  =  -2.0 -1.5*f23
      z2  =       +0.5*f23   
      z3  = +(2.0 +    f23)
c
      qk0(j,i,1,4) = z1*q(j,kdim1,i,1)+z2*q(j,k2,i,1)+z3*qk0(j,i,1,3)
      qk0(j,i,2,4) = z1*q(j,kdim1,i,2)+z2*q(j,k2,i,2)+z3*qk0(j,i,2,3)
      qk0(j,i,3,4) = z1*q(j,kdim1,i,3)+z2*q(j,k2,i,3)+z3*qk0(j,i,3,3)
      qk0(j,i,4,4) = z1*q(j,kdim1,i,4)+z2*q(j,k2,i,4)+z3*qk0(j,i,4,3)
      qk0(j,i,5,4) = z1*q(j,kdim1,i,5)+z2*q(j,k2,i,5)+z3*qk0(j,i,5,3)
c 
c     Turbulence
c     only need to do advanced model turbulence B.C.s on finest grid
c
      if (ivmx.ge.2.and.level.ge.lglobal.and.ntime.ne.0) then 
        if (abs(cqt).lt.1e-12) then
          vk0(j,i,1,3) = 0.
          vk0(j,i,1,4) = 0.
          if (ivmx.eq.5) then 
            tk0(j,i,1,3) = -tursav(j,kdim-1,i,1)
            tk0(j,i,2,3) = -tursav(j,kdim-1,i,2)
            tk0(j,i,1,4) = 2.*tk0(j,i,1,3)-tursav(j,kdim-1,i,1)
            tk0(j,i,2,4) = 2.*tk0(j,i,2,3)-tursav(j,kdim-1,i,2)
          elseif(ivmx.eq.6.or.ivmx.eq.7) then 
            c2b=cbar/tinf
            c2bp=c2b+1.0
            re=reue/xmach
            beta1=.075
            k=kdim-1
            tt=gamma*qk0(j,i,5,3)/qk0(j,i,1,3)
            fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
            dist=ccabs(smin(j,k,i))
            tk0(j,i,1,3) = 
     =        2.*(60.*fnu/(re**2*q(j,k,i,1)*beta1*dist**2))-
     +        tursav(j,k,i,1)
            tk0(j,i,2,3) = -tursav(j,k,i,2)
            tk0(j,i,1,4) = 2.*tk0(j,i,1,3)-tursav(j,kdim-1,i,1)
            tk0(j,i,2,4) = 2.*tk0(j,i,2,3)-tursav(j,kdim-1,i,2)
          endif
        else
          vk0(j,i,1,3) = vist3d(j,kdim-1,i)
          vk0(j,i,1,4) = 0.0
          do l=1,nummem
            tk0(j,i,l,3) = ccmax(tur10(l), turwj(l))
            tk0(j,i,l,4) = tk0(j,i,l,3)
          enddo
        endif
      endif
 1500 continue
 1600 continue
c
      end if
c
c******************************************************************************
c      i=1 boundary        viscous wall with T & cq specified       bctype 2036
c******************************************************************************
      if (nface.eq.1) then
c 
      do 2000 k=ksta,kend1
      kk = k-ksta+1
c
      do 1900 j=jsta,jend1
      jj = j-jsta+1
c     read data
      cq0           = bcdata(jj,kk,ip,2)
      cqa           = bcdata(jj,kk,ip,3)
      sjetx         = bcdata(jj,kk,ip,4)
      sjety         = bcdata(jj,kk,ip,5)
      sjetz         = bcdata(jj,kk,ip,6)
      t_cyc         = bcdata(jj,kk,ip,7)
      t_s99         = bcdata(jj,kk,ip,8)
      d_cyc         = bcdata(jj,kk,ip,9)
      
      call bc2036_check(t_cyc,t_s99,d_cyc,nou,bou,ibufdim,nbuf)

c
c     surface velocities
c
      uub = 0.
      vvb = 0.
      wwb = 0.
c
c     for dynamic mesh, set velocity at wall to grid velocity at wall
c     if irelv > 0; otherwise, set to zero
c
      if (iuns.gt.0 .and. irelv.gt.0) then
      uub = xtbi(j,k,1,1)
      vvb = xtbi(j,k,2,1)
      wwb = xtbi(j,k,3,1)
      end if

c
c     compute total (steady and unsteady) surface jet
c
      cqt = cq2036(time, cq0, cqa, t_cyc, t_s99, d_cyc)

c     compute direction cosines for surface jet
c
      sjet   = sqrt(sjetx*sjetx+sjety*sjety+sjetz*sjetz)
      snx    = si(j,k,1,1)
      sny    = si(j,k,1,2)
      snz    = si(j,k,1,3)
      if (sjet.ne.0.0) then
        sjetx = sjetx/sjet
        sjety = sjety/sjet
        sjetz = sjetz/sjet
      else
        sjetx = snx
        sjety = sny
        sjetz = snz
      endif
c     wall normal component
c     > 0 means inflow
      cqn = cqt*(sjetx*snx+sjety*sny+sjetz*snz)
      pb            = q(j,k,1,5)
      c2            = gamma*q(j,k,1,5)/q(j,k,1,1)
      
      if (real(cqn).lt.1e-12) then ! without jet
        dpb           = q(j,k,2,5)-q(j,k,1,5)
        pb            = pb - dpb/2.0
        if (real(pb).le.0.0)  pb = q(j,k,1,5)
        if (real(bcdata(jj,kk,ip,1)) .gt. 0.) then
          c2  = bcdata(jj,kk,ip,1)
        else if (real(bcdata(jj,kk,ip,1)) .lt. 0.) then
          c2=1.e0+gm1*0.5e0*xmach*xmach
        else
          if (iuns.gt.0 .and. irelv.gt.0) then
          xm2 = (q(j,k,1,2)-xtbi(j,k,1,1))**2+
     +          (q(j,k,1,3)-xtbi(j,k,2,1))**2+
     +          (q(j,k,1,4)-xtbi(j,k,3,1))**2
          else
          xm2 = q(j,k,1,2)**2+q(j,k,1,3)**2+q(j,k,1,4)**2
          end if
          xm2 = xm2/c2
          c2  = c2*(1.+0.5*gm1*xm2)
        end if
      end if

c     add directional components of surface jet
c
      uub = uub + sjetx*xmach*cqt*c2/(gamma*pb)
      vvb = vvb + sjety*xmach*cqt*c2/(gamma*pb)
      wwb = wwb + sjetz*xmach*cqt*c2/(gamma*pb)
c
      qi0(j,k,1,1) = gamma*pb/c2
      qi0(j,k,2,1) = uub
      qi0(j,k,3,1) = vvb
      qi0(j,k,4,1) = wwb
      qi0(j,k,5,1) = pb
      bci(j,k,1)   = 1.0
c
c     f23 = 0.0  -  2-point extrapolation
c           1.0  -  3-point extrapolation
c
      f23 = 0.0
c
      i2 = min(2,idim1)
      if (i2.eq.1) f23 = 0.0
c
      z1 =   2.0 +1.5*f23
      z2 =       -0.5*f23   
      z3 = -(2.0 +    f23)
c
      qi0(j,k,1,2) = z1*q(j,k,1,1) + z2*q(j,k,i2,1) + z3*qi0(j,k,1,1)
      qi0(j,k,2,2) = z1*q(j,k,1,2) + z2*q(j,k,i2,2) + z3*qi0(j,k,2,1)
      qi0(j,k,3,2) = z1*q(j,k,1,3) + z2*q(j,k,i2,3) + z3*qi0(j,k,3,1)
      qi0(j,k,4,2) = z1*q(j,k,1,4) + z2*q(j,k,i2,4) + z3*qi0(j,k,4,1)
      qi0(j,k,5,2) = z1*q(j,k,1,5) + z2*q(j,k,i2,5) + z3*qi0(j,k,5,1)
c 
c     Turbulence
c     only need to do advanced model turbulence B.C.s on finest grid
c
      if (ivmx.ge.2.and.level.ge.lglobal.and.ntime.ne.0) then 
        if (abs(cqt).lt.1e-12) then
          vi0(j,k,1,1) = 0.
          vi0(j,k,1,2) = 0.
          if (ivmx.eq.5) then 
            ti0(j,k,1,1) = -tursav(j,k,1,1)
            ti0(j,k,2,1) = -tursav(j,k,1,2)
            ti0(j,k,1,2) = 2.*ti0(j,k,1,1)-tursav(j,k,1,1)
            ti0(j,k,2,2) = 2.*ti0(j,k,2,1)-tursav(j,k,1,2)
          elseif(ivmx.eq.6.or.ivmx.eq.7) then 
            c2b=cbar/tinf
            c2bp=c2b+1.0
            re=reue/xmach
            beta1=.075
            i=1
            tt=gamma*qi0(j,k,5,1)/qi0(j,k,1,1)
            fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
            dist=ccabs(smin(j,k,i))
            ti0(j,k,1,1) = 
     =        2.*(60.*fnu/(re**2*q(j,k,i,1)*beta1*dist**2))-
     +        tursav(j,k,i,1)
            ti0(j,k,2,1) = -tursav(j,k,i,2)
            ti0(j,k,1,2) = 2.*ti0(j,k,1,1)-tursav(j,k,1,1)
            ti0(j,k,2,2) = 2.*ti0(j,k,2,1)-tursav(j,k,1,2)
          endif
        else
          vi0(j,k,1,1) = vist3d(j,k,1)
          vi0(j,k,1,2) = 0.
          do l=1,nummem
            ti0(j,k,l,1) = ccmax(tur10(l), turwj(l))
            ti0(j,k,l,2) = ti0(j,k,l,1)
          enddo
        endif
      endif
 1900 continue
 2000 continue
c
      end if
c
c******************************************************************************
c      i=idim boundary     viscous wall with T & cq specified       bctype 2036
c******************************************************************************
      if (nface.eq.2) then
c
      do 2400 k=ksta,kend1
      kk = k-ksta+1
c
      do 2300 j=jsta,jend1
      jj = j-jsta+1
c     read data
      cq0           = bcdata(jj,kk,ip,2)
      cqa           = bcdata(jj,kk,ip,3)
      sjetx         = bcdata(jj,kk,ip,4)
      sjety         = bcdata(jj,kk,ip,5)
      sjetz         = bcdata(jj,kk,ip,6)
      t_cyc         = bcdata(jj,kk,ip,7)
      t_s99         = bcdata(jj,kk,ip,8)
      d_cyc         = bcdata(jj,kk,ip,9)
      
      call bc2036_check(t_cyc,t_s99,d_cyc,nou,bou,ibufdim,nbuf)

c
c     surface velocities
c
      uub = 0.
      vvb = 0.
      wwb = 0.
c
c     for dynamic mesh, set velocity at wall to grid velocity at wall
c     if irelv > 0; otherwise, set to zero
c
      if (iuns.gt.0 .and. irelv.gt.0) then
        uub = xtbi(j,k,1,2)
        vvb = xtbi(j,k,2,2)
        wwb = xtbi(j,k,3,2)
      end if

c
c     compute total (steady and unsteady) surface jet
c
      cqt = cq2036(time, cq0, cqa, t_cyc, t_s99, d_cyc)

c     compute direction cosines for surface jet
c
      sjet = sqrt(sjetx*sjetx+sjety*sjety+sjetz*sjetz)
      snx  = -si(j,k,idim,1)
      sny  = -si(j,k,idim,2)
      snz  = -si(j,k,idim,3)
      if (sjet.ne.0.0) then
        sjetx = sjetx/sjet
        sjety = sjety/sjet
        sjetz = sjetz/sjet
      else
        sjetx = snx
        sjety = sny
        sjetz = snz
      endif
      cqn = cqt*(sjetx*snx+sjety*sny+sjetz*snz)
      pb  = q(j,k,idim1,5)
      c2  = gamma*q(j,k,idim1,5)/q(j,k,idim1,1)
      if (real(cqn).lt.1e-12) then ! without jet
        dpb           = q(j,k,idim1,5)-q(j,k,idim1-1,5)
        pb            = pb + dpb/2.0
        if (real(pb).le.0.0)  pb = q(j,k,idim1,5)
        if (real(bcdata(jj,kk,ip,1)) .gt. 0.) then
          c2  = bcdata(jj,kk,ip,1)
        else if (real(bcdata(jj,kk,ip,1)) .lt. 0.) then
          c2=1.e0+gm1*0.5e0*xmach*xmach
        else
          if (iuns.gt.0 .and. irelv.gt.0) then
          xm2 = (q(j,k,idim1,2)-xtbi(j,k,1,2))**2+
     +          (q(j,k,idim1,3)-xtbi(j,k,2,2))**2+
     +          (q(j,k,idim1,4)-xtbi(j,k,3,2))**2
          else
          xm2 = q(j,k,idim1,2)**2+q(j,k,idim1,3)**2+q(j,k,idim1,4)**2
          end if
          xm2 = xm2/c2
          c2  = c2*(1.+0.5*gm1*xm2)
        end if
      end if
c     add directional components of surface jet
c
      uub = uub + sjetx*xmach*cqt*c2/(gamma*pb)
      vvb = vvb + sjety*xmach*cqt*c2/(gamma*pb)
      wwb = wwb + sjetz*xmach*cqt*c2/(gamma*pb)
c
      qi0(j,k,1,3) = gamma*pb/c2
      qi0(j,k,2,3) = uub
      qi0(j,k,3,3) = vvb
      qi0(j,k,4,3) = wwb
      qi0(j,k,5,3) = pb
      bci(j,k,2)   = 1.0
c
c     f23 = 0.0  -  2-point extrapolation
c           1.0  -  3-point extrapolation
c
      f23 = 0.0
c
      i2 = max(1,idim-2)
      if (i2.eq.1) f23 = 0.0
c
      z1 =  -2.0 -1.5*f23
      z2 =       +0.5*f23   
      z3 = +(2.0 +    f23)
c
      qi0(j,k,1,4) = z1*q(j,k,idim1,1)+z2*q(j,k,i2,1)+z3*qi0(j,k,1,3)
      qi0(j,k,2,4) = z1*q(j,k,idim1,2)+z2*q(j,k,i2,2)+z3*qi0(j,k,2,3)
      qi0(j,k,3,4) = z1*q(j,k,idim1,3)+z2*q(j,k,i2,3)+z3*qi0(j,k,3,3)
      qi0(j,k,4,4) = z1*q(j,k,idim1,4)+z2*q(j,k,i2,4)+z3*qi0(j,k,4,3)
      qi0(j,k,5,4) = z1*q(j,k,idim1,5)+z2*q(j,k,i2,5)+z3*qi0(j,k,5,3)
c 
c     Turbulence
c     only need to do advanced model turbulence B.C.s on finest grid
c
      if (ivmx.ge.2.and.level.ge.lglobal.and.ntime.ne.0) then 
        if (abs(cqt).lt.1e-12) then
          vi0(j,k,1,3) = 0.
          vi0(j,k,1,4) = 0.
          if (ivmx.eq.5) then 
            ti0(j,k,1,3) = -tursav(j,k,idim-1,1)
            ti0(j,k,2,3) = -tursav(j,k,idim-1,2)
            ti0(j,k,1,4) = 2.*ti0(j,k,1,3)-tursav(j,k,idim-1,1)
            ti0(j,k,2,4) = 2.*ti0(j,k,2,3)-tursav(j,k,idim-1,2)
          elseif(ivmx.eq.6.or.ivmx.eq.7) then
            c2b=cbar/tinf
            c2bp=c2b+1.0
            re=reue/xmach
            beta1=.075
            i=idim-1
            tt=gamma*qi0(j,k,5,3)/qi0(j,k,1,3)
            fnu=c2bp*tt*sqrt(tt)/(c2b+tt)
            dist=ccabs(smin(j,k,i))
            ti0(j,k,1,3) = 
     =        2.*(60.*fnu/(re**2*q(j,k,i,1)*beta1*dist**2))-
     +        tursav(j,k,i,1)
            ti0(j,k,2,3) = -tursav(j,k,i,2)
            ti0(j,k,1,4) = 2.*ti0(j,k,1,3)-tursav(j,k,idim-1,1)
            ti0(j,k,2,4) = 2.*ti0(j,k,2,3)-tursav(j,k,idim-1,2)
          endif
        else
          vi0(j,k,1,3) = vist3d(j,k,idim-1)
          vi0(j,k,1,4) = 0.
          do l=1,nummem
            ti0(j,k,l,3) = ccmax(tur10(l), turwj(l))
            ti0(j,k,l,4) = ti0(j,k,l,3)
          enddo
        endif
      endif
 2300 continue
 2400 continue
c
      end if
c
      return
      end

c     d_cyc > 0 pulsed jet
c     d_cyc = 0 steady jet
c     d_cyc < 0 one-minus-cosine
c  a) pulsed jet (d_cyc>0):
c       d_cyc = duty cycle of the pulse between (0, 1)
c       t_cyc = period of unsteady jet
c       t_s99 = time of the pulsed jet achieving its maximum value
c         t_s99 = 0 => step function
c         t_s99 > 0 => tanh function
c         t_s99 < 0 => linear function
c  b) one-shot jet (d_cyc = 0):
c       t_cyc = jet start time
c       t_s99 = time of the pulsed jet achieving its maximum value
c         t_s99 = 0 => step function
c         t_s99 > 0 => tanh function
c         t_s99 < 0 => linear function
c  c) one-minus-cosine jet (d_cyc < 0):
c       d_cyc = -1 => continuous
c       d_cyc = -2 => single wave
c       t_s99 = time of the jet achieving its maximum value (half period)
c       t_cyc = jet start time
      function cq2036(time, cq0, cqu, t_cyc, t_s99, d_cyc) result(cqt)
#ifdef CMPLX
      implicit complex(a-h,o-z)
#endif
        parameter pi   =3.141592653589793
        parameter twopi=6.283185307179586
        cqt = cq0
        if (d_cyc.eq.0) then  ! one shot
          t = time-t_cyc
          if (t.lt.0) then 
            cqt = cq0
          elseif(t.ge.0.and.t.lt.abs(t_s99)) then 
            c = t/abs(t_s99)
            if (t_s99.gt.0) c = 0.5+0.5*tanh((c-0.5)*twopi)
            cqt = cq0 + c*cqu
          else 
            cqt = cq0 + cqu
          endif
        elseif(d_cyc.lt.0) then  ! one-minus-cosine
          t = time-t_cyc
          if (t.le.0) then 
            cqt = cq0
          elseif (int(d_cyc).eq.-1.or.t.le.2*t_s99) then
            cqt = cq0 + 0.5*cqu*(1.0-cos(pi*t/t_s99))
          else
            cqt = cq0
          endif
        else ! pulsed jet
          t = mod(time, t_cyc)
          t00 = 0.5*((1.0-d_cyc)*t_cyc-abs(t_s99)) 
          t01 = t00 + abs(t_s99)
          t11 = t_cyc - t00
          t10 = t11 - abs(t_s99)
          if (t.le.t00) then 
            cqt=cq0
          elseif (t.ge.t00.and.t.lt.t01) then 
            dt = t-t00
            c = dt/abs(t_s99)
            if (t_s99.gt.0) c = 0.5+0.5*tanh((c-0.5)*twopi)
            cqt = cq0 + c*cqu
          elseif (t.ge.t01.and.t.lt.t10) then 
            cqt = cq0+cqu
          elseif (t.ge.t10.and.t.lt.t11) then 
            dt = t11-t
            c = dt/abs(t_s99)
            if (t_s99.gt.0) c = 0.5+0.5*tanh((c-0.5)*twopi)
            cqt = cq0 + c*cqu
          else
            cqt = cq0
          endif
        endif
      end function cq2036

      subroutine bc2036_check(t_cyc,t_s99,d_cyc,nou,bou,ibufdim,nbuf)
#ifdef CMPLX
      implicit complex(a-h,o-z)
#endif
        character*120 bou(ibufdim,nbuf)
        dimension nou(nbuf)
        if (d_cyc .lt. 0) then ! one-minus-cosine
          if (t_s99.lt.1e-20) then 
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' stopping in bc2036...',
     .        ' one-minus-cosine jet, invalid half period t_s99=',
     .        t_s99
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
          endif
          if (t_cyc.lt.0) then 
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' stopping in bc2036...',
     .        ' one-minus-cosine jet, invalid start time t_cyc=',
     .        t_cyc
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
          endif
          idc = int(d_cyc)
          if (idc.ne.-1 .and. idc.ne.-2) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' stopping in bc2036...',
     .        ' one-minus-cosine jet, invalid d_cyc=', d_cyc
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
          endif
        elseif (d_cyc.eq.0) then ! one-shot
          if (t_cyc.lt.0) then 
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' stopping in bc2036...',
     .        ' one-shot jet, invalid start time t_cyc=',
     .        t_cyc
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
          endif
        else ! pulsed
          if (t_cyc.lt.0) then 
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)' stopping in bc2036...',
     .        ' pulsed jet, invalid period t_cyc=',
     .        t_cyc
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
          endif
          ! check duty cycle
          if (d_cyc.ge.1.0) then 
            write(bou(nou(1),1),*)' stopping in bc2036...',
     .        ' pulsed jet, invalid duty cycle d_cyc=',
     .        d_cyc
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
          endif
          ts99max = ccmin(d_cyc, 1.0-d_cyc)*t_cyc
          if (abs(t_s99).gt.ts99max) then 
            write(bou(nou(1),1),'(A,2(A,E10.3),A)')
     .      ' stopping in bc2036...',
     .      ' pulsed jet, invalid t_s99=', t_s99, 
     .      ' must not be greater than ', ts99max
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
          endif
        endif
      end subroutine bc2036_check
