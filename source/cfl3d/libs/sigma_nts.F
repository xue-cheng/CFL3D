      subroutine sigma_nts(jdim,kdim,idim,sj,sk,si,vol,
     .  turre,fnu,nou,bou,nbuf,ibufdim,iex,iex2,iex3,ux,sigma)
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      common /des/ cdes,ides,cddes,ides_sgs
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
      common /twod/ i2d
      common /zero/ iexp
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /axisym/ iaxi2plane,iaxi2planeturb,istrongturbdis,iforcev0
c
      dimension sj(jdim,kdim,idim-1,5),sk(jdim,kdim,idim-1,5),
     .          si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1)
      dimension turre(0-iex:jdim+iex,0-iex:kdim+iex,0-iex2:idim+iex2),
     .          fnu(0:jdim,0:kdim,0-iex3:idim+iex3)
      dimension ux(jdim-1,kdim-1,idim-1,9) 
      dimension sigma(jdim-1,kdim-1,idim-1)

      parameter cmu=0.09,cv1=7.1,ch1=3.0,ch2=1.0,ch3=2.0,small=1e-20
      
      re=reue/xmach
      uref  = sqrt(u0*u0+v0*v0+w0*w0)

      do i=1,idim-1
        do j=1,jdim-1
          do k=1,kdim-1
            deltaj = 2.*vol(j,k,i)/(sj(j,k,i,4)+sj(j+1,k,i,4))
            deltak = 2.*vol(j,k,i)/(sk(j,k,i,4)+sk(j,k+1,i,4))
            deltai = 2.*vol(j,k,i)/(si(j,k,i,4)+si(j,k,i+1,4))
            delta = ccmax(deltaj, deltak)
            if( i2d .ne. 1 .and. iaxi2planeturb .ne. 1 ) then
              delta = ccmax(delta,deltai)
            end if
            vort    = sqrt((ux(j,k,i,8) - ux(j,k,i,6))**2 +
     +                     (ux(j,k,i,3) - ux(j,k,i,7))**2 +
     +                     (ux(j,k,i,4) - ux(j,k,i,2))**2 ) 
            str     = sqrt((ux(j,k,i,8) + ux(j,k,i,6))**2 +
     +                     (ux(j,k,i,3) + ux(j,k,i,7))**2 +
     +                     (ux(j,k,i,4) + ux(j,k,i,2))**2 +
     +                     2.0*(ux(j,k,i,1)**2 + 
     +                          ux(j,k,i,5)**2 +
     +                          ux(j,k,i,9)**2) ) 
            velterm = 0.5*(vort**2+str**2)
c
            fK      = ccmax(sqrt(velterm),0.1*uref)
            
            chi     = turre(j,k,i)/fnu(j,k,i)
            fv1     = chi**3/(chi**3+cv1**3)   
            lturb   = sqrt((fv1*turre(j,k,i)+fnu(j,k,i))/
     *                (re*fK*(cmu**1.5)))
            fB      = ch3*vort*ccmax(str,vort)/ccmax(velterm,small)
            g       = tanh(fB**4)
            rl      = cdes*delta/ccmax(lturb,small)
            fA      = ch2*ccmaxcr(rl/ccmax(g,small)-0.5,0.0)
            sigma(j,k,i) = tanh(fA**ch1)
          enddo
        enddo
      enddo
      end subroutine sigma_nts