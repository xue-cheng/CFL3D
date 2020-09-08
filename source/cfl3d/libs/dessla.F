      subroutine dessla(delta,j,k,i,jdim,kdim,idim,iex,iex2,iex3,
     *                  x,y,z,vol,sj,sk,si,vor,ux,turre,fnu)
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
      common /twod/ i2d
      common /zero/ iexp
      common /des/ cdes,ides,cddes,ides_sla
      common /slasgs/ fkh_min, fkh_max, fkh_a1, fkh_a2
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension sj(jdim,kdim,idim-1,5)
      dimension sk(jdim,kdim,idim-1,5)
      dimension si(jdim,kdim,idim,5)
      dimension vol(jdim,kdim,idim-1),vor(jdim-1,kdim-1,idim-1)
      dimension ux(jdim-1,kdim-1,idim-1,9)
      dimension vvor(3), uvor(3), rn(3,8), cg(3), rl(3,8), sdw(3)
      dimension fnu(0:jdim,0:kdim,0-iex3:idim+iex3)
      dimension turre(0-iex:jdim+iex,0-iex:kdim+iex,0-iex2:idim+iex2)
      parameter osqrt3 = 0.5773502691896258
      parameter sqrt6  = 2.449489742783178
      xminn=10.**(-iexp+1)
      if (i2d.ne.0.and.ides_sla.gt.0.and.vor(j,k,i).gt.xminn) then
        vvor(1) = ux(j,k,i,8) - ux(j,k,i,6)
        vvor(2) = ux(j,k,i,3) - ux(j,k,i,7)
        vvor(3) = ux(j,k,i,4) - ux(j,k,i,2)
        uvor = vvor / vor(j,k,i) ! normalized
        rn(1,1) = x(j  ,k  ,i  )
        rn(1,2) = x(j+1,k  ,i  )
        rn(1,3) = x(j  ,k+1,i  )
        rn(1,4) = x(j+1,k+1,i  )
        rn(1,5) = x(j  ,k  ,i+1)
        rn(1,6) = x(j+1,k  ,i+1)
        rn(1,7) = x(j  ,k+1,i+1)
        rn(1,8) = x(j+1,k+1,i+1)
        rn(2,1) = y(j  ,k  ,i  )
        rn(2,2) = y(j+1,k  ,i  )
        rn(2,3) = y(j  ,k+1,i  )
        rn(2,4) = y(j+1,k+1,i  )
        rn(2,5) = y(j  ,k  ,i+1)
        rn(2,6) = y(j+1,k  ,i+1)
        rn(2,7) = y(j  ,k+1,i+1)
        rn(2,8) = y(j+1,k+1,i+1)
        rn(3,1) = z(j  ,k  ,i  )
        rn(3,2) = z(j+1,k  ,i  )
        rn(3,3) = z(j  ,k+1,i  )
        rn(3,4) = z(j+1,k+1,i  )
        rn(3,5) = z(j  ,k  ,i+1)
        rn(3,6) = z(j+1,k  ,i+1)
        rn(3,7) = z(j  ,k+1,i+1)
        rn(3,8) = z(j+1,k+1,i+1)

        cg = 0.
        do ir =1,8
          cg = cg + rn(:,ir)
        enddo
        cg = cg/8
        do ir =1,8
          rn(:,ir) = rn(:,ir) - cg
          call cross_prod(uvor, rn(:,ir), rl(:,ir))
        enddo
        call max_cross_distance(dmax, rl, 8)
        delta = osqrt3*dmax
        if (ides_sla.eq.2) then
c  Shear Layer Adapt                  
                  ! VTM
c               shat:
c                      ux      (vx+uy)/2   (wx+uz)/2
c                  (vx+uy)/2       vy      (wy+vz)/2
c                  (wx+uz)/2   (wy+vz)/2       wz
          sdw(1) =  ux(j,k,i,1)                   *vvor(1) + 
     +              0.5*(ux(j,k,i,2)+ux(j,k,i,4)) *vvor(2) + 
     +              0.5*(ux(j,k,i,3)+ux(j,k,i,7)) *vvor(3)
          sdw(2) =  0.5*(ux(j,k,i,2)+ux(j,k,i,4)) *vvor(1) +
     +              ux(j,k,i,5)                   *vvor(2) +
     +              0.5*(ux(j,k,i,6)+ux(j,k,i,8)) *vvor(3)
          sdw(3) =  0.5*(ux(j,k,i,3)+ux(j,k,i,7)) *vvor(1) +
     +              0.5*(ux(j,k,i,6)+ux(j,k,i,8)) *vvor(2) +
     +              ux(j,k,i,9)                   *vvor(3)
          call cross_prod(sdw, vvor, uvor)
          VTM = sqrt6*norm2(uvor)
          trs2=
     =      ux(j,k,i,1)**2+ux(j,k,i,5)**2+ux(j,k,i,9)**2 +
     +      0.5*(ux(j,k,i,2)+ux(j,k,i,4))**2 +
     +      0.5*(ux(j,k,i,3)+ux(j,k,i,7))**2 + 
     +      0.5*(ux(j,k,i,6)+ux(j,k,i,8))**2
          trs = ux(j,k,i,1)+ux(j,k,i,5)+ux(j,k,i,9)
          VTM = VTM/(dot3(vvor,vvor)*sqrt(3*trs2-trs**2)) 
          VTM = VTM*ccmaxcr(0.2*fnu(j,k,i)/turre(j,k,i),1.0)
          fkh = fkh_min+(fkh_max-fkh_min)/(fkh_a2-fkh_a1)*(VTM-fkh_a1)
          fkh = ccmax(fkh_min,ccmax(fkh_max,fkh))
          delta = delta*fkh

        endif
      endif

      end subroutine dessla
      
      subroutine cross_prod(v1,v2,v3)
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
        dimension v1(3), v2(3), v3(3)
        
        v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
        v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
        v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
      end subroutine cross_prod

      function dot3(v1,v2) result(v3)
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
        dimension v1(3), v2(3)
        v3 = v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
      end function dot3

      subroutine max_cross_distance(dmax, rl, n)
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
        dimension rl(3,n)
        
        dm2 = 0.0
        do im = 1,n-1
          do in = im+1,n
            dx = rl(1,im)-rl(1,in)
            dy = rl(2,im)-rl(2,in)
            dz = rl(3,im)-rl(3,in)
            d2 = dx*dx+dy*dy+dz*dz
            if (d2.gt.dm2) dm2 = d2
          enddo
        enddo
        dmax = sqrt(dm2)
      end subroutine max_cross_distance
