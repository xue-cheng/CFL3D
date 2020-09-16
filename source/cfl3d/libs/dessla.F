      subroutine dessla(delta,rn,vor,fnu,turre,
     .                  ux,uy,uz,vx,vy,vz,wx,wy,wz)
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c INPUTS
        dimension rn(3,8)
        dimension vvor(3) ! vorticity vector
        dimension uvor(3) ! normalized `vvor`, later used in VTM calc
c
        dimension rl(3,8), sdw(3)
c
        common /twod/ i2d
        common /zero/ iexp
        common /des/ cdes,ides,cddes,ides_sla
        common /slasgs/ fkh_min, fkh_max, fkh_a1, fkh_a2
        parameter osqrt3 = 0.5773502691896258
        parameter sqrt6  = 2.449489742783178
        xminn=10.**(-iexp+1)
c     only in cells which have non-zero vorticity
        if (i2d.ne.0.and.ides_sla.gt.0.and.vor.gt.xminn) then
          vvor(1) = wy - vz
          vvor(2) = uz - wx
          vvor(3) = vx - uy
          uvor = vvor / vor ! normalized
          do ir =1,8
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
            shat11 = ux
            shat12 = 0.5*(uy+vx)
            shat13 = 0.5*(wx+uz)
            shat22 = vy
            shat23 = 0.5*(wy+vz)
            shat33 = wz
            sdw(1) = shat11*vvor(1)+shat12*vvor(2)+sha13*vvor(3)
            sdw(2) = shat12*vvor(1)+shat22*vvor(2)+sha23*vvor(3)
            sdw(3) = shat13*vvor(1)+shat23*vvor(2)+sha33*vvor(3)
            call cross_prod(sdw, vvor, uvor)
            VTM = sqrt6*norm2(uvor)
            trs2=   shat11**2 +   shat22**2 +   shat33**2 +
     +            2*shat12**2 + 2*shat13**2 + 2*shat23**2
            trs = shat11+shat22+shat33
            VTM = VTM/(vor*vor*sqrt(3*trs2-trs**2)) 
            VTM = VTM*ccmaxcr(0.2*fnu/turre,1.0)
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

      subroutine max_cross_distance(dmax, rl, d, n)
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
        dimension rl(d,n)
        
        dm2 = 0.0
        do im = 1,n-1
          do in = im+1,n
            d2 = 0.0
            do id = 1,d
              d2 = d2+(rl(id,im)-rl(id,in))**2
            enddo
            if (d2.gt.dm2) dm2 = d2
          enddo
        enddo
        dmax = sqrt(dm2)
      end subroutine max_cross_distance
