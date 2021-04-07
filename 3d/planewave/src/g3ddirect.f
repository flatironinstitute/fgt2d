c
c**********************************************************************
      subroutine g3d_directcp_vec(nd,delta,dmax,sources,ns,charge,
     $           targ,ntarg,pot)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT
c     at the target point TARGET, due to a vector of 
c     charges at SOURCE(3,ns). 
c
c     pot(ii)  = \sum_j exp(-r^2/delta)*charge(ii,j)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)       : potential is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd
      real *8 sources(3,ns),targ(3,ntarg),xdiff,ydiff,zdiff,rr,r
      real *8 rtmp,delta
      real *8 pot(nd,ntarg)
      real *8 charge(nd,ns)
c
      do itarg=1,ntarg
      do i = 1,ns
         xdiff=targ(1,itarg)-sources(1,i)
         ydiff=targ(2,itarg)-sources(2,i)
         zdiff=targ(3,itarg)-sources(3,i)
         rr=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            do ii=1,nd
               pot(ii,itarg) = pot(ii,itarg) + rtmp*charge(ii,i)
            enddo
         endif
      enddo
      enddo
      
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine g3d_directcg_vec(nd,delta,dmax,sources,ns,charge,
     1             targ,ntarg,pot,grad)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT and gradients GRAD
c     at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c
c     pot(ii)  = \sum_j exp(-r^2/delta)*charge(ii,j)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     targ          :   location of the target

c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)       : potential is incremented
c     grad(nd,2)      : gradient is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd
      real *8 sources(3,ns),targ(3,ntarg)
      real *8 xdiff,ydiff,zdiff,rr,r,rtmp,delta
      real *8 pot(nd,ntarg),grad(nd,3,ntarg)
      real *8 dx,dy,dz,rinc
      real *8 charge(nd,ns)
c
      rfac=-2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         xdiff=targ(1,itarg)-sources(1,i)
         ydiff=targ(2,itarg)-sources(2,i)
         zdiff=targ(3,itarg)-sources(3,i)
         rr=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            dy = rfac*ydiff
            dz = rfac*zdiff
            do ii = 1,nd
               rinc = rtmp*charge(ii,i)
               pot(ii,itarg) = pot(ii,itarg) + rinc
               grad(ii,1,itarg) = grad(ii,1,itarg) + dx*rinc
               grad(ii,2,itarg) = grad(ii,2,itarg) + dy*rinc
               grad(ii,3,itarg) = grad(ii,3,itarg) + dz*rinc
            enddo
         endif
      enddo
      enddo
      return
      end
c
c
c
c
c
c
c
c
c**********************************************************************
      subroutine g3d_directch_vec(nd,delta,dmax,sources,ns,charge,
     1           targ,ntarg,pot,grad,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT, gradients GRAD
c     and Hessians at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     
c     pot(ii)  = \sum_j exp(-r^2/delta)*charge(ii,j)
c
c     grad = grad
c     hess = (dxx,dxy,dyy)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)       : potential is incremented
c     grad(nd,2)    : gradient is incremented
c     hess(nd,3)    : Hessian is incremented
c                     using ordering (dxx,dyy,dzz,dxy,dxz,dyz)
      
c---------------------------------------------------------------------
      integer i,ns,ifexpon,ii,nd
      real *8 sources(3,ns),targ(3,ntarg)
      real *8 xdiff,ydiff,zdiff,rr,r,delta
      real *8 pot(nd,ntarg),grad(nd,3,ntarg),hess(nd,6,ntarg)
      real *8 rtmp,dx,dy,dz,rinc,dxx,dyy,dzz
      real *8 dxy,dxz,dyz
      real *8 charge(nd,ns)
c
      rfac = -2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         xdiff=targ(1,itarg)-sources(1,i)
         ydiff=targ(2,itarg)-sources(2,i)
         zdiff=targ(3,itarg)-sources(3,i)
         rr=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            dy = rfac*ydiff
            dz = rfac*zdiff
            dxx = rfac + dx*dx
            dyy = rfac + dy*dy
            dzz = rfac + dz*dz
            dxy = dx*dy
            dxz = dx*dz
            dyz = dy*dz
c
            do ii = 1,nd
               rinc = rtmp*charge(ii,i)
               pot(ii,itarg) = pot(ii,itarg) + rinc
            
               grad(ii,1,itarg) = grad(ii,1,itarg) + dx*rinc
               grad(ii,2,itarg) = grad(ii,2,itarg) + dy*rinc
               grad(ii,3,itarg) = grad(ii,3,itarg) + dz*rinc

               hess(ii,1,itarg) = hess(ii,1,itarg) + dxx*rinc
               hess(ii,2,itarg) = hess(ii,2,itarg) + dyy*rinc
               hess(ii,3,itarg) = hess(ii,3,itarg) + dzz*rinc
               hess(ii,4,itarg) = hess(ii,4,itarg) + dxy*rinc
               hess(ii,5,itarg) = hess(ii,5,itarg) + dxz*rinc
               hess(ii,6,itarg) = hess(ii,6,itarg) + dyz*rinc
            enddo
         endif
      enddo
      enddo
      
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine g3d_directdp_vec(nd,delta,dmax,sources,ns,
     $           rnormal,dipstr,targ,ntarg,pot)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT
c     at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c
c     pot(ii)  = \sum_j \grad_s(exp(-r^2/delta)*charge(ii,j))*N(s)
c     
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     dipstr(nd,ns) :   dipole strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)   (complex *16)      : potential is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd
      real *8 sources(3,ns),rnormal(3,ns),targ(3,ntarg),xdiff,ydiff,zdiff
      real *8 dx,dy,dz,rtmp,rinc,delta,rr,r
      real *8 pot(nd,ntarg)
      real *8 dipstr(nd,ns)
c
      rfac = 2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         xdiff=targ(1,itarg)-sources(1,i)
         ydiff=targ(2,itarg)-sources(2,i)
         zdiff=targ(3,itarg)-sources(3,i)
         rr=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            dy = rfac*ydiff
            dz = rfac*zdiff
            
            rinc = (dx*rnormal(1,i)+dy*rnormal(2,i)
     1            +dz*rnormal(3,i))*rtmp
            do ii=1,nd
               pot(ii,itarg) = pot(ii,itarg) + rinc*dipstr(ii,i)
            enddo
         endif
      enddo
      enddo
      
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine g3d_directdg_vec(nd,delta,dmax,sources,ns,
     1             rnormal,dipstr,targ,ntarg,pot,grad)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT and gradients GRAD
c     at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     
c     pot(ii)  = \sum_j \grad_s(exp(-r^2/delta)*charge(ii,j))*N(s)
c
c     grad  = d(pot)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     rnormal(2,ns) :   normals at sources
c     dipstr(nd,ns) :   dipole strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)       : potential is incremented
c     grad(nd,2)    : gradient is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd
      real *8 sources(3,ns),rnormal(3,ns),targ(3,ntarg)
      real *8 xdiff,ydiff,zdiff,rr,r,rtmp,rinc,dx,dy,dz
      real *8 rincx,rincy,rincz,delta
      real *8 pot(nd,ntarg),grad(nd,3,ntarg)
      real *8 dipstr(nd,ns)
c
      rfac = 2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         xdiff=targ(1,itarg)-sources(1,i)
         ydiff=targ(2,itarg)-sources(2,i)
         zdiff=targ(3,itarg)-sources(3,i)
         rr=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
         if (rr .lt. dmax) then         
            dtmp = dexp(-rr/delta)

            dx = rfac*xdiff
            dy = rfac*ydiff
            dz = rfac*zdiff
            rinc = (dx*rnormal(1,i)+dy*rnormal(2,i)+dz*rnormal(3,i))
            rincx = rnormal(1,i)*rfac - rinc*dx
            rincy = rnormal(2,i)*rfac - rinc*dy
            rincz = rnormal(3,i)*rfac - rinc*dz
c     
            do ii = 1,nd
               rtmp = dtmp*dipstr(ii,i)
               pot(ii,itarg) = pot(ii,itarg) + rinc*rtmp
               grad(ii,1,itarg) = grad(ii,1,itarg) + rincx*rtmp
               grad(ii,2,itarg) = grad(ii,2,itarg) + rincy*rtmp
               grad(ii,3,itarg) = grad(ii,3,itarg) + rincz*rtmp
            enddo
         endif
      enddo
      enddo
      
      return
      end
c
c
c
c
c
c
c
c
c**********************************************************************
      subroutine g3d_directdh_vec(nd,delta,dmax,sources,ns,
     1           rnormal,dipstr,targ,ntarg,pot,grad,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT, gradients GRAD
c     and Hessians at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     
c     pot(ii)  = \sum_j \grad_s(exp(-r^2/delta)*charge(ii,j))*N(s)
c
c     grad  = d(pot)
c     hess = dxx,dxy,dyy
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     rnormal(2,ns) :   normals at sources
c     dipstr(nd,ns) :   dipole strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)     (complex *16)      : potential is incremented
c     grad(nd)  (complex *16)      : gradient is incremented
c     hess(nd)  (complex *16)      : Hessian is incremented
c---------------------------------------------------------------------
      integer i,ns,ifexpon,ii,nd
      real *8 sources(3,ns),rnormal(3,ns),targ(3,ntarg)
      real *8 xdiff,ydiff,zdiff,rr,r,rtmp,rinc,dx,dy,dz
      real *8 rincx,rincy,rincz,delta
      real *8 rincxx,rincyy,rinczz
      real *8 rincxy,rincxz,rincyz
      real *8 pot(nd,ntarg),grad(nd,3,ntarg),hess(nd,6,ntarg)
      real *8 dipstr(nd,ns)

c
      rfac = 2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         xdiff=targ(1,itarg)-sources(1,i)
         ydiff=targ(2,itarg)-sources(2,i)
         zdiff=targ(3,itarg)-sources(3,i)
         rr=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
         if (rr .lt. dmax) then         
            dtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            dy = rfac*ydiff
            dz = rfac*zdiff
            
            rinc = (dx*rnormal(1,i)+dy*rnormal(2,i)+dz*rnormal(3,i))
            rincx = rnormal(1,i)*rfac - rinc*dx
            rincy = rnormal(2,i)*rfac - rinc*dy
            rincz = rnormal(3,i)*rfac - rinc*dz
            
            rincxx = -rfac*(dx*rnormal(1,i)+rinc) - rincx*dx
            rincyy = -rfac*(dy*rnormal(2,i)+rinc) - rincy*dy
            rinczz = -rfac*(dz*rnormal(3,i)+rinc) - rincz*dz
            
            rincxy = -rfac*dx*rnormal(2,i) - rincx*dy
            rincxz = -rfac*dx*rnormal(3,i) - rincx*dz
            rincyz = -rfac*dy*rnormal(3,i) - rincy*dz
c
            do ii = 1,nd
               rtmp = dtmp*dipstr(ii,i)
               
               pot(ii,itarg) = pot(ii,itarg) + rinc*rtmp
               grad(ii,1,itarg) = grad(ii,1,itarg) + rincx*rtmp
               grad(ii,2,itarg) = grad(ii,2,itarg) + rincy*rtmp
               grad(ii,3,itarg) = grad(ii,3,itarg) + rincz*rtmp
               
               hess(ii,1,itarg) = hess(ii,1,itarg) + rincxx*rtmp
               hess(ii,2,itarg) = hess(ii,2,itarg) + rincyy*rtmp
               hess(ii,3,itarg) = hess(ii,3,itarg) + rinczz*rtmp

               hess(ii,4,itarg) = hess(ii,4,itarg) + rincxy*rtmp
               hess(ii,5,itarg) = hess(ii,5,itarg) + rincxz*rtmp
               hess(ii,6,itarg) = hess(ii,6,itarg) + rincyz*rtmp
            enddo
         endif
      enddo
      enddo

      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine g3d_directcdp_vec(nd,delta,dmax,sources,ns,
     $           charge,rnormal,dipstr,targ,ntarg,pot)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT
c     at the target point TARGET, due to a vector of 
c     charges and dipoles at SOURCE(2,ns). 
c
c     pot(ii)  = \sum_j exp(-r^2/delta)*charge(ii,j)
c                + \sum_j \grad_s(exp(-r^2/delta)*charge(ii,j))*N(s)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     rnormal(2,ns) :   normals at sources
c     dipstr(nd,ns) :   dipole strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)       : potential is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd
      real *8 sources(3,ns),rnormal(3,ns),targ(3,ntarg)
      real *8 xdiff,ydiff,zdiff,rr,r,rtmp,rinc,dx,dy,dz,delta
      real *8 pot(nd,ntarg)
      real *8 charge(nd,ns),dipstr(nd,ns)

c
      rfac = 2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         xdiff=targ(1,itarg)-sources(1,i)
         ydiff=targ(2,itarg)-sources(2,i)
         zdiff=targ(3,itarg)-sources(3,i)
         rr=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
         if (rr .lt. dmax) then         
            rtmp = dexp(-rr/delta)
            do ii=1,nd
               pot(ii,itarg) = pot(ii,itarg) + rtmp*charge(ii,i)
            enddo
c
            dx = rfac*xdiff
            dy = rfac*ydiff
            dz = rfac*zdiff
            rinc = (dx*rnormal(1,i)+dy*rnormal(2,i)
     1          +dz*rnormal(3,i))*rtmp
            do ii=1,nd
               pot(ii,itarg) = pot(ii,itarg) + rinc*dipstr(ii,i)
            enddo
         endif
      enddo
      enddo
      
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine g3d_directcdg_vec(nd,delta,dmax,sources,ns,
     1             charge,rnormal,dipstr,targ,ntarg,pot,grad)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT and gradients GRAD
c     at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     
c     pot(ii)  = \sum_j exp(-r^2/delta)*charge(ii,j)
c                + \sum_j \grad_s(exp(-r^2/delta)*charge(ii,j))*N(s)
c
c     grad  = d(pot)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     rnormal(2,ns) :   normals at sources
c     dipstr(nd,ns) :   charge strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)     (complex *16)      : potential is incremented
c     grad(nd)  (complex *16)      : gradient is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd
      real *8 sources(3,ns),rnormal(3,ns),targ(3,ntarg),delta
      real *8 xdiff,ydiff,zdiff,rr,r,rtmp,rinc,dx,dy,dz
      real *8 rincx,rincy,rincz
      real *8 pot(nd,ntarg),grad(nd,3,ntarg)
      real *8 charge(nd,ns),dipstr(nd,ns)

c
      rfac = -2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         xdiff=targ(1,itarg)-sources(1,i)
         ydiff=targ(2,itarg)-sources(2,i)
         zdiff=targ(3,itarg)-sources(3,i)
         rr=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
         if (rr .lt. dmax) then         
            dtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            dy = rfac*ydiff
            dz = rfac*zdiff
            do ii = 1,nd
               rinc = dtmp*charge(ii,i)
               pot(ii,itarg) = pot(ii,itarg) + rinc
               grad(ii,1,itarg) = grad(ii,1,itarg) + dx*rinc
               grad(ii,2,itarg) = grad(ii,2,itarg) + dy*rinc
               grad(ii,3,itarg) = grad(ii,3,itarg) + dz*rinc
            enddo
c
            rinc = -dx*rnormal(1,i)-dy*rnormal(2,i)-dz*rnormal(3,i)
            rincx = rinc*dx-rnormal(1,i)*rfac
            rincy = rinc*dy-rnormal(2,i)*rfac
            rincz = rinc*dz-rnormal(3,i)*rfac
c
            do ii = 1,nd
               rtmp = dtmp*dipstr(ii,i)
               pot(ii,itarg) = pot(ii,itarg) + rinc*rtmp
               grad(ii,1,itarg) = grad(ii,1,itarg) + rincx*rtmp
               grad(ii,2,itarg) = grad(ii,2,itarg) + rincy*rtmp
               grad(ii,3,itarg) = grad(ii,3,itarg) + rincz*rtmp
            enddo
         endif
      enddo
      enddo
      
      return
      end
c
c
c
c
c
c
c
c
c**********************************************************************
      subroutine g3d_directcdh_vec(nd,delta,dmax,sources,ns,
     1           charge,rnormal,dipstr,targ,ntarg,pot,grad,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT, gradients GRAD
c     and Hessians at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     
c
c     pot(ii)  = \sum_j exp(-r^2/delta)*charge(ii,j)
c                + \sum_j \grad_s(exp(-r^2/delta)*charge(ii,j))*N(s)
c
c     grad  = d(pot)
c     hess = dxx,dxy,dyy
c
c---------------------------------------------------------------------
c     INPUT:
c
c     delta         :   Gaussian variance
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     rnormal(2,ns) :   normals at sources
c     dipstr(nd,ns) :   dipole strengths
c     targ          :   location of the target
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)       : potential is incremented
c     grad(nd)      : gradient is incremented
c     hess(nd)      : Hessian is incremented
c---------------------------------------------------------------------
      integer i,ns,ifexpon,ii,nd
      real *8 sources(3,ns),targ(3,ntarg),rnormal(3,ns)
      real *8 xdiff,ydiff,zdiff,rr,r,rtmp,rinc,dx,dy,dz
      real *8 rincx,rincy,rincz,delta
      real *8 rincxx,rincyy,rinczz,rincxy,rincxz,rincyz
      real *8 dxx,dyy,dzz,dxy,dxz,dyz
      real *8 pot(nd,ntarg),grad(nd,3,ntarg),hess(nd,6,ntarg)
      real *8 charge(nd,ns),dipstr(nd,ns)


c
      rfac = -2.0d0/delta
      do itarg=1,ntarg
      do i = 1,ns
         xdiff=targ(1,itarg)-sources(1,i)
         ydiff=targ(2,itarg)-sources(2,i)
         zdiff=targ(3,itarg)-sources(3,i)
         rr=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
         if (rr .lt. dmax) then
            dtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            dy = rfac*ydiff
            dz = rfac*zdiff
            
            dxx = rfac + dx*dx
            dyy = rfac + dy*dy
            dzz = rfac + dz*dz

            dxy = dx*dy
            dxz = dx*dz
            dyz = dy*dz
c
            do ii = 1,nd
               rinc = dtmp*charge(ii,i)
               pot(ii,itarg) = pot(ii,itarg) + rinc
            
               grad(ii,1,itarg) = grad(ii,1,itarg) + dx*rinc
               grad(ii,2,itarg) = grad(ii,2,itarg) + dy*rinc
               grad(ii,3,itarg) = grad(ii,3,itarg) + dz*rinc

               hess(ii,1,itarg) = hess(ii,1,itarg) + dxx*rinc
               hess(ii,2,itarg) = hess(ii,2,itarg) + dyy*rinc
               hess(ii,3,itarg) = hess(ii,3,itarg) + dzz*rinc
               
               hess(ii,4,itarg) = hess(ii,4,itarg) + dxy*rinc
               hess(ii,5,itarg) = hess(ii,5,itarg) + dxz*rinc
               hess(ii,6,itarg) = hess(ii,6,itarg) + dyz*rinc
            enddo

ccc         rr=xdiff*xdiff+ydiff*ydiff
ccc         rtmp = dexp(-rr/delta)
            rinc = -dx*rnormal(1,i)-dy*rnormal(2,i)-dz*rnormal(3,i)
            rincx = rinc*dx-rnormal(1,i)*rfac
            rincy = rinc*dy-rnormal(2,i)*rfac
            rincz = rinc*dz-rnormal(3,i)*rfac
            
            rincxx = rfac*(rinc-dx*rnormal(1,i)) + rincx*dx
            rincyy = rfac*(rinc-dy*rnormal(2,i)) + rincy*dy
            rinczz = rfac*(rinc-dz*rnormal(3,i)) + rincz*dz

            rincxy = rincx*dy - rfac*dx*rnormal(2,i) 
            rincxz = rincx*dz - rfac*dx*rnormal(3,i) 
            rincyz = rincy*dz - rfac*dy*rnormal(3,i) 
            
c
            do ii = 1,nd
               rtmp = dtmp*dipstr(ii,i)
               pot(ii,itarg) = pot(ii,itarg) + rinc*rtmp
               grad(ii,1,itarg) = grad(ii,1,itarg) + rincx*rtmp
               grad(ii,2,itarg) = grad(ii,2,itarg) + rincy*rtmp
               grad(ii,3,itarg) = grad(ii,3,itarg) + rincz*rtmp
            
               hess(ii,1,itarg) = hess(ii,1,itarg) + rincxx*rtmp
               hess(ii,2,itarg) = hess(ii,2,itarg) + rincyy*rtmp
               hess(ii,3,itarg) = hess(ii,3,itarg) + rinczz*rtmp
               
               hess(ii,4,itarg) = hess(ii,4,itarg) + rincxy*rtmp
               hess(ii,5,itarg) = hess(ii,5,itarg) + rincxz*rtmp
               hess(ii,6,itarg) = hess(ii,6,itarg) + rincyz*rtmp
            enddo
         endif
      enddo
      enddo
      return
      end
c
c
c
