c     This file contains the following subroutines
c
c     g2d_directcp_vec : direct evaluation of the potential
c                     due to sources of strength charge()
c
c     g2d_directcg_vec : direct evaluation of the potential and gradient
c                     due to sources of strength charge()
c
c     g2d_directch_vec : direct evaluation of the potential and gradient and hessian
c                     due to sources of strength charge()
c
c     g2d_directdp_vec : direct evaluation of the potential
c                     due to dipoles
c
c     g2d_directdg_vec : direct evaluation of the potential and gradient
c                     due to dipoles
c
c     g2d_directdh_vec : direct evaluation of the potential and gradient and hessian
c                     due to dipoles
c
c     g2d_directcdp_vec : direct evaluation of the potential
c                     due to charges and dipoles
c
c     g2d_directcdg_vec : direct evaluation of the potential and gradient
c                     due to charges and dipoles
c
c     g2d_directcdh_vec : direct evaluation of the potential and gradient and hessian
c                     due to chrages and dipoles
c
c
c**********************************************************************
      subroutine g2d_directcp_vec(nd,delta,eps,sources,ns,charge,
     $           targ,pot)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT
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
c---------------------------------------------------------------------
      integer i,ns,ii,nd
      real *8 sources(2,ns),targ(2),xdiff,ydiff,rr,r
      real *8 rtmp,delta,eps,dmax
      real *8 pot(nd)
      real *8 charge(nd,ns)
c
      dmax=log(1/eps)*delta
      
      do i = 1,ns
         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            do ii=1,nd
               pot(ii) = pot(ii) + rtmp*charge(ii,i)
            enddo
         endif
      enddo
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine g2d_directcg_vec(nd,delta,eps,sources,ns,charge,
     1             targ,pot,grad)
      implicit none
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
      real *8 sources(2,ns),targ(2)
      real *8 xdiff,ydiff,rr,r,rtmp,delta,eps,dmax
      real *8 pot(nd),grad(nd,2)
      real *8 dx,dy,rinc,rfac
      real *8 charge(nd,ns)
c
      dmax=log(1/eps)*delta
      rfac=-2.0d0/delta
      
      do i = 1,ns
         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            dy = rfac*ydiff
            do ii = 1,nd
               rinc = rtmp*charge(ii,i)
               pot(ii) = pot(ii) + rinc
               grad(ii,1) = grad(ii,1) + dx*rinc
               grad(ii,2) = grad(ii,2) + dy*rinc
            enddo
         endif
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
      subroutine g2d_directch_vec(nd,delta,eps,sources,ns,charge,targ,
     1           pot,grad,hess)
      implicit none
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
c---------------------------------------------------------------------
      integer i,ns,ifexpon,ii,nd
      real *8 sources(2,ns),targ(2)
      real *8 xdiff,ydiff,rr,r,delta,eps,dmax
      real *8 pot(nd),grad(nd,2),hess(nd,3)
      real *8 rtmp,dx,dy,rinc,dxx,dxy,dyy,rfac
      real *8 charge(nd,ns)
c
      dmax=log(1/eps)*delta
      rfac = -2.0d0/delta
      
      do i = 1,ns
         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            dy = rfac*ydiff
            dxx = rfac + dx*dx
            dyy = rfac + dy*dy
            dxy = dx*dy
c
            do ii = 1,nd
               rinc = rtmp*charge(ii,i)
               pot(ii) = pot(ii) + rinc
               grad(ii,1) = grad(ii,1) + dx*rinc
               grad(ii,2) = grad(ii,2) + dy*rinc
               hess(ii,1) = hess(ii,1) + dxx*rinc
               hess(ii,2) = hess(ii,2) + dxy*rinc
               hess(ii,3) = hess(ii,3) + dyy*rinc
            enddo
         endif
      enddo

      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine g2d_directdp_vec(nd,delta,eps,sources,ns,
     $           rnormal,dipstr,targ,pot)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT
c     at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     We use the unscaled version of log
c     response: i.e., log|z|
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
      real *8 sources(2,ns),rnormal(2,ns),targ(2),xdiff,ydiff,rr,r
      real *8 dx,dy,rtmp,rinc,delta,rfac,eps,dmax
      real *8 pot(nd)
      real *8 dipstr(nd,ns)
c
      dmax=log(1/eps)*delta
      rfac = 2.0d0/delta
      
      do i = 1,ns
         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            dy = rfac*ydiff
            rinc = (dx*rnormal(1,i)+dy*rnormal(2,i))*rtmp
            do ii=1,nd
               pot(ii) = pot(ii) + rinc*dipstr(ii,i)
            enddo
         endif
      enddo
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine g2d_directdg_vec(nd,delta,eps,sources,ns,
     1           rnormal,dipstr,targ,pot,grad)
      implicit none
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
      real *8 sources(2,ns),rnormal(2,ns),targ(2)
      real *8 xdiff,ydiff,rr,r,rtmp,rinc,dx,dy,rincx,rincy,delta,rfac
      real *8 pot(nd),grad(nd,2),eps,dmax
      real *8 dipstr(nd,ns)
c
      dmax=log(1/eps)*delta
      rfac = 2.0d0/delta
      
      do i = 1,ns
         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            dy = rfac*ydiff
            rinc = (dx*rnormal(1,i)+dy*rnormal(2,i))
            rincx = rnormal(1,i)*rfac - rinc*dx
            rincy = rnormal(2,i)*rfac - rinc*dy
            rinc = rinc*rtmp
            rincx = rincx*rtmp
            rincy = rincy*rtmp
c
            do ii = 1,nd
               pot(ii) = pot(ii) + rinc*dipstr(ii,i)
               grad(ii,1) = grad(ii,1) + rincx*dipstr(ii,i)
               grad(ii,2) = grad(ii,2) + rincy*dipstr(ii,i)
            enddo
         endif
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
      subroutine g2d_directdh_vec(nd,delta,eps,sources,ns,
     1           rnormal,dipstr,targ,pot,grad,hess)
      implicit none
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
      real *8 sources(2,ns),rnormal(2,ns),targ(2)
      real *8 xdiff,ydiff,rr,r,rtmp,rinc,dx,dy,rincx,rincy,delta
      real *8 rincxx,rincxy,rincyy,rfac,eps,dmax
      real *8 pot(nd),grad(nd,2),hess(nd,3)
      real *8 dipstr(nd,ns)

c
      dmax=log(1/eps)*delta
      rfac = 2.0d0/delta
      
      do i = 1,ns
         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            dy = rfac*ydiff
            rinc = (dx*rnormal(1,i)+dy*rnormal(2,i))
            rincx = rnormal(1,i)*rfac - rinc*dx
            rincy = rnormal(2,i)*rfac - rinc*dy
            rincxx = -rfac*(dx*rnormal(1,i)+rinc) - rincx*dx
            rincyy = -rfac*(dy*rnormal(2,i)+rinc) - rincy*dy
            rincxy = -rfac*dx*rnormal(2,i) - rincx*dy
c     
            do ii = 1,nd
               rtmp = rtmp*dipstr(ii,i)
               pot(ii) = pot(ii) + rinc*rtmp
               grad(ii,1) = grad(ii,1) + rincx*rtmp
               grad(ii,2) = grad(ii,2) + rincy*rtmp
               hess(ii,1) = hess(ii,1) + rincxx*rtmp
               hess(ii,2) = hess(ii,2) + rincxy*rtmp
               hess(ii,3) = hess(ii,3) + rincyy*rtmp
            enddo
         endif
      enddo
      

      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine g2d_directcdp_vec(nd,delta,eps,sources,ns,charge,
     $           rnormal,dipstr,targ,pot)
      implicit none
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
      real *8 sources(2,ns),rnormal(2,ns),targ(2)
      real *8 xdiff,ydiff,rr,r,rtmp,dx,dy,delta,rfac,rinc
      real *8 pot(nd),eps,dmax
      real *8 charge(nd,ns),dipstr(nd,ns)

c
      dmax=log(1/eps)*delta
      rfac = 2.0d0/delta
      
      do i = 1,ns
         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if (rr .lt. dmax) then
         
            dx = rfac*xdiff
            dy = rfac*ydiff
         
            rtmp = dexp(-rr/delta)
            do ii=1,nd
               pot(ii) = pot(ii) + rtmp*charge(ii,i)
            enddo

            rinc = (dx*rnormal(1,i)+dy*rnormal(2,i))*rtmp
            do ii = 1,nd
               pot(ii) = pot(ii) + rinc*dipstr(ii,i)
            enddo
         endif
      enddo
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine g2d_directcdg_vec(nd,delta,eps,sources,ns,charge,
     1             rnormal,dipstr,targ,pot,grad)
      implicit none
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
      real *8 sources(2,ns),rnormal(2,ns),targ(2),delta
      real *8 xdiff,ydiff,rr,r,rtmp,rinc,dx,dy,rincx,rincy,rfac
      real *8 pot(nd),grad(nd,2),eps,dmax
      real *8 charge(nd,ns),dipstr(nd,ns)

c
      dmax=log(1/eps)*delta
      rfac = -2.0d0/delta
      
      do i = 1,ns
         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            dy = rfac*ydiff
            do ii = 1,nd
               rinc = rtmp*charge(ii,i)
               pot(ii) = pot(ii) + rinc
               grad(ii,1) = grad(ii,1) + dx*rinc
               grad(ii,2) = grad(ii,2) + dy*rinc
            enddo
c
            rinc = -dx*rnormal(1,i)-dy*rnormal(2,i)
            rincx = rinc*dx - rfac*rnormal(1,i) 
            rincy = rinc*dy - rfac*rnormal(2,i)
c
            do ii = 1,nd
               rtmp = rtmp*dipstr(ii,i)
               pot(ii) = pot(ii) + rinc*rtmp
               grad(ii,1) = grad(ii,1) + rincx*rtmp
               grad(ii,2) = grad(ii,2) + rincy*rtmp
            enddo
         endif
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
      subroutine g2d_directcdh_vec(nd,delta,eps,sources,ns,charge,
     1           rnormal,dipstr,targ,pot,grad,hess)
      implicit none
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
      real *8 sources(2,ns),targ(2),rnormal(2,ns)
      real *8 xdiff,ydiff,rr,r,rtmp,rinc,dx,dy,rincx,rincy,delta
      real *8 rincxx,rincxy,rincyy,dxx,dxy,dyy,rfac,eps,dmax
      real *8 pot(nd),grad(nd,2),hess(nd,3)
      real *8 charge(nd,ns),dipstr(nd,ns)


c
      dmax=log(1/eps)*delta
      rfac = -2.0d0/delta

      do i = 1,ns
         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if (rr .lt. dmax) then
            rtmp = dexp(-rr/delta)
            dx = rfac*xdiff
            dy = rfac*ydiff
            dxx = rfac + dx*dx
            dyy = rfac + dy*dy
            dxy = dx*dy
c
            do ii = 1,nd
               rinc = rtmp*charge(ii,i)
               pot(ii) = pot(ii) + rinc
               grad(ii,1) = grad(ii,1) + dx*rinc
               grad(ii,2) = grad(ii,2) + dy*rinc
               hess(ii,1) = hess(ii,1) + dxx*rinc
               hess(ii,2) = hess(ii,2) + dxy*rinc
               hess(ii,3) = hess(ii,3) + dyy*rinc
            enddo

ccc         rr=xdiff*xdiff+ydiff*ydiff
ccc         rtmp = dexp(-rr/delta)
            rinc = -dx*rnormal(1,i)-dy*rnormal(2,i)
            rincx = rinc*dx - rnormal(1,i)*rfac  
            rincy = rinc*dy - rnormal(2,i)*rfac 
            rincxx = rfac*(rinc-dx*rnormal(1,i)) + rincx*dx
            rincyy = rfac*(rinc-dy*rnormal(2,i)) + rincy*dy
            rincxy = rincx*dy - rfac*dx*rnormal(2,i) 
c
            do ii = 1,nd
               rtmp = rtmp*dipstr(ii,i)
               pot(ii) = pot(ii) + rinc*rtmp
               grad(ii,1) = grad(ii,1) + rincx*rtmp
               grad(ii,2) = grad(ii,2) + rincy*rtmp
               hess(ii,1) = hess(ii,1) + rincxx*rtmp
               hess(ii,2) = hess(ii,2) + rincxy*rtmp
               hess(ii,3) = hess(ii,3) + rincyy*rtmp
            enddo
         endif
      enddo
      
      return
      end
c
c
c
