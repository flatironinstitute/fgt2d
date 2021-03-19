c     This file contains the following subroutines
c
c     g2dformhc_vec : computes the Hermite expansion about the center CENT
c                     due to sources of strength charge()
c
c     g2dformhd_vec : computes the Hermite expansion about the center CENT
c                     due to dipoles
c
c     g2dformhcd_vec : computes the Hermite expansion about the center CENT
c                     due to charges and dipoles
c
C     g2dhevalp_vec : evaluates the Hermite expansion 
C                         potential only.
c      
C     g2dhevalg_vec : evaluates the Hermite expansion 
C                         pot/grad
c      
C     g2dhevalh_vec : evaluates the Hermite expansion 
C                         pot/grad/Hessian
      
c     g2dformlc_vec : computes the Taylor expansion about the center CENT
c                     due to sources of strength charge()
c
c     g2dformld_vec : computes the Taylor expansion about the center CENT
c                     due to dipoles
c
c     g2dformlcd_vec : computes the Taylor expansion about the center CENT
c                     due to charges and dipoles
c
C     g2dlevalp_vec : evaluates the local Taylor expansion 
C                         potential only.
c      
C     g2dlevalg_vec : evaluates the local Taylor expansion 
C                         pot/grad
c      
C     g2dlevalh_vec : evaluates the local Taylor expansion 
C                         pot/grad/Hessian
c      
C     g2dh2sx_h2s_vec : converts the Hermite expansion to the four hybrid
C                       SOE/X expansions and four SOE expansions
c
c     g2dherm2local_vec : converts the Hermite expansion to the Taylor expansion
c      
c     g2dsoe2local_vec : converts four SOE expansions to the local Taylor expansion
c
c     g2dsx2local_vec : converts four SOE/X expansions to the local Taylor expansion
c
c     g2dh2lmat : returns the matrix converting the Hermite expansion
C                 to the local Taylor expansion (1D only)
c      
c     g2dh2smat : returns the matrix converting the Hermite expansion
C                  to the SOE+ expansions (1D only)
C     
c     g2dh2xmat : returns the matrix converting the Hermite expansion 
C                 to the planewave expansion (1D only)
C     
C     g2ds2lmat : returns the matrix converting the SOE+ expansion
C                 to the local Taylor expansion (1D only)
C
C     g2dx2lmat : returns the matrix converting the planewave expansion
C                 to the local Taylor expansion (1D only)
C
c
c*********************************************************************
C
C form Hermite expansions (charge, dipole, charge & dipole)
C
C*********************************************************************
      subroutine g2dformhc_vec(nd,delta,sources,ns,charge,cent,
     1            nterms,ffexp)
C
C     This subroutine computes the Hermite expansion about
C     the center CENT due to the sources at locations
C     SOURCES  of strength CHARGE.
C
C     INPUT:
C
c     nd            = vector length (parallel)
C     delta         = Gaussian variance
C     sources   = source locations
C     ns        = number of sources
C     charge    = strengths of sources
C     cent      = center of the expansion
C     nterms    = number of terms in expansion
C
C     OUTPUT:
C
C     ffexp     = Hermite expansion
C---------------------------------------------------------------------
C
      implicit real *8 (a-h,o-z)
      integer ns,nterms
      real *8 cent(2),sources(2,ns),charge(nd,ns)
      real *8 ffexp(0:nterms,0:nterms,nd)
      real *8 x,y,chg
      real *8, allocatable ::  xp(:),yp(:)
C
C     initialize coefficients to zero.
C
      allocate(xp(0:nterms))
      allocate(yp(0:nterms))
c
      dsq = 1.0D0/dsqrt(delta)
C
C     accumulate expansion due to each source.
C
      do i=1,ns
         x = (sources(1,i) - cent(1))*dsq
         y = (sources(2,i) - cent(2))*dsq
         xp(0) = 1.0D0
         yp(0) = 1.0D0
         do k = 1,nterms
            tmp = dsqrt(dble(k))
            xp(k) = xp(k-1)*x/tmp
            yp(k) = yp(k-1)*y/tmp
         enddo
c
         do ind=1,nd
            chg = charge(ind,i)
            do k=0,nterms
               ytmp = chg*yp(k)
               do j=0,nterms
                  ffexp(j,k,ind) = ffexp(j,k,ind) + ytmp*xp(j)
               enddo
            enddo
         enddo
      enddo
      return
      end
c
C
C
C*********************************************************************C
      subroutine g2dformhd_vec(nd,delta,sources,ns,rnormal,dipstr,
     1            center,nterms,ffexp)
C
C     This subroutine computes the Hermite expansion about
C     the center CENT due to the dipole sources at locations
C     SOURCES  of strength CHARGE.
C
C     INPUT:
C
c     nd        = vector length (parallel)
C     delta     = Gaussian variance
C     sources   = source locations
C     ns        = number of sources
C     rnormal   = dipole directions
C     dipstr    = dipole strengths 
C     center    = center of the expansion
C     nterms    = number of terms in expansion
C
C     OUTPUT:
C
C     ffexp     = Hermite expansion
C---------------------------------------------------------------------
C
      implicit none
      integer ns,nterms,nd,i,j,ind,k,j1
      real *8 center(2),sources(2,ns)
      real *8 rnormal(2,ns),dipstr(nd,ns)
      real *8 ffexp(0:nterms,0:nterms,nd)
      real *8 x,y,delta,dsq,tmp,d1,d2
      real *8, allocatable ::  hx(:),hy(:)
      real *8, allocatable ::  dxhx(:),dyhy(:)
C
C     initialize coefficients to zero.
C
      allocate(hx(0:nterms))
      allocate(hy(0:nterms))
      allocate(dxhx(0:nterms),dyhy(0:nterms))
c
      dsq = 1.0D0/dsqrt(delta)
C
C     accumulate expansion due to each source.
C
      do i=1,ns
         x = (sources(1,i) - center(1))*dsq
         y = (sources(2,i) - center(2))*dsq
         hx(0) = 1.0D0
         hy(0) = 1.0D0
         dxhx(0) = 0.0D0
         dyhy(0) = 0.0D0
         do j1 = 1,nterms
            tmp = dsqrt(dble(j1))
            hy(j1)=hy(j1-1)*y/tmp
            hx(j1)=hx(j1-1)*x/tmp
            dxhx(j1) = j1*dsq*hx(j1-1)/tmp
            dyhy(j1) = j1*dsq*hy(j1-1)/tmp
         enddo
c
         do ind = 1,nd
            do k=0,nterms
               d1 = hy(k)*rnormal(1,i)*dipstr(ind,i)
               d2 = dyhy(k)*rnormal(2,i)*dipstr(ind,i)
               do j=0,nterms
                  ffexp(j,k,ind) = ffexp(j,k,ind) + 
     1                dxhx(j)*d1+hx(j)*d2
               enddo
            enddo
         enddo
      enddo
      return
      end
C
C
C
C*********************************************************************C
      subroutine g2dformhcd_vec(nd,delta,sources,ns,charge,rnormal,
     1            dipstr,center,nterms,ffexp)
C
C     This subroutine computes the Hermite expansion about
C     the center CENT due to the dipole sources at locations
C     SOURCES  of strength CHARGE.
C
C     INPUT:
C
c     nd        = vector length (parallel)
C     delta     = Gaussian variance
C     sources   = source locations
C     ns        = number of sources
C     charge    = charge strengths
C     rnormal   = dipole directions
C     dipstr    = dipole strengths 
C     center    = center of the expansion
C     nterms    = number of terms in expansion
C
C     OUTPUT:
C
C     ffexp     = Hermite expansion
C---------------------------------------------------------------------
C
      implicit none
      integer ns,nterms,nd,i,j,ind,k,j1
      real *8 center(2),sources(2,ns)
      real *8 rnormal(2,ns),dipstr(nd,ns),charge(nd,ns)
      real *8 ffexp(0:nterms,0:nterms,nd)
      real *8 x,y,delta,dsq,tmp,chg,ytmp,d1,d2
      real *8, allocatable ::  hx(:),hy(:)
      real *8, allocatable ::  dxhx(:),dyhy(:)
C
C     initialize coefficients to zero.
C
      allocate(hx(0:nterms))
      allocate(hy(0:nterms))
      allocate(dxhx(0:nterms),dyhy(0:nterms))
c
      dsq = 1.0D0/dsqrt(delta)
C
C     accumulate expansion due to each source.
C
      do i=1,ns
         x = (sources(1,i) - center(1))*dsq
         y = (sources(2,i) - center(2))*dsq
         hx(0) = 1.0D0
         hy(0) = 1.0D0
         dxhx(0) = 0.0D0
         dyhy(0) = 0.0D0
         do j1 = 1,nterms
            tmp = dsqrt(dble(j1))
            hy(j1)=hy(j1-1)*y/tmp
            hx(j1)=hx(j1-1)*x/tmp
            dxhx(j1) = j1*dsq*hx(j1-1)/tmp
            dyhy(j1) = j1*dsq*hy(j1-1)/tmp
         enddo
c
         do ind = 1,nd
            chg = charge(ind,i)
            do k=0,nterms
               ytmp = chg*hy(k)
               d1 = hy(k)*rnormal(1,i)*dipstr(ind,i)
               d2 = dyhy(k)*rnormal(2,i)*dipstr(ind,i)
               do j=0,nterms
                  ffexp(j,k,ind) = ffexp(j,k,ind) + 
     1                dxhx(j)*d1+hx(j)*d2+hx(j)*ytmp
               enddo
            enddo
         enddo
      enddo
      return
      end
c
c
c
c
c*********************************************************************
C
C eval Hermite expansions (pot, pot+grad, pot+grad+hess)
C
C*********************************************************************
      subroutine g2dhevalp_vec(nd,delta,cent,nterms,ffexp
     1    targ,ntarg,pot)
C
C     This subroutine evaluates the far field expansion FFEXP about
C     CENT at location TARG.
C
C     INPUT:
C
c     nd       = vector length (parallel)
C     delta    = Gaussian variance
C     ffexp    = coefficients of far field expansion
C     nterms   = number of terms in expansion
C     cent     = center of the expansion
C     targ     = target location
C
C     OUTPUT:
C
C     POT = evaluated potential
C---------------------------------------------------------------------
C
      implicit real *8 (a-h,o-z)
      integer nterms
      real *8 cent(2),pot(nd,ntarg)
      real *8 ffexp(0:nterms,0:nterms,nd)
      real *8 x,y,targ(2,ntarg)
      real *8, allocatable ::  hexpx(:),hexpy(:)
C
C     initialize coefficients to zero.
C
      allocate(hexpx(0:nterms))
      allocate(hexpy(0:nterms))
C
C     tabulate Hermite polynomials in HEXPX and HEXPY.
C
      dsq = 1.0d0/dsqrt(delta)
      do itarg=1,ntarg
         x = (targ(1,itarg) - cent(1))*dsq
         y = (targ(2,itarg) - cent(2))*dsq
         facx = dexp(-x*x)
         facy = dexp(-y*y)
         hexpx(0) = 1.0d0*facx
         hexpy(0) = 1.0d0*facy
         hexpx(1) = 2.0d0*x*facx
         hexpy(1) = 2.0d0*y*facy
         do i = 1,nterms-1
            hexpx(i+1)=2.0d0*(x/dsqrt((i+1)*1.0d0)*hexpx(i)
     1          -dsqrt(1.0d0*i/(i+1))*hexpx(i-1))
            hexpy(i+1)=2.0d0*(y/dsqrt((i+1)*1.0d0)*hexpy(i)
     1          -dsqrt(1.0d0*i/(i+1))*hexpy(i-1))
         enddo
c
C---- EVALUATE EXPANSION
c
         do ind = 1,nd
            dd=0
            do k = 0,nterms
               d2=0
               do j = 0,nterms
                  d2=d2+hexpx(j)*ffexp(j,k,ind)
               enddo
               dd=dd+d2*hexpy(k)
            enddo
            pot(ind,itarg) = pot(ind,itarg) + dd
         enddo
      enddo
      
      return
      end
C
c
C***********************************************************************
      subroutine g2dhevalg_vec(nd,delta,cent,nterms,ffexp
     1           targ,ntarg,pot,grad)

C
C     This subroutine evaluates the far field expansion FFEXP about
C     CENT at location (XP,YP).
C
C     INPUT:
C
c     nd       = vector length (parallel)
C     delta    = Gaussian variance
C     ffexp    = coefficients of far field expansion
C     nterms   = number of terms in expansion
C     cent     = center of the expansion
C     targ     = target location
C
C     OUTPUT:
C
C     pot      = evaluated potential
C     grad     = evaluated gradient
C---------------------------------------------------------------------
C
      implicit real *8 (a-h,o-z)
      integer nterms
      real *8 cent(2),pot(nd,ntarg),grad(nd,2,ntarg)
      real *8 ffexp(0:nterms,0:nterms,nd)
      real *8 x,y,targ(2,ntarg)
      real *8, allocatable ::  hexpx(:),hexpy(:)
      real *8, allocatable ::  dhexpx(:),dhexpy(:)
C
C     initialize coefficients to zero.
C
      allocate(hexpx(0:nterms+1))
      allocate(dhexpx(0:nterms))
      allocate(hexpy(0:nterms+1))
      allocate(dhexpy(0:nterms))
C
C     tabulate Hermite polynomials in HEXPX and HEXPY.
C
      dsq = 1.0d0/dsqrt(delta)
      do itarg=1,ntarg
         x = (targ(1,itarg) - cent(1))*dsq
         y = (targ(2,itarg) - cent(2))*dsq
         facx = dexp(-x*x)
         facy = dexp(-y*y)
         hexpx(0) = 1.0d0*facx
         hexpy(0) = 1.0d0*facy
         hexpx(1) = 2.0d0*x*facx
         hexpy(1) = 2.0d0*y*facy
         do i = 1,nterms
            hexpx(i+1)=2.0d0*(x/dsqrt((i+1)*1.0d0)*hexpx(i)
     1          -dsqrt(1.0d0*i/(i+1))*hexpx(i-1))
            hexpy(i+1)=2.0d0*(y/dsqrt((i+1)*1.0d0)*hexpy(i)
     1          -dsqrt(1.0d0*i/(i+1))*hexpy(i-1))
         enddo
         do i = 0,nterms
            dhexpx(i)=-hexpx(i+1)*dsqrt(i+1.0d0)
            dhexpy(i)=-hexpy(i+1)*dsqrt(i+1.0d0)
         enddo
c
C---- EVALUATE EXPANSION
c
         do ind = 1,nd
            dd=0
            gx = 0.0d0
            gy = 0.0d0
            do k = 0,nterms
               d2=0
               g1=0
               g2=0
               do j = 0,nterms
                  d2=d2+hexpx(j)*ffexp(j,k,ind)
                  g1 = g1 + dhexpx(j)*ffexp(j,k,ind)
                  g2 = g2 + hexpx(j)*ffexp(j,k,ind)
               enddo
               dd=dd+d2*hexpy(k)
               gx=gx+g1*hexpy(k)
               gy=gy+g2*dhexpy(k)
            enddo
            pot(ind,itarg) = pot(ind,itarg) + dd
            grad(ind,1,itarg) = grad(ind,1,itarg) + gx*dsq
            grad(ind,2,itarg) = grad(ind,2,itarg) + gy*dsq
         enddo
      enddo
      
      return
      end
C
C
c
C***********************************************************************
      subroutine g2dhevalh_vec(nd,delta,cent,nterms,ffexp,
     1           targ,ntarg,pot,grad,hess)

C
C     This subroutine evaluates the far field expansion FFEXP about
C     CENT at location (XP,YP).
C
C     INPUT:
C
c     nd       = vector length (parallel)
C     delta    = Gaussian variance
C     ffexp    = coefficients of far field expansion
C     nterms   = number of terms in expansion
C     cent     = center of the expansion
C     targ     = target location
C
C     OUTPUT:
C
C     pot      = evaluated potential
C     grad     = evaluated gradient
C     hess     = evaluated Hessian
C---------------------------------------------------------------------
C
      implicit real *8 (a-h,o-z)
      integer nterms
      real *8 cent(2),pot(nd,ntarg),grad(nd,2,ntarg),hess(nd,3,ntarg)
      real *8 ffexp(0:nterms,0:nterms,nd)
      real *8 x,y,targ(2,ntarg)
      real *8, allocatable ::  hexpx(:),hexpy(:)
      real *8, allocatable ::  dhexpx(:),dhexpy(:)
      real *8, allocatable ::  ddhexpx(:),ddhexpy(:)
C
C     initialize coefficients to zero.
C
      allocate(hexpx(0:nterms+2))
      allocate(dhexpx(0:nterms))
      allocate(ddhexpx(0:nterms))
      allocate(hexpy(0:nterms+2))
      allocate(dhexpy(0:nterms))
      allocate(ddhexpy(0:nterms))
C
C     tabulate Hermite polynomials in HEXPX and HEXPY.
C
      dsq = 1.0d0/dsqrt(delta)
      do itarg=1,ntarg
         x = (targ(1,itarg) - cent(1))*dsq
         y = (targ(2,itarg) - cent(2))*dsq
         facx = dexp(-x*x)
         facy = dexp(-y*y)
         hexpx(0) = 1.0d0*facx
         hexpy(0) = 1.0d0*facy
         hexpx(1) = 2.0d0*x*facx
         hexpy(1) = 2.0d0*y*facy
         do i = 1,nterms+1
            hexpx(i+1)=2.0d0*(x/dsqrt((i+1)*1.0d0)*hexpx(i)
     1          -dsqrt(1.0d0*i/(i+1))*hexpx(i-1))
            hexpy(i+1)=2.0d0*(y/dsqrt((i+1)*1.0d0)*hexpy(i)
     1          -dsqrt(1.0d0*i/(i+1))*hexpy(i-1))
         enddo
         do i = 0,nterms
            dhexpx(i)=-hexpx(i+1)*dsqrt(i+1.0d0)
            ddhexpx(i)=hexpx(i+2)*dsqrt(i+1.0d0)*dsqrt(i+2.0d0)
            dhexpy(i)=-hexpy(i+1)*dsqrt(i+1.0d0)
            ddhexpy(i)=hexpy(i+2)*dsqrt(i+1.0d0)*dsqrt(i+2.0d0)
         enddo
c
c
C---- EVALUATE EXPANSION
c
         do ind = 1,nd
            dd=0
            gx = 0.0d0
            gy = 0.0d0
            hxx = 0.0d0
            hyy = 0.0d0
            hxy = 0.0d0
            do k = 0,nterms
               d2=0
               g1=0
               g2=0
               h11=0
               h22=0
               h12=0
               do j = 0,nterms
                  d2  = d2  +   hexpx(j)*ffexp(j,k,ind)
                  g1  = g1  +  dhexpx(j)*ffexp(j,k,ind)
                  g2  = g2  +   hexpx(j)*ffexp(j,k,ind)               
                  h11 = h11 + ddhexpx(j)*ffexp(j,k,ind)
                  h22 = h22 +   hexpx(j)*ffexp(j,k,ind)
                  h12 = h12 +  dhexpx(j)*ffexp(j,k,ind)
               enddo
               dd  = dd  + d2   *hexpy(k)
               gx  = gx  + g1   *hexpy(k)
               gy  = gy  + g2  *dhexpy(k)
               hxx = hxx + h11  *hexpy(k)
               hyy = hyy + h22*ddhexpy(k)
               hxy = hxy + h12 *dhexpy(k)
            enddo
         
            pot(ind,itarg)    = pot(ind,itarg) + dd
            grad(ind,1,itarg) = grad(ind,1,itarg) + gx*dsq
            grad(ind,2,itarg) = grad(ind,2,itarg) + gy*dsq
            hess(ind,1,itarg) = hess(ind,1,itarg) + hxx*dsq*dsq
            hess(ind,2,itarg) = hess(ind,2,itarg) + hxy*dsq*dsq
            hess(ind,3,itarg) = hess(ind,3,itarg) + hyy*dsq*dsq
         enddo
      enddo
      
      return
      end
C
C************************************************************************
C
C     form local subroutines (charge, dipole, charge+dipole)
C
C***********************************************************************
      subroutine g2dformlc_vec(nd,delta,sources,ns,charge,center,
     1            nlocal,local)
C
C     This subroutine computes the Taylor expansions about
C     the center CENT due to sources of strength charge().
C
C     INPUT
C
c     nd            = vector length (parallel)
C     delta         = Gaussian variance
C     sources(2,ns) = coordinates of sources
C     ns            = number of sources
C     charge        = strength of sources
C     center        = center of the expansion
C     nlocal        = number of terms in local exp
C
C     OUTPUT:
C
C     local         = local expansions 
C
      implicit real*8 (a-h,o-z)
      real *8 center(2),sources(2,ns),charge(nd,ns)
      real *8 local(0:nlocal,0:nlocal,nd)
      real *8 dd, d2
      real *8, allocatable :: hexpx(:),hexpy(:)
c
      allocate(hexpx(0:nlocal),hexpy(0:nlocal))
c
      dsq = 1.0D0/dsqrt(delta)
C
      do i=1,ns
         x = (sources(1,i) - center(1))*dsq
         y = (sources(2,i) - center(2))*dsq
c

         facx=dexp(-x*x)
         facy=dexp(-y*y)

         hexpx(0)=1.0d0*facx
         hexpy(0)=1.0d0*facy
         hexpx(1)=2.0d0*x*facx
         hexpy(1)=2.0d0*y*facy

         do j=1,nlocal-1
            hexpx(j+1)=2.0d0*(x*hexpx(j)-hexpx(j-1))/(j+1)
            hexpy(j+1)=2.0d0*(y*hexpy(j)-hexpy(j-1))/(j+1)
         enddo

         do ind=1,nd
            chg=charge(ind,i)
            do j1=0,nlocal
               ytmp = chg*hexpy(j1) 
               do j2=0,nlocal
                  local(j2,j1,ind) = local(j2,j1,ind)+ytmp*hexpx(j2)
               enddo
            enddo
            
         enddo
      enddo
      return
      end
C
C
C
      subroutine g2dformld_vec(nd,delta,sources,ns,rnormal,dipstr,
     1            center,nlocal,local)
C
C     This subroutine computes the Taylor expansions about
C     the center CENT due to sources of strength charge().
C
C     INPUT
C
c     nd            = vector length (parallel)
C     delta         = Gaussian variance
C     sources(2,ns) = coordinates of sources
C     ns            = number of sources
C     rnormal   = dipole directions
C     dipstr    = dipole strengths 
C     center        = center of the expansion
C     nlocal        = number of terms in local exp
C
C     OUTPUT:
C
C     local         = local expansions 
C
      implicit real*8 (a-h,o-z)
      real *8 center(2),sources(2,ns),dipstr(nd,ns)
      real *8 rnormal(2,ns)
      real *8 local(0:nlocal,0:nlocal,nd)
      real *8 dd, d2
      real *8, allocatable :: hexpx(:),hexpy(:)
      real *8, allocatable :: dhexpx(:),dhexpy(:)
c
      allocate(hexpx(0:nlocal+1),hexpy(0:nlocal+1))
      allocate(dhexpx(0:nlocal),dhexpy(0:nlocal))
c
      dsq = 1.0D0/dsqrt(delta)
C
      do i=1,ns
         x = (sources(1,i) - center(1))*dsq
         y = (sources(2,i) - center(2))*dsq
c

         facx=dexp(-x*x)
         facy=dexp(-y*y)

         hexpx(0)=1.0d0*facx
         hexpy(0)=1.0d0*facy
         hexpx(1)=2.0d0*x*facx
         hexpy(1)=2.0d0*y*facy

         do j=1,nlocal
            hexpx(j+1)=2.0d0*(x*hexpx(j)-hexpx(j-1))/(j+1)
            hexpy(j+1)=2.0d0*(y*hexpy(j)-hexpy(j-1))/(j+1)
         enddo
         do j = 0,nlocal
            dhexpx(j)=-hexpx(j+1)*dsq*(j+1.0d0)
            dhexpy(j)=-hexpy(j+1)*dsq*(j+1.0d0)
         enddo

         do ind=1,nd
            do j1=0,nlocal
               do j2=0,nlocal
                  local(j2,j1,ind) = local(j2,j1,ind)+
     1          (dhexpx(j2)*hexpy(j1)*rnormal(1,i)+
     2           hexpx(j2)*dhexpy(j1)*rnormal(2,i))*dipstr(ind,i) 
               enddo
            enddo
            
         enddo
      enddo
      return
      end
C
C
C
C
C
      subroutine g2dformlcd_vec(nd,delta,sources,ns,charge,rnormal,
     1            dipstr,center,nlocal,local)
C
C     This subroutine computes the Taylor expansions about
C     the center CENT due to sources of strength charge().
C
C     INPUT
C
c     nd            = vector length (parallel)
C     delta         = Gaussian variance
C     sources(2,ns) = coordinates of sources
C     ns            = number of sources
C     rnormal   = dipole directions
C     dipstr    = dipole strengths 
C     center        = center of the expansion
C     nlocal        = number of terms in local exp
C
C     OUTPUT:
C
C     local         = local expansions 
C
      implicit real*8 (a-h,o-z)
      real *8 center(2),sources(2,ns),dipstr(nd,ns)
      real *8 rnormal(2,ns), charge(nd,ns)
      real *8 local(0:nlocal,0:nlocal,nd)
      real *8 dd, d2
      real *8, allocatable :: hexpx(:),hexpy(:)
      real *8, allocatable :: dhexpx(:),dhexpy(:)
c
      allocate(hexpx(0:nlocal+1),hexpy(0:nlocal+1))
      allocate(dhexpx(0:nlocal),dhexpy(0:nlocal))
c
      dsq = 1.0D0/dsqrt(delta)
C
      do i=1,ns
         x = (sources(1,i) - center(1))*dsq
         y = (sources(2,i) - center(2))*dsq
c

         facx=dexp(-x*x)
         facy=dexp(-y*y)

         hexpx(0)=1.0d0*facx
         hexpy(0)=1.0d0*facy
         hexpx(1)=2.0d0*x*facx
         hexpy(1)=2.0d0*y*facy

         do j=1,nlocal
            hexpx(j+1)=2.0d0*(x*hexpx(j)-hexpx(j-1))/(j+1)
            hexpy(j+1)=2.0d0*(y*hexpy(j)-hexpy(j-1))/(j+1)
         enddo
         do j = 0,nlocal
            dhexpx(j)=-hexpx(j+1)*dsq*(j+1.0d0)
            dhexpy(j)=-hexpy(j+1)*dsq*(j+1.0d0)
         enddo

         do ind=1,nd
            chg=charge(ind,i)
            do j1=0,nlocal
               ytmp=chg*hexpy(j1)
               do j2=0,nlocal
                  local(j2,j1,ind) = local(j2,j1,ind)+
     1          (dhexpx(j2)*hexpy(j1)*rnormal(1,i)+
     2           hexpx(j2)*dhexpy(j1)*rnormal(2,i))*dipstr(ind,i) 
                  local(j2,j1,ind) = local(j2,j1,ind)
     1                +ytmp*hexpx(j2)
               enddo
            enddo
            
         enddo
      enddo
      return
      end
C
C
C
C
C************************************************************************
C
C     eval local subroutines (pot, pot+grad, pot+grad+hess)
C
C***********************************************************************
C
      subroutine g2dlevalp_vec(nd,delta,center,nlocal,local,
     1              targ,ntarg,pot)
C
C     This subroutine evaluates the local expansion about
C     CENT at location TARG.
C
C     INPUT
C
c     nd            = vector length (for vector input)
C     delta         = Gaussian variance
C     center        = center of the expansion
C     nlocal        = number of terms in local expansions
C     local         = local Taylor expansion
c
C     targ          = target location
C
C     OUTPUT:
C
C     pot           = potential (or vectorized potentials) incremented
C
      implicit real*8 (a-h,o-z)
      real *8 local(0:nlocal,0:nlocal,nd)
      real *8 center(2),targ(2,ntarg)
      real *8 pot(nd,ntarg),xp(0:100),yp(0:100)
C
      dsq = 1.0D0/dsqrt(delta)
C
      do itarg = 1,ntarg
c
         x = (targ(1,itarg) - center(1))*dsq
         y = (targ(2,itarg) - center(2))*dsq

         xp(0)=1.0d0
         yp(0)=1.0d0
         do i=1,nlocal
            xp(i)=x*xp(i-1)
            yp(i)=y*yp(i-1)
         enddo
c
         do ind=1,nd
            dd=0
            do j1=0,nlocal
               d2=0
               do j2=0,nlocal
                  d2=d2+local(j2,j1,ind)*xp(j2)
               enddo
               dd=dd+d2*yp(j1)
            enddo
            pot(ind,itarg) = pot(ind,itarg)+dd
         enddo
      enddo

      return
      end
C
C
C
C
      subroutine g2dlevalg_vec(nd,delta,center,nlocal,local,
     1              targ,ntarg,pot,grad)
C
C     This subroutine evaluates the local expansion about
C     CENT at location TARG.
C
C     INPUT
C
c     nd            = vector length (for vector input)
C     delta         = Gaussian variance
C     center        = center of the expansion
C     nlocal        = number of terms in local expansions
C     local         = local Taylor expansion
c
C     targ          = target location
C
C     OUTPUT:
C
C     pot           = potential (or vectorized potentials) incremented
C
      implicit real*8 (a-h,o-z)
      real *8 local(0:nlocal,0:nlocal,nd)
      real *8 center(2),targ(2,ntarg)
      real *8 pot(nd,ntarg),xp(0:100),yp(0:100)
      real *8 grad(nd,2,ntarg)
      real *8, allocatable :: xpx(:),ypy(:)
C
      allocate(xpx(0:nlocal))
      allocate(ypy(0:nlocal))
C
      dsq = 1.0D0/dsqrt(delta)
C
      do itarg = 1,ntarg
c
         x = (targ(1,itarg) - center(1))*dsq
         y = (targ(2,itarg) - center(2))*dsq

         xp(0)=1.0d0
         yp(0)=1.0d0
         do i=1,nlocal
            xp(i)=x*xp(i-1)
            yp(i)=y*yp(i-1)
         enddo
c
         xpx(0) = 0.0d0
         ypy(0) = 0.0d0
         do i = 1,nlocal
            xpx(i) = i*xp(i-1)*dsq
            ypy(i) = i*yp(i-1)*dsq
         enddo
c
         do ind=1,nd
            dd=0
            g1=0
            g2=0
            do j1=0,nlocal
               d2=0
               d2x=0
               do j2=0,nlocal
                  d2=d2+local(j2,j1,ind)*xp(j2)
                  d2x=d2x+local(j2,j1,ind)*xpx(j2)
               enddo
               dd=dd+d2*yp(j1)
               g1=g1+d2x*yp(j1)
               g2=g2+d2*ypy(j1)
            enddo
            pot(ind,itarg) = pot(ind,itarg)+dd
            grad(ind,1,itarg) = grad(ind,1,itarg)+g1
            grad(ind,2,itarg) = grad(ind,2,itarg)+g2
         enddo
      enddo

      return
      end
C
C
C
C
      subroutine g2dlevalh_vec(nd,delta,center,nlocal,local,
     1              targ,ntarg,pot,grad,hess)
C
C     This subroutine evaluates the local expansion about
C     CENT at location TARG.
C
C     INPUT
C
c     nd            = vector length (for vector input)
C     delta         = Gaussian variance
C     center        = center of the expansion
C     nlocal        = number of terms in local expansions
C     local         = local Taylor expansion
c
C     targ          = target location
C
C     OUTPUT:
C
C     pot           = potential (or vectorized potentials) incremented
C
      implicit real*8 (a-h,o-z)
      real *8 local(0:nlocal,0:nlocal,nd)
      real *8 center(2),targ(2,ntarg)
      real *8 pot(nd,ntarg),xp(0:100),yp(0:100)
      real *8 grad(nd,2,ntarg)
      real *8 hess(nd,3,ntarg)
      real *8, allocatable :: xpx(:),ypy(:)
      real *8, allocatable :: xpxx(:),ypyy(:)
C
      allocate(xpx(0:nlocal))
      allocate(ypy(0:nlocal))
      allocate(xpxx(0:nlocal))
      allocate(ypyy(0:nlocal))
C
      dsq = 1.0D0/dsqrt(delta)
C
      do itarg = 1,ntarg
c
         x = (targ(1,itarg) - center(1))*dsq
         y = (targ(2,itarg) - center(2))*dsq

         xp(0)=1.0d0
         yp(0)=1.0d0
         do i=1,nlocal
            xp(i)=x*xp(i-1)
            yp(i)=y*yp(i-1)
         enddo
c
         xpx(0) = 0.0d0
         ypy(0) = 0.0d0
         do i = 1,nlocal
            xpx(i) = i*xp(i-1)*dsq
            ypy(i) = i*yp(i-1)*dsq
         enddo
c
         xpxx(0) = 0.0d0
         ypyy(0) = 0.0d0
         xpxx(1) = 0.0d0
         ypyy(1) = 0.0d0
         do i = 2,nlocal
            xpxx(i) = i*(i-1)*xp(i-2)*dsq*dsq
            ypyy(i) = i*(i-1)*yp(i-2)*dsq*dsq
         enddo

         do ind=1,nd
            dd=0
            g1=0
            g2=0
            h1=0
            h2=0
            h3=0
            do j1=0,nlocal
               d2=0
               d2x=0
               d2xx=0
               do j2=0,nlocal
                  d2=d2+local(j2,j1,ind)*xp(j2)
                  d2x=d2x+local(j2,j1,ind)*xpx(j2)
                  d2xx=d2xx+local(j2,j1,ind)*xpxx(j2)
               enddo
               dd=dd+d2*yp(j1)
               g1=g1+d2x*yp(j1)
               g2=g2+d2*ypy(j1)
               
               h1=h1+d2xx*yp(j1)
               h2=h2+d2x*ypy(j1)
               h3=h3+d2*ypyy(j1)
            enddo
            pot(ind,itarg) = pot(ind,itarg)+dd
            grad(ind,1,itarg) = grad(ind,1,itarg)+g1
            grad(ind,2,itarg) = grad(ind,2,itarg)+g2
            hess(ind,1,itarg) = hess(ind,1,itarg)+h1
            hess(ind,2,itarg) = hess(ind,2,itarg)+h2
            hess(ind,3,itarg) = hess(ind,3,itarg)+h3
         enddo
      enddo

      return
      end
C
C***********************************************************************
C
C 2D translations (h2l, h2sx_h2s_real, h2sx_h2s)
C
c***********************************************************************      
      subroutine g2dherm2local_vec(nd,nlocal,ntermsh,h2l,
     1    hexp,local)
C
C     This subroutine converts the Hermite expansion to 
C     local Taylor expansion
C
C     INPUT
C
c     nd            = vector length (parallel)
C     nlocal        = number of terms in local exp
C     ntermsh       = number of terms in Hermite exp
C     h2l           = matrix converting the Hermite expansion to the local expansion
C     hexp          = Hermite expansions
C      
C     OUTPUT:
C
C     local         = local Taylor expansions
C      
      implicit real*8 (a-h,o-z)
      real *8 h2l(0:ntermsh,0:nlocal)
      real *8 hexp(0:ntermsh,0:ntermsh,nd)
      real *8 local(0:nlocal,0:nlocal,nd)
      
      real *8, allocatable :: lxhy(:,:)
c      
c
      allocate(lxhy(0:nlocal,0:ntermsh))

      do ind=1,nd

cccc         do k1=0,nlocal
cccc            do j2=0,ntermsh
cccc               cd=0
cccc               do j1=0,ntermsh
cccc                  cd=cd+h2l(j1,k1)*hexp(j1,j2,ind)
cccc               enddo
cccc               lxhy(k1,j2)=cd
cccc            enddo
cccc         enddo

         do k1=0,nlocal,2
            do j2=0,ntermsh
               cd=0
               do j1=0,ntermsh,2
                  cd=cd+h2l(j1,k1)*hexp(j1,j2,ind)
               enddo
               lxhy(k1,j2)=cd
            enddo
         enddo

         do k1=1,nlocal,2
            do j2=0,ntermsh
               cd=0
               do j1=1,ntermsh,2
                  cd=cd+h2l(j1,k1)*hexp(j1,j2,ind)
               enddo
               lxhy(k1,j2)=cd
            enddo
         enddo

cccc         do k1=0,nlocal
cccc            do k2=0,nlocal
cccc               cd=0
cccc               do j2=0,ntermsh
cccc                  cd=cd+h2l(j2,k1)*lxhy(k2,j2)
cccc               enddo
cccc               local(k2,k1,ind)=local(k2,k1,ind)+cd
cccc            enddo
cccc         enddo

         do k1=0,nlocal,2
            do k2=0,nlocal
               cd=0
               do j2=0,ntermsh,2
                  cd=cd+h2l(j2,k1)*lxhy(k2,j2)
               enddo
               local(k2,k1,ind)=local(k2,k1,ind)+cd
            enddo
         enddo
         
         do k1=1,nlocal,2
            do k2=0,nlocal
               cd=0
               do j2=1,ntermsh,2
                  cd=cd+h2l(j2,k1)*lxhy(k2,j2)
               enddo
               local(k2,k1,ind)=local(k2,k1,ind)+cd
            enddo
         enddo
         
      enddo
         
      return
      end
C
C
C
C
C
      subroutine g2dh2sx_h2s_real_vec(nd,ntermsh,npw,nsoe,h2s,h2x,
     1    hexp,sxall,soeall)
C
C     This subroutine converts the Hermite expansion to the four hybrid
C     SOE/X expansions  and four SOE expansions
C
C     INPUT
C
c     nd            = vector length (parallel)
C     delta         = Gaussian variance
C     center        = center of the expansion
C     ntermsh       = number of terms in Hermite exp
C     nsoe          = number of exponentials in SOE 
C     h2s           = matrix converting the Hermite expansion to the SOE+ expansion (1D)
C     h2x           = matrix converting the Hermite expansion to the X expansion (1D)
C     hexp          = Hermite expansion
C      
C     OUTPUT:
C
C     sxall         = hybrid SOE/X expansions in order px,mx,xp,xm
C     soeall        = all four (unrolled) SOE expansions are incremented
c                     in the order: soepp,soepm,soemp,soemm
C      
      implicit real*8 (a-h,o-z)
      complex *16 h2s(0:ntermsh,nsoe)
      real *8 h2x(0:ntermsh,-npw:npw)
      real *8 hexp(0:ntermsh,0:ntermsh,nd)
      complex *16 sxall(4*(2*npw+1)*nsoe/2,nd)
      complex *16 soeall(4*nsoe*nsoe/2,nd)

      complex *16 pe,me,po,mo
      complex *16 pxe,mxe,xpe,xme,pxo,mxo,xpo,xmo
      complex *16, allocatable :: ph(:,:),mh(:,:),hp(:,:),hm(:,:)
      complex *16 px(0:3),mx(0:3),xp(0:3),xm(0:3)
      complex *16 eye
      
      eye = dcmplx(0,1)
c
      allocate(ph(0:ntermsh,nsoe/2))
      allocate(mh(0:ntermsh,nsoe/2))
      allocate(hp(0:ntermsh,nsoe))
      allocate(hm(0:ntermsh,nsoe))

c      
      dsq = 1.0D0/dsqrt(delta)
C
      nsxexp = (2*npw+1)*nsoe/2
      nexp = nsoe*nsoe/2
      
      do ind=1,nd
         do j1=1,nsoe/2
            do j2=0,ntermsh
               pe=0
               me=0
               po=0
               mo=0
               
               do k=0,ntermsh,2
                  pe=pe+h2s(k,j1)*hexp(k,j2,ind)
                  me=me+h2s(k,j1)*hexp(j2,k,ind)                  
               enddo
               do k=1,ntermsh,2
                  po=po+h2s(k,j1)*hexp(k,j2,ind)
                  mo=mo+h2s(k,j1)*hexp(j2,k,ind)                  
               enddo
               ph(j2,j1)=pe+po
               mh(j2,j1)=pe-po
               hp(j2,j1)=me+mo
               hm(j2,j1)=me-mo
            enddo
         enddo

         
         do j1=1,nsoe/2
            do j2=0,ntermsh
               hp(j2,j1+nsoe/2)=conjg(hp(j2,j1))
               hm(j2,j1+nsoe/2)=conjg(hm(j2,j1))
            enddo
         enddo

         j=0
         do j1=1,nsoe/2
            do j2=1,nsoe
               j=j+1
               pe=0
               me=0
               po=0
               mo=0
               
               do k=0,ntermsh,2
                  pe=pe+h2s(k,j1)*hp(k,j2)
                  me=me+h2s(k,j1)*hm(k,j2)
               enddo
               do k=1,ntermsh,2
                  po=po+h2s(k,j1)*hp(k,j2)
                  mo=mo+h2s(k,j1)*hm(k,j2)                  
               enddo
               soeall(j,ind) = soeall(j,ind) + pe+po
               soeall(nexp+j,ind) = soeall(nexp+j,ind)+me+mo
               soeall(2*nexp+j,ind) = soeall(2*nexp+j,ind)+pe-po
               soeall(3*nexp+j,ind) = soeall(3*nexp+j,ind)+me-mo
            enddo
         enddo
         
         j = 0
         do j1=1,nsoe/2
            do j2=-npw,0
               j = j+1

               do k=0,3
                  px(k)=0
                  mx(k)=0
                  xp(k)=0
                  xm(k)=0
               enddo
               
c     use real h2x translation matrix
c     performance actually degraded due to slower memory access 
               do k=0,ntermsh,4
                  px(0)=px(0)+h2x(k,j2)*ph(k,j1)
                  mx(0)=mx(0)+h2x(k,j2)*mh(k,j1)
                  xp(0)=xp(0)+h2x(k,j2)*hp(k,j1)
                  xm(0)=xm(0)+h2x(k,j2)*hm(k,j1)
               enddo
               do k=1,ntermsh,4
                  px(1)=px(1)+h2x(k,j2)*ph(k,j1)
                  mx(1)=mx(1)+h2x(k,j2)*mh(k,j1)
                  xp(1)=xp(1)+h2x(k,j2)*hp(k,j1)
                  xm(1)=xm(1)+h2x(k,j2)*hm(k,j1)
               enddo
               do k=2,ntermsh,4
                  px(2)=px(2)+h2x(k,j2)*ph(k,j1)
                  mx(2)=mx(2)+h2x(k,j2)*mh(k,j1)
                  xp(2)=xp(2)+h2x(k,j2)*hp(k,j1)
                  xm(2)=xm(2)+h2x(k,j2)*hm(k,j1)
               enddo
               do k=3,ntermsh,4
                  px(3)=px(3)+h2x(k,j2)*ph(k,j1)
                  mx(3)=mx(3)+h2x(k,j2)*mh(k,j1)
                  xp(3)=xp(3)+h2x(k,j2)*hp(k,j1)
                  xm(3)=xm(3)+h2x(k,j2)*hm(k,j1)
               enddo

               pxe=px(0)-px(2)
               mxe=mx(0)-mx(2)
               xpe=xp(0)-xp(2)
               xme=xm(0)-xm(2)

               pxo=eye*(px(1)-px(3))
               mxo=eye*(mx(1)-mx(3))
               xpo=eye*(xp(1)-xp(3))
               xmo=eye*(xm(1)-xm(3))
               
               sxall(         j,ind) = sxall(         j,ind)
     1                               + pxe+pxo
               sxall(  nsxexp+j,ind) = sxall(  nsxexp+j,ind)
     1                               + mxe+mxo
               sxall(2*nsxexp+j,ind) = sxall(2*nsxexp+j,ind)
     1                               + xpe+xpo
               sxall(3*nsxexp+j,ind) = sxall(3*nsxexp+j,ind)
     1                               + xme+xm0

               if (j2 .lt. 0) then
                  jp=j-2*j2
                  
                  sxall(         jp,ind) = sxall(         jp,ind)
     1                                   + pxe-pxo
                  sxall(  nsxexp+jp,ind) = sxall(  nsxexp+jp,ind)
     1                                   + mxe-mxo
                  sxall(2*nsxexp+jp,ind) = sxall(2*nsxexp+jp,ind)
     1                                   + xpe-xpo
                  sxall(3*nsxexp+jp,ind) = sxall(3*nsxexp+jp,ind)
     1                                   + xme-xmo
               endif                  
            enddo

            j=j+npw
         enddo
      enddo
         
      return
      end
C
C
C
C
C
      subroutine g2dh2sx_h2s_vec(nd,ntermsh,npw,nsoe,h2s,h2x,
     1    hexp,sxall,soeall)
C
C     This subroutine converts the Hermite expansion to the four hybrid
C     SOE/X expansions  and four SOE expansions
C
C     INPUT
C
c     nd            = vector length (parallel)
C     delta         = Gaussian variance
C     center        = center of the expansion
C     ntermsh       = number of terms in Hermite exp
C     nsoe          = number of exponentials in SOE 
C     h2s           = matrix converting the Hermite expansion to the SOE+ expansion (1D)
C     h2x           = matrix converting the Hermite expansion to the X expansion (1D)
C     hexp          = Hermite expansion
C      
C     OUTPUT:
C
C     sxall         = hybrid SOE/X expansions in order px,mx,xp,xm
C     soeall        = all four (unrolled) SOE expansions are incremented
c                     in the order: soepp,soepm,soemp,soemm
C      
      implicit real*8 (a-h,o-z)
      complex *16 h2s(0:ntermsh,nsoe)
      complex *16 h2x(0:ntermsh,-npw:npw)
      real *8 hexp(0:ntermsh,0:ntermsh,nd)
      complex *16 sxall(4*(2*npw+1)*nsoe/2,nd)
      complex *16 soeall(4*nsoe*nsoe/2,nd)

      complex *16 pe,me,po,mo
      complex *16 px,mx,xp,xm
      complex *16, allocatable :: ph(:,:),mh(:,:),hp(:,:),hm(:,:)
      complex *16, allocatable :: pxe(:),pxo(:),mxe(:),mxo(:)
      complex *16, allocatable :: xpe(:),xpo(:),xme(:),xmo(:)
c
      allocate(ph(0:ntermsh,nsoe/2))
      allocate(mh(0:ntermsh,nsoe/2))
      allocate(hp(0:ntermsh,nsoe))
      allocate(hm(0:ntermsh,nsoe))

      allocate(pxe(-npw:0))
      allocate(pxo(-npw:0))
      allocate(mxe(-npw:0))
      allocate(mxo(-npw:0))
      allocate(xpe(-npw:0))
      allocate(xpo(-npw:0))
      allocate(xme(-npw:0))
      allocate(xmo(-npw:0))
c      
      dsq = 1.0D0/dsqrt(delta)
C
      nsxexp = (2*npw+1)*nsoe/2
      nexp = nsoe*nsoe/2
      
      do ind=1,nd
         do j1=1,nsoe/2
            do j2=0,ntermsh
               pe=0
               me=0
               po=0
               mo=0
               
               do k=0,ntermsh,2
                  pe=pe+h2s(k,j1)*hexp(k,j2,ind)
                  me=me+h2s(k,j1)*hexp(j2,k,ind)                  
               enddo
               do k=1,ntermsh,2
                  po=po+h2s(k,j1)*hexp(k,j2,ind)
                  mo=mo+h2s(k,j1)*hexp(j2,k,ind)                  
               enddo
               ph(j2,j1)=pe+po
               mh(j2,j1)=pe-po
               hp(j2,j1)=me+mo
               hm(j2,j1)=me-mo
            enddo
         enddo

         
         do j1=1,nsoe/2
            do j2=0,ntermsh
               hp(j2,j1+nsoe/2)=conjg(hp(j2,j1))
               hm(j2,j1+nsoe/2)=conjg(hm(j2,j1))
            enddo
         enddo

         j=0
         do j1=1,nsoe/2
            do j2=1,nsoe
               j=j+1
               pe=0
               me=0
               po=0
               mo=0
               
               do k=0,ntermsh,2
                  pe=pe+h2s(k,j1)*hp(k,j2)
                  me=me+h2s(k,j1)*hm(k,j2)
               enddo
               do k=1,ntermsh,2
                  po=po+h2s(k,j1)*hp(k,j2)
                  mo=mo+h2s(k,j1)*hm(k,j2)                  
               enddo
               soeall(j,ind) = soeall(j,ind) + pe+po
               soeall(nexp+j,ind) = soeall(nexp+j,ind)+me+mo
               soeall(2*nexp+j,ind) = soeall(2*nexp+j,ind)+pe-po
               soeall(3*nexp+j,ind) = soeall(3*nexp+j,ind)+me-mo
            enddo
         enddo
         
         j = 0
         do j1=1,nsoe/2
cccc            do j2=-npw,npw
cccc               j = j+1
cccc               px=0
cccc               mx=0
cccc               xp=0
cccc               xm=0
cccc               
cccc               do k=0,ntermsh
cccc                  px=px+h2x(k,j2)*ph(k,j1)
cccc                  mx=mx+h2x(k,j2)*mh(k,j1)
cccc                  xp=xp+h2x(k,j2)*hp(k,j1)
cccc                  xm=xm+h2x(k,j2)*hm(k,j1)
cccc               enddo
cccc
cccc               sxall(j,ind) = sxall(j,ind)+ px
cccc               sxall(nsxexp+j,ind) = sxall(nsxexp+j,ind)+ mx
cccc               sxall(2*nsxexp+j,ind) = sxall(2*nsxexp+j,ind)+ xp
cccc               sxall(3*nsxexp+j,ind) = sxall(3*nsxexp+j,ind) + xm
cccc            enddo
            do j2=-npw,0
               j = j+1
               pxe(j2)=0
               pxo(j2)=0
               mxe(j2)=0
               mxo(j2)=0
               xpe(j2)=0
               xpo(j2)=0
               xme(j2)=0
               xmo(j2)=0
c              reduce one complex multiplication using the symmetry in h2x
               do k=0,ntermsh,2
                  pxe(j2)=pxe(j2)+h2x(k,j2)*ph(k,j1)
                  mxe(j2)=mxe(j2)+h2x(k,j2)*mh(k,j1)
                  xpe(j2)=xpe(j2)+h2x(k,j2)*hp(k,j1)
                  xme(j2)=xme(j2)+h2x(k,j2)*hm(k,j1)
               enddo
               
               do k=1,ntermsh,2
                  pxo(j2)=pxo(j2)+h2x(k,j2)*ph(k,j1)
                  mxo(j2)=mxo(j2)+h2x(k,j2)*mh(k,j1)
                  xpo(j2)=xpo(j2)+h2x(k,j2)*hp(k,j1)
                  xmo(j2)=xmo(j2)+h2x(k,j2)*hm(k,j1)
               enddo

               sxall(j,ind) = sxall(j,ind)
     1             + pxe(j2)+pxo(j2)
               sxall(nsxexp+j,ind) = sxall(nsxexp+j,ind)
     1             + mxe(j2)+mxo(j2)
               sxall(2*nsxexp+j,ind) = sxall(2*nsxexp+j,ind)
     1             + xpe(j2)+xpo(j2)
               sxall(3*nsxexp+j,ind) = sxall(3*nsxexp+j,ind)
     1             + xme(j2)+xmo(j2)
            enddo
            
            do j2=1,npw
               j=j+1
               sxall(j,ind) = sxall(j,ind)
     1             + pxe(-j2)-pxo(-j2)
               sxall(nsxexp+j,ind) = sxall(nsxexp+j,ind)
     1             + mxe(-j2)-mxo(-j2)
               sxall(2*nsxexp+j,ind) = sxall(2*nsxexp+j,ind)
     1             + xpe(-j2)-xpo(-j2)
               sxall(3*nsxexp+j,ind) = sxall(3*nsxexp+j,ind)
     1             + xme(-j2)-xmo(-j2)
            enddo
         enddo
      enddo
         
      return
      end
C
C
C
C
C
C***********************************************************************
C
C 2D translations (soe2local, sx2local)
C
c***********************************************************************      
      subroutine g2dsoe2local_vec(nd,nlocal,nsoe,s2l,
     1    soeall,local)
C
C     This subroutine converts four SOE expansions to local
C     expansions
C
C     INPUT
C
c     nd            = vector length (parallel)
C     nlocal        = number of terms in local exp
C     nsoe          = number of exponentials in SOE 
C     s2l           = matrix converting SOE expansions to the local expansion
C     soeall        = SOE expansions
C      
C     OUTPUT:
C
C     local         = local Taylor expansions
C      
      implicit real*8 (a-h,o-z)
      complex *16 s2l(nsoe,0:nlocal)
      real *8 local(0:nlocal,0:nlocal,nd)
      complex *16 soeall(4*nsoe*nsoe/2,nd)
      complex *16, allocatable :: sply(:,:),smly(:,:)
      complex *16 cp,cm,cd,yp,ym
c      
c
      allocate(sply(nsoe/2,0:nlocal))
      allocate(smly(nsoe/2,0:nlocal))

      nexp = nsoe*nsoe/2
      do ind=1,nd
         do k2=0,nlocal,2
            j=0
            do j1=1,nsoe/2
               cp=0
               cm=0
               do j2=1,nsoe
                  j=j+1
                  yp=s2l(j2,k2)
                  
                  cp=cp+yp*(soeall(j,ind)+soeall(nexp+j,ind)) 
                  cm=cm+yp*(soeall(2*nexp+j,ind)+soeall(3*nexp+j,ind))
               enddo
               sply(j1,k2)=cp
               smly(j1,k2)=cm
            enddo
         enddo
         
         do k2=1,nlocal,2
            j=0
            do j1=1,nsoe/2
               cp=0
               cm=0
               do j2=1,nsoe
                  j=j+1
                  yp=s2l(j2,k2)
                  
                  cp=cp+yp*(soeall(j,ind)-soeall(nexp+j,ind)) 
                  cm=cm+yp*(soeall(2*nexp+j,ind)-soeall(3*nexp+j,ind))
               enddo
               sply(j1,k2)=cp
               smly(j1,k2)=cm
            enddo
         enddo

         do k1=0,nlocal
ccc         k2 even            
            do k2=0,nlocal,2
               cd=0
               do j1=1,nsoe/2
                  cd=cd+s2l(j1,k2)*(sply(j1,k1)+smly(j1,k1))
               enddo
               local(k2,k1,ind)=local(k2,k1,ind)+dble(cd)/2
            enddo
ccc         k2 odd            
            do k2=1,nlocal,2
               cd=0
               do j1=1,nsoe/2
                  cd=cd+s2l(j1,k2)*(sply(j1,k1)-smly(j1,k1))
               enddo
               local(k2,k1,ind)=local(k2,k1,ind)+dble(cd)/2
            enddo            
         enddo
      enddo
         
      return
      end
C
C
C
C
C
      subroutine g2dsx2local_vec(nd,nlocal,npw,nsoe,x2l,s2l,
     1    sxall,local)
C
C     This subroutine converts four hybrid SOE/X expansions to 
C     local expansions
C
C     INPUT
C
c     nd            = vector length (parallel)
C     nlocal        = number of terms in local exp
C     npw           = number of terms in planewave exp
C     nsoe          = number of exponentials in SOE 
C     x2l           = matrix converting the planewave expansion to the local expansion
C     s2l     = matrices converting the SOE expansions to the local expansion
C     sxall         = SOE/X expansions
C      
C     OUTPUT:
C
C     local         = local Taylor expansions
C      
      implicit real*8 (a-h,o-z)
      complex *16 x2l(-npw:npw,0:nlocal)
      complex *16 s2l(nsoe,0:nlocal)
      complex *16 sxall(2*(2*npw+1)*nsoe,nd)
      real *8 local(0:nlocal,0:nlocal,nd)
      complex *16, allocatable :: sply(:,:),smly(:,:)
      complex *16, allocatable :: lxsp(:,:),lxsm(:,:)
      complex *16 cd,cd1,cd2,cd3,cd4
c      
c
      allocate(sply(nsoe/2,0:nlocal))
      allocate(smly(nsoe/2,0:nlocal))
      allocate(lxsp(nsoe/2,0:nlocal))
      allocate(lxsm(nsoe/2,0:nlocal))

      np0= npw+1

      nexp = (2*npw+1)*nsoe/2
      do ind=1,nd

cccc         do k=0,nlocal
cccc            j=0
cccc            do j1=1,nsoe/2
cccc               cd1=0
cccc               cd2=0
cccc               cd3=0
cccc               cd4=0
cccc               do j2=-npw,npw
cccc                  j=j+1
cccc                  
cccc                  cd1=cd1+sxall(j,ind)*x2l(j2,k)  
cccc                  cd2=cd2+sxall(nexp+j,ind)*x2l(j2,k)
cccc                  cd3=cd3+sxall(2*nexp+j,ind)*x2l(j2,k)
cccc                  cd4=cd4+sxall(3*nexp+j,ind)*x2l(j2,k)
cccc               enddo
cccc               sply(j1,k)=cd1
cccc               smly(j1,k)=cd2
cccc               lxsp(j1,k)=cd3
cccc               lxsm(j1,k)=cd4
cccc            enddo
cccc         enddo
         do k=0,nlocal,2
            j=0
            do j1=1,nsoe/2
               cd1=0
               cd2=0
               cd3=0
               cd4=0
               if (k.eq. 0) then
                  jp=j+np0
                  cd1=cd1+sxall(jp,ind)
                  cd2=cd2+sxall(nexp+jp,ind)
                  cd3=cd3+sxall(2*nexp+jp,ind)
                  cd4=cd4+sxall(3*nexp+jp,ind)
               endif
cccc           reduce one complex multiplication using the symmetry in x2l
c              but it could be a wash due to memory access
               do j2=-npw,-1
                  j=j+1
                  jp=j-2*j2
                  cd1=cd1+( sxall(j,ind)
     1                     +sxall(jp,ind) )*x2l(j2,k)  
                  cd2=cd2+( sxall(nexp+j,ind)
     1                     +sxall(nexp+jp,ind) )*x2l(j2,k)
                  cd3=cd3+( sxall(2*nexp+j,ind)
     1                     +sxall(2*nexp+jp,ind) )*x2l(j2,k)
                  cd4=cd4+( sxall(3*nexp+j,ind)
     1                     +sxall(3*nexp+jp,ind) )*x2l(j2,k)
               enddo

               j=j+np0
               
               sply(j1,k)=cd1
               smly(j1,k)=cd2
               lxsp(j1,k)=cd3
               lxsm(j1,k)=cd4
            enddo
         enddo
         
         do k=1,nlocal,2
            j=0
            do j1=1,nsoe/2
               cd1=0
               cd2=0
               cd3=0
               cd4=0
               do j2=-npw,-1
                  j=j+1
                  jp=j-2*j2
                  cd1=cd1+( sxall(j,ind)
     1                     -sxall(jp,ind) )*x2l(j2,k)  
                  cd2=cd2+( sxall(nexp+j,ind)
     1                     -sxall(nexp+jp,ind) )*x2l(j2,k)
                  cd3=cd3+( sxall(2*nexp+j,ind)
     1                     -sxall(2*nexp+jp,ind) )*x2l(j2,k)
                  cd4=cd4+( sxall(3*nexp+j,ind)
     1                     -sxall(3*nexp+jp,ind) )*x2l(j2,k)
               enddo

               j=j+np0
               
               sply(j1,k)=cd1
               smly(j1,k)=cd2
               lxsp(j1,k)=cd3
               lxsm(j1,k)=cd4
            enddo
         enddo         
c
c        k1 even
c         
         do k1=0,nlocal,2
c           k2 even
            do k2=0,nlocal,2
               cd=0
               do j1=1,nsoe/2
                  cd=cd+(sply(j1,k1)+smly(j1,k1))*s2l(j1,k2)
     1                 +(lxsp(j1,k2)+lxsm(j1,k2))*s2l(j1,k1)
               enddo
               local(k2,k1,ind)=local(k2,k1,ind)+dble(cd)
            enddo
c           k2 odd            
            do k2=1,nlocal,2
               cd=0
               do j1=1,nsoe/2
                  cd=cd+(sply(j1,k1)-smly(j1,k1))*s2l(j1,k2)
     1                 +(lxsp(j1,k2)+lxsm(j1,k2))*s2l(j1,k1)
               enddo
               local(k2,k1,ind)=local(k2,k1,ind)+dble(cd)
            enddo
         enddo
c
c        k1 odd
c         
         do k1=1,nlocal,2
c           k2 even            
            do k2=0,nlocal,2
               cd=0
               do j1=1,nsoe/2
                  cd=cd+(sply(j1,k1)+smly(j1,k1))*s2l(j1,k2)
     1                 +(lxsp(j1,k2)-lxsm(j1,k2))*s2l(j1,k1)
               enddo
               local(k2,k1,ind)=local(k2,k1,ind)+dble(cd)
            enddo
c           k2 odd
            do k2=1,nlocal,2
               cd=0
               do j1=1,nsoe/2
                  cd=cd+(sply(j1,k1)-smly(j1,k1))*s2l(j1,k2)
     1                 +(lxsp(j1,k2)-lxsm(j1,k2))*s2l(j1,k1)
               enddo
               local(k2,k1,ind)=local(k2,k1,ind)+dble(cd)
            enddo            
         enddo
      enddo
         
      return
      end
C
C
C
C
C
C************************************************************************
C
C 1d translation matrices
C
C**********************************************************************
      subroutine g2dh2lmat(nlocal,ntermsh,h2l)
C
C     This subroutine returns the matrix converting the 
C     Hermite expansion to the local expansion
C     
C
C     INPUT
C
C     nlocal       = number of terms in local exp
C     ntermsh       = number of terms in the Hermite exp
C
C     OUTPUT:
C
C     h2l        = translation matrix
C
      implicit real*8 (a-h,o-z)
      real *8 h2l(0:ntermsh,0:nlocal)
      real *8 fac(0:100),hexpx(0:100)
      real *8, allocatable :: sqc(:,:)

      ntot = ntermsh+nlocal
      allocate(sqc(0:ntot,0:ntot))
      call getsqrtbinomialcoeffs(ntot,sqc)
      
      fac(0)=1.0d0
      do i=1,100
         fac(i)=fac(i-1)*sqrt(i*1.0d0)
      enddo

      x=0
      hexpx(0)=1.0d0
      hexpx(1)=2*x
      
      do i=1,99
         hexpx(i+1)=2.0d0*(x/dsqrt((i+1)*1.0d0)*hexpx(i)
     1          -dsqrt(1.0d0*i/(i+1))*hexpx(i-1))
      enddo

      do i=0,ntermsh
         do j=0,nlocal
            h2l(i,j)=(-1)**j*hexpx(i+j)*sqc(i+j,i)/fac(j)
         enddo
      enddo
      
      return
      end
C
C
C
C
      subroutine g2dh2smat(ntermsh,nsoe,ws,ts,h2s)
C
C     This subroutine returns a matrix converting the Hermite expansion
C     to the SOE+ expansions
C
C     INPUT
C
C     ntermsh       = number of terms in Hermite exp
C     nsoe          = number of exponentials in SOE 
c
C     ws,ts         = SOE weights and nodes
C
C     OUTPUT:
C
C     h2s           = translation matrix
C
      implicit real*8 (a-h,o-z)
      complex *16 ws(*), ts(*)
      complex *16 h2s(0:ntermsh,nsoe)
      real *8 fac(0:100)

      fac(0)=1.0d0
      do i=1,100
         fac(i)=fac(i-1)/sqrt(dble(i))
      enddo

      
      do i=1,nsoe
         do j=0,ntermsh
            h2s(j,i)=ws(i)*ts(i)**j*fac(j)
         enddo
      enddo

      return
      end
C
C
C
C
      subroutine g2drh2xmat(ntermsh,npw,ws,ts,rh2x)
C
C     This subroutine returns the real matrix converting the Hermite expansion
C     to the planewave expansions, i.e., i**j factor is ignored.
C
C     INPUT
C
C     ntermsh       = number of terms in Hermite exp
C     npw          = number of exponentials in the planewave exp
c
C     ws,ts         = planewave exp weights and nodes
C
C     OUTPUT:
C
C     rh2x       = real translation matrix, i.e., i^j factor is left out.
C
      implicit real*8 (a-h,o-z)
      real *8 ws(-npw:npw), ts(-npw:npw)
      real *8 rh2x(0:ntermsh,-npw:npw)
      real *8 fac(0:100)

      fac(0)=1.0d0
      do i=1,100
         fac(i)=fac(i-1)/sqrt(dble(i))
      enddo
      
      do i=-npw,npw
         do j=0,ntermsh
            rh2x(j,i)=ws(i)*(-ts(i))**j*fac(j)
         enddo
      enddo

      return
      end
C
C
C
C
      subroutine g2dh2xmat(ntermsh,npw,ws,ts,h2x)
C
C     This subroutine returns the matrix converting the Hermite expansion
C     to the planewave expansions
C
C     INPUT
C
C     ntermsh       = number of terms in Hermite exp
C     npw          = number of exponentials in the planewave exp
c
C     ws,ts         = planewave exp weights and nodes
C
C     OUTPUT:
C
C     h2x       = translation matrices
C
      implicit real*8 (a-h,o-z)
      real *8 ws(-npw:npw), ts(-npw:npw)
      complex *16 h2x(0:ntermsh,-npw:npw)
      real *8 fac(0:100)
      complex *16 eye
      
      eye = dcmplx(0,1)

      fac(0)=1.0d0
      do i=1,100
         fac(i)=fac(i-1)/sqrt(dble(i))
      enddo
      
      do i=-npw,npw
         do j=0,ntermsh
            h2x(j,i)=ws(i)*(-eye*ts(i))**j*fac(j)
         enddo
      enddo

      return
      end
C
C
C
C      
      subroutine g2ds2lmat(nlocal,nsoe,ws,ts,s2l)
C
C     This subroutine returns the matrices converting the SOE+ expansions
C     to the local expansion
C
C     INPUT
C
C     nlocal       = number of terms in local exp
C     nsoe          = number of exponentials in SOE 
C     ws,ts         = SOE weights and nodes
C
C     OUTPUT:
C
C     s2l        = translation matrix
C
      implicit real*8 (a-h,o-z)
      complex *16 ws(*), ts(*)
      complex *16 s2l(nsoe,0:nlocal)
      real *8 fac(0:100)

      fac(0)=1.0d0
      do i=1,100
         fac(i)=fac(i-1)*i
      enddo
      
      do i=1,nsoe
         do j=0,nlocal
            s2l(i,j)=(-ts(i))**j/fac(j)
         enddo
      enddo

      return
      end
C
C
C
C
      subroutine g2dx2lmat(nlocal,npw,ws,ts,x2l)
C
C     This subroutine returns the matrices converting the planewave expansions
C     to the local expansion
C
C     INPUT
C
C     nlocal       = number of terms in local exp
C     npw          = number of terms in the planvewave exp 
c
C     ws,ts         = planewave exp weights and nodes
C
C     OUTPUT:
C
C     x2l        = translation matrix
C
      implicit real*8 (a-h,o-z)
      real *8 ws(-npw:npw), ts(-npw:npw)
      complex *16 x2l(-npw:npw,0:nlocal)
      complex *16 eye
      real *8 fac(0:100)
      
      eye = dcmplx(0,1)
      
      fac(0)=1.0d0
      do i=1,100
         fac(i)=fac(i-1)*i
      enddo
      
      do j=0,nlocal
         do i=-npw,npw
            x2l(i,j)=(eye*ts(i))**j/fac(j)
         enddo
      enddo

      return
      end
C
C
C
C***********************************************************************
      subroutine g2dhermzero_vec(nd,hexp,ntermsh)
      implicit none
C***********************************************************************
c
c     This subroutine sets a vector multipole expansion to zero.
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :   vector length (number of mpole expansions)
c     ntermsh :   order of multipole expansion
C---------------------------------------------------------------------
c     OUTPUT:
c
c     hexp  :   coeffs for the expansion set to zero.
C---------------------------------------------------------------------
      integer n,ntermsh,nd,ii
      real *8 hexp((ntermsh+1)*(ntermsh+1),nd)
c
      do ii=1,nd
      do n=1,(ntermsh+1)*(ntermsh+1)
         hexp(n,ii)=0.0d0
      enddo
      enddo
      return
      end
c      
c      
c      
c      
C***********************************************************************
      subroutine g2dlocalzero_vec(nd,local,nlocal)
      implicit none
C***********************************************************************
c
c     This subroutine sets a vector local expansion to zero.
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :   vector length (number of mpole expansions)
c     nlocal :   order of local expansion
C---------------------------------------------------------------------
c     OUTPUT:
c
c     local  :   coeffs for the expansion set to zero.
C---------------------------------------------------------------------
      integer n,nlocal,nd,ii
      real *8 local((nlocal+1)*(nlocal+1),nd)
c
      do ii=1,nd
      do n=1,(nlocal+1)*(nlocal+1)
         local(n,ii)=0.0d0
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
      subroutine getsqrtbinomialcoeffs(n,dc)
C***********************************************************************
c
c     computes the sqrt of the binomial coefficients, from fmmcommon.f
c-----------------------------------------------------------------------
      implicit none
      integer i,j,n
      real *8 dc(0:n,0:n)
      real *8, allocatable :: d(:,:)
      allocate(d(0:n,0:n))

      do i=0,n
        do j=0,n
          d(i,j) = 0
          dc(i,j) = 0
        enddo
      enddo

      do i=0,n
        d(i,0) = 1
        dc(i,0) = 1
      enddo
      do i=1,n
        d(i,i) = 1
        dc(i,i) = 1
        do j=i+1,n
          d(j,i) =d(j-1,i)+d(j-1,i-1)
          dc(j,i) = sqrt(d(j,i))
        enddo
      enddo

      return
      end


