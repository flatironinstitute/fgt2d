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
c     g2dformpwc_vec : computes the PW expansion due to charges
c                      type 1 NUFFT is called once
c      
c     g2dformpwd_vec : computes the PW expansion due to dipoles
c                      type 1 NUFFT is called three times using finufft3d1many
c      
c     g2dformpwcd_vec : computes the PW expansion due to charges and dipoles
c                      type 1 NUFFT is called four times using finufft3d1many
c      
c     g2dpwevalp_vec : evaluates the PW expansion
c                      potential only, type 2 NUFFT is called once
c                      using finufft3d2many
c      
c     g2dpwevalg_vec : evaluates the PW expansion
c                      potential + gradient
c                      4 type 2 NUFFTs are called using finufft3d2many 
c      
c     g2dpwevalh_vec : evaluates the PW expansion
c                      potential + gradient + hessian
c                      10 type 2 NUFFTs are called using finufft3d2many
c      
c     get_pwnodes : returns PW weights and nodes, midpoint rule is used
c                   so the number is always an even number!
c
c     pw_translation_matrices : returns precomputed translation matrices
c                   for PW mp to loc translations
c
c     g2dshiftpw_vec : translates PW expansion via precomputed translation 
c                      matrices
c
c     g2dcopypwexp_vec : copy PW expansion
c
c
c
c
c     Note: the plane wave expansion uses different convention here.
c
c     In the directory "planewave", the plane wave expansion is stored
c     in x, y in sequal, but only half of x are stored
c
c     Here, only half of y are stored, for the convenience of zero
c     padding for NUFFT calls. The order is still x, then y.
c
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
      do ind=1,nd
         do k=0,nterms
            do j=0,nterms
               ffexp(j,k,ind)=0
            enddo
         enddo
      enddo
      
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
                  ffexp(j,k,ind) = ffexp(j,k,ind)+ ytmp*xp(j)
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
      do ind=1,nd
         do k=0,nterms
            do j=0,nterms
               ffexp(j,k,ind)=0
            enddo
         enddo
      enddo
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
                  ffexp(j,k,ind) = ffexp(j,k,ind)+dxhx(j)*d1+hx(j)*d2
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
      do ind=1,nd
         do k=0,nterms
            do j=0,nterms
               ffexp(j,k,ind)=0
            enddo
         enddo
      enddo

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
                  ffexp(j,k,ind) = ffexp(j,k,ind)+
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
            r1 = rnormal(1,i)*dipstr(ind,i)
            r2 = rnormal(2,i)*dipstr(ind,i)
            do j1=0,nlocal
               ry1 = hexpy(j1)*r1
               ry2 = dhexpy(j1)*r2
               do j2=0,nlocal
                  local(j2,j1,ind) = local(j2,j1,ind)+
     1                dhexpx(j2)*ry1+hexpx(j2)*ry2
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
            r1 = rnormal(1,i)*dipstr(ind,i)
            r2 = rnormal(2,i)*dipstr(ind,i)
            
            do j1=0,nlocal
               ytmp=chg*hexpy(j1)
               dy1 = hexpy(j1)*r1
               dy2 = dhexpy(j1)*r2
               do j2=0,nlocal
                  local(j2,j1,ind) = local(j2,j1,ind)+
     1                dhexpx(j2)*dy1+hexpx(j2)*dy2+ytmp*hexpx(j2)
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
      real *8 pot(nd,ntarg),xp(0:200),yp(0:200)
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
      real *8 pot(nd,ntarg)
      real *8 grad(nd,2,ntarg)
      real *8, allocatable :: xp(:),yp(:)
      real *8, allocatable :: xpx(:),ypy(:)
C
      allocate(xp(0:nlocal))
      allocate(yp(0:nlocal))
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
      real *8 pot(nd,ntarg),xp(0:200),yp(0:200)
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
c*********************************************************************
C
C form PW expansions (charge, dipole, charge & dipole)
C
C*********************************************************************
      subroutine g2dformpwc_vec(nd,delta,eps,sources,ns,charge,
     1            cent,npw,ws,ts,nexp,wnufft,ffexp)
C
C     This subroutine computes the PW expansion about
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
C     npw    = number of terms in PW expansion
C
C     OUTPUT:
C
C     ffexp     = PW expansion
C---------------------------------------------------------------------
C
      implicit real *8 (a-h,o-z)
      integer ns,npw
      real *8 cent(2),sources(2,ns),charge(nd,ns)
      real *8 ws(-npw/2:npw/2-1), ts(-npw/2:npw/2-1)
      real *8 wnufft(nexp)
      complex *16 ffexp(nexp,nd)
      real *8 x,y,chg
      real *8, allocatable ::  xj(:),yj(:)
      complex *16, allocatable :: cj(:,:),fk(:,:,:),wj(:)
c     this (since unallocated) used to pass a NULL ptr to FINUFFT...
      integer*8, allocatable :: null
      integer*8 ns8, npw8
      complex *16 eye
      eye = dcmplx(0,1)

      npw2=npw/2
      
      allocate(xj(ns))
      allocate(yj(ns))
      allocate(cj(ns,nd))
      allocate(wj(ns))
      allocate(fk(-npw2:npw2-1,-npw2:npw2-1,nd))
c
      dsq = 2*ts(0)/dsqrt(delta)
C
C     accumulate expansion due to each source.
C
      do j=1,ns
         xj(j) = (sources(1,j) - cent(1))*dsq
         yj(j) = (sources(2,j) - cent(2))*dsq
         wj(j) = exp(-0.5d0*eye*(xj(j)+yj(j)))
      enddo

      do ind = 1,nd
         do j=1,ns
            cj(j,ind) = charge(ind,j)*wj(j)
         enddo
      enddo
      
      ns8=ns
      npw8=npw
      iflag = -1
      call finufft2d1many(nd,ns8,xj,yj,cj,iflag,eps,npw8,npw8,
     1       fk,null,ier)
cccc  print *, 'ier=', ier
      do ind=1,nd
         j=0
         do k2=-npw2,-1
            do k1=-npw2,npw2-1
               j=j+1
               ffexp(j,ind)=fk(k1,k2,ind)*wnufft(j)
            enddo
         enddo
      enddo

      return
      end
c
C
C
C
      subroutine g2dformpwd_vec(nd,delta,eps,sources,ns,rnormal,
     1    dipstr,cent,npw,ws,ts,nexp,wnufft,ffexp)
C
C     This subroutine computes the PW expansion about
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
C     npw    = number of terms in PW expansion
C
C     OUTPUT:
C
C     ffexp     = PW expansion
C---------------------------------------------------------------------
C
      implicit real *8 (a-h,o-z)
      integer ns,npw
      real *8 cent(2),sources(2,ns),dipstr(nd,ns)
      real *8 rnormal(2,ns)
      real *8 ws(-npw/2:npw/2-1), ts(-npw/2:npw/2-1)
      real *8 wnufft(nexp)
      complex *16 ffexp(nexp,nd)
      real *8 x,y,chg
      real *8, allocatable ::  xj(:),yj(:),zj(:)
      complex *16, allocatable :: cj(:,:,:),fk(:,:,:,:),wj(:)
c     this (since unallocated) used to pass a NULL ptr to FINUFFT...
      integer*8, allocatable :: null
      integer*8 ns8, npw8
      complex *16 eye,ztmp
      eye = dcmplx(0,1)

      npw2=npw/2
      
      allocate(xj(ns))
      allocate(yj(ns))
      allocate(cj(ns,2,nd))
      allocate(wj(ns))
      allocate(fk(-npw2:npw2-1,-npw2:npw2-1,2,nd))
c
      dsq0 = 1/dsqrt(delta)
      dsq = 2*ts(0)/dsqrt(delta)
C
C     accumulate expansion due to each source.
C
      do j=1,ns
         xj(j) = (sources(1,j) - cent(1))*dsq
         yj(j) = (sources(2,j) - cent(2))*dsq
         wj(j) = exp(-0.5d0*eye*(xj(j)+yj(j)))
      enddo
      
      do ind = 1,nd
         do j=1,ns
            ztmp = -eye*wj(j)*(dipstr(ind,j)*dsq0)
            cj(j,1,ind) = ztmp*rnormal(1,j)
            cj(j,2,ind) = ztmp*rnormal(2,j)
         enddo
      enddo
      
      ns8=ns
      npw8=npw
      iflag = -1
c     total number of NUFFTs
      ntrans = 2*nd
      call finufft2d1many(ntrans,ns8,xj,yj,cj,iflag,eps,
     1       npw8,npw8,fk,null,ier)
cccc  print *, 'ier=', ier
      do ind=1,nd
         j=0
         do k2=-npw2,-1
            do k1=-npw2,npw2-1
               j=j+1
               ffexp(j,ind)=wnufft(j)*(fk(k1,k2,1,ind)*ts(k1)+
     1             fk(k1,k2,2,ind)*ts(k2))
            enddo
         enddo
      enddo

      return
      end
c
C
C
C
      subroutine g2dformpwcd_vec(nd,delta,eps,sources,ns,charge,
     1    rnormal,dipstr,cent,npw,ws,ts,nexp,wnufft,ffexp)
C
C     This subroutine computes the PW expansion about
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
C     npw    = number of terms in PW expansion
C
C     OUTPUT:
C
C     ffexp     = PW expansion
C---------------------------------------------------------------------
C
      implicit real *8 (a-h,o-z)
      integer ns,npw
      real *8 cent(2),sources(2,ns),dipstr(nd,ns)
      real *8 rnormal(2,ns),charge(nd,ns)
      real *8 ws(-npw/2:npw/2-1), ts(-npw/2:npw/2-1)
      real *8 wnufft(nexp)
      complex *16 ffexp(nexp,nd)
      real *8 x,y,chg
      real *8, allocatable ::  xj(:),yj(:)
      complex *16, allocatable :: cj(:,:,:),fk(:,:,:,:),wj(:)
c     this (since unallocated) used to pass a NULL ptr to FINUFFT...
      integer*8, allocatable :: null
      integer*8 ns8, npw8
      complex *16 eye,ztmp
      eye = dcmplx(0,1)

      npw2=npw/2
      
      allocate(xj(ns))
      allocate(yj(ns))
      allocate(cj(ns,0:2,nd))
      allocate(wj(ns))
      allocate(fk(-npw2:npw2-1,-npw2:npw2-1,0:2,nd))
c
      dsq0 = 1/dsqrt(delta)
      dsq = 2*ts(0)/dsqrt(delta)
C
C     accumulate expansion due to each source.
C
      do j=1,ns
         xj(j) = (sources(1,j) - cent(1))*dsq
         yj(j) = (sources(2,j) - cent(2))*dsq
         wj(j) = exp(-0.5d0*eye*(xj(j)+yj(j)))
      enddo
      
      do ind = 1,nd
         do j=1,ns
            cj(j,0,ind) = charge(ind,j)*wj(j)
            ztmp = -eye*wj(j)*(dipstr(ind,j)*dsq0)
            cj(j,1,ind) = ztmp*rnormal(1,j)
            cj(j,2,ind) = ztmp*rnormal(2,j)
         enddo
      enddo
      
      ns8=ns
      npw8=npw
      iflag = -1
c     total number of NUFFTs
      ntrans = 3*nd
      call finufft2d1many(ntrans,ns8,xj,yj,cj,iflag,eps,
     1       npw8,npw8,fk,null,ier)
cccc  print *, 'ier=', ier
      do ind=1,nd
         j=0
         do k2=-npw2,-1
            do k1=-npw2,npw2-1
               j=j+1
               ffexp(j,ind)=wnufft(j)*(fk(k1,k2,0,ind)+
     1             fk(k1,k2,1,ind)*ts(k1)+
     1             fk(k1,k2,2,ind)*ts(k2))
            enddo
         enddo
      enddo

      return
      end
c
C
C
C      
c*********************************************************************
C
C evaluate PW expansions (potential, pot + grad, pot + grad + hess)
C
C*********************************************************************
C
C
c
C
C
C
      subroutine g2dpwevalp_vec(nd,delta,eps,center,npw,ws,ts,
     1              nexp,pwexp,targ,nt,pot)
C
C     This subroutine evaluates the plane wave 
C     expansions about CENTER at location TARG
C     potential only
C
C     INPUT
C
c     nd            = vector length (for multiple charges at same locations)
C     delta         = Gaussian variance
C     charge        = strength of sources
C     center        = center of the expansion
C     npw           = number of Fourier plane waves
C     ws,ts         = planewave weights and nodes
C     pwexp         = pw expansions 
C     targ          = target
C
C     OUTPUT:
C     pot           = potential (or vectorized potentials) incremented
C
      implicit real *8 (a-h,o-z)
      integer nd,npw,ntarg
      real *8 delta,center(2),targ(2,nt)
      real *8 pot(nd,nt)
      
      real *8 ws(-npw/2:npw/2-1),ts(-npw/2:npw/2-1)
      complex *16 pwexp(nexp,nd)
      integer*8 nt8, npw8

      integer i,ind,j1,j2,j,npw2,itarg,k
      real *8 x,y,dsq
      complex *16 eye
      complex *16 qqx,qqy,qq1,qq2
      
      complex *16 z,cd
      real *8, allocatable ::  xj(:),yj(:)
      complex *16, allocatable :: cj(:,:),fk(:,:),wj(:)
c     this (since unallocated) used to pass a NULL ptr to FINUFFT...
      integer*8, allocatable :: null
C
      eye = dcmplx(0,1)
      
      allocate(xj(nt))
      allocate(yj(nt))
      allocate(cj(nt,nd))
      allocate(wj(nt))
      allocate(fk(2*nexp,nd))
      
      dsq = 2*ts(0)/dsqrt(delta)
C
      do j=1,nt
         xj(j) = (targ(1,j) - center(1))*dsq
         yj(j) = (targ(2,j) - center(2))*dsq
         wj(j) = exp(0.5d0*eye*(xj(j)+yj(j)))
      enddo
      
      do ind = 1,nd
         do j=1,nexp
            fk(j,ind)=pwexp(j,ind)
         enddo
         do j=nexp+1,2*nexp
            fk(j,ind)=0
         enddo
      enddo

      nt8=nt
      npw8=npw
      iflag = 1
      call finufft2d2many(nd,nt8,xj,yj,cj,iflag,eps,npw8,npw8,
     1    fk,null,ier)
      
      do ind=1,nd
         do j=1,nt
            pot(ind,j)=pot(ind,j)+dreal(cj(j,ind)*wj(j))*2
         enddo
      enddo
c
      return
      end
C
C
c
C
C
C
      subroutine g2dpwevalg_vec(nd,delta,eps,center,npw,ws,ts,
     1              nexp,pwexp,targ,nt,pot,grad)
C
C     This subroutine evaluates the plane wave 
C     expansions about CENTER at location TARG
C     potential + gradient
C
C     INPUT
C
c     nd            = vector length (for multiple charges at same locations)
C     delta         = Gaussian variance
C     charge        = strength of sources
C     center        = center of the expansion
C     npw           = number of Fourier plane waves
C     ws,ts         = planewave weights and nodes
C     pwexp         = pw expansions 
C     targ          = target
C
C     OUTPUT:
C     pot           = potential (or vectorized potentials) incremented
C     grad          = gradient (or vectorized gradients) incremented
C
      implicit real *8 (a-h,o-z)
      integer nd,npw,nt
      real *8 delta,center(2),targ(2,nt)
      real *8 pot(nd,nt),grad(nd,2,nt)
      
      real *8 ws(-npw/2:npw/2-1),ts(-npw/2:npw/2-1)
      complex *16 pwexp(nexp,nd)
      integer*8 nt8, npw8

      integer i,ind,j1,j2,j,npw2,itarg,k
      real *8 x,y,dsq
      complex *16 eye
      complex *16 qqx,qqy,qq1,qq2
      
      complex *16 z,cd
      real *8, allocatable ::  xj(:),yj(:)
      complex *16, allocatable :: cj(:,:,:),fk(:,:,:),wj(:)
c     this (since unallocated) used to pass a NULL ptr to FINUFFT...
      integer*8, allocatable :: null
C
      eye = dcmplx(0,1)
      
      allocate(xj(nt))
      allocate(yj(nt))
      allocate(cj(nt,0:2,nd))
      allocate(wj(nt))
      allocate(fk(2*nexp,0:2,nd))

      dsq = 2*ts(0)/dsqrt(delta)
C
      do j=1,nt
         xj(j) = (targ(1,j) - center(1))*dsq
         yj(j) = (targ(2,j) - center(2))*dsq
         wj(j) = exp(0.5d0*eye*(xj(j)+yj(j)))
      enddo

      z = eye/dsqrt(delta)
      
      do ind = 1,nd
         do j=1,nexp
            fk(j,0,ind)=pwexp(j,ind)
         enddo

         j=0
         do j2=-npw/2,-1
            do j1=-npw/2,npw/2-1
               j=j+1
               fk(j,1,ind)=pwexp(j,ind)*z*ts(j1)
               fk(j,2,ind)=pwexp(j,ind)*z*ts(j2)
            enddo
         enddo
         
         do k=0,2
            do j=nexp+1,2*nexp
               fk(j,k,ind)=0
            enddo
         enddo
      enddo

      nt8=nt
      npw8=npw
      ntrans=3*nd
      iflag = 1
      call finufft2d2many(ntrans,nt8,xj,yj,cj,iflag,eps,
     1    npw8,npw8,fk,null,ier)
      
      do ind=1,nd
         do j=1,nt
            pot(ind,j)=pot(ind,j)+dreal(cj(j,0,ind)*wj(j))*2
            grad(ind,1,j)=grad(ind,1,j)+dreal(cj(j,1,ind)*wj(j))*2
            grad(ind,2,j)=grad(ind,2,j)+dreal(cj(j,2,ind)*wj(j))*2
         enddo
      enddo
c
      return
      end
C
C
c
C
C
C
      subroutine g2dpwevalh_vec(nd,delta,eps,center,npw,ws,ts,
     1              nexp,pwexp,targ,nt,pot,grad,hess)
C
C     This subroutine evaluates the plane wave 
C     expansions about CENTER at location TARG
C     potential + gradient + hess
C
C     INPUT
C
c     nd            = vector length (for multiple charges at same locations)
C     delta         = Gaussian variance
C     charge        = strength of sources
C     center        = center of the expansion
C     npw           = number of Fourier plane waves
C     ws,ts         = planewave weights and nodes
C     pwexp         = pw expansions 
C     targ          = target
C
C     OUTPUT:
C     pot           = potential (or vectorized potentials) incremented
C     grad          = gradient (or vectorized gradients) incremented
C     hess          = hessian (or vectorized hessians) incremented
C
      implicit real *8 (a-h,o-z)
      integer nd,npw,ntarg
      real *8 delta,center(2),targ(2,nt)
      real *8 pot(nd,nt),grad(nd,2,nt),hess(nd,3,nt)
      
      real *8 ws(-npw/2:npw/2-1),ts(-npw/2:npw/2-1)
      complex *16 pwexp(nexp,nd)
      integer*8 nt8, npw8

      integer i,ind,j1,j2,j,npw2,itarg,k
      real *8 x,y,dsq
      complex *16 eye
      complex *16 qqx,qqy,qq1,qq2
      
      complex *16 z,cd
      real *8, allocatable ::  xj(:),yj(:)
      complex *16, allocatable :: cj(:,:,:),fk(:,:,:),wj(:)
      complex *16, allocatable :: cp(:)
      
c     this (since unallocated) used to pass a NULL ptr to FINUFFT...
      integer*8, allocatable :: null
C
      eye = dcmplx(0,1)
      
      allocate(xj(nt))
      allocate(yj(nt))

      allocate(cj(nt,0:5,nd))
      allocate(wj(nt))
      allocate(fk(2*nexp,0:5,nd))
      allocate(cp(-npw/2:npw/2-1))
      
      dsq = 2*ts(0)/dsqrt(delta)
C
      do j=1,nt
         xj(j) = (targ(1,j) - center(1))*dsq
         yj(j) = (targ(2,j) - center(2))*dsq
         wj(j) = exp(0.5d0*eye*(xj(j)+yj(j)))
      enddo

      z = eye/dsqrt(delta)
      do j=-npw/2,npw/2-1
         cp(j) = z*ts(j)
      enddo
      
      do ind = 1,nd
         do j=1,nexp
            fk(j,0,ind)=pwexp(j,ind)
         enddo

         j=0
         do j2=-npw/2,-1
            do j1=-npw/2,npw/2-1
               j=j+1
               fk(j,1,ind)=pwexp(j,ind)*cp(j1)
               fk(j,2,ind)=pwexp(j,ind)*cp(j2)
                  
               fk(j,3,ind)=pwexp(j,ind)*cp(j1)*cp(j1)
               fk(j,4,ind)=pwexp(j,ind)*cp(j1)*cp(j2)
               fk(j,5,ind)=pwexp(j,ind)*cp(j2)*cp(j2)
            enddo
         enddo
         
         do k=0,5
            do j=nexp+1,2*nexp
               fk(j,k,ind)=0
            enddo
         enddo
      enddo

      nt8=nt
      npw8=npw
      ntrans=6*nd
      iflag = 1
      
      call finufft2d2many(ntrans,nt8,xj,yj,cj,iflag,eps,
     1    npw8,npw8,fk,null,ier)
      
      do ind=1,nd
         do j=1,nt
            pot(ind,j)=pot(ind,j)+dreal(cj(j,0,ind)*wj(j))*2
            
            grad(ind,1,j)=grad(ind,1,j)+dreal(cj(j,1,ind)*wj(j))*2
            grad(ind,2,j)=grad(ind,2,j)+dreal(cj(j,2,ind)*wj(j))*2
            
            hess(ind,1,j)=hess(ind,1,j)+dreal(cj(j,3,ind)*wj(j))*2
            hess(ind,2,j)=hess(ind,2,j)+dreal(cj(j,4,ind)*wj(j))*2
            hess(ind,3,j)=hess(ind,3,j)+dreal(cj(j,5,ind)*wj(j))*2
         enddo
      enddo
c
      return
      end
C
C
c
C
C
C
c*********************************************************************
C
C get plane wave approximation nodes and weights
C
C*********************************************************************
      subroutine get_pwnodes(pmax,npw,ws,ts)
C
C     Get planewave exp weights,nodes
C
      implicit real *8 (a-h,o-z)
      real *8 ws(-npw/2:npw/2-1),ts(-npw/2:npw/2-1)

      pi = 4.0d0*datan(1.0d0)
      npw2=npw/2
      h = pmax/npw2
      w = h/(2.0d0*dsqrt(pi))

      do j =-npw2,npw2-1
         ts(j) = (j+0.5d0)*h
         ws(j) = w*dexp(-ts(j)*ts(j)/4)
      enddo
c
      return
      end
C
C
c
C
c*********************************************************************
C
C shift PW expansions (mp to loc at the cutoff level)
C
C*********************************************************************
      subroutine nufft_weights(npw,ws,ts,
     1    nexp,wnufft)
C
C     This subroutine precomputes all translation matrices for all SOE/X
C     expansions from child to parent or vice versa.
C
c     used in mp to mp or loc to loc stage
c      
C     INPUT
C
c     nd      = vector length (for multiple charges at same locations)
C     npw     = number of terms in plane wave exp
C     nmax    = number of different translation lengths in the whole scheme 
C     ts      = pw nodes
c      
C     OUTPUT:
C
C     wshift  = table of translation matrices for PW shift 
C
      implicit real *8 (a-h,o-z)
      real *8 ws(-npw/2:npw/2-1),ts(-npw/2:npw/2-1)
      
      real *8 wnufft(nexp)
      
      complex *16 eye,ztmp
C
      eye = dcmplx(0,1)
      
      j=0
      do j2=-npw/2,-1
         do j1=-npw/2,npw/2-1
            j=j+1
            wnufft(j) = ws(j2)*ws(j1)
         enddo
      enddo
c
      return
      end
c
c
c     
c
      subroutine pw_translation_matrices(xmin,npw,ts,nmax,
     1    nexp,wshift)
C
C     This subroutine precomputes all translation matrices for all SOE/X
C     expansions from child to parent or vice versa.
C
c     used in mp to mp or loc to loc stage
c      
C     INPUT
C
c     nd      = vector length (for multiple charges at same locations)
C     npw     = number of terms in plane wave exp
C     nmax    = number of different translation lengths in the whole scheme 
C     ts      = pw nodes
c      
C     OUTPUT:
C
C     wshift  = table of translation matrices for PW shift 
C
      implicit real *8 (a-h,o-z)
      real *8 ts(-npw/2:npw/2-1)
      
      complex *16 wshift(nexp,-nmax:nmax,-nmax:nmax)
      
      complex *16,allocatable:: ww(:,:)

      complex *16 eye,ztmp
C
      eye = dcmplx(0,1)
      
      allocate(ww(-npw/2:npw/2-1,-nmax:nmax))


      do j1=-npw/2,npw/2-1
         do k1=-nmax,nmax
            ww(j1,k1) = exp(eye*ts(j1)*xmin*k1)
         enddo
      enddo
      
      do k1=-nmax,nmax
      do k2=-nmax,nmax
         j=0
         do j1=-npw/2,-1
         do j2=-npw/2,npw/2-1
            j=j+1
            wshift(j,k2,k1) = ww(j2,k2)*ww(j1,k1)
         enddo
         enddo
      enddo
      enddo
c
      return
      end
c
c
c     
c
      subroutine g2dshiftpw_vec(nd,nexp,pwexp1,
     1              pwexp2,wshift)
C
C     This subroutine converts the PW expansion (pwexp1) about
C     the center (CENT1) into an PW expansion (pwexp2) about 
C     (CENT2) using precomputed translation matrix wshift.
C
C     INPUT
C
c     nd      = vector length (for vector input)
C     delta   = Gaussian variance
C     nn      = number of terms in PW expansion
C     pwexp1  = original expansion 
C     wshift  = precomputed PW exp translation matrix 
C
C     OUTPUT:
C
C     pwexp2 = shifted expansion 
C
      implicit none
      integer nd,j,ind,nexp
      complex *16 pwexp1(nexp,nd)
      complex *16 pwexp2(nexp,nd)
      complex *16 wshift(nexp)

C
      do ind=1,nd
         do j=1,nexp
            pwexp2(j,ind) = pwexp2(j,ind)
     1          +pwexp1(j,ind)*wshift(j)
         enddo
      enddo
c
      return
      end
c
C
c
C
      subroutine g2dcopypwexp_vec(nd,nexp,pwexp1,
     1              pwexp2)
C
C     This subroutine copy one PW expansion (pwexp1) 
C     to an PW expansion (pwexp2).
C
C     INPUT
C
c     nd      = vector length (for vector input)
C     nn      = number of terms in PW expansion
C     pwexp1  = original expansion 
C
C     OUTPUT:
C
C     pwexp2 = copied expansion 
C
      implicit none
      integer nd,nn,j,ind,nexp
      complex *16 pwexp1(nexp,nd)
      complex *16 pwexp2(nexp,nd)

C
      do ind=1,nd
         do j=1,nexp
            pwexp2(j,ind) = pwexp1(j,ind)
         enddo
      enddo
c
      return
      end
c
C
c
C
