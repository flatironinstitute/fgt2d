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
c     g2dpwevalp_vec : evaluates the PW expansion
c                      potential only
c      
c     g2dpwevalg_vec : evaluates the PW expansion
c                      potential + gradient
c      
c     g2dpwevalh_vec : evaluates the PW expansion
c                      potential + gradient + hessian
c      
C     g2dh2pw_vec : converts the Hermite expansion to the PW expansion
c      
c     g2dh2pw_real_vec : same as above, use REAL 1D h2pw matrix
c
c     g2dherm2local_vec : converts the Hermite expansion to the Taylor expansion
c      
c     g2dpw2local_vec : converts PW expansions to the local Taylor expansion
c
c     g2dpw2local_real_vec : converts PW expansion to local expansion
c                            USE REAL 1D pw2l matrix
c      
c     g2dh2lmat : returns the matrix converting the Hermite expansion
C                 to the local Taylor expansion (1D only)
c      
c     g2dh2xmat : returns the matrix converting the Hermite expansion 
C                 to the planewave expansion (1D only)
C     
c     g2drh2xmat : returns REAL matrix converting the Hermite expansion 
C                 to the planewave expansion (1D only)
C                 That is, i^j factor is ignored
C     
C     g2dx2lmat : returns the matrix converting the planewave expansion
C                 to the local Taylor expansion (1D only)
C
C     g2drx2lmat : returns REAL matrix converting the planewave expansion
C                 to the local Taylor expansion (1D only),
C                 that is, i^j factor is ignored.
c
c     get_pwnodes : returns PW weights and nodes, midpoint rule is used
c                   so the number is always an even number!
c
c     pw_translation_matrices : returns precomputed translation matrices
c                   for PW mp to loc translations
c
c     merge_split_pw_matrices : returns precomputed tranlation matrices
c                   for PW mp2mp or loc2loc translations
c
c     g2dshiftpw_vec : translates PW expansion via precomputed translation 
c                      matrices
c
c     g2dcopypwexp_vec : copy PW expansion
c
c     g2dshiftpw0_vec : translates PW expansion from one center to another
c                       without any precomputation
c
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
      subroutine g2dh2pw_vec(nd,ntermsh,npw,h2x,
     1    hexp,pwexp)
C
C     This subroutine converts the Hermite expansion to the PW
C     expansion
C
C     INPUT
C
c     nd            = vector length (parallel)
C     delta         = Gaussian variance
C     center        = center of the expansion
C     ntermsh       = number of terms in Hermite exp
C     npw           = number of exponentials in PW exp 
C     h2x           = matrix converting the Hermite expansion to the X expansion (1D)
C     hexp          = Hermite expansion
C      
C     OUTPUT:
C
C     hwexp         = plane wave expansion
C      
      implicit real*8 (a-h,o-z)
      complex *16 h2x(npw,0:ntermsh)
      real *8 hexp(0:ntermsh,0:ntermsh,nd)
      complex *16 pwexp(npw/2,npw,nd)

      complex *16, allocatable :: pxhy(:,:)
      complex *16 eye,cd
      
      eye = dcmplx(0,1)
c
      allocate(pxhy(npw/2,0:ntermsh))

c      
      do ind=1,nd
         do k2=1,npw/2
            do j2=0,ntermsh
               cd=0
               do j1=0,ntermsh
                  cd=cd+h2x(k2,j1)*hexp(j1,j2,ind)
               enddo
               pxhy(k2,j2)=cd
            enddo
         enddo

         do k1=1,npw
            do k2=1,npw/2
               cd=0
               do j2=0,ntermsh
                  cd=cd+h2x(k1,j2)*pxhy(k2,j2)
               enddo
               pwexp(k2,k1,ind) = pwexp(k2,k1,ind) + cd
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
      subroutine g2dh2pw_real_vec(nd,ntermsh,npw,h2x,
     1    hexp,pwexp)
C
C     This subroutine converts the Hermite expansion to the PW
C     expansion
C
C     INPUT
C
c     nd            = vector length (parallel)
C     delta         = Gaussian variance
C     center        = center of the expansion
C     ntermsh       = number of terms in Hermite exp
C     npw           = number of exponentials in PW exp 
C     h2x           = REAL matrix converting the Hermite expansion to the X expansion (1D)
C     hexp          = Hermite expansion
C      
C     OUTPUT:
C
C     hwexp         = plane wave expansion
C      
      implicit real*8 (a-h,o-z)
      real *8 h2x(npw,0:ntermsh)
      real *8 hexp(0:ntermsh,0:ntermsh,nd)
      complex *16 pwexp(npw/2,npw,nd)

      complex *16, allocatable :: pxhy(:,:)
      complex *16 eye,c0,c1,c2,c3,cd
      
      eye = dcmplx(0,1)
c
cccc      allocate(pxhy(npw/2,0:ntermsh))
      allocate(pxhy(0:ntermsh,npw/2))

c      
      do ind=1,nd
         do k2=1,npw/2
            do j2=0,ntermsh
               c0=0
               c1=0
               c2=0
               c3=0
               do j1=0,ntermsh,4
                  c0=c0+h2x(k2,j1)*hexp(j1,j2,ind)
               enddo
               do j1=1,ntermsh,4
                  c1=c1+h2x(k2,j1)*hexp(j1,j2,ind)
               enddo
               do j1=2,ntermsh,4
                  c2=c2+h2x(k2,j1)*hexp(j1,j2,ind)
               enddo
               do j1=3,ntermsh,4
                  c3=c3+h2x(k2,j1)*hexp(j1,j2,ind)
               enddo
cccc               pxhy(k2,j2)= c0-c2+eye*(c1-c3)
               pxhy(j2,k2)= c0-c2+eye*(c1-c3)
            enddo
         enddo

         do k1=1,npw
            do k2=1,npw/2
               c0=0
               c1=0
               c2=0
               c3=0
               do j2=0,ntermsh,4
cccc                  c0=c0+h2x(k1,j2)*pxhy(k2,j2)
                  c0=c0+h2x(k1,j2)*pxhy(j2,k2)
               enddo
               do j2=1,ntermsh,4
cccc                  c1=c1+h2x(k1,j2)*pxhy(k2,j2)
                  c1=c1+h2x(k1,j2)*pxhy(j2,k2)
               enddo
               do j2=2,ntermsh,4
cccc                  c2=c2+h2x(k1,j2)*pxhy(k2,j2)
                  c2=c2+h2x(k1,j2)*pxhy(j2,k2)
               enddo
               do j2=3,ntermsh,4
cccc                  c3=c3+h2x(k1,j2)*pxhy(k2,j2)
                  c3=c3+h2x(k1,j2)*pxhy(j2,k2)
               enddo
               pwexp(k2,k1,ind) = pwexp(k2,k1,ind) +
     1             + c0-c2+eye*(c1-c3)
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
C 2D translations pw2local
C
c***********************************************************************      
      subroutine g2dpw2local_vec(nd,nlocal,npw,pw2l,
     1    pwexp,local)
C
C     This subroutine converts four SOE expansions to local
C     expansions
C
C     INPUT
C
c     nd            = vector length (parallel)
C     nlocal        = number of terms in local exp
C     npw          = number of exponentials in PW exp
C     pw2l           = matrix converting PW expansion to the local expansion
C     pwexp        = PW expansions
C      
C     OUTPUT:
C
C     local         = local Taylor expansions
C      
      implicit real*8 (a-h,o-z)
      complex *16 pw2l(npw,0:nlocal)
      real *8 local(0:nlocal,0:nlocal,nd)
      complex *16 pwexp(npw/2,npw,nd)
      complex *16, allocatable :: lxpy(:,:)
      complex *16 cd
c      
c
      allocate(lxpy(0:nlocal,npw))

      npw2=npw/2

      do ind=1,nd
         do k2=0,nlocal
            do j1=1,npw
               cd=0
               do j2=1,npw/2
                  cd=cd+pw2l(j2,k2)*pwexp(j2,j1,ind)
               enddo
               lxpy(k2,j1)=cd
            enddo
         enddo

cccc         do k1=0,nlocal
cccc            do k2=0,nlocal
cccc               cd=0
cccc               do j1=1,npw
cccc                  cd=cd+pw2l(j1,k1)*lxpy(k2,j1)
cccc               enddo
cccc               local(k2,k1,ind)=local(k2,k1,ind)+dble(cd)*2
cccc            enddo
cccc         enddo
         do k1=0,nlocal,2
            do k2=0,nlocal
               cd=0
               do j1=1,npw2
                  cd=cd+pw2l(j1,k1)*(lxpy(k2,j1)+lxpy(k2,j1+npw2))
               enddo
               local(k2,k1,ind)=local(k2,k1,ind)+dble(cd)*2
            enddo
         enddo
         
         do k1=1,nlocal,2
            do k2=0,nlocal
               cd=0
               do j1=1,npw2
                  cd=cd+pw2l(j1,k1)*(lxpy(k2,j1)-lxpy(k2,j1+npw2))
               enddo
               local(k2,k1,ind)=local(k2,k1,ind)+dble(cd)*2
            enddo
         enddo
      enddo
         
      return
      end
C
C
C
C
      subroutine g2dpw2local_real_vec(nd,nlocal,npw,pw2l,
     1    pwexp,local)
C
C     This subroutine converts four SOE expansions to local
C     expansions
C
C     INPUT
C
c     nd            = vector length (parallel)
C     nlocal        = number of terms in local exp
C     npw          = number of exponentials in PW exp
C     pw2l         = REAL matrix converting PW expansion to the local expansion
C     pwexp        = PW expansions
C      
C     OUTPUT:
C
C     local         = local Taylor expansions
C      
      implicit real*8 (a-h,o-z)
      real *8 pw2l(npw,0:nlocal)
      real *8 local(0:nlocal,0:nlocal,nd)
      complex *16 pwexp(npw/2,npw,nd)
      complex *16, allocatable :: lxpy(:,:)
      complex *16 eye,cd
c      
c
      eye = dcmplx(0,1)
      
      allocate(lxpy(0:nlocal,npw))

      npw2=npw/2

      do ind=1,nd
         do k2=0,nlocal,4
            do j1=1,npw
               cd=0
               do j2=1,npw/2
                  cd=cd+pw2l(j2,k2)*pwexp(j2,j1,ind)
               enddo
               lxpy(k2,j1)=cd
            enddo
         enddo
         do k2=1,nlocal,4
            do j1=1,npw
               cd=0
               do j2=1,npw/2
                  cd=cd+pw2l(j2,k2)*pwexp(j2,j1,ind)
               enddo
               lxpy(k2,j1)=cd*eye
            enddo
         enddo
         do k2=2,nlocal,4
            do j1=1,npw
               cd=0
               do j2=1,npw/2
                  cd=cd+pw2l(j2,k2)*pwexp(j2,j1,ind)
               enddo
               lxpy(k2,j1)=-cd
            enddo
         enddo
         do k2=3,nlocal,4
            do j1=1,npw
               cd=0
               do j2=1,npw/2
                  cd=cd+pw2l(j2,k2)*pwexp(j2,j1,ind)
               enddo
               lxpy(k2,j1)=-eye*cd
            enddo
         enddo

cccc         do k1=0,nlocal
cccc            do k2=0,nlocal
cccc               cd=0
cccc               do j1=1,npw
cccc                  cd=cd+pw2l(j1,k1)*lxpy(k2,j1)
cccc               enddo
cccc               local(k2,k1,ind)=local(k2,k1,ind)+dble(cd)*2
cccc            enddo
cccc         enddo
         do k1=0,nlocal,4
            do k2=0,nlocal
               cd=0
               do j1=1,npw2
                  cd=cd+pw2l(j1,k1)*(lxpy(k2,j1)+lxpy(k2,j1+npw2))
               enddo
               local(k2,k1,ind)=local(k2,k1,ind)+dble(cd)*2
            enddo
         enddo
         do k1=1,nlocal,4
            do k2=0,nlocal
               cd=0
               do j1=1,npw2
                  cd=cd+pw2l(j1,k1)*(lxpy(k2,j1)-lxpy(k2,j1+npw2))
               enddo
               local(k2,k1,ind)=local(k2,k1,ind)-aimag(cd)*2
            enddo
         enddo
         do k1=2,nlocal,4
            do k2=0,nlocal
               cd=0
               do j1=1,npw2
                  cd=cd+pw2l(j1,k1)*(lxpy(k2,j1)+lxpy(k2,j1+npw2))
               enddo
               local(k2,k1,ind)=local(k2,k1,ind)-dble(cd)*2
            enddo
         enddo
         do k1=3,nlocal,4
            do k2=0,nlocal
               cd=0
               do j1=1,npw2
                  cd=cd+pw2l(j1,k1)*(lxpy(k2,j1)-lxpy(k2,j1+npw2))
               enddo
               local(k2,k1,ind)=local(k2,k1,ind)+aimag(cd)*2
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
      real *8, allocatable :: sqc(:,:),fac(:),hexpx(:)

      ntot = ntermsh+nlocal
      allocate(sqc(0:ntot,0:ntot),fac(0:nlocal),hexpx(0:ntot))
      
      call getsqrtbinomialcoeffs(ntot,sqc)
      
      fac(0)=1.0d0
      do i=1,nlocal
         fac(i)=fac(i-1)*sqrt(i*1.0d0)
      enddo

      x=0
      hexpx(0)=1.0d0
      hexpx(1)=2*x
      
      do i=1,ntermsh+nlocal-1
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
      real *8 ws(npw), ts(npw)
      complex *16 h2x(npw,0:ntermsh)
      real *8, allocatable:: fac(:)
      complex *16 eye
      
      eye = dcmplx(0,1)
      
      allocate(fac(0:ntermsh))

      fac(0)=1.0d0
      do i=1,ntermsh
         fac(i)=fac(i-1)/sqrt(dble(i))
cccc         fac(i)=fac(i-1)/i
      enddo
      
      do j=0,ntermsh
         do i=1,npw
            h2x(i,j)=ws(i)*(-eye*ts(i))**j*fac(j)
         enddo
      enddo
      
      return
      end
C
C
C
C
      subroutine g2drh2xmat(ntermsh,npw,ws,ts,h2x)
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
C     h2x       = real translation matrices
C
      implicit real*8 (a-h,o-z)
      real *8 ws(npw), ts(npw)
      real *8 h2x(npw,0:ntermsh)
      real *8, allocatable:: fac(:)
      complex *16 eye
      
      eye = dcmplx(0,1)
      
      allocate(fac(0:ntermsh))

      fac(0)=1.0d0
      do i=1,ntermsh
         fac(i)=fac(i-1)/sqrt(dble(i))
cccc         fac(i)=fac(i-1)/i
      enddo
      
      do j=0,ntermsh
         do i=1,npw
            h2x(i,j)=ws(i)*(-ts(i))**j*fac(j)
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
      real *8 ws(npw), ts(npw)
      complex *16 x2l(npw,0:nlocal)
      complex *16 eye
      real *8, allocatable:: fac(:)

      allocate(fac(0:nlocal))
      
      eye = dcmplx(0,1)
      
      fac(0)=1.0d0
      do i=1,nlocal
         fac(i)=fac(i-1)*i
      enddo
      
      do j=0,nlocal
         do i=1,npw
            x2l(i,j)=(eye*ts(i))**j/fac(j)
         enddo
      enddo

      return
      end
C
C
C
      subroutine g2drx2lmat(nlocal,npw,ws,ts,x2l)
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
C     x2l        = real translation matrix
C
      implicit real*8 (a-h,o-z)
      real *8 ws(npw), ts(npw)
      real *8 x2l(npw,0:nlocal)
      complex *16 eye
      real *8, allocatable:: fac(:)

      allocate(fac(0:nlocal))
      
      eye = dcmplx(0,1)
      
      fac(0)=1.0d0
      do i=1,nlocal
         fac(i)=fac(i-1)*i
      enddo
      
      do j=0,nlocal
         do i=1,npw
            x2l(i,j)=(ts(i))**j/fac(j)
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
c
c      
c      
c      
c      
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
      subroutine g2dpwevalp_vec(nd,delta,center,npw,wx,tx,
     1              pwexp,targ,ntarg,pot)
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
C     wx,tx         = planewave weights and nodes
C     pwexp         = pw expansions 
C     targ          = target
C
C     OUTPUT:
C     pot           = potential (or vectorized potentials) incremented
C
      implicit none
      integer nd,npw,ntarg
      real *8 delta,center(2),targ(2,ntarg)
      real *8 pot(nd,ntarg)
      
      real *8 wx(npw),tx(npw)
      complex *16 pwexp(npw/2,npw,nd)

      integer i,ind,j1,j2,j,npw2,itarg,k
      real *8 x,y,dsq
      complex *16 eye
      complex *16 qqx,qqy,qq1,qq2
      
      complex *16 ww1(100)
      complex *16 ww2(100)

      complex *16 z,cd
C
      eye = dcmplx(0,1)
      dsq = 1.0D0/dsqrt(delta)
C
      npw2=npw/2
      do itarg=1,ntarg
         x = (targ(1,itarg) - center(1))*dsq
         y = (targ(2,itarg) - center(2))*dsq

         qqx = cdexp(eye*tx(1)*x)
         qqy = cdexp(eye*tx(1)*y)
         qq1 = qqx
         qq2 = qqy
         qqx = qqx*qqx
         qqy = qqy*qqy
         

         do j1=1,npw2
            ww1(j1) = qq1
            ww2(j1) = qq2            
            qq1 = qq1*qqx
            qq2 = qq2*qqy
            ww1(j1+npw2) = dconjg(ww1(j1))
            ww2(j1+npw2) = dconjg(ww2(j1))
         enddo
c
         do ind = 1,nd
            z=0
            do j2=1,npw
               cd=0
               do j1=1,npw/2
                  cd=cd+pwexp(j1,j2,ind)*ww1(j1)
               enddo
               z=z+cd*ww2(j2)
            enddo
            
            pot(ind,itarg) = pot(ind,itarg)+dreal(z)*2
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
      subroutine g2dpwevalg_vec(nd,delta,center,npw,wx,tx,
     1              pwexp,targ,ntarg,pot,grad)
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
C     wx,tx         = planewave weights and nodes
C     pwexp         = pw expansions 
C     targ          = target
C
C     OUTPUT:
C     pot           = potential (or vectorized potentials) incremented
C     grad          = gradient (or vectorized gradients) incremented
C
      implicit none
      integer nd,npw,ntarg
      real *8 delta,center(2),targ(2,ntarg)
      real *8 pot(nd,ntarg),grad(nd,2,ntarg)
      
      real *8 wx(npw),tx(npw)
      complex *16 pwexp(npw/2,npw,nd)

      integer i,ind,j1,j2,j,npw2,itarg,k
      real *8 x,y,dsq
      complex *16 eye
      complex *16 qqx,qqy,qq1,qq2
      
      complex *16 ww1(100),ww1x(100)
      complex *16 ww2(100),ww2y(100)

      complex *16 z,zx,zy,c,cx,cy
C
      eye = dcmplx(0,1)
      dsq = 1.0D0/dsqrt(delta)
C
      npw2=npw/2
      do itarg=1,ntarg
         x = (targ(1,itarg) - center(1))*dsq
         y = (targ(2,itarg) - center(2))*dsq

         qqx = cdexp(eye*tx(1)*x)
         qqy = cdexp(eye*tx(1)*y)
         qq1 = qqx
         qq2 = qqy
         qqx = qqx*qqx
         qqy = qqy*qqy
         

         do j1=1,npw2
            ww1(j1) = qq1
            ww2(j1) = qq2            
            ww1(j1+npw2) = dconjg(ww1(j1))
            ww2(j1+npw2) = dconjg(ww2(j1))

            ww1x(j1)= eye*tx(j1)*dsq*ww1(j1)
            ww2y(j1)= eye*tx(j1)*dsq*ww2(j1)
            ww1x(j1+npw2)=dconjg(ww1x(j1))
            ww2y(j1+npw2)=dconjg(ww2y(j1))

            qq1 = qq1*qqx
            qq2 = qq2*qqy
         enddo
c
         do ind = 1,nd
            z=0
            zx=0
            zy=0
            do j2=1,npw
               c=0
               cx=0
               cy=0
               do j1=1,npw/2
                  c = c+pwexp(j1,j2,ind)*ww1(j1)
                  cx = cx+pwexp(j1,j2,ind)*ww1x(j1)
cccc                  cy=cy+pwexp(j1,j2,ind)*ww1(j1)
               enddo
               z = z+c*ww2(j2)
               zx = zx+cx*ww2(j2)
               zy = zy+c*ww2y(j2)
            enddo
            
            pot(ind,itarg) = pot(ind,itarg)+dreal(z)*2
            grad(ind,1,itarg) = grad(ind,1,itarg)+dreal(zx)*2
            grad(ind,2,itarg) = grad(ind,2,itarg)+dreal(zy)*2
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
      subroutine g2dpwevalh_vec(nd,delta,center,npw,wx,tx,
     1              pwexp,targ,ntarg,pot,grad,hess)
C
C     This subroutine evaluates the plane wave 
C     expansions about CENTER at location TARG
C     potential + gradient + hessian
C
C     INPUT
C
c     nd            = vector length (for multiple charges at same locations)
C     delta         = Gaussian variance
C     charge        = strength of sources
C     center        = center of the expansion
C     npw           = number of Fourier plane waves
C     wx,tx         = planewave weights and nodes
C     pwexp         = pw expansions 
C     targ          = target
C
C     OUTPUT:
C     pot           = potential (or vectorized potentials) incremented
C     grad          = gradient (or vectorized gradients) incremented
C     hess          = Hessian (or vectorized Hessians) incremented
C
      implicit none
      integer nd,npw,ntarg
      real *8 delta,center(2),targ(2,ntarg)
      real *8 pot(nd,ntarg),grad(nd,2,ntarg),hess(nd,3,ntarg)
      
      real *8 wx(npw),tx(npw)
      complex *16 pwexp(npw/2,npw,nd)

      integer i,ind,j1,j2,j,npw2,itarg,k
      real *8 x,y,dsq
      complex *16 eye
      complex *16 qqx,qqy,qq1,qq2
      
      complex *16 ww1(100),ww1x(100),ww1xx(100)
      complex *16 ww2(100),ww2y(100),ww2yy(100)

      complex *16 z,zx,zy,zxx,zxy,zyy
      complex *16 c,cx,cy,cxx,cxy,cyy

C
      eye = dcmplx(0,1)
      dsq = 1.0D0/dsqrt(delta)
C
      npw2=npw/2
      do itarg=1,ntarg
         x = (targ(1,itarg) - center(1))*dsq
         y = (targ(2,itarg) - center(2))*dsq

         qqx = cdexp(eye*tx(1)*x)
         qqy = cdexp(eye*tx(1)*y)
         qq1 = qqx
         qq2 = qqy
         qqx = qqx*qqx
         qqy = qqy*qqy
         

         do j1=1,npw2
            ww1(j1) = qq1
            ww2(j1) = qq2            
            ww1(j1+npw2) = dconjg(ww1(j1))
            ww2(j1+npw2) = dconjg(ww2(j1))

            ww1x(j1)= eye*tx(j1)*dsq*ww1(j1)
            ww2y(j1)= eye*tx(j1)*dsq*ww2(j1)
            ww1x(j1+npw2)=dconjg(ww1x(j1))
            ww2y(j1+npw2)=dconjg(ww2y(j1))

            ww1xx(j1)= eye*tx(j1)*dsq*ww1x(j1)
            ww2yy(j1)= eye*tx(j1)*dsq*ww2y(j1)
            ww1xx(j1+npw2)=dconjg(ww1xx(j1))
            ww2yy(j1+npw2)=dconjg(ww2yy(j1))

            qq1 = qq1*qqx
            qq2 = qq2*qqy
         enddo
c
         do ind = 1,nd
            z=0
            zx=0
            zy=0
            zxx=0
            zxy=0
            zyy=0
            do j2=1,npw
               c=0
               cx=0
               cy=0
               cxx=0
               cxy=0
               cyy=0               
               do j1=1,npw/2
                  c=c+pwexp(j1,j2,ind)*ww1(j1)
                  cx=cx+pwexp(j1,j2,ind)*ww1x(j1)
                  cxx=cxx+pwexp(j1,j2,ind)*ww1xx(j1)                  
               enddo
               z=z+c*ww2(j2)
               zx=zx+cx*ww2(j2)
               zy=zy+c*ww2y(j2)
               
               zxx=zxx+cxx*ww2(j2)
               zxy=zxy+cx*ww2y(j2)
               zyy=zyy+c*ww2yy(j2)
            enddo
            
            pot(ind,itarg) = pot(ind,itarg)+dreal(z)*2
            grad(ind,1,itarg) = grad(ind,1,itarg)+dreal(zx)*2
            grad(ind,2,itarg) = grad(ind,2,itarg)+dreal(zy)*2
            hess(ind,1,itarg) = hess(ind,1,itarg)+dreal(zxx)*2
            hess(ind,2,itarg) = hess(ind,2,itarg)+dreal(zxy)*2
            hess(ind,3,itarg) = hess(ind,3,itarg)+dreal(zyy)*2
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
      real *8 ws(npw),ts(npw)

      pi = 4.0d0*datan(1.0d0)
      npw2=npw/2
      h = pmax/npw2
      w = h/(2.0d0*dsqrt(pi))

      do j =1,npw2
         ts(j) = (j-0.5d0)*h
         ws(j) = w*dexp(-ts(j)*ts(j)/4)
         ts(j+npw2)=-ts(j)
         ws(j+npw2)=ws(j)
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
C shift PW expansions (mp to loc, mp to mp, loc to loc)
C
C*********************************************************************
      subroutine pw_translation_matrices(xmin,npw,nmax,
     1              wshift,ts)
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
      real *8 ts(npw)
      
      complex *16 wshift(npw/2,npw,-nmax:nmax,-nmax:nmax)
      
      complex *16,allocatable:: ww(:,:)

      complex *16 eye,ztmp
C
      eye = dcmplx(0,1)
      
      allocate(ww(npw,-nmax:nmax))


      do j1=1,npw
         ztmp = exp(eye*ts(j1)*xmin)
         ww(j1,0)=1
         do k1=1,nmax
            ww(j1,k1) = ztmp
            ztmp = ztmp*ztmp
            ww(j1,-k1) = dconjg(ww(j1,k1))
         enddo
      enddo
      
      do k1=-nmax,nmax
         do k2=-nmax,nmax
            do j1=1,npw
               do j2=1,npw/2
                  wshift(j2,j1,k2,k1) = ww(j2,k2)*ww(j1,k1)
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
      subroutine merge_split_pw_matrices(xmin,npw,nmax,
     1           wshift,ts)
C
C     This subroutine precomputes all translation matrices for 
c     PW translations from child to parent or vice versa.
C
C     INPUT
C
c     xmin     = half of the scaled (by 1/sqrt(delta) size of the box 
c                at the finest level
C     npw     = number of terms in PW expansion
C     nmax     = number of different translation lengths in the whole scheme 
C     ws,ts    = PW weights and nodes
C
C     OUTPUT:
C
C     wshift   = table of translation matrices for PW  shift used in
c                merge mp and split loc stage
C
      implicit real *8 (a-h,o-z)
      real *8 xmin
      complex *16 wshift(npw/2,npw,4,nmax)
      real *8 ts(npw)
      complex *16 ztmp
      complex *16 eye
      complex *16, allocatable:: ww(:,:)
C
      eye =dcmplx(0,1)
      
      allocate(ww(npw,nmax))


      do j1=1,npw
         ztmp = exp(eye*ts(j1)*xmin)
         do k1=1,nmax
            ww(j1,k1) = ztmp
            ztmp = ztmp*ztmp
         enddo
      enddo
      
      do k1=1,nmax
         do j1=1,npw
            do j2=1,npw/2
c              pp               
               wshift(j2,j1,1,k1) = ww(j2,k1)*ww(j1,k1)
c              pm
               wshift(j2,j1,2,k1) = ww(j2,k1)*conjg(ww(j1,k1))
c              mp
               wshift(j2,j1,3,k1) = conjg(ww(j2,k1))*ww(j1,k1)
c              mm
               wshift(j2,j1,4,k1) = conjg(ww(j2,k1))*conjg(ww(j1,k1))
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
      subroutine g2dshiftpw_vec(nd,nn,pwexp1,
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
      integer nd,nn,j,j1,j2,ind,nexp
cccc      complex *16 pwexp1(nn/2,nn,nd)
cccc      complex *16 pwexp2(nn/2,nn,nd)
cccc      complex *16 wshift(nn/2,nn)
      
      complex *16 pwexp1(nn*nn/2,nd)
      complex *16 pwexp2(nn*nn/2,nd)
      complex *16 wshift(nn*nn/2)

C
      nexp = nn*nn/2
      
      do ind=1,nd
cccc         do j1=1,nn
cccc            do j2=1,nn/2
cccc               pwexp2(j2,j1,ind) = pwexp2(j2,j1,ind)
cccc     1             +pwexp1(j2,j1,ind)*wshift(j2,j1)
cccc            enddo
cccc         enddo
         do j=1,nexp
               pwexp2(j,ind) = pwexp2(j,ind)
     1             +pwexp1(j,ind)*wshift(j)
         enddo
      enddo
c
      return
      end
c
C
c
C
      subroutine g2dcopypwexp_vec(nd,nn,pwexp1,
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
      integer nd,nn,j,j1,j2,ind,nexp
      complex *16 pwexp1(nn/2,nn,nd)
      complex *16 pwexp2(nn/2,nn,nd)

C
      do ind=1,nd
         do j1=1,nn
            do j2=1,nn/2
               pwexp2(j2,j1,ind) = pwexp1(j2,j1,ind)
            enddo
         enddo
      enddo
c
      return
      end
c
C
c
C
      subroutine g2dshiftpw0_vec(nd,delta,nn,pwexp1,
     1              cent1,pwexp2,cent2,ws,ts)
C
C     This subroutine converts the PW expansion (pwexp1) about
C     the center (CENT1) into an PW  expansion (pwexp2) about 
C     (CENT2).
C
C     INPUT
C
c     nd      = vector length (for vector input)
C     delta   = Gaussian variance
C     nn      = number of terms in PW expansion
C     pwexp1 = original expansion 
C     cent1   = center of soeall1
C     cent2   = center soeall2
C     ws,ts         = PW exp weights and nodes
C
C     OUTPUT:
C
C     pwexp2 = shifted expansion 
C
      implicit none
      integer nd,nn,j,j1,j2,ind,nexp
      real *8 cent1(2),cent2(2),x,y,delta,dsq
      complex *16 pwexp1(nn/2,nn,nd)
      complex *16 pwexp2(nn/2,nn,nd)
      complex *16, allocatable :: wshift(:,:)
      complex *16 ww1(100),ww2(100),eye
      real *8 ws(nn),ts(nn)
C
      allocate(wshift(nn/2,nn))
      eye = dcmplx(0,1)
      dsq = 1.0D0/dsqrt(delta)
      x = (cent2(1) - cent1(1))*dsq
      y = (cent2(2) - cent1(2))*dsq
      
      do j1=1,nn/2
         ww1(j1) = cdexp(eye*ts(j1)*x)
         ww2(j1) = cdexp(eye*ts(j1)*y)
         ww2(nn/2+j1) = dconjg(ww2(j1))
      enddo
c
      do j1=1,nn
         do j2=1,nn/2
            wshift(j2,j1) = ww1(j2)*ww2(j1)
         enddo
      enddo
      
      do ind=1,nd
         do j1=1,nn
            do j2=1,nn/2
               pwexp2(j2,j1,ind) = pwexp2(j2,j1,ind)
     1             +pwexp1(j2,j1,ind)*wshift(j2,j1)
            enddo
         enddo
      enddo
c
      return
      end
c
C
c
C
