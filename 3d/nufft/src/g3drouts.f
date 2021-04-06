c     This file contains the following subroutines
c
c     g3dformhc_vec : computes the Hermite expansion about the center CENT
c                     due to sources of strength charge()
c
c     g3dformhd_vec : computes the Hermite expansion about the center CENT
c                     due to dipoles
c
c     g3dformhcd_vec : computes the Hermite expansion about the center CENT
c                     due to charges and dipoles
c
C     g3dhevalp_vec : evaluates the Hermite expansion 
C                         potential only.
c      
C     g3dhevalg_vec : evaluates the Hermite expansion 
C                         pot/grad
c      
C     g3dhevalh_vec : evaluates the Hermite expansion 
C                         pot/grad/Hessian
      
c     g3dformlc_vec : computes the Taylor expansion about the center CENT
c                     due to sources of strength charge()
c
c     g3dformld_vec : computes the Taylor expansion about the center CENT
c                     due to dipoles
c
c     g3dformlcd_vec : computes the Taylor expansion about the center CENT
c                     due to charges and dipoles
c
C     g3dlevalp_vec : evaluates the local Taylor expansion 
C                         potential only.
c      
C     g3dlevalg_vec : evaluates the local Taylor expansion 
C                         pot/grad
c      
C     g3dlevalh_vec : evaluates the local Taylor expansion 
C                         pot/grad/Hessian
c
c     g3dformpwc_vec : computes the PW expansion due to charges
c                      type 1 NUFFT is called once
c      
c     g3dformpwd_vec : computes the PW expansion due to dipoles
c                      type 1 NUFFT is called three times using finufft3d1many
c      
c     g3dformpwcd_vec : computes the PW expansion due to charges and dipoles
c                      type 1 NUFFT is called four times using finufft3d1many
c      
c     g3dpwevalp_vec : evaluates the PW expansion
c                      potential only, type 2 NUFFT is called once
c                      using finufft3d2many
c      
c     g3dpwevalg_vec : evaluates the PW expansion
c                      potential + gradient
c                      4 type 2 NUFFTs are called using finufft3d2many 
c      
c     g3dpwevalh_vec : evaluates the PW expansion
c                      potential + gradient + hessian
c                      10 type 2 NUFFTs are called using finufft3d2many
c      
c     get_pwnodes : returns PW weights and nodes, midpoint rule is used
c                   so the number is always an even number!
c
c     pw_translation_matrices : returns precomputed translation matrices
c                   for PW mp to loc translations
c
c     g3dshiftpw_vec : translates PW expansion via precomputed translation 
c                      matrices
c
c     g3dcopypwexp_vec : copy PW expansion
c
c
c
c
c     Note: the plane wave expansion uses different convention here.
c
c     In the directory "planewave", the plane wave expansion is stored
c     in x, y, z in sequal, but only half of x are stored
c
c     Here, only half of z are stored, for the convenience of zero
c     padding for NUFFT calls. The order is still x, y, then z.
c
c*********************************************************************
C
C form Hermite expansions (charge, dipole, charge & dipole)
C
C*********************************************************************
      subroutine g3dformhc_vec(nd,delta,sources,ns,charge,cent,
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
      real *8 cent(3),sources(3,ns),charge(nd,ns)
      real *8 ffexp(0:nterms,0:nterms,0:nterms,nd)
      real *8 x,y,chg
      real *8, allocatable ::  xp(:),yp(:),zp(:)
C
C     initialize coefficients to zero.
C
      do ind=1,nd
         do i=0,nterms
            do k=0,nterms
               do j=0,nterms
                  ffexp(j,k,i,ind)=0
               enddo
            enddo
         enddo
      enddo
      
      allocate(xp(0:nterms))
      allocate(yp(0:nterms))
      allocate(zp(0:nterms))
c
      dsq = 1.0D0/dsqrt(delta)
C
C     accumulate expansion due to each source.
C
      do i=1,ns
         x = (sources(1,i) - cent(1))*dsq
         y = (sources(2,i) - cent(2))*dsq
         z = (sources(3,i) - cent(3))*dsq
         xp(0) = 1.0D0
         yp(0) = 1.0D0
         zp(0) = 1.0D0
         do k = 1,nterms
            tmp = dsqrt(dble(k))
            xp(k) = xp(k-1)*x/tmp
            yp(k) = yp(k-1)*y/tmp
            zp(k) = zp(k-1)*z/tmp
         enddo
c
         do ind=1,nd
            chg = charge(ind,i)
            do l=0,nterms
               ztmp = chg*zp(l)
               do k=0,nterms
                  rtmp = ztmp*yp(k)
                  do j=0,nterms
                     ffexp(j,k,l,ind) = ffexp(j,k,l,ind)+ rtmp*xp(j)
cccc                     ffexp(j,k,l,ind) = ffexp(j,k,l,ind)
cccc     1                   + charge(ind,i)*xp(j)*yp(k)*zp(l)
                  enddo
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
      subroutine g3dformhd_vec(nd,delta,sources,ns,rnormal,dipstr,
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
      implicit real *8 (a-h,o-z)
      integer ns,nterms,nd,i,j,ind,k,j1
      real *8 center(3),sources(3,ns)
      real *8 rnormal(3,ns),dipstr(nd,ns)
      real *8 ffexp(0:nterms,0:nterms,0:nterms,nd)
      real *8 x,y,delta,dsq,tmp,d1,d2
      real *8, allocatable ::  hx(:),hy(:),hz(:)
      real *8, allocatable ::  dxhx(:),dyhy(:),dzhz(:)
C
C     initialize coefficients to zero.
C
      do ind=1,nd
         do l=0,nterms
         do k=0,nterms
            do j=0,nterms
               ffexp(j,k,l,ind)=0
            enddo
         enddo
         enddo
      enddo
      allocate(hx(0:nterms))
      allocate(hy(0:nterms))
      allocate(hz(0:nterms))
      allocate(dxhx(0:nterms),dyhy(0:nterms),dzhz(0:nterms))
c
      dsq = 1.0D0/dsqrt(delta)
C
C     accumulate expansion due to each source.
C
      do i=1,ns
         x = (sources(1,i) - center(1))*dsq
         y = (sources(2,i) - center(2))*dsq
         z = (sources(3,i) - center(3))*dsq
         hx(0) = 1.0D0
         hy(0) = 1.0D0
         hz(0) = 1.0D0
         dxhx(0) = 0.0D0
         dyhy(0) = 0.0D0
         dzhz(0) = 0.0D0
         do j1 = 1,nterms
            tmp = dsqrt(dble(j1))
            hy(j1)=hy(j1-1)*y/tmp
            hx(j1)=hx(j1-1)*x/tmp
            hz(j1)=hz(j1-1)*z/tmp
            dxhx(j1) = j1*dsq*hx(j1-1)/tmp
            dyhy(j1) = j1*dsq*hy(j1-1)/tmp
            dzhz(j1) = j1*dsq*hz(j1-1)/tmp
         enddo
c
         do ind = 1,nd
            r1 = rnormal(1,i)*dipstr(ind,i)
            r2 = rnormal(2,i)*dipstr(ind,i)
            r3 = rnormal(3,i)*dipstr(ind,i)
            do l=0,nterms
               d1 = hz(l)*r1
               d2 = hz(l)*r2
               d3 = dzhz(l)*r3
               do k=0,nterms
                  c1 = hy(k)*d1
                  c2 = dyhy(k)*d2
                  c3 = hy(k)*d3
                  do j=0,nterms
                     ffexp(j,k,l,ind) = ffexp(j,k,l,ind)
     1                   +dxhx(j)*c1+hx(j)*(c2+c3)
               enddo
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
      subroutine g3dformhcd_vec(nd,delta,sources,ns,charge,rnormal,
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
      implicit real *8 (a-h,o-z)
      integer ns,nterms,nd,i,j,ind,k,j1
      real *8 center(3),sources(3,ns)
      real *8 rnormal(3,ns),dipstr(nd,ns),charge(nd,ns)
      real *8 ffexp(0:nterms,0:nterms,0:nterms,nd)
      real *8 x,y,delta,dsq,tmp,chg,ytmp,d1,d2
      real *8, allocatable ::  hx(:),hy(:),hz(:)
      real *8, allocatable ::  dxhx(:),dyhy(:),dzhz(:)
C
C     initialize coefficients to zero.
C
      do ind=1,nd
         do l=0,nterms
         do k=0,nterms
            do j=0,nterms
               ffexp(j,k,l,ind)=0
            enddo
         enddo
         enddo
      enddo

      allocate(hx(0:nterms))
      allocate(hy(0:nterms))
      allocate(hz(0:nterms))
      allocate(dxhx(0:nterms),dyhy(0:nterms),dzhz(0:nterms))
c
      dsq = 1.0D0/dsqrt(delta)
C
C     accumulate expansion due to each source.
C
      do i=1,ns
         x = (sources(1,i) - center(1))*dsq
         y = (sources(2,i) - center(2))*dsq
         z = (sources(3,i) - center(3))*dsq
         hx(0) = 1.0D0
         hy(0) = 1.0D0
         hz(0) = 1.0D0
         dxhx(0) = 0.0D0
         dyhy(0) = 0.0D0
         dzhz(0) = 0.0D0
         do j1 = 1,nterms
            tmp = dsqrt(dble(j1))
            hx(j1)=hx(j1-1)*x/tmp
            hy(j1)=hy(j1-1)*y/tmp
            hz(j1)=hz(j1-1)*z/tmp
            dxhx(j1) = j1*dsq*hx(j1-1)/tmp
            dyhy(j1) = j1*dsq*hy(j1-1)/tmp
            dzhz(j1) = j1*dsq*hz(j1-1)/tmp
         enddo
c
         do ind = 1,nd
            chg = charge(ind,i)
            r1 = rnormal(1,i)*dipstr(ind,i)
            r2 = rnormal(2,i)*dipstr(ind,i)
            r3 = rnormal(3,i)*dipstr(ind,i)
            
            do l=0,nterms
               ztmp=chg*hz(l)

               d1 = hz(l)*r1
               d2 = hz(l)*r2
               d3 = dzhz(l)*r3
               
               do k=0,nterms
                  rtmp=ztmp*hy(k)

                  c1 = hy(k)*d1
                  c2 = dyhy(k)*d2
                  c3 = hy(k)*d3
                  do j=0,nterms
                     ffexp(j,k,l,ind) = ffexp(j,k,l,ind)+
     1                   dxhx(j)*c1+hx(j)*(c2+c3+rtmp)
               enddo
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
      subroutine g3dhevalp_vec(nd,delta,cent,nterms,ffexp,
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
      real *8 cent(3),pot(nd,ntarg)
      real *8 ffexp(0:nterms,0:nterms,0:nterms,nd)
      real *8 x,y,targ(3,ntarg)
      real *8, allocatable ::  hexpx(:),hexpy(:),hexpz(:)
C
C     initialize coefficients to zero.
C
      allocate(hexpx(0:nterms))
      allocate(hexpy(0:nterms))
      allocate(hexpz(0:nterms))
C
C     tabulate Hermite polynomials in HEXPX and HEXPY.
C
      dsq = 1.0d0/dsqrt(delta)
      do itarg=1,ntarg
         x = (targ(1,itarg) - cent(1))*dsq
         y = (targ(2,itarg) - cent(2))*dsq
         z = (targ(3,itarg) - cent(3))*dsq
         facx = dexp(-x*x)
         facy = dexp(-y*y)
         facz = dexp(-z*z)
         hexpx(0) = 1.0d0*facx
         hexpy(0) = 1.0d0*facy
         hexpz(0) = 1.0d0*facz
         hexpx(1) = 2.0d0*x*facx
         hexpy(1) = 2.0d0*y*facy
         hexpz(1) = 2.0d0*z*facz
         do i = 1,nterms-1
            hexpx(i+1)=2.0d0*(x/dsqrt((i+1)*1.0d0)*hexpx(i)
     1          -dsqrt(1.0d0*i/(i+1))*hexpx(i-1))
            hexpy(i+1)=2.0d0*(y/dsqrt((i+1)*1.0d0)*hexpy(i)
     1          -dsqrt(1.0d0*i/(i+1))*hexpy(i-1))
            hexpz(i+1)=2.0d0*(z/dsqrt((i+1)*1.0d0)*hexpz(i)
     1          -dsqrt(1.0d0*i/(i+1))*hexpz(i-1))
         enddo
c
C---- EVALUATE EXPANSION
c
         do ind = 1,nd
            dd=0
            do l = 0,nterms
               d1=0
               do k = 0,nterms
                  d2=0
                  do j = 0,nterms
                     d2=d2+hexpx(j)*ffexp(j,k,l,ind)
                  enddo
                  d1=d1+d2*hexpy(k)
               enddo
               dd = dd + d1*hexpz(l)
            enddo
            pot(ind,itarg) = pot(ind,itarg) + dd
         enddo
      enddo
      
      return
      end
C
c
C***********************************************************************
      subroutine g3dhevalg_vec(nd,delta,cent,nterms,ffexp,
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
      real *8 cent(3),pot(nd,ntarg),grad(nd,3,ntarg)
      real *8 ffexp(0:nterms,0:nterms,0:nterms,nd)
      real *8 x,y,z,targ(3,ntarg)
      real *8, allocatable ::  hexpx(:),hexpy(:),hexpz(:)
      real *8, allocatable ::  dhexpx(:),dhexpy(:),dhexpz(:)
C
C     initialize coefficients to zero.
C
      allocate(hexpx(0:nterms+1))
      allocate(dhexpx(0:nterms))
      allocate(hexpy(0:nterms+1))
      allocate(dhexpy(0:nterms))
      allocate(hexpz(0:nterms+1))
      allocate(dhexpz(0:nterms))
C
C     tabulate Hermite polynomials in HEXPX and HEXPY.
C
      dsq = 1.0d0/dsqrt(delta)
      do itarg=1,ntarg
         x = (targ(1,itarg) - cent(1))*dsq
         y = (targ(2,itarg) - cent(2))*dsq
         z = (targ(3,itarg) - cent(3))*dsq
         facx = dexp(-x*x)
         facy = dexp(-y*y)
         facz = dexp(-z*z)
         hexpx(0) = 1.0d0*facx
         hexpy(0) = 1.0d0*facy
         hexpz(0) = 1.0d0*facz
         hexpx(1) = 2.0d0*x*facx
         hexpy(1) = 2.0d0*y*facy
         hexpz(1) = 2.0d0*z*facz
         do i = 1,nterms
            hexpx(i+1)=2.0d0*(x/dsqrt((i+1)*1.0d0)*hexpx(i)
     1          -dsqrt(1.0d0*i/(i+1))*hexpx(i-1))
            hexpy(i+1)=2.0d0*(y/dsqrt((i+1)*1.0d0)*hexpy(i)
     1          -dsqrt(1.0d0*i/(i+1))*hexpy(i-1))
            hexpz(i+1)=2.0d0*(z/dsqrt((i+1)*1.0d0)*hexpz(i)
     1          -dsqrt(1.0d0*i/(i+1))*hexpz(i-1))
         enddo
         do i = 0,nterms
            dhexpx(i)=-hexpx(i+1)*dsqrt(i+1.0d0)
            dhexpy(i)=-hexpy(i+1)*dsqrt(i+1.0d0)
            dhexpz(i)=-hexpz(i+1)*dsqrt(i+1.0d0)
         enddo
c
C---- EVALUATE EXPANSION
c
         do ind = 1,nd
            g = 0
            gx = 0.0d0
            gy = 0.0d0
            gz = 0.0d0
            do l = 0,nterms
               d = 0
               dx = 0
               dy = 0
               dz = 0
               do k = 0,nterms
                  c=0
                  cx=0
                  cy=0
                  cz=0
                  do j = 0,nterms
                     c=c+hexpx(j)*ffexp(j,k,l,ind)
                     cx = cx + dhexpx(j)*ffexp(j,k,l,ind)
                     cy = cy + hexpx(j)*ffexp(j,k,l,ind)
                     cz = cz + hexpx(j)*ffexp(j,k,l,ind)
                  enddo
                  d=d+c*hexpy(k)
                  dx=dx+cx*hexpy(k)
                  dy=dy+cy*dhexpy(k)
                  dz=dz+cz*hexpy(k)
               enddo
               g=g+d*hexpz(l)
               gx=gx+dx*hexpz(l)
               gy=gy+dy*hexpz(l)
               gz=gz+dz*dhexpz(l)
            enddo
            pot(ind,itarg) = pot(ind,itarg) + g
            grad(ind,1,itarg) = grad(ind,1,itarg) + gx*dsq
            grad(ind,2,itarg) = grad(ind,2,itarg) + gy*dsq
            grad(ind,3,itarg) = grad(ind,3,itarg) + gz*dsq
         enddo
      enddo
      
      return
      end
C
C
c
C***********************************************************************
      subroutine g3dhevalh_vec(nd,delta,cent,nterms,ffexp,
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
      real *8 cent(3),pot(nd,ntarg),grad(nd,3,ntarg),hess(nd,6,ntarg)
      real *8 ffexp(0:nterms,0:nterms,0:nterms,nd)
      real *8 x,y,targ(3,ntarg)
      real *8, allocatable ::  hexpx(:),hexpy(:),hexpz(:)
      real *8, allocatable ::  dhexpx(:),dhexpy(:),dhexpz(:)
      real *8, allocatable ::  ddhexpx(:),ddhexpy(:),ddhexpz(:)
C
C     initialize coefficients to zero.
C
      allocate(hexpx(0:nterms+2))
      allocate(dhexpx(0:nterms))
      allocate(ddhexpx(0:nterms))
      allocate(hexpy(0:nterms+2))
      allocate(dhexpy(0:nterms))
      allocate(ddhexpy(0:nterms))
      allocate(hexpz(0:nterms+2))
      allocate(dhexpz(0:nterms))
      allocate(ddhexpz(0:nterms))
C
C     tabulate Hermite polynomials in HEXPX and HEXPY.
C
      dsq = 1.0d0/dsqrt(delta)
      do itarg=1,ntarg
         x = (targ(1,itarg) - cent(1))*dsq
         y = (targ(2,itarg) - cent(2))*dsq
         y = (targ(3,itarg) - cent(3))*dsq
         facx = dexp(-x*x)
         facy = dexp(-y*y)
         facz = dexp(-z*z)
         hexpx(0) = 1.0d0*facx
         hexpy(0) = 1.0d0*facy
         hexpz(0) = 1.0d0*facz
         hexpx(1) = 2.0d0*x*facx
         hexpy(1) = 2.0d0*y*facy
         hexpz(1) = 2.0d0*z*facz
         do i = 1,nterms+1
            hexpx(i+1)=2.0d0*(x/dsqrt((i+1)*1.0d0)*hexpx(i)
     1          -dsqrt(1.0d0*i/(i+1))*hexpx(i-1))
            hexpy(i+1)=2.0d0*(y/dsqrt((i+1)*1.0d0)*hexpy(i)
     1          -dsqrt(1.0d0*i/(i+1))*hexpy(i-1))
            hexpz(i+1)=2.0d0*(z/dsqrt((i+1)*1.0d0)*hexpz(i)
     1          -dsqrt(1.0d0*i/(i+1))*hexpz(i-1))
         enddo
         do i = 0,nterms
            dhexpx(i)=-hexpx(i+1)*dsqrt(i+1.0d0)
            ddhexpx(i)=hexpx(i+2)*dsqrt(i+1.0d0)*dsqrt(i+2.0d0)
            dhexpy(i)=-hexpy(i+1)*dsqrt(i+1.0d0)
            ddhexpy(i)=hexpy(i+2)*dsqrt(i+1.0d0)*dsqrt(i+2.0d0)
            dhexpz(i)=-hexpz(i+1)*dsqrt(i+1.0d0)
            ddhexpz(i)=hexpz(i+2)*dsqrt(i+1.0d0)*dsqrt(i+2.0d0)
         enddo
c
c
C---- EVALUATE EXPANSION
c
         do ind = 1,nd
            g=0
            gx = 0.0d0
            gy = 0.0d0
            gz = 0.0d0
            hxx = 0.0d0
            hyy = 0.0d0
            hzz = 0.0d0
            hxy = 0.0d0
            hxz = 0.0d0
            hyz = 0.0d0
            do l=0,nterms
               d = 0
               dx = 0
               dy = 0
               dz = 0
               dxx = 0
               dyy = 0
               dzz = 0
               dxy = 0
               dxz = 0
               dzz = 0
               do k = 0,nterms
                  c=0
                  cx=0
                  cy=0
                  cz=0
                  cxx=0
                  cyy=0
                  czz=0
                  cxy=0
                  cxz=0
                  czz=0
                  do j = 0,nterms
                     c  = c  +   hexpx(j)*ffexp(j,k,l,ind)
                     cx  = cx  +  dhexpx(j)*ffexp(j,k,l,ind)
                     cy  = cy  +   hexpx(j)*ffexp(j,k,l,ind)               
                     cz  = cz  +   hexpx(j)*ffexp(j,k,l,ind)
                     
                     cxx = cxx + ddhexpx(j)*ffexp(j,k,l,ind)
                     cyy = cyy +   hexpx(j)*ffexp(j,k,l,ind)
                     czz = czz +   hexpx(j)*ffexp(j,k,l,ind)
                     cxy = cxy +  dhexpx(j)*ffexp(j,k,l,ind)
                     cxz = cxz +  dhexpx(j)*ffexp(j,k,l,ind)
                     cyz = cyz +   hexpx(j)*ffexp(j,k,l,ind)
                  enddo
                  d   = d   + c    *hexpy(k)
                  dx  = dx  + cx   *hexpy(k)
                  dy  = dy  + cy  *dhexpy(k)
                  dz  = dz  + cz   *hexpy(k)
                  dxx = dxx + cxx  *hexpy(k)
                  dyy = dyy + cyy*ddhexpy(k)
                  dzz = dzz + czz  *hexpy(k)
                  dxy = dxy + cxy *dhexpy(k)
                  dxz = dxz + cxz  *hexpy(k)
                  dyz = dyz + cyz *dhexpy(k)
               enddo
               g = g+d*hexpz(l)
               gx=gz+dx*hexpz(l)
               gy=gy+dy*hexpz(l)
               gz=gz+dz*hexpz(l)

               hxx=hxx+dxx*hexpz(l)
               hyy=hyy+dyy*hexpz(l)
               hzz=hzz+dzz*ddhexpz(l)

               hxy=hxy+dxy*hexpz(l)
               hxz=hxz+dxz*dhexpz(l)
               hyz=hyz+dzz*dhexpz(l)
            enddo
         
            pot(ind,itarg)    = pot(ind,itarg) + g
            grad(ind,1,itarg) = grad(ind,1,itarg) + gx*dsq
            grad(ind,2,itarg) = grad(ind,2,itarg) + gy*dsq
            grad(ind,3,itarg) = grad(ind,3,itarg) + gz*dsq
            
            hess(ind,1,itarg) = hess(ind,1,itarg) + hxx*dsq*dsq
            hess(ind,2,itarg) = hess(ind,2,itarg) + hxy*dsq*dsq
            hess(ind,3,itarg) = hess(ind,3,itarg) + hzz*dsq*dsq
            
            hess(ind,4,itarg) = hess(ind,4,itarg) + hxy*dsq*dsq
            hess(ind,5,itarg) = hess(ind,5,itarg) + hxz*dsq*dsq
            hess(ind,6,itarg) = hess(ind,6,itarg) + hyz*dsq*dsq
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
      subroutine g3dformlc_vec(nd,delta,sources,ns,charge,center,
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
      real *8 center(3),sources(3,ns),charge(nd,ns)
      real *8 local(0:nlocal,0:nlocal,0:nlocal,nd)
      real *8 dd, d2
      real *8, allocatable :: hexpx(:),hexpy(:),hexpz(:)
c
      allocate(hexpx(0:nlocal),hexpy(0:nlocal),hexpz(0:nlocal))
c
      dsq = 1.0D0/dsqrt(delta)
C
      do i=1,ns
         x = (sources(1,i) - center(1))*dsq
         y = (sources(2,i) - center(2))*dsq
         y = (sources(3,i) - center(3))*dsq
c

         facx=dexp(-x*x)
         facy=dexp(-y*y)
         facz=dexp(-z*z)

         hexpx(0)=1.0d0*facx
         hexpy(0)=1.0d0*facy
         hexpz(0)=1.0d0*facz
         hexpx(1)=2.0d0*x*facx
         hexpy(1)=2.0d0*y*facy
         hexpz(1)=2.0d0*z*facz

         do j=1,nlocal-1
            hexpx(j+1)=2.0d0*(x*hexpx(j)-hexpx(j-1))/(j+1)
            hexpy(j+1)=2.0d0*(y*hexpy(j)-hexpy(j-1))/(j+1)
            hexpz(j+1)=2.0d0*(z*hexpz(j)-hexpz(j-1))/(j+1)
         enddo

         do ind=1,nd
            chg=charge(ind,i)
            do j1=0,nlocal
               ztmp = chg*hexpz(j1) 
               do j2=0,nlocal
                  ytmp=ztmp*hexpy(j2)
                  do j3=0,nlocal
                     local(j3,j2,j1,ind) = local(j3,j2,j1,ind)
     1                   +ytmp*hexpx(j3)
                  enddo
               enddo
            enddo
            
         enddo
      enddo
      return
      end
C
C
C
      subroutine g3dformld_vec(nd,delta,sources,ns,rnormal,dipstr,
     1            center,nlocal,local)
C
C     This subroutine computes the Taylor expansions about
C     the center CENT due to sources of strength charge().
C
C     INPUT
C
c     nd            = vector length (parallel)
C     delta         = Gaussian variance
C     sources(3,ns) = coordinates of sources
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
      real *8 center(3),sources(3,ns),dipstr(nd,ns)
      real *8 rnormal(3,ns)
      real *8 local(0:nlocal,0:nlocal,0:nlocal,nd)
      real *8 dd, d2
      real *8, allocatable :: hexpx(:),hexpy(:),hexpz(:)
      real *8, allocatable :: dhexpx(:),dhexpy(:),dhexpz(:)
c
      allocate(hexpx(0:nlocal+1),hexpy(0:nlocal+1),hexpz(0:nlocal+1))
      allocate(dhexpx(0:nlocal),dhexpy(0:nlocal),dhexpz(0:nlocal))
c
      dsq = 1.0D0/dsqrt(delta)
C
      do i=1,ns
         x = (sources(1,i) - center(1))*dsq
         y = (sources(2,i) - center(2))*dsq
         z = (sources(3,i) - center(3))*dsq
c

         facx=dexp(-x*x)
         facy=dexp(-y*y)
         facz=dexp(-z*z)

         hexpx(0)=1.0d0*facx
         hexpy(0)=1.0d0*facy
         hexpz(0)=1.0d0*facz
         hexpx(1)=2.0d0*x*facx
         hexpy(1)=2.0d0*y*facy
         hexpz(1)=2.0d0*z*facz

         do j=1,nlocal
            hexpx(j+1)=2.0d0*(x*hexpx(j)-hexpx(j-1))/(j+1)
            hexpy(j+1)=2.0d0*(y*hexpy(j)-hexpy(j-1))/(j+1)
            hexpz(j+1)=2.0d0*(z*hexpz(j)-hexpz(j-1))/(j+1)
         enddo
         do j = 0,nlocal
            dhexpx(j)=-hexpx(j+1)*dsq*(j+1.0d0)
            dhexpy(j)=-hexpy(j+1)*dsq*(j+1.0d0)
            dhexpz(j)=-hexpz(j+1)*dsq*(j+1.0d0)
         enddo

         do ind=1,nd
            r1 = rnormal(1,i)*dipstr(ind,i)
            r2 = rnormal(2,i)*dipstr(ind,i)
            r3 = rnormal(3,i)*dipstr(ind,i)
            do j1=0,nlocal
               rz1 = hexpz(j1)*r1
               rz2 = hexpz(j1)*r2
               rz3 = dhexpz(j1)*r3
               do j2=0,nlocal
                  ry1 = rz1*hexpy(j2)
                  ry2 = rz2*dhexpy(j2)
                  ry3 = rz3*hexpy(j2)
                  do j3=0,nlocal
                     local(j3,j2,j1,ind) = local(j3,j2,j1,ind)+
     1                   dhexpx(j2)*ry1+hexpx(j2)*(ry2+ry3)
                  enddo
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
      subroutine g3dformlcd_vec(nd,delta,sources,ns,charge,rnormal,
     1            dipstr,center,nlocal,local)
C
C     This subroutine computes the Taylor expansions about
C     the center CENT due to sources of strength charge().
C
C     INPUT
C
c     nd            = vector length (parallel)
C     delta         = Gaussian variance
C     sources(3,ns) = coordinates of sources
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
      real *8 center(3),sources(3,ns),dipstr(nd,ns)
      real *8 rnormal(3,ns), charge(nd,ns)
      real *8 local(0:nlocal,0:nlocal,0:nlocal,nd)
      real *8 dd, d2
      real *8, allocatable :: hexpx(:),hexpy(:),hexpz(:)
      real *8, allocatable :: dhexpx(:),dhexpy(:),dhexpz(:)
c
      allocate(hexpx(0:nlocal+1),hexpy(0:nlocal+1),hexpz(0:nlocal+1))
      allocate(dhexpx(0:nlocal),dhexpy(0:nlocal),dhexpz(0:nlocal))
c
      dsq = 1.0D0/dsqrt(delta)
C
      do i=1,ns
         x = (sources(1,i) - center(1))*dsq
         y = (sources(2,i) - center(2))*dsq
         z = (sources(3,i) - center(3))*dsq
c

         facx=dexp(-x*x)
         facy=dexp(-y*y)
         facz=dexp(-z*z)

         hexpx(0)=1.0d0*facx
         hexpy(0)=1.0d0*facy
         hexpz(0)=1.0d0*facz
         hexpx(1)=2.0d0*x*facx
         hexpy(1)=2.0d0*y*facy
         hexpz(1)=2.0d0*z*facz

         do j=1,nlocal
            hexpx(j+1)=2.0d0*(x*hexpx(j)-hexpx(j-1))/(j+1)
            hexpy(j+1)=2.0d0*(y*hexpy(j)-hexpy(j-1))/(j+1)
            hexpz(j+1)=2.0d0*(z*hexpz(j)-hexpz(j-1))/(j+1)
         enddo
         do j = 0,nlocal
            dhexpx(j)=-hexpx(j+1)*dsq*(j+1.0d0)
            dhexpy(j)=-hexpy(j+1)*dsq*(j+1.0d0)
            dhexpz(j)=-hexpz(j+1)*dsq*(j+1.0d0)
         enddo

         do ind=1,nd
            chg=charge(ind,i)
            r1 = rnormal(1,i)*dipstr(ind,i)
            r2 = rnormal(2,i)*dipstr(ind,i)
            r3 = rnormal(3,i)*dipstr(ind,i)
            
            do j1=0,nlocal
               ztmp=chg*hexpz(j1)
               dz1 = hexpz(j1)*r1
               dz2 = hexpz(j1)*r2
               dz3 = dhexpz(j1)*r3
               do j2=0,nlocal
                  ytmp=ztmp*hexpy(j2)
                  dy1=dz1*hexpy(j2)
                  dy2=dz2*dhexpy(j2)
                  dy3=dz3*hexpy(j2)
                  do j3=0,nlocal
                     local(j3,j2,j1,ind) = local(j3,j2,j1,ind)+
     1                   dhexpx(j2)*dy1+hexpx(j2)*(dy2+dy3+ytmp)
               enddo
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
      subroutine g3dlevalp_vec(nd,delta,center,nlocal,local,
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
      real *8 local(0:nlocal,0:nlocal,0:nlocal,nd)
      real *8 center(3),targ(3,ntarg)
      real *8 pot(nd,ntarg),xp(0:200),yp(0:200),zp(0:200)
C
      dsq = 1.0D0/dsqrt(delta)
C
      do itarg = 1,ntarg
c
         x = (targ(1,itarg) - center(1))*dsq
         y = (targ(2,itarg) - center(2))*dsq
         z = (targ(3,itarg) - center(3))*dsq

         xp(0)=1.0d0
         yp(0)=1.0d0
         zp(0)=1.0d0
         do i=1,nlocal
            xp(i)=x*xp(i-1)
            yp(i)=y*yp(i-1)
            zp(i)=z*zp(i-1)
         enddo
c
         do ind=1,nd
            dd=0
            do j1=0,nlocal
               d1=0
               do j2=0,nlocal
                  d2=0
                  do j3=0,nlocal
                     d2=d2+local(j3,j2,j1,ind)*xp(j3)
cccc                     dd=dd+local(j3,j2,j1,ind)*xp(j3)
cccc     1                   *yp(j2)*zp(j1)
                  enddo
                  d1=d1+d2*yp(j2)
               enddo
               dd=dd+d1*zp(j1)
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
      subroutine g3dlevalg_vec(nd,delta,center,nlocal,local,
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
      real *8 local(0:nlocal,0:nlocal,0:nlocal,nd)
      real *8 center(3),targ(3,ntarg)
      real *8 pot(nd,ntarg)
      real *8 grad(nd,3,ntarg)
      real *8, allocatable :: xp(:),yp(:),zp(:)
      real *8, allocatable :: xpx(:),ypy(:),zpz(:)
C
      allocate(xp(0:nlocal))
      allocate(yp(0:nlocal))
      allocate(zp(0:nlocal))
      allocate(xpx(0:nlocal))
      allocate(ypy(0:nlocal))
      allocate(zpz(0:nlocal))
C
      dsq = 1.0D0/dsqrt(delta)
C
      do itarg = 1,ntarg
c
         x = (targ(1,itarg) - center(1))*dsq
         y = (targ(2,itarg) - center(2))*dsq
         z = (targ(3,itarg) - center(3))*dsq

         xp(0)=1.0d0
         yp(0)=1.0d0
         zp(0)=1.0d0
         do i=1,nlocal
            xp(i)=x*xp(i-1)
            yp(i)=y*yp(i-1)
            zp(i)=z*zp(i-1)
         enddo
c
         xpx(0) = 0.0d0
         ypy(0) = 0.0d0
         zpz(0) = 0.0d0
         do i = 1,nlocal
            xpx(i) = i*xp(i-1)*dsq
            ypy(i) = i*yp(i-1)*dsq
            zpz(i) = i*zp(i-1)*dsq
         enddo
c
         do ind=1,nd
            dd=0
            g1=0
            g2=0
            g3=0
            do j1=0,nlocal
               d1=0
               d1x=0
               d1y=0
               do j2=0,nlocal
                  d2=0
                  d2x=0
                  do j3=0,nlocal
                     d2=d2+local(j3,j2,j1,ind)*xp(j3)
                     d2x=d2x+local(j3,j2,j1,ind)*xpx(j3)
                  enddo
                  d1=d1+d2*yp(j2)
                  d1x=d1x+d2x*yp(j2)
                  d1y=d1y+d2*ypy(j2)
               enddo
               dd=dd+d1*zp(j1)
               g1=g1+d1x*zp(j1)
               g2=g2+d1y*zp(j1)
               g3=g3+d1*zpz(j1)
            enddo
            pot(ind,itarg) = pot(ind,itarg)+dd
            grad(ind,1,itarg) = grad(ind,1,itarg)+g1
            grad(ind,2,itarg) = grad(ind,2,itarg)+g2
            grad(ind,3,itarg) = grad(ind,3,itarg)+g3
         enddo
      enddo

      return
      end
C
C
C
C
      subroutine g3dlevalh_vec(nd,delta,center,nlocal,local,
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
      real *8 local(0:nlocal,0:nlocal,0:nlocal,nd)
      real *8 center(3),targ(3,ntarg)
      real *8 pot(nd,ntarg)
      real *8 grad(nd,3,ntarg)
      real *8 hess(nd,6,ntarg)
      real *8, allocatable :: xp(:),yp(:),zp(:)
      real *8, allocatable :: xpx(:),ypy(:),zpz(:)
      real *8, allocatable :: xpxx(:),ypyy(:),zpzz(:)
C
      allocate(xp(0:nlocal))
      allocate(yp(0:nlocal))
      allocate(zp(0:nlocal))
      
      allocate(xpx(0:nlocal))
      allocate(ypy(0:nlocal))
      allocate(zpz(0:nlocal))

      allocate(xpxx(0:nlocal))
      allocate(ypyy(0:nlocal))
      allocate(zpzz(0:nlocal))
C
      dsq = 1.0D0/dsqrt(delta)
C
      do itarg = 1,ntarg
c
         x = (targ(1,itarg) - center(1))*dsq
         y = (targ(2,itarg) - center(2))*dsq
         z = (targ(3,itarg) - center(3))*dsq

         xp(0)=1.0d0
         yp(0)=1.0d0
         zp(0)=1.0d0
         do i=1,nlocal
            xp(i)=x*xp(i-1)
            yp(i)=y*yp(i-1)
            zp(i)=z*zp(i-1)
         enddo
c
         xpx(0) = 0.0d0
         ypy(0) = 0.0d0
         zpz(0) = 0.0d0
         do i = 1,nlocal
            xpx(i) = i*xp(i-1)*dsq
            ypy(i) = i*yp(i-1)*dsq
            zpz(i) = i*zp(i-1)*dsq
         enddo
c
         xpxx(0) = 0.0d0
         ypyy(0) = 0.0d0
         zpzz(0) = 0.0d0
         xpxx(1) = 0.0d0
         ypyy(1) = 0.0d0
         zpzz(1) = 0.0d0
         do i = 2,nlocal
            xpxx(i) = i*(i-1)*xp(i-2)*dsq*dsq
            ypyy(i) = i*(i-1)*yp(i-2)*dsq*dsq
            zpzz(i) = i*(i-1)*zp(i-2)*dsq*dsq
         enddo

         do ind=1,nd
            dd=0
            g1=0
            g2=0
            g3=0
            h11=0
            h22=0
            h33=0
            h12=0
            h13=0
            h23=0
            do j1=0,nlocal
               d1=0
               d1x=0
               d1y=0

               d1xx=0
               d1yy=0
               d1xy=0
               do j2=0,nlocal
                  d2=0
                  d2x=0
                  d2xx=0
                  do j3=0,nlocal
                     d2=d2+local(j3,j2,j1,ind)*xp(j3)
                     d2x=d2x+local(j3,j2,j1,ind)*xpx(j3)
                     d2xx=d2xx+local(j3,j2,j1,ind)*xpxx(j3)
                  enddo
                  d1=d1+d2*yp(j2)
                  
                  d1x=d1x+d2x*yp(j2)
                  d1y=d1y+d2*ypy(j2)
                  
                  d1xy=d1xy+d2x*ypy(j2)
                  
                  d1xx=d1xx+d2xx*yp(j2)
                  d1yy=d1yy+d2*ypyy(j2)
               enddo
               dd=dd+d1*zp(j1)
               
               g1=g1+d1x*zp(j1)
               g2=g2+d1y*zp(j1)
               g3=g3+d1*zpz(j1)
               
               h11=h11+d1xx*zp(j1)
               h22=h22+d1yy*zp(j1)
               h33=h33+d1*zpzz(j1)
               
               h12=h12+d1xy*zp(j1)
               h13=h13+d1x*zpz(j1)
               h23=h23+d1y*zpz(j1)
            enddo
            pot(ind,itarg) = pot(ind,itarg)+dd
            
            grad(ind,1,itarg) = grad(ind,1,itarg)+g1
            grad(ind,2,itarg) = grad(ind,2,itarg)+g2
            grad(ind,3,itarg) = grad(ind,3,itarg)+g3

            hess(ind,1,itarg) = hess(ind,1,itarg)+h11
            hess(ind,2,itarg) = hess(ind,2,itarg)+h22
            hess(ind,3,itarg) = hess(ind,3,itarg)+h33

            hess(ind,4,itarg) = hess(ind,4,itarg)+h12
            hess(ind,5,itarg) = hess(ind,5,itarg)+h13
            hess(ind,6,itarg) = hess(ind,6,itarg)+h23
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
      subroutine g3dherm2local_vec(nd,nlocal,ntermsh,h2l,
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
      real *8 hexp(0:ntermsh,0:ntermsh,0:ntermsh,nd)
      real *8 local(0:nlocal,0:nlocal,0:nlocal,nd)
      
      real *8, allocatable :: lxhyz(:,:,:),lxyhz(:,:,:)
c      
c
      allocate(lxhyz(0:nlocal,0:ntermsh,0:ntermsh))
      allocate(lxyhz(0:nlocal,0:nlocal,0:ntermsh))

      do ind=1,nd
c        transform in x
         do k1=0,nlocal,2
            do j3=0,ntermsh
            do j2=0,ntermsh
               cd=0
               do j1=0,ntermsh,2
                  cd=cd+h2l(j1,k1)*hexp(j1,j2,j3,ind)
               enddo
               lxhyz(k1,j2,j3)=cd
            enddo
            enddo
         enddo

         do k1=1,nlocal,2
            do j3=0,ntermsh
            do j2=0,ntermsh
               cd=0
               do j1=1,ntermsh,2
                  cd=cd+h2l(j1,k1)*hexp(j1,j2,j3,ind)
               enddo
               lxhyz(k1,j2,j3)=cd
            enddo
            enddo
         enddo

c        transform in y         
         do k1=0,nlocal,2
            do k2=0,nlocal
               do j3=0,ntermsh
                  cd=0
                  do j2=0,ntermsh,2
                     cd=cd+h2l(j2,k1)*lxhyz(k2,j2,j3)
                  enddo
                  lxyhz(k2,k1,j3)=cd
               enddo
            enddo
         enddo
         do k1=1,nlocal,2
            do k2=0,nlocal
               do j3=0,ntermsh
                  cd=0
                  do j2=1,ntermsh,2
                     cd=cd+h2l(j2,k1)*lxhyz(k2,j2,j3)
                  enddo
                  lxyhz(k2,k1,j3)=cd
               enddo
            enddo
         enddo
         

c        transfrom in z         
         do k1=0,nlocal,2
            do k2=0,nlocal
               do k3=0,nlocal
                  cd=0
                  do j2=0,ntermsh,2
                     cd=cd+h2l(j2,k1)*lxyhz(k3,k2,j2)
                  enddo
                  local(k3,k2,k1,ind)=local(k3,k2,k1,ind)+cd
               enddo
            enddo
         enddo
         do k1=1,nlocal,2
            do k2=0,nlocal
               do k3=0,nlocal
                  cd=0
                  do j2=1,ntermsh,2
                     cd=cd+h2l(j2,k1)*lxyhz(k3,k2,j2)
                  enddo
                  local(k3,k2,k1,ind)=local(k3,k2,k1,ind)+cd
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
C***********************************************************************
      subroutine g3dhermzero_vec(nd,hexp,ntermsh)
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
      real *8 hexp((ntermsh+1)*(ntermsh+1)*(ntermsh+1),nd)
c
      do ii=1,nd
      do n=1,(ntermsh+1)*(ntermsh+1)*(ntermsh+1)
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
      subroutine g3dlocalzero_vec(nd,local,nlocal)
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
      real *8 local((nlocal+1)*(nlocal+1)*(nlocal+1),nd)
c
      do ii=1,nd
      do n=1,(nlocal+1)*(nlocal+1)*(nlocal+1)
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
      subroutine g3dformpwc_vec(nd,delta,eps,sources,ns,charge,
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
      real *8 cent(3),sources(3,ns),charge(nd,ns)
      real *8 ws(-npw/2:npw/2-1), ts(-npw/2:npw/2-1)
      real *8 wnufft(nexp)
      complex *16 ffexp(nexp,nd)
      real *8 x,y,chg
      real *8, allocatable ::  xj(:),yj(:),zj(:)
      complex *16, allocatable :: cj(:,:),fk(:,:,:,:),wj(:)
c     this (since unallocated) used to pass a NULL ptr to FINUFFT...
      integer*8, allocatable :: null
      integer*8 ns8, npw8
      complex *16 eye
      eye = dcmplx(0,1)

      npw2=npw/2
      
      allocate(xj(ns))
      allocate(yj(ns))
      allocate(zj(ns))
      allocate(cj(ns,nd))
      allocate(wj(ns))
      allocate(fk(-npw2:npw2-1,-npw2:npw2-1,-npw2:npw2-1,nd))
c
      dsq = 2*ts(0)/dsqrt(delta)
C
C     accumulate expansion due to each source.
C
      do j=1,ns
         xj(j) = (sources(1,j) - cent(1))*dsq
         yj(j) = (sources(2,j) - cent(2))*dsq
         zj(j) = (sources(3,j) - cent(3))*dsq
         wj(j) = exp(-0.5d0*eye*(xj(j)+yj(j)+zj(j)))
      enddo

      do ind = 1,nd
         do j=1,ns
            cj(j,ind) = charge(ind,j)*wj(j)
         enddo
      enddo
      
      ns8=ns
      npw8=npw
      iflag = -1
      call finufft3d1many(nd,ns8,xj,yj,zj,cj,iflag,eps,npw8,npw8,npw8,
     1       fk,null,ier)
cccc  print *, 'ier=', ier
      do ind=1,nd
         j=0
         do k3=-npw2,-1
            do k2=-npw2,npw2-1
               do k1=-npw2,npw2-1
                  j=j+1
                  ffexp(j,ind)=fk(k1,k2,k3,ind)*wnufft(j)
               enddo
            enddo
         enddo
      enddo

      return
      end
c
C
C
C
      subroutine g3dformpwd_vec(nd,delta,eps,sources,ns,rnormal,
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
      real *8 cent(3),sources(3,ns),dipstr(nd,ns)
      real *8 rnormal(3,ns)
      real *8 ws(-npw/2:npw/2-1), ts(-npw/2:npw/2-1)
      real *8 wnufft(nexp)
      complex *16 ffexp(nexp,nd)
      real *8 x,y,chg
      real *8, allocatable ::  xj(:),yj(:),zj(:)
      complex *16, allocatable :: cj(:,:,:),fk(:,:,:,:,:),wj(:)
c     this (since unallocated) used to pass a NULL ptr to FINUFFT...
      integer*8, allocatable :: null
      integer*8 ns8, npw8
      complex *16 eye,ztmp
      eye = dcmplx(0,1)

      npw2=npw/2
      
      allocate(xj(ns))
      allocate(yj(ns))
      allocate(zj(ns))
      allocate(cj(ns,3,nd))
      allocate(wj(ns))
      allocate(fk(-npw2:npw2-1,-npw2:npw2-1,-npw2:npw2-1,3,nd))
c
      dsq0 = 1/dsqrt(delta)
      dsq = 2*ts(0)/dsqrt(delta)
C
C     accumulate expansion due to each source.
C
      do j=1,ns
         xj(j) = (sources(1,j) - cent(1))*dsq
         yj(j) = (sources(2,j) - cent(2))*dsq
         zj(j) = (sources(3,j) - cent(3))*dsq
         wj(j) = exp(-0.5d0*eye*(xj(j)+yj(j)+zj(j)))
      enddo
      
      do ind = 1,nd
         do j=1,ns
            ztmp = -eye*wj(j)*(dipstr(ind,j)*dsq0)
            cj(j,1,ind) = ztmp*rnormal(1,j)
            cj(j,2,ind) = ztmp*rnormal(2,j)
            cj(j,3,ind) = ztmp*rnormal(3,j)
         enddo
      enddo
      
      ns8=ns
      npw8=npw
      iflag = -1
c     total number of NUFFTs
      ntrans = 3*nd
      call finufft3d1many(ntrans,ns8,xj,yj,zj,cj,iflag,eps,
     1       npw8,npw8,npw8,fk,null,ier)
cccc  print *, 'ier=', ier
      do ind=1,nd
         j=0
         do k3=-npw2,-1
            do k2=-npw2,npw2-1
               do k1=-npw2,npw2-1
                  j=j+1
                  ffexp(j,ind)=wnufft(j)*(fk(k1,k2,k3,1,ind)*ts(k1)+
     1                fk(k1,k2,k3,2,ind)*ts(k2)+
     2                fk(k1,k2,k3,3,ind)*ts(k3))
               enddo
            enddo
         enddo
      enddo

      return
      end
c
C
C
C
      subroutine g3dformpwcd_vec(nd,delta,eps,sources,ns,charge,
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
      real *8 cent(3),sources(3,ns),dipstr(nd,ns)
      real *8 rnormal(3,ns),charge(nd,ns)
      real *8 ws(-npw/2:npw/2-1), ts(-npw/2:npw/2-1)
      real *8 wnufft(nexp)
      complex *16 ffexp(nexp,nd)
      real *8 x,y,chg
      real *8, allocatable ::  xj(:),yj(:),zj(:)
      complex *16, allocatable :: cj(:,:,:),fk(:,:,:,:,:),wj(:)
c     this (since unallocated) used to pass a NULL ptr to FINUFFT...
      integer*8, allocatable :: null
      integer*8 ns8, npw8
      complex *16 eye,ztmp
      eye = dcmplx(0,1)

      npw2=npw/2
      
      allocate(xj(ns))
      allocate(yj(ns))
      allocate(zj(ns))
      allocate(cj(ns,0:3,nd))
      allocate(wj(ns))
      allocate(fk(-npw2:npw2-1,-npw2:npw2-1,-npw2:npw2-1,0:3,nd))
c
      dsq0 = 1/dsqrt(delta)
      dsq = 2*ts(0)/dsqrt(delta)
C
C     accumulate expansion due to each source.
C
      do j=1,ns
         xj(j) = (sources(1,j) - cent(1))*dsq
         yj(j) = (sources(2,j) - cent(2))*dsq
         zj(j) = (sources(3,j) - cent(3))*dsq
         wj(j) = exp(-0.5d0*eye*(xj(j)+yj(j)+zj(j)))
      enddo
      
      do ind = 1,nd
         do j=1,ns
            cj(j,0,ind) = charge(ind,j)*wj(j)
            ztmp = -eye*wj(j)*(dipstr(ind,j)*dsq0)
            cj(j,1,ind) = ztmp*rnormal(1,j)
            cj(j,2,ind) = ztmp*rnormal(2,j)
            cj(j,3,ind) = ztmp*rnormal(3,j)
         enddo
      enddo
      
      ns8=ns
      npw8=npw
      iflag = -1
c     total number of NUFFTs
      ntrans = 4*nd
      call finufft3d1many(ntrans,ns8,xj,yj,zj,cj,iflag,eps,
     1       npw8,npw8,npw8,fk,null,ier)
cccc  print *, 'ier=', ier
      do ind=1,nd
         j=0
         do k3=-npw2,-1
            do k2=-npw2,npw2-1
               do k1=-npw2,npw2-1
                  j=j+1
                  ffexp(j,ind)=wnufft(j)*(fk(k1,k2,k3,0,ind)+
     1                fk(k1,k2,k3,1,ind)*ts(k1)+
     1                fk(k1,k2,k3,2,ind)*ts(k2)+
     2                fk(k1,k2,k3,3,ind)*ts(k3))
               enddo
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
      subroutine g3dpwevalp_vec(nd,delta,eps,center,npw,ws,ts,
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
      real *8 delta,center(3),targ(3,nt)
      real *8 pot(nd,nt)
      
      real *8 ws(-npw/2:npw/2-1),ts(-npw/2:npw/2-1)
      complex *16 pwexp(nexp,nd)
      integer*8 nt8, npw8

      integer i,ind,j1,j2,j,npw2,itarg,k
      real *8 x,y,dsq
      complex *16 eye
      complex *16 qqx,qqy,qq1,qq2
      
      complex *16 z,cd
      real *8, allocatable ::  xj(:),yj(:),zj(:)
      complex *16, allocatable :: cj(:,:),fk(:,:),wj(:)
c     this (since unallocated) used to pass a NULL ptr to FINUFFT...
      integer*8, allocatable :: null
C
      eye = dcmplx(0,1)
      
      allocate(xj(nt))
      allocate(yj(nt))
      allocate(zj(nt))
      allocate(cj(nt,nd))
      allocate(wj(nt))
      allocate(fk(2*nexp,nd))
      
      dsq = 2*ts(0)/dsqrt(delta)
C
      do j=1,nt
         xj(j) = (targ(1,j) - center(1))*dsq
         yj(j) = (targ(2,j) - center(2))*dsq
         zj(j) = (targ(3,j) - center(3))*dsq
         wj(j) = exp(0.5d0*eye*(xj(j)+yj(j)+zj(j)))
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
      call finufft3d2many(nd,nt8,xj,yj,zj,cj,iflag,eps,npw8,npw8,npw8,
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
      subroutine g3dpwevalg_vec(nd,delta,eps,center,npw,ws,ts,
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
      real *8 delta,center(3),targ(3,nt)
      real *8 pot(nd,nt),grad(nd,3,nt)
      
      real *8 ws(-npw/2:npw/2-1),ts(-npw/2:npw/2-1)
      complex *16 pwexp(nexp,nd)
      integer*8 nt8, npw8

      integer i,ind,j1,j2,j,npw2,itarg,k
      real *8 x,y,dsq
      complex *16 eye
      complex *16 qqx,qqy,qq1,qq2
      
      complex *16 z,cd
      real *8, allocatable ::  xj(:),yj(:),zj(:)
      complex *16, allocatable :: cj(:,:,:),fk(:,:,:),wj(:)
c     this (since unallocated) used to pass a NULL ptr to FINUFFT...
      integer*8, allocatable :: null
C
      eye = dcmplx(0,1)
      
      allocate(xj(nt))
      allocate(yj(nt))
      allocate(zj(nt))
      allocate(cj(nt,0:3,nd))
      allocate(wj(nt))
      allocate(fk(2*nexp,0:3,nd))

      dsq = 2*ts(0)/dsqrt(delta)
C
      do j=1,nt
         xj(j) = (targ(1,j) - center(1))*dsq
         yj(j) = (targ(2,j) - center(2))*dsq
         zj(j) = (targ(3,j) - center(3))*dsq
         wj(j) = exp(0.5d0*eye*(xj(j)+yj(j)+zj(j)))
      enddo

      z = eye/dsqrt(delta)
      
      do ind = 1,nd
         do j=1,nexp
            fk(j,0,ind)=pwexp(j,ind)
         enddo

         j=0
         do j3=-npw/2,-1
            do j2=-npw/2,npw/2-1
               do j1=-npw/2,npw/2-1
                  j=j+1
                  fk(j,1,ind)=pwexp(j,ind)*z*ts(j1)
                  fk(j,2,ind)=pwexp(j,ind)*z*ts(j2)
                  fk(j,3,ind)=pwexp(j,ind)*z*ts(j3)
               enddo
            enddo
         enddo
         
         do k=0,3
            do j=nexp+1,2*nexp
               fk(j,k,ind)=0
            enddo
         enddo
      enddo

      nt8=nt
      npw8=npw
      ntrans=4*nd
      iflag = 1
      call finufft3d2many(ntrans,nt8,xj,yj,zj,cj,iflag,eps,
     1    npw8,npw8,npw8,fk,null,ier)
      
      do ind=1,nd
         do j=1,nt
            pot(ind,j)=pot(ind,j)+dreal(cj(j,0,ind)*wj(j))*2
            grad(ind,1,j)=grad(ind,1,j)+dreal(cj(j,1,ind)*wj(j))*2
            grad(ind,2,j)=grad(ind,2,j)+dreal(cj(j,2,ind)*wj(j))*2
            grad(ind,3,j)=grad(ind,3,j)+dreal(cj(j,3,ind)*wj(j))*2
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
      subroutine g3dpwevalh_vec(nd,delta,eps,center,npw,ws,ts,
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
      real *8 delta,center(3),targ(3,nt)
      real *8 pot(nd,nt),grad(nd,3,nt),hess(nd,6,nt)
      
      real *8 ws(-npw/2:npw/2-1),ts(-npw/2:npw/2-1)
      complex *16 pwexp(nexp,nd)
      integer*8 nt8, npw8

      integer i,ind,j1,j2,j,npw2,itarg,k
      real *8 x,y,dsq
      complex *16 eye
      complex *16 qqx,qqy,qq1,qq2
      
      complex *16 z,cd
      real *8, allocatable ::  xj(:),yj(:),zj(:)
      complex *16, allocatable :: cj(:,:,:),fk(:,:,:),wj(:)
      complex *16, allocatable :: cp(:)
      
c     this (since unallocated) used to pass a NULL ptr to FINUFFT...
      integer*8, allocatable :: null
C
      eye = dcmplx(0,1)
      
      allocate(xj(nt))
      allocate(yj(nt))
      allocate(zj(nt))
      allocate(cj(nt,0:9,nd))
      allocate(wj(nt))
      allocate(fk(2*nexp,0:9,nd))
      allocate(cp(-npw/2:npw/2-1))
      
      dsq = 2*ts(0)/dsqrt(delta)
C
      do j=1,nt
         xj(j) = (targ(1,j) - center(1))*dsq
         yj(j) = (targ(2,j) - center(2))*dsq
         zj(j) = (targ(3,j) - center(3))*dsq
         wj(j) = exp(0.5d0*eye*(xj(j)+yj(j)+zj(j)))
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
         do j3=-npw/2,-1
            do j2=-npw/2,npw/2-1
               do j1=-npw/2,npw/2-1
                  j=j+1
                  fk(j,1,ind)=pwexp(j,ind)*cp(j1)
                  fk(j,2,ind)=pwexp(j,ind)*cp(j2)
                  fk(j,3,ind)=pwexp(j,ind)*cp(j3)
                  
                  fk(j,4,ind)=pwexp(j,ind)*cp(j1)*cp(j1)
                  fk(j,5,ind)=pwexp(j,ind)*cp(j2)*cp(j2)
                  fk(j,6,ind)=pwexp(j,ind)*cp(j3)*cp(j3)

                  fk(j,7,ind)=pwexp(j,ind)*cp(j1)*cp(j2)
                  fk(j,8,ind)=pwexp(j,ind)*cp(j1)*cp(j3)
                  fk(j,9,ind)=pwexp(j,ind)*cp(j2)*cp(j3)
               enddo
            enddo
         enddo
         
         do k=0,9
            do j=nexp+1,2*nexp
               fk(j,k,ind)=0
            enddo
         enddo
      enddo

      nt8=nt
      npw8=npw
      ntrans=10*nd
      iflag = 1
      
      call finufft3d2many(ntrans,nt8,xj,yj,zj,cj,iflag,eps,
     1    npw8,npw8,npw8,fk,null,ier)
      
      do ind=1,nd
         do j=1,nt
            pot(ind,j)=pot(ind,j)+dreal(cj(j,0,ind)*wj(j))*2
            
            grad(ind,1,j)=grad(ind,1,j)+dreal(cj(j,1,ind)*wj(j))*2
            grad(ind,2,j)=grad(ind,2,j)+dreal(cj(j,2,ind)*wj(j))*2
            grad(ind,3,j)=grad(ind,3,j)+dreal(cj(j,3,ind)*wj(j))*2
            
            hess(ind,1,j)=hess(ind,1,j)+dreal(cj(j,4,ind)*wj(j))*2
            hess(ind,2,j)=hess(ind,2,j)+dreal(cj(j,5,ind)*wj(j))*2
            hess(ind,3,j)=hess(ind,3,j)+dreal(cj(j,6,ind)*wj(j))*2
            hess(ind,4,j)=hess(ind,4,j)+dreal(cj(j,7,ind)*wj(j))*2
            hess(ind,5,j)=hess(ind,5,j)+dreal(cj(j,8,ind)*wj(j))*2
            hess(ind,6,j)=hess(ind,6,j)+dreal(cj(j,9,ind)*wj(j))*2
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
      do j1=-npw/2,-1
         do j2=-npw/2,npw/2-1
            do j3=-npw/2,npw/2-1
               j=j+1
               wnufft(j) = ws(j3)*ws(j2)*ws(j1)
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
      
      complex *16 wshift(nexp,-nmax:nmax,-nmax:nmax,-nmax:nmax)
      
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
      do k3=-nmax,nmax
         j=0
         do j1=-npw/2,-1
         do j2=-npw/2,npw/2-1
         do j3=-npw/2,npw/2-1
            j=j+1
            wshift(j,k3,k2,k1) = ww(j3,k3)*ww(j2,k2)*ww(j1,k1)
         enddo
         enddo
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
      subroutine g3dshiftpw_vec(nd,nexp,pwexp1,
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
      subroutine g3dcopypwexp_vec(nd,nexp,pwexp1,
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
