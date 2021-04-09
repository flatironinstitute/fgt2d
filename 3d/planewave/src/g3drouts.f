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
c                      
c      
c     g3dformpwd_vec : computes the PW expansion due to dipoles
c                      
c      
c     g3dformpwcd_vec : computes the PW expansion due to charges and dipoles
c                      
c      
c     g3dpwevalp_vec : evaluates the PW expansion
c                      potential only
c      
c     g3dpwevalg_vec : evaluates the PW expansion
c                      potential + gradient
c      
c     g3dpwevalh_vec : evaluates the PW expansion
c                      potential + gradient + hessian
c      
c     g3dherm2local_vec : converts the Hermite expansion to the Taylor expansion
c      
C     g3dh2pw_vec : converts the Hermite expansion to the PW expansion
c      
c     g3dh2pw_real_vec : same as above, use REAL 1D h2pw matrix
c
c     g3dpw2local_vec : converts PW expansions to the local Taylor expansion
c
c     g3dpw2local_real_vec : converts PW expansion to local expansion
c                            USE REAL 1D pw2l matrix
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
c     g3dshiftpw_vec : translates PW expansion via precomputed translation 
c                      matrices
c
c     g3dcopypwexp_vec : copy PW expansion
c
c     g3dshiftpw0_vec : translates PW expansion from one center to another
c                       without any precomputation
c
c
c     g3dh2lmat : returns the matrix converting the Hermite expansion
C                 to the local Taylor expansion (1D only)
c      
c     g3dh2xmat : returns the matrix converting the Hermite expansion 
C                 to the planewave expansion (1D only)
C     
c     g3drh2xmat : returns REAL matrix converting the Hermite expansion 
C                 to the planewave expansion (1D only)
C                 That is, i^j factor is ignored
C     
C     g3dx2lmat : returns the matrix converting the planewave expansion
C                 to the local Taylor expansion (1D only)
C
C     g3drx2lmat : returns REAL matrix converting the planewave expansion
C                 to the local Taylor expansion (1D only),
C                 that is, i^j factor is ignored.
c
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
      real *8, allocatable :: xp(:),yp(:),zp(:),sqk(:)
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
      allocate(sqk(nterms))
c
      dsq = 1.0D0/dsqrt(delta)
C
C     accumulate expansion due to each source.
C
      do k=1,nterms
         sqk(k)=dsqrt(dble(k))
      enddo
      
      do i=1,ns
         x = (sources(1,i) - cent(1))*dsq
         y = (sources(2,i) - cent(2))*dsq
         z = (sources(3,i) - cent(3))*dsq
         xp(0) = 1.0D0
         yp(0) = 1.0D0
         zp(0) = 1.0D0
         do k = 1,nterms
            xp(k) = xp(k-1)*x/sqk(k)
            yp(k) = yp(k-1)*y/sqk(k)
            zp(k) = zp(k-1)*z/sqk(k)
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
      real *8, allocatable :: sqk(:)
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
      allocate(sqk(nterms))
c
      dsq = 1.0D0/dsqrt(delta)
C
C     accumulate expansion due to each source.
C
      do k=1,nterms
         sqk(k)=dsqrt(dble(k))
      enddo
c
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
            hy(j1)=hy(j1-1)*y/sqk(j1)
            hx(j1)=hx(j1-1)*x/sqk(j1)
            hz(j1)=hz(j1-1)*z/sqk(j1)
            dxhx(j1) = j1*dsq*hx(j1-1)/sqk(j1)
            dyhy(j1) = j1*dsq*hy(j1-1)/sqk(j1)
            dzhz(j1) = j1*dsq*hz(j1-1)/sqk(j1)
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
                  c2 = dyhy(k)*d2+hy(k)*d3
                  do j=0,nterms
                     ffexp(j,k,l,ind) = ffexp(j,k,l,ind)
     1                   +dxhx(j)*c1+hx(j)*c2
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
                  c1 = hy(k)*d1
                  c2 = dyhy(k)*d2+hy(k)*(d3+ztmp)
                  do j=0,nterms
                     ffexp(j,k,l,ind) = ffexp(j,k,l,ind)+
     1                   dxhx(j)*c1+hx(j)*c2
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
C     ntarg    = number of targets
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
C     ntarg    = number of targets
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
C     ntarg    = number of targets
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
C     sources(3,ns) = coordinates of sources
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
C     rnormal       = dipole directions
C     dipstr        = dipole strengths 
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
C     rnormal       = dipole directions
C     dipstr        = dipole strengths 
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
C     targ          = target location
C     ntarg         = number of targets
C
C     OUTPUT:
C
C     pot           = potential (or vectorized potentials) incremented
C
      implicit real*8 (a-h,o-z)
      real *8 local(0:nlocal,0:nlocal,0:nlocal,nd)
      real *8 center(3),targ(3,ntarg)
      real *8 pot(nd,ntarg),xp(0:100),yp(0:100),zp(0:100)
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
C     targ          = target location
C     ntarg         = number of targets
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
C     targ          = target location
C     ntarg         = number of targets
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
                  d1xx=d1xx+d2xx*yp(j2)
                  
                  d1y=d1y+d2*ypy(j2)
                  d1xy=d1xy+d2x*ypy(j2)
                  
                  d1yy=d1yy+d2*ypyy(j2)
               enddo
               dd=dd+d1*zp(j1)
               g1=g1+d1x*zp(j1)
               g2=g2+d1y*zp(j1)

               h11=h11+d1xx*zp(j1)
               h22=h22+d1yy*zp(j1)
               h12=h12+d1xy*zp(j1)
               
               g3=g3+d1*zpz(j1)
               h13=h13+d1x*zpz(j1)
               h23=h23+d1y*zpz(j1)
               
               h33=h33+d1*zpzz(j1)
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
C 3D translations (h2l, h2pw, h2pw_real)
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
      subroutine g3dh2pw_vec(nd,ntermsh,npw,h2x,
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
cccc      complex *16 h2x(npw,0:ntermsh)
      complex *16 h2x(0:ntermsh,npw)
      real *8 hexp(0:ntermsh,0:ntermsh,0:ntermsh,nd)
      complex *16 pwexp(npw,npw,npw/2,nd)

      complex *16, allocatable :: pxhyz(:,:,:),pxyhz(:,:,:)
      complex *16 eye,cd
      
      eye = dcmplx(0,1)
c
      allocate(pxhyz(npw,0:ntermsh,0:ntermsh))
      allocate(pxyhz(npw,npw,0:ntermsh))

c      
      do ind=1,nd
c        transform in x
         do j1=0,ntermsh
            do j2=0,ntermsh
               do k3=1,npw/2
                  cd=0
                  do j3=0,ntermsh
cccc                     cd=cd+h2x(k3,j3)*hexp(j3,j2,j1,ind)
                     cd=cd+h2x(j3,k3)*hexp(j3,j2,j1,ind)
                  enddo
                  pxhyz(k3,j2,j1)=cd
                  pxhyz(k3+npw/2,j2,j1)=conjg(cd)
               enddo
            enddo
         enddo

c        transform in y
         do j1=0,ntermsh
            do k2=1,npw
               do k3=1,npw
                  cd=0
                  do j2=0,ntermsh
cccc                     cd=cd+h2x(k2,j2)*pxhyz(k3,j2,j1)
                     cd=cd+h2x(j2,k2)*pxhyz(k3,j2,j1)
                  enddo
                  pxyhz(k3,k2,j1)=cd
               enddo
            enddo
         enddo

c        transform in z
         do k1=1,npw/2
            do k2=1,npw
               do k3=1,npw
                  cd=0
                  do j1=0,ntermsh
cccc                     cd=cd+h2x(k1,j1)*pxyhz(k3,k2,j1)
                     cd=cd+h2x(j1,k1)*pxyhz(k3,k2,j1)
                  enddo
                  pwexp(k3,k2,k1,ind) = pwexp(k3,k2,k1,ind) + cd
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
      subroutine g3dh2pw_real_vec(nd,ntermsh,npw,h2x,
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
      real *8 h2x(0:ntermsh,npw)
      real *8 hexp(0:ntermsh,0:ntermsh,0:ntermsh,nd)
      complex *16 pwexp(npw,npw,npw/2,nd)

      complex *16, allocatable :: pxhyz(:,:,:),pxyhz(:,:,:)
      complex *16 eye,cd,c0,c1,c2,c3,ce,co
      
      eye = dcmplx(0,1)
c
      allocate(pxhyz(npw,0:ntermsh,0:ntermsh))
      allocate(pxyhz(npw,npw,0:ntermsh))
c
      npw2=npw/2
      
      do ind=1,nd
c        transform in x
         do j1=0,ntermsh
            do j2=0,ntermsh
               do k3=1,npw/2
                  c0=0
                  c1=0
                  c2=0
                  c3=0
                  do j3=0,ntermsh,4
                     c0=c0+h2x(j3,k3)*hexp(j3,j2,j1,ind)
                  enddo
                  do j3=1,ntermsh,4
                     c1=c1+h2x(j3,k3)*hexp(j3,j2,j1,ind)
                  enddo
                  do j3=2,ntermsh,4
                     c2=c2+h2x(j3,k3)*hexp(j3,j2,j1,ind)
                  enddo
                  do j3=3,ntermsh,4
                     c3=c3+h2x(j3,k3)*hexp(j3,j2,j1,ind)
                  enddo
                  ce=c0-c2
                  co=eye*(c1-c3)
                  pxhyz(k3,j2,j1)=ce+co
                  pxhyz(k3+npw2,j2,j1)=ce-co
               enddo
            enddo
         enddo

c        transform in y
         do j1=0,ntermsh
            do k2=1,npw2
               do k3=1,npw
                  c0=0
                  c1=0
                  c2=0
                  c3=0                  
                  do j2=0,ntermsh,4
                     c0=c0+h2x(j2,k2)*pxhyz(k3,j2,j1)
                  enddo
                  do j2=1,ntermsh,4
                     c1=c1+h2x(j2,k2)*pxhyz(k3,j2,j1)
                  enddo
                  do j2=2,ntermsh,4
                     c2=c2+h2x(j2,k2)*pxhyz(k3,j2,j1)
                  enddo
                  do j2=3,ntermsh,4
                     c3=c3+h2x(j2,k2)*pxhyz(k3,j2,j1)
                  enddo
                  ce=c0-c2
                  co=eye*(c1-c3)
                  pxyhz(k3,k2,j1)=ce+co
                  pxyhz(k3,k2+npw2,j1)=ce-co
               enddo
            enddo
         enddo
         
c        transform in z
         do k1=1,npw/2
            do k2=1,npw
               do k3=1,npw
                  c0=0
                  c1=0
                  c2=0
                  c3=0
                  do j1=0,ntermsh,4
                     c0=c0+h2x(j1,k1)*pxyhz(k3,k2,j1)
                  enddo
                  do j1=1,ntermsh,4
                     c1=c1+h2x(j1,k1)*pxyhz(k3,k2,j1)
                  enddo
                  do j1=2,ntermsh,4
                     c2=c2+h2x(j1,k1)*pxyhz(k3,k2,j1)
                  enddo
                  do j1=3,ntermsh,4
                     c3=c3+h2x(j1,k1)*pxyhz(k3,k2,j1)
                  enddo
                  pwexp(k3,k2,k1,ind) = pwexp(k3,k2,k1,ind) +
     1                c0-c2+eye*(c1-c3)
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
C
C 3D translations pw2local, pw2local_real
C
c***********************************************************************      
      subroutine g3dpw2local_vec(nd,nlocal,npw,pw2l,
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
      real *8 local(0:nlocal,0:nlocal,0:nlocal,nd)
      complex *16 pwexp(npw,npw,npw/2,nd)
      complex *16, allocatable :: lxpyz(:,:,:),lxypz(:,:,:)
      complex *16 cd
c      
c
      allocate(lxpyz(0:nlocal,npw,npw/2))
      allocate(lxypz(0:nlocal,0:nlocal,npw/2))

      npw2=npw/2

      do ind=1,nd
c        transform in x
         do j1=1,npw/2
            do j2=1,npw
               do k3=0,nlocal
                  cd=0
                  do j3=1,npw
                     cd=cd+pw2l(j3,k3)*pwexp(j3,j2,j1,ind)
                  enddo
                  lxpyz(k3,j2,j1)=cd
               enddo
            enddo
         enddo

c        transform in y
         do j1=1,npw/2
            do k2=0,nlocal
               do k3=0,nlocal
                  cd=0
                  do j2=1,npw
                     cd=cd+pw2l(j2,k2)*lxpyz(k3,j2,j1)
                  enddo
                  lxypz(k3,k2,j1)=cd
               enddo
            enddo
         enddo

c        transform in z
         do k1=0,nlocal
            do k2=0,nlocal
               do k3=0,nlocal
                  cd=0
                  do j1=1,npw/2
                     cd=cd+pw2l(j1,k1)*lxypz(k3,k2,j1)
                  enddo
                  local(k3,k2,k1,ind)=dble(cd)*2
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
      subroutine g3dpw2local_real_vec(nd,nlocal,npw,pw2l,
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
      real *8 local(0:nlocal,0:nlocal,0:nlocal,nd)
      complex *16 pwexp(npw,npw,npw/2,nd)
      complex *16, allocatable :: lxpyz(:,:,:),lxypz(:,:,:)
      complex *16 cd,eye
c      
c
      allocate(lxpyz(0:nlocal,npw,npw/2))
      allocate(lxypz(0:nlocal,0:nlocal,npw/2))
c      
c
      eye = dcmplx(0,1)
      
      npw2=npw/2

      do ind=1,nd
c        transform in x
         do j1=1,npw/2
            do j2=1,npw
               do k3=0,nlocal,4
                  cd=0
                  do j3=1,npw2
                     cd=cd+pw2l(j3,k3)*(pwexp(j3,j2,j1,ind)+
     1                   pwexp(j3+npw2,j2,j1,ind))
                  enddo
                  lxpyz(k3,j2,j1)=cd
               enddo
               do k3=1,nlocal,4
                  cd=0
                  do j3=1,npw2
                     cd=cd+pw2l(j3,k3)*(pwexp(j3,j2,j1,ind)-
     1                   pwexp(j3+npw2,j2,j1,ind))
                  enddo
                  lxpyz(k3,j2,j1)=cd*eye
               enddo
               do k3=2,nlocal,4
                  cd=0
                  do j3=1,npw2
                     cd=cd+pw2l(j3,k3)*(pwexp(j3,j2,j1,ind)+
     1                   pwexp(j3+npw2,j2,j1,ind))
                  enddo
                  lxpyz(k3,j2,j1)=-cd
               enddo
               do k3=3,nlocal,4
                  cd=0
                  do j3=1,npw2
                     cd=cd+pw2l(j3,k3)*(pwexp(j3,j2,j1,ind)-
     1                   pwexp(j3+npw2,j2,j1,ind))
                  enddo
                  lxpyz(k3,j2,j1)=-eye*cd
               enddo
            enddo
         enddo
c        transform in y
         do j1=1,npw/2
            do k2=0,nlocal,4
               do k3=0,nlocal
                  cd=0
                  do j2=1,npw2
                     cd=cd+pw2l(j2,k2)*(lxpyz(k3,j2,j1)+
     1                   lxpyz(k3,j2+npw2,j1))
                  enddo
                  lxypz(k3,k2,j1)=cd
               enddo
            enddo
            do k2=1,nlocal,4
               do k3=0,nlocal
                  cd=0
                  do j2=1,npw2
                     cd=cd+pw2l(j2,k2)*(lxpyz(k3,j2,j1)-
     1                   lxpyz(k3,j2+npw2,j1))
                  enddo
                  lxypz(k3,k2,j1)=cd*eye
               enddo
            enddo
            do k2=2,nlocal,4
               do k3=0,nlocal
                  cd=0
                  do j2=1,npw/2
                     cd=cd+pw2l(j2,k2)*(lxpyz(k3,j2,j1)+
     1                   lxpyz(k3,j2+npw2,j1))
                  enddo
                  lxypz(k3,k2,j1)=-cd
               enddo
            enddo
            do k2=3,nlocal,4
               do k3=0,nlocal
                  cd=0
                  do j2=1,npw/2
                     cd=cd+pw2l(j2,k2)*(lxpyz(k3,j2,j1)-
     1                   lxpyz(k3,j2+npw2,j1))
                  enddo
                  lxypz(k3,k2,j1)=-eye*cd
               enddo
            enddo
         enddo
         
c        transform in z
         do k1=0,nlocal,4
            do k2=0,nlocal
               do k3=0,nlocal
                  cd=0
                  do j1=1,npw/2
                     cd=cd+pw2l(j1,k1)*lxypz(k3,k2,j1)
                  enddo
                  local(k3,k2,k1,ind)=dble(cd)*2
               enddo
            enddo
         enddo
         do k1=1,nlocal,4
            do k2=0,nlocal
               do k3=0,nlocal
                  cd=0
                  do j1=1,npw/2
                     cd=cd+pw2l(j1,k1)*lxypz(k3,k2,j1)
                  enddo
                  local(k3,k2,k1,ind)=-aimag(cd)*2
               enddo
            enddo
         enddo
         do k1=2,nlocal,4
            do k2=0,nlocal
               do k3=0,nlocal
                  cd=0
                  do j1=1,npw/2
                     cd=cd+pw2l(j1,k1)*lxypz(k3,k2,j1)
                  enddo
                  local(k3,k2,k1,ind)=-dble(cd)*2
               enddo
            enddo
         enddo
         do k1=3,nlocal,4
            do k2=0,nlocal
               do k3=0,nlocal
                  cd=0
                  do j1=1,npw/2
                     cd=cd+pw2l(j1,k1)*lxypz(k3,k2,j1)
                  enddo
                  local(k3,k2,k1,ind)=aimag(cd)*2
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
C
C************************************************************************
C
C     form PW exp subroutines (charge, dipole, charge+dipole)
C
C***********************************************************************
      subroutine g3dformpwc_vec(nd,delta,sources,ns,charge,center,
     1    npw,ws,ts,pwexp)
C
C     This subroutine computes the PW expansions about
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
c     npw           = number of terms in 1d PW exp
c     ws,ts         = weights and nodes in 1d PW exp
C
C     OUTPUT:
C
C     pwexp         = 3d PW exp in the order of x, y, z (only one half in z)
C                     incremented
c      
      implicit real*8 (a-h,o-z)
      real *8 center(3),sources(3,ns),charge(nd,ns)
      real *8 ws(npw),ts(npw)
      complex *16 pwexp(npw,npw,npw/2,nd)

      integer i,ind,j1,j2,j3,j,npw2,itarg,k
      real *8 x,y,z,dsq
      complex *16 eye,ytmp,ztmp
      complex *16 qqx,qqy,qqz,qq1,qq2,qq3
      
      complex *16 ww1(100)
      complex *16 ww2(100)
      complex *16 ww3(100)
C
      eye = dcmplx(0,1)
      dsq = 1.0D0/dsqrt(delta)
C
      npw2=npw/2
C
      do i=1,ns
         x = (sources(1,i) - center(1))*dsq
         y = (sources(2,i) - center(2))*dsq
         z = (sources(3,i) - center(3))*dsq
c
         qqx = cdexp(-eye*ts(1)*x)
         qqy = cdexp(-eye*ts(1)*y)
         qqz = cdexp(-eye*ts(1)*z)
         
         qq1 = qqx
         qq2 = qqy
         qq3 = qqz
         
         qqx = qqx*qqx
         qqy = qqy*qqy
         qqz = qqz*qqz
         

         do j1=1,npw2
            ww1(j1) = qq1*ws(j1)
            ww2(j1) = qq2*ws(j1)            
            ww3(j1) = qq3*ws(j1)            
            qq1 = qq1*qqx
            qq2 = qq2*qqy
            qq3 = qq3*qqz
            
            ww1(j1+npw2) = dconjg(ww1(j1))
            ww2(j1+npw2) = dconjg(ww2(j1))
            ww3(j1+npw2) = dconjg(ww3(j1))
         enddo
c
         do ind = 1,nd
            chg=charge(ind,i)
            do j3=1,npw/2
               ztmp=chg*ww3(j3)
               do j2=1,npw
                  ytmp=ztmp*ww2(j2)
                  do j1=1,npw
                     pwexp(j1,j2,j3,ind)=pwexp(j1,j2,j3,ind)+
     1                   ytmp*ww1(j1)
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
      subroutine g3dformpwd_vec(nd,delta,sources,ns,rnormal,dipstr,
     1    center,npw,ws,ts,pwexp)
C
C     This subroutine computes the PW expansions about
C     the center CENT due to sources of dipoles.
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
c     npw           = number of terms in 1d PW exp
c     ws,ts         = weights and nodes in 1d PW exp
C
C     OUTPUT:
C
C     pwexp         = 3d PW exp in the order of x, y, z (only one half in z)
C                     incremented
c      
      implicit real*8 (a-h,o-z)
      real *8 center(3),sources(3,ns)
      real *8 rnormal(3,ns),dipstr(nd,ns)
      real *8 ws(npw),ts(npw)
      complex *16 pwexp(npw,npw,npw/2,nd)

      integer i,ind,j1,j2,j3,j,npw2,itarg,k
      real *8 x,y,z,dsq
      complex *16 eye,r1,r2,r3,rz1,rz2,rz3,ry1,ry2,ry3
      complex *16 qqx,qqy,qqz,qq1,qq2,qq3,ztmp
      
      complex *16 ww1(100),ww1x(100)
      complex *16 ww2(100),ww2y(100)
      complex *16 ww3(100),ww3z(100)
C
      eye = dcmplx(0,1)
      dsq = 1.0D0/dsqrt(delta)
C
      npw2=npw/2
C
      do i=1,ns
         x = (sources(1,i) - center(1))*dsq
         y = (sources(2,i) - center(2))*dsq
         z = (sources(3,i) - center(3))*dsq
c
         qqx = cdexp(-eye*ts(1)*x)
         qqy = cdexp(-eye*ts(1)*y)
         qqz = cdexp(-eye*ts(1)*z)
         
         qq1 = qqx
         qq2 = qqy
         qq3 = qqz
         
         qqx = qqx*qqx
         qqy = qqy*qqy
         qqz = qqz*qqz
         

         do j1=1,npw2
            ww1(j1) = qq1*ws(j1)
            ww2(j1) = qq2*ws(j1)            
            ww3(j1) = qq3*ws(j1)            
            qq1 = qq1*qqx
            qq2 = qq2*qqy
            qq3 = qq3*qqz
            
            ww1(j1+npw2) = dconjg(ww1(j1))
            ww2(j1+npw2) = dconjg(ww2(j1))
            ww3(j1+npw2) = dconjg(ww3(j1))
         enddo

         do j1=1,npw
            ztmp = -eye*ts(j1)*dsq
            ww1x(j1) = ww1(j1)*ztmp
            ww2y(j1) = ww2(j1)*ztmp
            ww3z(j1) = ww3(j1)*ztmp
         enddo
c
         do ind = 1,nd
            r1 = rnormal(1,i)*dipstr(ind,i)
            r2 = rnormal(2,i)*dipstr(ind,i)
            r3 = rnormal(3,i)*dipstr(ind,i)            
            do j3=1,npw/2
               rz1 = ww3(j3)*r1
               rz2 = ww3(j3)*r2
               rz3 = ww3z(j3)*r3               
               do j2=1,npw
                  ry1 = rz1*ww2(j2)
                  ry2 = rz2*ww2y(j2)+rz3*ww2(j2)
                  do j1=1,npw
                     pwexp(j1,j2,j3,ind)=pwexp(j1,j2,j3,ind)+
     1                   ww1x(j1)*ry1+ww1(j1)*ry2
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
      subroutine g3dformpwcd_vec(nd,delta,sources,ns,charge,
     1    rnormal,dipstr,center,npw,ws,ts,pwexp)
C
C     This subroutine computes the PW expansions about
C     the center CENT due to charges and dipoles.
C
C     INPUT
C
c     nd            = vector length (parallel)
C     delta         = Gaussian variance
C     sources(2,ns) = coordinates of sources
C     ns            = number of sources
C     charge        = strength of sources
C     rnormal       = dipole directions
C     dipstr        = dipole strengths 
C     center        = center of the expansion
c     npw           = number of terms in 1d PW exp
c     ws,ts         = weights and nodes in 1d PW exp
C
C     OUTPUT:
C
C     pwexp         = 3d PW exp in the order of x, y, z (only one half in z)
C                     incremented
c      
      implicit real*8 (a-h,o-z)
      real *8 center(3),sources(3,ns),charge(nd,ns)
      real *8 rnormal(3,ns),dipstr(nd,ns)
      real *8 ws(npw),ts(npw)
      complex *16 pwexp(npw,npw,npw/2,nd)

      integer i,ind,j1,j2,j3,j,npw2,itarg,k
      real *8 x,y,z,dsq
      complex *16 eye,r1,r2,r3,rz1,rz2,rz3,ry1,ry2,ry3
      complex *16 qqx,qqy,qqz,qq1,qq2,qq3,ztmp
      
      complex *16 ww1(100),ww1x(100)
      complex *16 ww2(100),ww2y(100)
      complex *16 ww3(100),ww3z(100)
C
      eye = dcmplx(0,1)
      dsq = 1.0D0/dsqrt(delta)
C
      npw2=npw/2
C
      do i=1,ns
         x = (sources(1,i) - center(1))*dsq
         y = (sources(2,i) - center(2))*dsq
         z = (sources(3,i) - center(3))*dsq
c
         qqx = cdexp(-eye*ts(1)*x)
         qqy = cdexp(-eye*ts(1)*y)
         qqz = cdexp(-eye*ts(1)*z)
         
         qq1 = qqx
         qq2 = qqy
         qq3 = qqz
         
         qqx = qqx*qqx
         qqy = qqy*qqy
         qqz = qqz*qqz
         

         do j1=1,npw2
            ww1(j1) = qq1*ws(j1)
            ww2(j1) = qq2*ws(j1)            
            ww3(j1) = qq3*ws(j1)            
            qq1 = qq1*qqx
            qq2 = qq2*qqy
            qq3 = qq3*qqz
            
            ww1(j1+npw2) = dconjg(ww1(j1))
            ww2(j1+npw2) = dconjg(ww2(j1))
            ww3(j1+npw2) = dconjg(ww3(j1))
         enddo

         do j1=1,npw
            ztmp = -eye*ts(j1)*dsq
            ww1x(j1) = ww1(j1)*ztmp
            ww2y(j1) = ww2(j1)*ztmp
            ww3z(j1) = ww3(j1)*ztmp
         enddo
c
         do ind = 1,nd
            chg = charge(ind,i)
            r1 = rnormal(1,i)*dipstr(ind,i)
            r2 = rnormal(2,i)*dipstr(ind,i)
            r3 = rnormal(3,i)*dipstr(ind,i)            
            do j3=1,npw/2
               ztmp = ww3(j3)*chg
               rz1 = ww3(j3)*r1
               rz2 = ww3(j3)*r2
               rz3 = ww3z(j3)*r3               
               do j2=1,npw
                  ry1 = rz1*ww2(j2)
                  ry2 = rz2*ww2y(j2)+(rz3+ztmp)*ww2(j2)
                  do j1=1,npw
                     pwexp(j1,j2,j3,ind)=pwexp(j1,j2,j3,ind)+
     1                   ww1x(j1)*ry1+ww1(j1)*ry2
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
      subroutine g3dpwevalp_vec(nd,delta,center,npw,wx,tx,
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
      real *8 delta,center(3),targ(3,ntarg)
      real *8 pot(nd,ntarg)
      
      real *8 wx(npw),tx(npw)
      complex *16 pwexp(npw,npw,npw/2,nd)

      integer i,ind,j1,j2,j3,j,npw2,itarg,k
      real *8 x,y,z,dsq
      complex *16 eye
      complex *16 qqx,qqy,qqz,qq1,qq2,qq3
      
      complex *16 ww1(100)
      complex *16 ww2(100)
      complex *16 ww3(100)

      complex *16 c1,c2,c3
C
      eye = dcmplx(0,1)
      dsq = 1.0D0/dsqrt(delta)
C
      npw2=npw/2
      do itarg=1,ntarg
         x = (targ(1,itarg) - center(1))*dsq
         y = (targ(2,itarg) - center(2))*dsq
         z = (targ(3,itarg) - center(3))*dsq

         qqx = cdexp(eye*tx(1)*x)
         qqy = cdexp(eye*tx(1)*y)
         qqz = cdexp(eye*tx(1)*z)
         
         qq1 = qqx
         qq2 = qqy
         qq3 = qqz
         
         qqx = qqx*qqx
         qqy = qqy*qqy
         qqz = qqz*qqz
         

         do j1=1,npw2
            ww1(j1) = qq1
            ww2(j1) = qq2            
            ww3(j1) = qq3            
            qq1 = qq1*qqx
            qq2 = qq2*qqy
            qq3 = qq3*qqz
            
            ww1(j1+npw2) = dconjg(ww1(j1))
            ww2(j1+npw2) = dconjg(ww2(j1))
            ww3(j1+npw2) = dconjg(ww3(j1))
         enddo
c
         do ind = 1,nd
            c3=0
            do j3=1,npw/2
               c2=0
               do j2=1,npw
                  c1=0
                  do j1=1,npw
                     c1=c1+pwexp(j1,j2,j3,ind)*ww1(j1)
                  enddo
                  c2=c2+c1*ww2(j2)
               enddo
               c3=c3+c2*ww3(j3)
            enddo
            
            pot(ind,itarg) = pot(ind,itarg)+dreal(c3)*2
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
      subroutine g3dpwevalg_vec(nd,delta,center,npw,wx,tx,
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
      implicit real*8 (a-h,o-z)
      integer nd,npw,ntarg
      real *8 delta,center(3),targ(3,ntarg)
      real *8 pot(nd,ntarg),grad(nd,3,ntarg)
      
      real *8 wx(npw),tx(npw)
      complex *16 pwexp(npw,npw,npw/2,nd)

      integer i,ind,j1,j2,j,npw2,itarg,k
      real *8 x,y,z,dsq
      complex *16 eye
      complex *16 qqx,qqy,qqz,qq1,qq2,qq3
      
      complex *16 ww1(100),ww1x(100)
      complex *16 ww2(100),ww2y(100)
      complex *16 ww3(100),ww3z(100)

      complex *16 cd,g1,g2,g3,d1,d1x,d1y,d2,d2x
C
      eye = dcmplx(0,1)
      dsq = 1.0D0/dsqrt(delta)
C
      npw2=npw/2
      do itarg=1,ntarg
         x = (targ(1,itarg) - center(1))*dsq
         y = (targ(2,itarg) - center(2))*dsq
         z = (targ(3,itarg) - center(3))*dsq

         qqx = cdexp(eye*tx(1)*x)
         qqy = cdexp(eye*tx(1)*y)
         qqz = cdexp(eye*tx(1)*z)
         qq1 = qqx
         qq2 = qqy
         qq3 = qqz
         qqx = qqx*qqx
         qqy = qqy*qqy
         qqz = qqz*qqz
         

         do j1=1,npw2
            ww1(j1) = qq1
            ww2(j1) = qq2            
            ww3(j1) = qq3            
            ww1(j1+npw2) = dconjg(ww1(j1))
            ww2(j1+npw2) = dconjg(ww2(j1))
            ww3(j1+npw2) = dconjg(ww3(j1))

            ww1x(j1)= eye*tx(j1)*dsq*ww1(j1)
            ww2y(j1)= eye*tx(j1)*dsq*ww2(j1)
            ww3z(j1)= eye*tx(j1)*dsq*ww3(j1)
            
            ww1x(j1+npw2)=dconjg(ww1x(j1))
            ww2y(j1+npw2)=dconjg(ww2y(j1))
            ww3z(j1+npw2)=dconjg(ww3z(j1))

            qq1 = qq1*qqx
            qq2 = qq2*qqy
            qq3 = qq3*qqz
         enddo
c
         do ind = 1,nd
            cd=0
            g1=0
            g2=0
            g3=0
            do j3=1,npw/2
               d1=0
               d1x=0
               d1y=0
               do j2=1,npw
                  d2=0
                  d2x=0
                  do j1=1,npw
                     d2 = d2+pwexp(j1,j2,j3,ind)*ww1(j1)
                     d2x = d2x+pwexp(j1,j2,j3,ind)*ww1x(j1)
                  enddo
                  d1 = d1+d2*ww2(j2)
                  d1x = d1x+d2x*ww2(j2)
                  d1y = d1y+d2*ww2y(j2)
               enddo
               cd=cd+d1*ww3(j3)
               g1=g1+d1x*ww3(j3)
               g2=g2+d1y*ww3(j3)
               g3=g3+d1*ww3z(j3)
            enddo
            
            pot(ind,itarg) = pot(ind,itarg)+dreal(cd)*2
            grad(ind,1,itarg) = grad(ind,1,itarg)+dreal(g1)*2
            grad(ind,2,itarg) = grad(ind,2,itarg)+dreal(g2)*2
            grad(ind,3,itarg) = grad(ind,3,itarg)+dreal(g3)*2
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
      subroutine g3dpwevalh_vec(nd,delta,center,npw,wx,tx,
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
      implicit real*8 (a-h,o-z)
      integer nd,npw,ntarg
      real *8 delta,center(3),targ(3,ntarg)
      real *8 pot(nd,ntarg),grad(nd,3,ntarg),hess(nd,6,ntarg)
      
      real *8 wx(npw),tx(npw)
      complex *16 pwexp(npw,npw,npw/2,nd)

      integer i,ind,j1,j2,j,npw2,itarg,k
      real *8 x,y,z,dsq
      complex *16 eye
      complex *16 qqx,qqy,qqz,qq1,qq2,qq3
      
      complex *16 ww1(100),ww1x(100),ww1xx(100)
      complex *16 ww2(100),ww2y(100),ww2yy(100)
      complex *16 ww3(100),ww3z(100),ww3zz(100)

      complex *16 dd,g1,g2,g3,d1,d1x,d1y,d2,d2x,d2xx
      complex *16 h11,h22,h33,h12,h13,h23
      complex *16 d1xx,d1yy,d1xy

C
      eye = dcmplx(0,1)
      dsq = 1.0D0/dsqrt(delta)
C
      npw2=npw/2
      do itarg=1,ntarg
         x = (targ(1,itarg) - center(1))*dsq
         y = (targ(2,itarg) - center(2))*dsq
         z = (targ(3,itarg) - center(3))*dsq

         qqx = cdexp(eye*tx(1)*x)
         qqy = cdexp(eye*tx(1)*y)
         qqz = cdexp(eye*tx(1)*z)
         qq1 = qqx
         qq2 = qqy
         qq3 = qqz
         qqx = qqx*qqx
         qqy = qqy*qqy
         qqz = qqz*qqz
         

         do j1=1,npw2
            ww1(j1) = qq1
            ww2(j1) = qq2            
            ww3(j1) = qq3            
            ww1(j1+npw2) = dconjg(ww1(j1))
            ww2(j1+npw2) = dconjg(ww2(j1))
            ww3(j1+npw2) = dconjg(ww3(j1))

            ww1x(j1)= eye*tx(j1)*dsq*ww1(j1)
            ww2y(j1)= eye*tx(j1)*dsq*ww2(j1)
            ww3z(j1)= eye*tx(j1)*dsq*ww3(j1)
            
            ww1x(j1+npw2)=dconjg(ww1x(j1))
            ww2y(j1+npw2)=dconjg(ww2y(j1))
            ww3z(j1+npw2)=dconjg(ww3z(j1))

            ww1xx(j1)= eye*tx(j1)*dsq*ww1x(j1)
            ww2yy(j1)= eye*tx(j1)*dsq*ww2y(j1)
            ww3zz(j1)= eye*tx(j1)*dsq*ww3z(j1)
            
            ww1xx(j1+npw2)=dconjg(ww1xx(j1))
            ww2yy(j1+npw2)=dconjg(ww2yy(j1))
            ww3zz(j1+npw2)=dconjg(ww3zz(j1))
            
            qq1 = qq1*qqx
            qq2 = qq2*qqy
            qq3 = qq3*qqz
         enddo
c
         do ind = 1,nd
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
            do j3=1,npw/2
               d1=0
               d1x=0
               d1y=0

               d1xx=0
               d1yy=0
               d1xy=0               
               do j2=1,npw
                  d2=0
                  d2x=0
                  d2xx=0
                  do j1=1,npw
                     d2=d2+pwexp(j1,j2,j3,ind)*ww1(j1)
                     d2x=d2x+pwexp(j1,j2,j3,ind)*ww1x(j1)
                     d2xx=d2xx+pwexp(j1,j2,j3,ind)*ww1xx(j1)                  
                  enddo
                  d1=d1+d2*ww2(j2)
                  d1x=d1x+d2x*ww2(j2)
                  d1xx=d1xx+d2xx*ww2(j2)
                  
                  d1y=d1y+d2*ww2y(j2)
                  d1xy=d1xy+d2x*ww2y(j2)
                  
                  d1yy=d1yy+d2*ww2yy(j2)
               enddo
               dd=dd+d1*ww3(j3)
               g1=g1+d1x*ww3(j3)
               g2=g2+d1y*ww3(j3)

               h11=h11+d1xx*ww3(j3)
               h22=h22+d1yy*ww3(j3)
               h12=h12+d1xy*ww3(j3)
               
               g3=g3+d1*ww3z(j3)
               h13=h13+d1x*ww3z(j3)
               h23=h23+d1y*ww3z(j3)
               
               h33=h33+d1*ww3zz(j3)
            enddo
            pot(ind,itarg) = pot(ind,itarg)+dreal(dd)*2
            
            grad(ind,1,itarg) = grad(ind,1,itarg)+dreal(g1)*2
            grad(ind,2,itarg) = grad(ind,2,itarg)+dreal(g2)*2
            grad(ind,3,itarg) = grad(ind,3,itarg)+dreal(g3)*2

            hess(ind,1,itarg) = hess(ind,1,itarg)+dreal(h11)*2
            hess(ind,2,itarg) = hess(ind,2,itarg)+dreal(h22)*2
            hess(ind,3,itarg) = hess(ind,3,itarg)+dreal(h33)*2

            hess(ind,4,itarg) = hess(ind,4,itarg)+dreal(h12)*2
            hess(ind,5,itarg) = hess(ind,5,itarg)+dreal(h13)*2
            hess(ind,6,itarg) = hess(ind,6,itarg)+dreal(h23)*2
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
c      
c*********************************************************************      
C
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
      subroutine pw_translation_matrices(xmin,npw,ts,nmax,
     1              wshift)
C
C     This subroutine precomputes all translation matrices for PW
C     expansions for mp to loc at the cutoff level.
c      
C     INPUT
C
c     xmin    = scaled (by 1/sqrt(delta) boxsize at the cutoff level
C     npw     = number of terms in 1d plane wave expansion
C     ts      = 1d pw expansion nodes
C     nmax    = number of different translation lengths in the whole scheme 
c      
C     OUTPUT:
C
C     wshift  = table of translation matrices for PW shift 
C
      implicit real *8 (a-h,o-z)
      real *8 ts(npw)
      
      complex *16 wshift(npw*npw*npw/2,-nmax:nmax,-nmax:nmax,-nmax:nmax)
      
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
      do k3=-nmax,nmax
         j=0   
         do j1=1,npw/2
         do j2=1,npw
         do j3=1,npw
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
      subroutine merge_split_pw_matrices(xmin,npw,ts,nmax,
     1           wshift)
C
C     This subroutine precomputes all translation matrices for 
c     PW translations from child to parent or vice versa.
C
C     INPUT
C
c     xmin     = half of the scaled (by 1/sqrt(delta) size of the box 
c                at the finest level
C     npw      = number of terms in 1d PW expansion
C     ws,ts    = real *8, 1d PW expansion weights and nodes
C     nmax     = number of different translation lengths in the whole scheme 
C
C     OUTPUT:
C
C     wshift   = table of translation matrices for PW  shift used in
c                merge mp and split loc stage
C
      implicit real *8 (a-h,o-z)
      real *8 xmin
      complex *16 wshift(npw*npw*npw/2,8,nmax)
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
         j=0
         do j1=1,npw/2
         do j2=1,npw
         do j3=1,npw
            j=j+1
c           pp p              
            wshift(j,1,k1) = ww(j3,k1)*ww(j2,k1)
     1          *ww(j1,k1)
c           pm p
            wshift(j,2,k1) = ww(j3,k1)*conjg(ww(j2,k1))
     1          *ww(j1,k1)
c           mp p
            wshift(j,3,k1) = conjg(ww(j3,k1))*ww(j2,k1)
     1          *ww(j1,k1)
c           mm p
            wshift(j,4,k1) = conjg(ww(j3,k1))*conjg(ww(j2,k1))
     1          *ww(j1,k1)
            
c           pp m              
            wshift(j,5,k1) = ww(j3,k1)*ww(j2,k1)
     1          *conjg(ww(j1,k1))
c           pm m
            wshift(j,6,k1) = ww(j3,k1)*conjg(ww(j2,k1))
     1          *conjg(ww(j1,k1))
c           mp m
            wshift(j,7,k1) = conjg(ww(j3,k1))*ww(j2,k1)
     1          *conjg(ww(j1,k1))
c           mm m
            wshift(j,8,k1) = conjg(ww(j3,k1))*conjg(ww(j2,k1))
     1          *conjg(ww(j1,k1))
            
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
C     nexp      = number of terms in PW expansion
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
C     This subroutine adds one PW expansion (pwexp1) 
C     to another PW expansion (pwexp2).
C
C     INPUT
C
c     nd        = vector length (for vector input)
C     nexp      = total number of terms in 3D PW expansion
C     pwexp1    = original expansion 
C
C     OUTPUT:
C
C     pwexp2 = the pw expansion is incremented 
C
      implicit none
      integer nd,j,ind,nexp
      complex *16 pwexp1(nexp,nd)
      complex *16 pwexp2(nexp,nd)

C
      do ind=1,nd
         do j=1,nexp
            pwexp2(j,ind) = pwexp2(j,ind) + pwexp1(j,ind)
         enddo
      enddo
c
      return
      end
c
C
c
C
      subroutine g3dshiftpw0_vec(nd,delta,nn,pwexp1,
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
      integer nd,nn,j,j1,j2,j3,ind,nexp
      real *8 cent1(3),cent2(3),x,y,z,delta,dsq
      complex *16 pwexp1(nn/2,nn,nn,nd)
      complex *16 pwexp2(nn/2,nn,nn,nd)
      complex *16, allocatable :: wshift(:,:,:)
      complex *16 ww1(100),ww2(100),ww3(100),eye
      real *8 ws(nn),ts(nn)
C
      allocate(wshift(nn/2,nn,nn))
      eye = dcmplx(0,1)
      dsq = 1.0D0/dsqrt(delta)
      x = (cent2(1) - cent1(1))*dsq
      y = (cent2(2) - cent1(2))*dsq
      z = (cent2(3) - cent1(3))*dsq
      
      do j1=1,nn/2
         ww1(j1) = cdexp(eye*ts(j1)*x)
         ww2(j1) = cdexp(eye*ts(j1)*y)
         ww3(j1) = cdexp(eye*ts(j1)*z)
         
         ww2(nn/2+j1) = dconjg(ww2(j1))
         ww3(nn/2+j1) = dconjg(ww3(j1))
      enddo
c
      do j3=1,nn
      do j2=1,nn
      do j1=1,nn/2
         wshift(j1,j2,j3) = ww1(j1)*ww2(j2)*ww3(j3)
      enddo
      enddo
      enddo
      
      do ind=1,nd
         do j3=1,nn
         do j2=1,nn
         do j1=1,nn/2
            pwexp2(j1,j2,j3,ind) = pwexp2(j1,j2,j3,ind)
     1          +pwexp1(j1,j2,j3,ind)*wshift(j1,j2,j3)
         enddo
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
C************************************************************************
C
C 1d translation matrices
C
C**********************************************************************
      subroutine g3dh2lmat(nlocal,ntermsh,h2l)
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
      subroutine g3dh2xmat(ntermsh,npw,ws,ts,h2x)
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
cccc      complex *16 h2x(npw,0:ntermsh)
      complex *16 h2x(0:ntermsh,npw)
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
cccc            h2x(i,j)=ws(i)*(-eye*ts(i))**j*fac(j)
            h2x(j,i)=ws(i)*(-eye*ts(i))**j*fac(j)
         enddo
      enddo
      
      return
      end
C
C
C
C
      subroutine g3drh2xmat(ntermsh,npw,ws,ts,h2x)
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
cccc      real *8 h2x(npw,0:ntermsh)
      real *8 h2x(0:ntermsh,npw)
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
cccc            h2x(i,j)=ws(i)*(-ts(i))**j*fac(j)
            h2x(j,i)=ws(i)*(-ts(i))**j*fac(j)
         enddo
      enddo
      
      return
      end
C
C
C
C            
      subroutine g3dx2lmat(nlocal,npw,ws,ts,x2l)
C
C     This subroutine returns the matrices converting the planewave expansions
C     to the local expansion
C
C     INPUT
C
C     nlocal       = number of terms in local exp
C     npw          = number of terms in the planvewave exp 
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
      subroutine g3drx2lmat(nlocal,npw,ws,ts,x2l)
C
C     This subroutine returns the matrices converting the planewave expansions
C     to the local expansion
C
C     INPUT
C
C     nlocal       = number of terms in local exp
C     npw          = number of terms in the planvewave exp 
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
