c     This file contains the following subroutines
c
c     g2dformsxc_vec : computes the SOE expansion about the center CENT
c                     due to sources of strength charge()
c
C     g2dsxevalp_vec : evaluates the SOE expansion 
C                         potential only.
c      
c     g2dshiftsx_vec : shift one SOE/X expansion
c                     
c     g2dshiftall_vec : shift all four SOE/X expansions
c
c     get_pwnodes : returns plane wave approximation nodes and weights, obtained by
C                   the trapezoidal rule
c
c*********************************************************************
C
C form SOE/X expansions (charge only)
C
C*********************************************************************
      subroutine g2dformsxc_vec(nd,delta,sources,ns,charge,center,pmax,
     1            npw,nns,ws,ts,spall)
C
C     This subroutine computes the four hybrid 
C     SOE/plane wave expansions about CENTER
C     due to sources of strength charge().
C
C     INPUT
C
c     nd            = vector length (for multiple charges at same locations)
C     delta         = Gaussian variance
C     sources(2,ns) = coordinates of sources
C     ns            = number of sources
C     charge        = strength of sources
C     center        = center of the expansion
C     pmax          = range of Fourier integral for plane waves
C     npw           = number of Fourier plane waves
C     nns           = number of SOE terms (half of
C                     full quadrature - no complex conjugates)
C     ws,ts         = SOE weights and nodes
C
C     OUTPUT:
C
C     spall    = hybrid expansions in order px,mx,xp,xm
C
      implicit none
      integer ninbox,npw,nd,ns,nns,i,ii,ind,j1,j2,j,nexp
      real *8 center(2),sources(2,ns),charge(nd,ns),pmax
      complex *16 spall(4*(2*npw+1)*nns,nd)
      real *8 x,xp,xq,delta,dsq,y,pi,w
      real *8, allocatable :: hexpx(:,:),hexpy(:,:)
      real *8 ww(-100:100)
      complex *16 qqq,eye
      complex *16 qqx,qqy,qq1,qq2
      complex *16 ws(nns),ts(nns)
      complex *16 ww1p(6),ww2p(6),ww1m(6),ww2m(6)
      complex *16, allocatable :: wtall(:,:)
      complex *16, allocatable :: ww1c(:,:),ww2c(:,:)
      complex *16 ww1(-100:100)
      complex *16 ww2(-100:100)
C
C     initialize coefficients to zero.
C
      allocate(wtall(4*(2*npw+1)*nns,nd))
c      
      pi = 4.0d0*datan(1.0d0)
      eye = dcmplx(0,1)
      w = pmax/(2*2.0d0*npw*dsqrt(pi))
      do j =-npw,npw
         xq = j*pmax/npw
         ww(j) = 2*w*dexp(-xq*xq/4)
      enddo
c
      dsq = 1.0D0/dsqrt(delta)
      allocate(ww1c(-npw:npw,nd))
      allocate(ww2c(-npw:npw,nd))
c
      do i=1,ns
         x = (sources(1,i) - center(1))*dsq
         y = (sources(2,i) - center(2))*dsq
         do j1=1,nns
            qqq = cdexp(ts(j1)*x)
            ww1p(j1) = ws(j1)*qqq
            ww1m(j1) = ws(j1)/qqq
            qqq = cdexp(ts(j1)*y)
            ww2p(j1) = ws(j1)*qqq
            ww2m(j1) = ws(j1)/qqq
         enddo
c
         qqx = cdexp(-eye*pmax*x/npw)
         qqy = cdexp(-eye*pmax*y/npw)
         qq1 = 1.0d0
         qq2 = 1.0d0
         ww1(0) = ww(0)
         ww2(0) = ww(0)
         do j1=1,npw
            qq1 = qq1*qqx
            qq2 = qq2*qqy
            ww1(j1) = ww(j1)*qq1
            ww2(j1) = ww(j1)*qq2
            ww1(-j1) = dconjg(ww1(j1))
            ww2(-j1) = dconjg(ww2(j1))
         enddo
         do ind=1,nd
            do j1=-npw,npw
             ww1c(j1,ind)=ww1(j1)*charge(ind,i)
             ww2c(j1,ind)=ww2(j1)*charge(ind,i)
            enddo
         enddo
c
         do ind=1,nd
            nexp = (2*npw+1)*nns
            j = 0
            do j1=1,nns
            do j2=-npw,npw
               j = j+1
               wtall(j,ind) = ww1p(j1)*ww2c(j2,ind)
               wtall(nexp+j,ind) = ww1m(j1)*ww2c(j2,ind)
               wtall(2*nexp+j,ind) = ww2p(j1)*ww1c(j2,ind)
               wtall(3*nexp+j,ind) = ww2m(j1)*ww1c(j2,ind)
            enddo
            enddo
         enddo
c
         do ind=1,nd
         do j=1,4*(2*npw+1)*nns
               spall(j,ind) = spall(j,ind) + wtall(j,ind)
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
C evaluate SOE/X expansions (potential only)
C
C*********************************************************************
      subroutine g2dsxevalp_vec(nd,delta,center,pmax,npw,nns,
     1              spall,targ,ntarg,ws,ts,pot)
C
C     This subroutine evaluates the hybrid SOE/plane wave 
C     expansions about CENTER at location TARG
C
C     INPUT
C
c     nd            = vector length (for multiple charges at same locations)
C     delta         = Gaussian variance
C     charge        = strength of sources
C     center        = center of the expansion
C     pmax          = range of Fourier integral for plane waves
C     npw           = number of Fourier plane waves
C     nns           = number of SOE terms (half of
C                     full quadrature - no complex conjugates)
C     spall         = hybrid expansions in order px,mx,xp,xm
C     targ          = target
C     ws,ts         = SOE weights and nodes
C
C     OUTPUT:
C
C     pot           = potential is incremented
C
      implicit none
      integer nd,nns,i,ind,j1,j2,j,npw,nexp,itarg,ntarg
      real *8 center(2),targ(2,ntarg)
      real *8 x,y,xpoint,pot(nd,ntarg),delta,pmax,dsq
      real *8 facx,facy,pi,w,xq
      complex *16 z,eye
      complex *16 spall(4*(2*npw+1)*nns,nd)
      complex *16 ws(nns),ts(nns)
      complex *16 ww1p(20),ww2p(20),ww1m(20),ww2m(20)
      complex *16, allocatable :: wtall(:)
      complex *16, allocatable :: sumz(:)
      complex *16 ww1(-100:100)
      complex *16 ww2(-100:100)
C
C     initialize coefficients to zero.
C
      allocate(wtall(4*(2*npw+1)*nns))
      allocate(sumz(nd))
c
      dsq = 1.0D0/dsqrt(delta)
C
c
      pi = 4.0d0*datan(1.0d0)
      eye = dcmplx(0,1)
      w = pmax/(2*2.0d0*npw*dsqrt(pi))
      do itarg=1,ntarg
         x = (targ(1,itarg) - center(1))*dsq
         y = (targ(2,itarg) - center(2))*dsq
         do j =-npw,0
            xq = j*pmax/npw
            ww1(j) = cdexp(eye*xq*x)
            ww2(j) = cdexp(eye*xq*y)
         enddo

         do j=1,npw
            ww1(j) = conjg(ww1(-j))
            ww2(j) = conjg(ww2(-j))
         enddo
         
         do j1=1,nns
            ww1p(j1) = cdexp(-ts(j1)*x)
            ww2p(j1) = cdexp(-ts(j1)*y)
            ww1m(j1) = 1/ww1p(j1)
            ww2m(j1) = 1/ww2p(j1)
         enddo
c
         j = 0
         nexp = (2*npw+1)*nns
         do j1=1,nns
            do j2=-npw,npw
               j = j+1
               wtall(j) = ww1p(j1)*ww2(j2)
               wtall(nexp+j) = ww1m(j1)*ww2(j2)
               wtall(2*nexp+j) = ww2p(j1)*ww1(j2)
               wtall(3*nexp+j) = ww2m(j1)*ww1(j2)
            enddo
         enddo
c
         do ind=1,nd
            sumz(ind) = 0.0d0
            do j=1,4*nexp
               sumz(ind) = sumz(ind) +
     1             wtall(j)*spall(j,ind)
            enddo
            pot(ind,itarg) = pot(ind,itarg) + dble(sumz(ind))
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
      real *8 ws(-npw:npw),ts(-npw:npw)

      pi = 4.0d0*datan(1.0d0)
      h = pmax/npw
      w = h/(2.0d0*dsqrt(pi))

      do j =-npw,npw
         ts(j) = j*h
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
C shift SOE/X expansions (mp to loc, mp to mp, loc to loc)
C
C*********************************************************************
c
C
C
c            
      subroutine shiftsx_translation_matrices(xmin,npw,nns,nmax,
     1              wshift,ts,tx)
C
C     This subroutine precomputes all translation matrices for all SOE/X
C     expansions from child to parent or vice versa.
C
c     used in merge mp or split loc stage
c      
C     INPUT
C
c     nd      = vector length (for multiple charges at same locations)
C     npw     = number of terms in plane wave exp
C     nns    = number of SOE terms
C     nmax    = number of different translation lengths in the whole scheme 
C     ts      = SOE nodes
C     tx      = X nodes
c      
C     OUTPUT:
C
C     wshift  = table of translation matrices for SOE/X shift used in
c                merge mp and split loc stage
C
      implicit real *8 (a-h,o-z)
      complex *16 ts(nns)
      real *8 tx(-npw:npw)
      
      complex *16 wshift((2*npw+1)*nns,0:nmax,-nmax:nmax)
      
      complex *16,allocatable:: wws(:,:)
      complex *16,allocatable:: wwx(:,:)

      complex *16 eye
C
      eye = dcmplx(0,1)
      
      allocate(wws(nns,0:nmax))
      allocate(wwx(-npw:npw,-nmax:nmax))
      
      do j2=0,nmax
         do j1=1,nns
            wws(j1,j2)=exp(-ts(j1)*xmin*j2)
         enddo
      enddo

      do j2=-nmax,nmax
         do j1=-npw,npw
            wwx(j1,j2)=exp(eye*tx(j1)*xmin*j2)
         enddo
      enddo

      nexp=(2*npw+1)*nns
      do m=0,nmax
      do n=-nmax,nmax
         j = 0
         do j1=1,nns
         do j2=-npw,npw
            j = j+1
            wshift(j,m,n) = wws(j1,m)*wwx(j2,n)
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
      subroutine g2dshiftsx_vec(nd,npw,nns,sxall1,
     1              sxall2,wshift,isx)
C
C     This subroutine converts the hybrid SOE/X expansion
C     (shall1) about the center (CENT1) into a hybrid SOE/X 
C     expansion about (CENT2) for px direction.
C
C     INPUT
C
c     nd      = vector length (for multiple charges at same locations)
C     delta   = Gaussian variance
C     npw     = number of terms in plane wave exp
C     nns     = number of SOE terms
C     sxall1  = original expansion FFEXP
C     wshift   = precomputed shift array
C     isx     = SOE/X expansion type index
c               0-px, 1-mx, 2-xp, 3-xm
C
C     OUTPUT:
C
C     sxall2  = shifted expansion 
C
      implicit real *8 (a-h,o-z)
      complex *16 sxall1(4*(2*npw+1)*nns,nd)
      complex *16 sxall2(4*(2*npw+1)*nns,nd)
      complex *16 wshift((2*npw+1)*nns)
C
      nexp = (2*npw+1)*nns
      do ind=1,nd
      do j=1,nexp
         sxall2(j+isx*nexp,ind) = sxall2(j+isx*nexp,ind)
     1       +sxall1(j+isx*nexp,ind)*wshift(j)
      enddo
      enddo
c
      return
      end
c
C
C
c
      subroutine g2dshiftsxall2_vec(nd,npw,nns,sxall,
     1              wshift)
C
C     This subroutine converts the hybrid SOE/X expansion
C     (shall1) about the center (CENT1) into a hybrid SOE/X 
C     expansion about (CENT2) for px direction.
C
C     INPUT
C
c     nd      = vector length (for multiple charges at same locations)
C     delta   = Gaussian variance
C     npw     = number of terms in plane wave exp
C     nns     = number of SOE terms
C     sxall  = original expansion FFEXP
C     wshift   = precomputed shift array
C
C     OUTPUT:
C
C     sxall2  = shifted expansion 
C
      implicit real *8 (a-h,o-z)
      complex *16 sxall(4*(2*npw+1)*nns,nd)
      complex *16 wshift((2*npw+1)*nns)
C
      nexp = (2*npw+1)*nns
      do ind=1,nd
        do isx=0,3
          do j=1,nexp
            sxall(j+isx*nexp,ind) = sxall(j+isx*nexp,ind)*wshift(j)
          enddo
        enddo
      enddo
c
      return
      end
c
C
C
c            
      subroutine g2dcopysx_vec(nd,npw,nns,sxall1,
     1              sxall2,isx)
C
C     This subroutine converts the hybrid SOE/X expansion
C     (shall1) about the center (CENT1) into a hybrid SOE/X 
C     expansion about (CENT2) for px direction.
C
C     INPUT
C
c     nd      = vector length (for multiple charges at same locations)
C     delta   = Gaussian variance
C     npw     = number of terms in plane wave exp
C     nns     = number of SOE terms
C     sxall1  = original expansion FFEXP
C
C     OUTPUT:
C
C     sxall2  = shifted expansion 
C
      implicit real *8 (a-h,o-z)
      complex *16 sxall1(4*(2*npw+1)*nns,nd)
      complex *16 sxall2(4*(2*npw+1)*nns,nd)
C
      nexp = (2*npw+1)*nns
      do ind=1,nd
         do j=1,nexp
            sxall2(j+isx*nexp,ind) = sxall1(j+isx*nexp,ind)
         enddo
      enddo
c
      return
      end
c
c
c     
c
C***********************************************************************
      subroutine g2dsxzero_vec(nd,sxall,npw,nsoehalf)
      implicit none
C***********************************************************************
c
c     This subroutine sets a vector multipole expansion to zero.
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :   vector length (number of mpole expansions)
c     nterms :   order of multipole expansion
C---------------------------------------------------------------------
c     OUTPUT:
c
c     mpole  :   coeffs for the expansion set to zero.
C---------------------------------------------------------------------
      integer n,nsoehalf,npw,nd,ii
      complex *16 sxall(4*(2*npw+1)*nsoehalf,nd)
c
      do ii=1,nd
      do n=1,4*(2*npw+1)*nsoehalf
         sxall(n,ii)=0.0d0
      enddo
      enddo
      return
      end
c
C

