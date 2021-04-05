ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date: 2021-03-22 (Monday, 22 March 2021) $
c    $Revision:  $
c
c
c-----------------------------------------------------------------------------
c
c      g2dhlterms - determine number of terms in Hermite/local expansions 
c           for boxes of size "bsize" with Gaussian variance "delta".
c
c      g2dsoeterms - determine number of terms in the SOE expansion 
c           for a given precision
c
c      g2dpwterms - determine number of terms in the plane wave expansions 
c           for boxes of size "bsize" with Gaussian variance "delta".
c
c-----------------------------------------------------------------------------
c
c
c
      subroutine g2dhlterms(ndim,bsize,delta,eps,nterms)
c      
c     Determine the number of terms in the Taylor/Hermite expansions 
c     for boxes with side length = bsize
c
c     The error formula is from Wan and Karniadakis 2006 JCP 219 pp. 7-12.
c     Strain derived a similar formula in SIAM J. Sci. Stat. Comput. 12 (5) 
c     (1991) 1131–1139, essentially with only "ndim" missing.
c
c
c
c
c     INPUT:
c     ndim : dimension of the underlying space
c     bsize: side length of the box
c     delta: Gaussian variance
c     eps: desired precision
c    
c     OUTPUT:
c     nterms: number of terms in Taylor/Hermite expansion
c      
      implicit real*8 (a-h,o-z)
      integer nterms,p0,p
      real*8 delta, bsize
c
      e1 = exp(1.0d0)

      r = bsize/sqrt(2.0d0*delta)

      p0 = int(e1*r*r) + 1

      
      rk = 1.09*(8*atan(1.0d0))**(-0.25d0)

      
      do p=p0,300
         rp = r*sqrt(e1/p)
         error=ndim*rk*((p*1.0d0)**(-0.25d0))*(rp**p)/(1-rp)
         if (error .lt. eps) then
            nterms = p
            exit
         endif
      enddo

cccc      print *, p0, r*sqrt(2.0d0), nterms
      
      return
      end 
c
c
c
c
      subroutine g2dsoeterms(eps,nsoe)
c      
c     Determine the number of terms in the SOE expansion 
c     for a given precision
c
c     The SOE approximation is valid globally, independent
c     of box size and delta.
c
c
c     INPUT:
c     eps: desired precision
c    
c     OUTPUT:
c     nsoe: number of terms in the SOE expansion
c      
      implicit real*8 (a-h,o-z)

      if (eps .gt. 4d-4) then
         nsoe = 4
      elseif (eps .gt. 7d-6) then
         nsoe = 6
      elseif (eps .gt. 7d-8) then
         nsoe = 8
      elseif (eps .gt. 2d-9) then
         nsoe = 10
      elseif (eps .gt. 4d-11) then
         nsoe = 12
      elseif (eps .gt. 6d-13) then
         nsoe = 14
      elseif (eps .gt. 5d-14) then
         nsoe = 18
      else
         nsoe = 18
         write(*,*) 'higher precision of SOE approx not implemented'
         write(*,*) 'use nsoe=18, with roughly 13-digit accuracy'
      endif
      
      return
      end 
c
c
c
c
      subroutine g2dpwterms(bsize,delta,eps,pmax,npw)
c      
c     Determine the number of terms in the plane wave expansions 
c     for boxes with side length = bsize
c
c
c     Method: 
c
c     Let G(x)=exp(-x^2), i.e., x is already scaled by sqrt(delta).
c     Then its Fourier transform is \hat{G}(k) = exp(-k^2/4).
c      
c     By the Poisson summation formula, we have
c      
c     \sum G(x+2\pi n/h) = \frac{h}{2sqrt{\pi} \sum \hat{G}(m h) exp(i m h x)
c
c     where both sums are from -\infty to +\infty.
c
c     (1) Truncation on the Fourier series on the RHS requires that
c
c            \hat{G}(pmax) < eps.
c
c     This leads to 
c
c            pmax = 2*sqrt(log(1/eps)).
c
c     (2) We also need to make sure the aliasing error on the LHS is small, 
c     that is,
c
c            G(2\pi/h - rmax)< eps.
c
c     This leads to
c
c           h = 2 \pi / (rmax + sqrt(log(1/eps))).
c
c     INPUT:
c     bsize: side length of the box
c     delta: Gaussian variance
c     eps: desired precision
c    
c     OUTPUT:
c     npw: number of terms in the plane wave expansion
c     pmax : maximum plane wave number
c
      implicit real*8 (a-h,o-z)
      real*8 delta, bsize
c
      pi = 4.0d0*atan(1.0d0)
      rmax = bsize/sqrt(delta)

      d = sqrt(log(1.0d0/eps))

      pmax = 2.0d0*d

      h = 2*pi/(rmax+d)

      npw = int(pmax/h)

cccc      print *, h, pmax, pmax/h, npw
      
      return
      end 
c
c
c
c

      
