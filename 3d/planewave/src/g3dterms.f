c-----------------------------------------------------------------------------
c
c      g3dhlterms - determine number of terms in Hermite/local expansions 
c           for boxes of size "bsize" with Gaussian variance "delta".
c
c      g3dpwterms - determine number of terms in the plane wave expansions 
c           for boxes of size "bsize" with Gaussian variance "delta".
c
c-----------------------------------------------------------------------------
c
c
c
      subroutine g3dhlterms(ndim,bsize,delta,eps,nterms)
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
      subroutine g3dpwterms(bsize,delta,eps,pmax,npw)
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

cccc      npw = int(pmax/h)+1
      npw = int(pmax/h)

      npw = 2*npw
cccc      print *, h, pmax, pmax/h, npw
      
      return
      end 
c
c
c
c

      
