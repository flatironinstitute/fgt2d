cc Copyright (C) 2021: Leslie Greengard, Shidong Jiang, Manas Rachh
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date: 2021-03-22 (Monday, 22 March 2021) $
c    $Revision:  $
c
c
c-----------------------------------------------------------------------------
c
c      g2dterms - determine number of terms in Hermite/local expansions 
c           for boxes of size "bsize" with Gaussian variance "delta".
c
c-----------------------------------------------------------------------------
c
c
c
      subroutine g2dterms(ndim,bsize,delta,eps,nterms)
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
c     eps: desired precision
c     delta: Gaussian variance
c     bsize: side length of the box
c    
c     OUTPUT:
c     nterms: number of terms in Taylor/Hermite expansion
c      
      implicit real*8 (a-h,o-z)
      integer nterms,p0
      real*8 delta, bsize
c
      e1 = exp(1.0d0)

      r = bsize/sqrt(2.0d0*delta)

      p0 = ceiling(e1*r*r) + 1

      
      rk = 1.09*(8*atan(1.0d0))**(-0.25d0)

      
      do p=p0,300
         rp = r*sqrt(e1/p)
         error=ndim*rk*((p*1.0d0)**(-0.25d0))*(rp**p)/(1-rp)
         if (error .lt. eps) then
            nterms = p
            exit
         endif
      enddo

      print *, p0, r*sqrt(2.0d0), nterms
      
      return
      end 
c
c
c
c

      
