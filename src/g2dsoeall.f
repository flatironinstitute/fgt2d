c     This file contains the following subroutines
c
c     g2dformsc_vec : computes the SOE expansion about the center CENT
c                     due to sources of strength charge()
c
c     g2dformsd_vec : computes the SOE expansion about the center CENT
c                     due to dipoles
c
c     g2dformscd_vec : computes the SOE expansion about the center CENT
c                     due to charges and dipoles
c
C     g2dsevalp_vec : evaluates the SOE expansion 
C                         potential only.
c      
C     g2dsevalg_vec : evaluates the SOE expansion 
C                         pot/grad
c      
C     g2dsevalh_vec : evaluates the SOE expansion 
C                         pot/grad/Hessian
c      
c     g2dshiftsoe_vec : shift one SOE expansion
c                     
c     g2dshiftsoeall_vec : shift all four SOE expansions
c
c     get_soe1 : returns SOE approximation nodes and weights, obtained by
C                Trefethen's code
C     
c     get_soe2 : returns SOE approximation nodes and weights, obtained by
C                least squares optimization with get_soe1 as the initial 
C                guess for nsoe <=14. 
C                
C                nsoe = 18 is obtained by optimal contour + balanced trunction
c                optimization, which has better Hessian accuracy as compared
c                with nsoe = 14.
c
c*********************************************************************
C
C form SOE expansions (charge, dipole, charge & dipole)
C
C*********************************************************************
      subroutine g2dformsc_vec(nd,delta,sources,ns,charge,center,
     1            nn,ws,ts,soeall)
C
C     This subroutine computes the four (SOE=Sum-Of-Exponential)
C     expansions about CENTER due to sources of strength charge().
C
C     INPUT
C
c     nd            = vector length (for vector input)
C     delta         = Gaussian variance
C     sources(2,ns) = coordinates of sources
C     ns            = number of sources
C     charge        = strength of sources
C     center        = center of the expansion
C     nn            = number of exponentials 
C     ws,ts         = SOE weights and nodes
C
C     OUTPUT:
C
C     soeall     = all four (unrolled) SOE expansions are incremented
c                     in the order: soepp,soepm,soemp,soemm  
C
      implicit none
      integer ninbox,nterms,nd,ns,nn,i,ii,ind,j1,j2,j,nexp
      real *8 center(2),sources(2,ns),charge(nd,ns),pmax
      complex *16 soeall(4*nn*nn/2,nd)
      real *8 x,xp,xq,delta,dsq,y
      complex *16 ws(nn),ts(nn)
      complex *16, allocatable :: ww1p(:,:),ww1m(:,:)
      complex *16 ww2p(12),ww2m(12)
      complex *16 qqq,ccc
      complex *16, allocatable :: wtall(:,:)
C
C     initialize coefficients to zero.
C
      allocate(wtall(4*72,nd))
      allocate(ww1p(12,nd))
      allocate(ww1m(12,nd))
      dsq = 1.0D0/dsqrt(delta)
C
      do i=1,ns
         x = (sources(1,i) - center(1))*dsq
         y = (sources(2,i) - center(2))*dsq
         do j1=1,nn/2
            qqq = cdexp(ts(j1)*x)
ccc            qqq = 1.0d0
            do ind = 1,nd
               ccc = charge(ind,i)*ws(j1)
               ww1p(j1,ind) = ccc*qqq
               ww1m(j1,ind) = ccc/qqq
            enddo
            qqq = cdexp(ts(j1)*y)
            ww2p(j1) = ws(j1)*qqq
            ww2p(nn/2+j1) = dconjg(ww2p(j1))
            ww2m(j1) = ws(j1)/qqq
            ww2m(nn/2+j1) = dconjg(ww2m(j1))
         enddo
ccc         write(6,*) ' ww1p ',ww1p
ccc         write(6,*) ' ww2p ',ww2p
c
         nexp = nn*nn/2
         do ind = 1,nd
            j = 0
            do j1=1,nn/2
            do j2=1,nn
               j = j+1
               wtall(j,ind) = ww1p(j1,ind)*ww2p(j2)
               wtall(nexp+j,ind) = ww1p(j1,ind)*ww2m(j2)
               wtall(2*nexp+j,ind) = ww1m(j1,ind)*ww2p(j2)
               wtall(3*nexp+j,ind) = ww1m(j1,ind)*ww2m(j2)
            enddo
            enddo
         enddo
c
         do ind=1,nd
         do j=1,4*nn*nn/2
            soeall(j,ind) = soeall(j,ind) +
     1         wtall(j,ind)
         enddo
         enddo
      enddo
      return
      end
C
C
C
C
      subroutine g2dformsd_vec(nd,delta,sources,ns,rnormal,dipstr,
     1            center,nn,ws,ts,soeall)
C
C     This subroutine computes the four SOE expansions about
C     CENTER due to dipole sources of strength dipstr().
C
C     INPUT
C
c     nd            = vector length (for vector input)
C     delta         = Gaussian variance
C     sources(2,ns) = coordinates of sources
C     ns            = number of sources
C     rnormal       = normals
C     dipstr        = dipole strengths
C     center        = center of the expansion
C     nn            = number of exponentials 
C     ws,ts         = SOE weights and nodes
C
C     OUTPUT:
C
C     soeall     = all four (unrolled) SOE expansions are incremented
c                     in the order: pp,pm,mp,mm  
C
      implicit none
      integer ninbox,nterms,nd,ns,nn,i,ii,ind,j1,j2,j,nexp
      real *8 center(2),sources(2,ns),dipstr(nd,ns),pmax
      real *8 rnormal(2,ns)
      complex *16 soeall(4*nn*nn/2,nd)
      real *8 x,xp,xq,delta,dsq,y
      complex *16 ws(nn),ts(nn)
      complex *16, allocatable :: ww1p(:,:),ww1m(:,:)
      complex *16, allocatable :: ww1px(:,:),ww1mx(:,:)
      complex *16 ww2p(12),ww2m(12)
      complex *16 ww2py(12),ww2my(12)
C
      complex *16 qqq,ccc
      complex *16, allocatable :: wtxall(:,:)
      complex *16, allocatable :: wtyall(:,:)
C
C     initialize coefficients to zero.
C
      allocate(wtxall(4*72,nd))
      allocate(wtyall(4*72,nd))
      allocate(ww1p(12,nd))
      allocate(ww1m(12,nd))
      allocate(ww1px(12,nd))
      allocate(ww1mx(12,nd))
      dsq = 1.0D0/dsqrt(delta)
C
      do i=1,ns
         x = (sources(1,i) - center(1))*dsq
         y = (sources(2,i) - center(2))*dsq
         do j1=1,nn/2
            qqq = cdexp(ts(j1)*x)
ccc            qqq = 1.0d0
            do ind = 1,nd
               ccc = dipstr(ind,i)*ws(j1)
               ww1p(j1,ind) = ccc*qqq
               ww1px(j1,ind) = ts(j1)*dsq*ww1p(j1,ind)
               ww1m(j1,ind) = ccc/qqq
               ww1mx(j1,ind) = -ts(j1)*dsq*ww1m(j1,ind)
            enddo
            qqq = cdexp(ts(j1)*y)
            ww2p(j1) = ws(j1)*qqq
            ww2py(j1) = ts(j1)*dsq*ww2p(j1)
            ww2p(nn/2+j1) = dconjg(ww2p(j1))
            ww2py(nn/2+j1) = dconjg(ww2py(j1))
            ww2m(j1) = ws(j1)/qqq
            ww2my(j1) = -ts(j1)*dsq*ww2m(j1)
            ww2m(nn/2+j1) = dconjg(ww2m(j1))
            ww2my(nn/2+j1) = dconjg(ww2my(j1))
         enddo
c
         nexp = nn*nn/2
         do ind = 1,nd
            j = 0
            do j1=1,nn/2
            do j2=1,nn
               j = j+1
               wtxall(j,ind) = ww1px(j1,ind)*ww2p(j2)
               wtxall(nexp+j,ind) = ww1px(j1,ind)*ww2m(j2)
               wtxall(2*nexp+j,ind) = ww1mx(j1,ind)*ww2p(j2)
               wtxall(3*nexp+j,ind) = ww1mx(j1,ind)*ww2m(j2)
               wtyall(j,ind) = ww1p(j1,ind)*ww2py(j2)
               wtyall(nexp+j,ind) = ww1p(j1,ind)*ww2my(j2)
               wtyall(2*nexp+j,ind) = ww1m(j1,ind)*ww2py(j2)
               wtyall(3*nexp+j,ind) = ww1m(j1,ind)*ww2my(j2)
            enddo
            enddo
         enddo
c
         do ind=1,nd
         do j=1,4*nn*nn/2
            soeall(j,ind) = soeall(j,ind) +
     1         wtxall(j,ind)*rnormal(1,i)+wtyall(j,ind)*rnormal(2,i)
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
      subroutine g2dformscd_vec(nd,delta,sources,ns,charge,
     1            rnormal,dipstr,center,nn,ws,ts,soeall)
C
C     This subroutine computes the four SOE expansions about
C     CENTER due to dipole sources of strength dipstr().
C
C     INPUT
C
c     nd            = vector length (for vector input)
C     delta         = Gaussian variance
C     sources(2,ns) = coordinates of sources
C     ns            = number of sources
C     charge        = charge strengths
C     rnormal       = normals
C     dipstr        = dipole strengths
C     center        = center of the expansion
C     nn            = number of exponentials 
C     ws,ts         = SOE weights and nodes
C
C     OUTPUT:
C
C     soeall     = all four (unrolled) SOE expansions are incremented
c                     in the order: pp,pm,mp,mm  
C
      implicit none
      integer ninbox,nterms,nd,ns,nn,i,ii,ind,j1,j2,j,nexp
      real *8 center(2),sources(2,ns),dipstr(nd,ns),pmax
      real *8 charge(nd,ns)
      real *8 rnormal(2,ns)
      complex *16 soeall(4*nn*nn/2,nd)
      real *8 x,xp,xq,delta,dsq,y
      complex *16 ws(nn),ts(nn)
      complex *16, allocatable :: ww1pc(:,:),ww1mc(:,:)
      complex *16, allocatable :: ww1p(:,:),ww1m(:,:)
      complex *16, allocatable :: ww1px(:,:),ww1mx(:,:)
      complex *16 ww2p(12),ww2m(12)
      complex *16 ww2py(12),ww2my(12)
C
      complex *16 qqq,ccc
      complex *16, allocatable :: wtall(:,:)
      complex *16, allocatable :: wtxall(:,:)
      complex *16, allocatable :: wtyall(:,:)
C
C     initialize coefficients to zero.
C
      allocate(wtall(4*72,nd))
      allocate(wtxall(4*72,nd))
      allocate(wtyall(4*72,nd))
      allocate(ww1p(12,nd))
      allocate(ww1m(12,nd))
      allocate(ww1pc(12,nd))
      allocate(ww1mc(12,nd))
      allocate(ww1px(12,nd))
      allocate(ww1mx(12,nd))
      dsq = 1.0D0/dsqrt(delta)
C
      do i=1,ns
         x = (sources(1,i) - center(1))*dsq
         y = (sources(2,i) - center(2))*dsq
         do j1=1,nn/2
            qqq = cdexp(ts(j1)*x)
ccc            qqq = 1.0d0
            do ind = 1,nd
               ccc = dipstr(ind,i)*ws(j1)
               ww1p(j1,ind) = ccc*qqq
               ww1px(j1,ind) = ts(j1)*dsq*ww1p(j1,ind)
               ww1m(j1,ind) = ccc/qqq
               ww1mx(j1,ind) = -ts(j1)*dsq*ww1m(j1,ind)
            enddo
            do ind = 1,nd
               ccc = charge(ind,i)*ws(j1)
               ww1pc(j1,ind) = ccc*qqq
               ww1mc(j1,ind) = ccc/qqq
            enddo
            qqq = cdexp(ts(j1)*y)
            ww2p(j1) = ws(j1)*qqq
            ww2py(j1) = ts(j1)*dsq*ww2p(j1)
            ww2p(nn/2+j1) = dconjg(ww2p(j1))
            ww2py(nn/2+j1) = dconjg(ww2py(j1))
            ww2m(j1) = ws(j1)/qqq
            ww2my(j1) = -ts(j1)*dsq*ww2m(j1)
            ww2m(nn/2+j1) = dconjg(ww2m(j1))
            ww2my(nn/2+j1) = dconjg(ww2my(j1))
         enddo
c
         nexp = nn*nn/2
         do ind = 1,nd
            j = 0
            do j1=1,nn/2
            do j2=1,nn
               j = j+1
               wtall(j,ind) = ww1pc(j1,ind)*ww2p(j2)
               wtall(nexp+j,ind) = ww1pc(j1,ind)*ww2m(j2)
               wtall(2*nexp+j,ind) = ww1mc(j1,ind)*ww2p(j2)
               wtall(3*nexp+j,ind) = ww1mc(j1,ind)*ww2m(j2)
c
               wtxall(j,ind) = ww1px(j1,ind)*ww2p(j2)
               wtxall(nexp+j,ind) = ww1px(j1,ind)*ww2m(j2)
               wtxall(2*nexp+j,ind) = ww1mx(j1,ind)*ww2p(j2)
               wtxall(3*nexp+j,ind) = ww1mx(j1,ind)*ww2m(j2)
               wtyall(j,ind) = ww1p(j1,ind)*ww2py(j2)
               wtyall(nexp+j,ind) = ww1p(j1,ind)*ww2my(j2)
               wtyall(2*nexp+j,ind) = ww1m(j1,ind)*ww2py(j2)
               wtyall(3*nexp+j,ind) = ww1m(j1,ind)*ww2my(j2)
            enddo
            enddo
         enddo
c
         do ind=1,nd
         do j=1,4*nn*nn/2
            soeall(j,ind) = soeall(j,ind) +
     1         wtxall(j,ind)*rnormal(1,i)+wtyall(j,ind)*rnormal(2,i)
            soeall(j,ind) = soeall(j,ind) + wtall(j,ind)
         enddo
         enddo
      enddo
      return
      end
C
C
C
C
c*********************************************************************
C
C eval SOE expansions (pot, pot+grad, pot+grad+hess)
C
C*********************************************************************
      subroutine g2dsevalp_vec(nd,delta,center,nn,ws,ts,soeall,
     1              targ,ntarg,pot)
C
C     This subroutine evaluates all SOE expansions about
C     CENT at location TARG.
C
C     INPUT
C
c     nd            = vector length (for vector input)
C     delta         = Gaussian variance
C     center        = center of the expansion
C     nn            = number of exponentials 
C     soeall        = all four (unrolled) SOE expansions are incremented
c                     in the order: pp,pm,mp,mm  
C     targ          = target location
C     ws,ts         = SOE weights and nodes
C
C     OUTPUT:
C
C     pot           = potential (or vectorized potentials) incremented
C
      implicit none
      integer nterms,nd,nn,i,ind,j1,j2,j,nexp,itarg,ntarg
      real *8 center(2),targ(2,ntarg)
      real *8 x,y,pot(nd,ntarg),delta,pmax,dsq
      complex *16, allocatable :: sumz(:)
      complex *16 soeall(4*nn*nn/2,nd)
      complex *16 ws(nn),ts(nn)
      complex *16 ww1p(12),ww2p(12),ww1m(12),ww2m(12)
      complex *16 wtall(4*72)
C
C     initialize coefficients to zero.
C
      allocate(sumz(nd))
      dsq = 1.0D0/dsqrt(delta)
C
      do itarg = 1,ntarg
         do ind = 1,nd
            sumz(ind) = 0.0d0
         enddo
c
         x = (targ(1,itarg) - center(1))*dsq
         y = (targ(2,itarg) - center(2))*dsq
         do j1=1,nn/2
            ww1p(j1) = cdexp(-ts(j1)*x)
            ww2p(j1) = cdexp(-ts(j1)*y)
            ww1m(j1) = 1/ww1p(j1)
            ww2m(j1) = 1/ww2p(j1)
            ww1p(nn/2+j1) = dconjg(ww1p(j1))
            ww2p(nn/2+j1) = dconjg(ww2p(j1))
            ww1m(nn/2+j1) = dconjg(ww1m(j1))
            ww2m(nn/2+j1) = dconjg(ww2m(j1))
         enddo
c
         j = 0
         nexp = nn*nn/2
         do j1=1,nn/2
         do j2=1,nn
            j = j+1
            wtall(j) = ww1p(j1)*ww2p(j2)
            wtall(nexp+j) = ww1p(j1)*ww2m(j2)
            wtall(2*nexp+j) = ww1m(j1)*ww2p(j2)
            wtall(3*nexp+j) = ww1m(j1)*ww2m(j2)
ccc           write(6,*) 'j1,j2,texp',j1,j2,texp(ind,j1,j2)
         enddo
         enddo
c
         do ind=1,nd
         do j=1,4*nn*nn/2
ccc         write(6,*) '3 in evaltp  nexp is ',nexp
ccc         write(6,*) 'talbotall ',talbotall(j,1)
ccc         write(6,*) 'wtall ',wtall(j)
               sumz(ind) = sumz(ind) +
     1         wtall(j)*soeall(j,ind)
         enddo
         enddo
c
         do ind = 1,nd
            pot(ind,itarg) = pot(ind,itarg)+dble(sumz(ind))/2
         enddo
      enddo
      return
      end
C
C
C
C
      subroutine g2dsevalg_vec(nd,delta,center,nn,ws,ts,soeall,
     1              targ,ntarg,pot,grad)
C
C     This subroutine evaluates all SOE expansions about
C     CENTER at location TARG.
C
C     INPUT
C
c     nd            = vector length (for vector input)
C     delta         = Gaussian variance
C     center        = center of the expansion
C     nn            = number of exponentials 
C     soeall     = all four (unrolled) SOE expansions are incremented
c                     in the order: pp,pm,mp,mm  
C     targ          = target location
C     ws,ts         = SOE weights and nodes
C
C     OUTPUT:
C
C     pot           = potential (or vectorized potentials) incremented
C     grad          = gradient (or vectorized gradients) incremented
C
      implicit none
      integer nterms,nd,nn,i,ind,j1,j2,j,nexp,itarg,ntarg
      real *8 center(2),targ(2,ntarg)
      real *8 x,y,pot(nd,ntarg),grad(nd,2,ntarg),delta,pmax,dsq
      complex *16, allocatable ::  sumz(:),sumzx(:),sumzy(:)
      complex *16 soeall(4*nn*nn/2,nd)
      complex *16 ws(nn),ts(nn)
      complex *16 ww1p(12),ww2p(12),ww1m(12),ww2m(12)
      complex *16 ww1px(12),ww2py(12),ww1mx(12),ww2my(12)
      complex *16 wtall(4*72)
      complex *16 wtxall(4*72)
      complex *16 wtyall(4*72)
C
C     initialize coefficients to zero.
C
      allocate(sumz(nd))
      allocate(sumzx(nd))
      allocate(sumzy(nd))
      dsq = 1.0D0/dsqrt(delta)
C
      do itarg = 1,ntarg
         do ind = 1,nd
            sumz(ind) = 0.0d0
            sumzx(ind) = 0.0d0
            sumzy(ind) = 0.0d0
         enddo
c
         x = (targ(1,itarg) - center(1))*dsq
         y = (targ(2,itarg) - center(2))*dsq
         do j1=1,nn/2
            ww1p(j1) = cdexp(-ts(j1)*x)
            ww2p(j1) = cdexp(-ts(j1)*y)
            ww1m(j1) = cdexp(ts(j1)*x)
            ww2m(j1) = cdexp(ts(j1)*y)
c
            ww1px(j1) = -ts(j1)*dsq*ww1p(j1)
            ww2py(j1) = -ts(j1)*dsq*ww2p(j1)
            ww1mx(j1) = ts(j1)*dsq*ww1m(j1) 
            ww2my(j1) = ts(j1)*dsq*ww2m(j1)
c
            ww1p(nn/2+j1) = dconjg(ww1p(j1))
            ww2p(nn/2+j1) = dconjg(ww2p(j1))
            ww1m(nn/2+j1) = dconjg(ww1m(j1))
            ww2m(nn/2+j1) = dconjg(ww2m(j1))
c
            ww1px(nn/2+j1) = dconjg(ww1px(j1))
            ww2py(nn/2+j1) = dconjg(ww2py(j1))
            ww1mx(nn/2+j1) = dconjg(ww1mx(j1))
            ww2my(nn/2+j1) = dconjg(ww2my(j1))
         enddo
c
         j = 0
         nexp = nn*nn/2
         do j1=1,nn/2
         do j2=1,nn
            j = j+1
            wtall(j) = ww1p(j1)*ww2p(j2)
            wtall(nexp+j) = ww1p(j1)*ww2m(j2)
            wtall(2*nexp+j) = ww1m(j1)*ww2p(j2)
            wtall(3*nexp+j) = ww1m(j1)*ww2m(j2)
c
            wtxall(j) = ww1px(j1)*ww2p(j2)
            wtxall(nexp+j) = ww1px(j1)*ww2m(j2)
            wtxall(2*nexp+j) = ww1mx(j1)*ww2p(j2)
            wtxall(3*nexp+j) = ww1mx(j1)*ww2m(j2)
c
            wtyall(j) = ww1p(j1)*ww2py(j2)
            wtyall(nexp+j) = ww1p(j1)*ww2my(j2)
            wtyall(2*nexp+j) = ww1m(j1)*ww2py(j2)
            wtyall(3*nexp+j) = ww1m(j1)*ww2my(j2)
ccc           write(6,*) 'j1,j2,texp',j1,j2,texp(ind,j1,j2)
         enddo
         enddo
c
         do ind=1,nd
         do j=1,4*nn*nn/2
            sumz(ind) = sumz(ind) +
     1      wtall(j)*soeall(j,ind)
            sumzx(ind) = sumzx(ind) +
     1      wtxall(j)*soeall(j,ind)
            sumzy(ind) = sumzy(ind) +
     1      wtyall(j)*soeall(j,ind)
         enddo
         enddo
c
         do ind = 1,nd
            pot(ind,itarg) = pot(ind,itarg)+dreal(sumz(ind))/2
            grad(ind,1,itarg) = grad(ind,1,itarg)+dreal(sumzx(ind))/2
            grad(ind,2,itarg) = grad(ind,2,itarg)+dreal(sumzy(ind))/2
         enddo
      enddo
      return
      end
C
C
C
C
C
      subroutine g2dsevalh_vec(nd,delta,center,nn,ws,ts,soeall,
     1              targ,ntarg,pot,grad,hess)
C
C     This subroutine evaluates all SOE expansions soeall about
C     CENTER at location TARG
C
C     INPUT
C
c     nd            = vector length (for vector input)
C     delta         = Gaussian variance
C     center        = center of the expansion
C     nn            = number of exponentials 
C     soeall     = all four (unrolled) SOE expansions are incremented
c                     in the order: pp,pm,mp,mm  
C     targ          = target location
C     ws,ts         = SOE weights and nodes
C
C     OUTPUT:
C
C     pot           = potential (or vectorized potentials) incremented
C     grad          = gradient (or vectorized gradients) incremented
C     hess          = Hessian (or vectorized Hessians) incremented
C
      implicit none
      integer nterms,nd,nn,i,ind,j1,j2,j,nexp,itarg,ntarg
      real *8 center(2),targ(2,ntarg)
      real *8 pot(nd,ntarg),grad(nd,2,ntarg),hess(nd,3,ntarg)
      real *8 x,y,delta,pmax,dsq
      complex *16, allocatable :: sumz(:),sumzx(:),sumzy(:)
      complex *16, allocatable ::  sumzxx(:),sumzxy(:),sumzyy(:)
      complex *16 soeall(4*nn*nn/2,nd)
      complex *16 ws(nn),ts(nn)
      complex *16 ww1p(12),ww2p(12),ww1m(12),ww2m(12)
      complex *16 ww1px(12),ww2py(12),ww1mx(12),ww2my(12)
      complex *16 ww1pxx(12),ww2pyy(12),ww1mxx(12),ww2myy(12)
      complex *16 wtall(4*72)
      complex *16 wtxall(4*72)
      complex *16 wtyall(4*72)
      complex *16 wtxxall(4*72)
      complex *16 wtxyall(4*72)
      complex *16 wtyyall(4*72)
C
C     initialize coefficients to zero.
C
      allocate(sumz(nd))
      allocate(sumzx(nd))
      allocate(sumzy(nd))
      allocate(sumzxx(nd))
      allocate(sumzxy(nd))
      allocate(sumzyy(nd))
      dsq = 1.0D0/dsqrt(delta)
C
      do itarg = 1,ntarg
         do ind = 1,nd
            sumz(ind) = 0.0d0
            sumzx(ind) = 0.0d0
            sumzy(ind) = 0.0d0
            sumzxx(ind) = 0.0d0
            sumzxy(ind) = 0.0d0
            sumzyy(ind) = 0.0d0
         enddo
c
         x = (targ(1,itarg) - center(1))*dsq
         y = (targ(2,itarg) - center(2))*dsq
         do j1=1,nn/2
c
            ww1p(j1) = cdexp(-ts(j1)*x)
            ww2p(j1) = cdexp(-ts(j1)*y)
            ww1m(j1) = cdexp(ts(j1)*x)
            ww2m(j1) = cdexp(ts(j1)*y)
c
            ww1px(j1) = -ts(j1)*dsq*ww1p(j1)
            ww2py(j1) = -ts(j1)*dsq*ww2p(j1)
            ww1mx(j1) = ts(j1)*dsq*ww1m(j1) 
            ww2my(j1) = ts(j1)*dsq*ww2m(j1)
c
            ww1p(nn/2+j1) = dconjg(ww1p(j1))
            ww2p(nn/2+j1) = dconjg(ww2p(j1))
            ww1m(nn/2+j1) = dconjg(ww1m(j1))
            ww2m(nn/2+j1) = dconjg(ww2m(j1))
c
            ww1px(nn/2+j1) = dconjg(ww1px(j1))
            ww2py(nn/2+j1) = dconjg(ww2py(j1))
            ww1mx(nn/2+j1) = dconjg(ww1mx(j1))
            ww2my(nn/2+j1) = dconjg(ww2my(j1))
c
            ww1pxx(j1) = -ts(j1)*dsq*ww1px(j1)
            ww2pyy(j1) = -ts(j1)*dsq*ww2py(j1)
            ww1mxx(j1) = ts(j1)*dsq*ww1mx(j1) 
            ww2myy(j1) = ts(j1)*dsq*ww2my(j1)
c
            ww1pxx(nn/2+j1) = dconjg(ww1pxx(j1))
            ww2pyy(nn/2+j1) = dconjg(ww2pyy(j1))
            ww1mxx(nn/2+j1) = dconjg(ww1mxx(j1))
            ww2myy(nn/2+j1) = dconjg(ww2myy(j1))
         enddo
c
         j = 0
         nexp = nn*nn/2
         do j1=1,nn/2
         do j2=1,nn
            j = j+1
            wtall(j) = ww1p(j1)*ww2p(j2)
            wtall(nexp+j) = ww1p(j1)*ww2m(j2)
            wtall(2*nexp+j) = ww1m(j1)*ww2p(j2)
            wtall(3*nexp+j) = ww1m(j1)*ww2m(j2)
c
            wtxall(j) = ww1px(j1)*ww2p(j2)
            wtxall(nexp+j) = ww1px(j1)*ww2m(j2)
            wtxall(2*nexp+j) = ww1mx(j1)*ww2p(j2)
            wtxall(3*nexp+j) = ww1mx(j1)*ww2m(j2)
c
            wtyall(j) = ww1p(j1)*ww2py(j2)
            wtyall(nexp+j) = ww1p(j1)*ww2my(j2)
            wtyall(2*nexp+j) = ww1m(j1)*ww2py(j2)
            wtyall(3*nexp+j) = ww1m(j1)*ww2my(j2)
c
            wtxxall(j) = ww1pxx(j1)*ww2p(j2)
            wtxxall(nexp+j) = ww1pxx(j1)*ww2m(j2)
            wtxxall(2*nexp+j) = ww1mxx(j1)*ww2p(j2)
            wtxxall(3*nexp+j) = ww1mxx(j1)*ww2m(j2)
c
            wtxyall(j) = ww1px(j1)*ww2py(j2)
            wtxyall(nexp+j) = ww1px(j1)*ww2my(j2)
            wtxyall(2*nexp+j) = ww1mx(j1)*ww2py(j2)
            wtxyall(3*nexp+j) = ww1mx(j1)*ww2my(j2)
c
            wtyyall(j) = ww1p(j1)*ww2pyy(j2)
            wtyyall(nexp+j) = ww1p(j1)*ww2myy(j2)
            wtyyall(2*nexp+j) = ww1m(j1)*ww2pyy(j2)
            wtyyall(3*nexp+j) = ww1m(j1)*ww2myy(j2)
         enddo
         enddo
c
         do ind=1,nd
         do j=1,4*nn*nn/2
            sumz(ind) = sumz(ind) +
     1      wtall(j)*soeall(j,ind)
            sumzx(ind) = sumzx(ind) +
     1      wtxall(j)*soeall(j,ind)
            sumzy(ind) = sumzy(ind) +
     1      wtyall(j)*soeall(j,ind)
            sumzxx(ind) = sumzxx(ind) +
     1      wtxxall(j)*soeall(j,ind)
            sumzxy(ind) = sumzxy(ind) +
     1      wtxyall(j)*soeall(j,ind)
            sumzyy(ind) = sumzyy(ind) +
     1      wtyyall(j)*soeall(j,ind)
         enddo
         enddo
c
         do ind = 1,nd
            pot(ind,itarg) = pot(ind,itarg)+dreal(sumz(ind))/2
            grad(ind,1,itarg) = grad(ind,1,itarg)+dreal(sumzx(ind))/2
            grad(ind,2,itarg) = grad(ind,2,itarg)+dreal(sumzy(ind))/2
            hess(ind,1,itarg) = hess(ind,1,itarg)+dreal(sumzxx(ind))/2
            hess(ind,2,itarg) = hess(ind,2,itarg)+dreal(sumzxy(ind))/2
            hess(ind,3,itarg) = hess(ind,3,itarg)+dreal(sumzyy(ind))/2
         enddo
      enddo
      return
      end
C
C
c
c
c
c*********************************************************************
C
C get SOE approximation nodes and weights
C
C*********************************************************************
      subroutine get_soe1(n,ws,ts)
c     SOE nodes and weights, n=4-12 uses Trefethen's code
c     n = 14, 18 from Shidong's code. n=16 not useful after testing.
c      
      implicit real *8 (a-h,o-z)
      complex *16 ts(*), ws(*)
      complex *16 ima
      data ima/(0.0d0,1.0d0)/
      
      dimension tr2(2),ti2(2),wr2(2),wi2(2)
      dimension tr3(3),ti3(3),wr3(3),wi3(3)
      dimension tr4(4),ti4(4),wr4(4),wi4(4)
      dimension tr5(5),ti5(5),wr5(5),wi5(5)
      dimension tr6(6),ti6(6),wr6(6),wi6(6)
      dimension tr7(7),ti7(7),wr7(7),wi7(7)
      
      dimension tr9(9),ti9(9),wr9(9),wi9(9)
c
c     soe approximation of the Gaussian with 4 exponentials
c     by symmetry, only 2 are needed
c     maximum absolute error = 3.568315082992379e-03 occured at x=0
c      
      data  tr2/
     1    1.323323306255279d+00,
     2    1.286223575135529d+00/
      data ti2/
     1    4.503154475908946d-01,
     2    1.422044092014818d+00/
      data wr2/
     1    1.704458483262654d+00,
     2   -7.080267983456461d-01/
      data wi2/
     1    4.523245320947962d+00,
     2   -4.574481278850558d-01/
c      
c     soe approximation of the Gaussian with 6 exponentials
c     by symmetry, only 3 are needed
c     maximum absolute error = 7.227866286507734e-05 occured at x=0
c      
      
      data tr3/
     1    1.593950878423775d+00,
     2    1.573920008119307d+00,
     3    1.527361492380405d+00/
      data ti3/
     1    3.742679039086306d-01,
     2    1.148334280704164d+00,
     3    2.028502256426840d+00/
      data wr3/
     1    3.342068172658133d+00,
     2   -2.531024297645422d+00,
     3    1.888838463244245d-01/
      data wi3/
     1    8.749200259833632d+00,
     2   -1.422348569452409d+00,
     3   -3.842837716207391d-03/
c      
c     soe approximation of the Gaussian with 8 exponentials
c     by symmetry, only 4 are needed
c     maximum absolute error = 1.256988759079647e-06 occured at x=0
c      
      data tr4/
     1    1.824272997628246d+00,
     2    1.811178322883317d+00,
     3    1.782903690503146d+00,
     4    1.732538970083639d+00/
      data ti4/
     1    3.271493923766529d-01,
     2    9.940411307226505d-01,
     3    1.705653711142857d+00,
     4    2.531843356999449d+00/
      data wr4/
     1    6.664582346450906d+00,
     2   -6.732794546758117d+00,
     3    1.097606962739752d+00,
     4   -2.939601942129933d-02/
      data wi4/
     1    1.731489069676853d+01,
     2   -3.579372008322773d+00,
     3   -1.674072276003963d-01,
     4    1.926961376278205d-02/
c      
c     soe approximation of the Gaussian with 10 exponentials
c     by symmetry, only 5 are needed
c     maximum absolute error = 2.004231713925719e-08 occured at x=0
c      
      data tr5/
     1    2.028383023345492d+00,
     2    2.018946454671973d+00,
     3    1.999156351285416d+00,
     4    1.966574319686306d+00,
     5    1.914430557121770d+00/
      data ti5/
     1    2.942876306371169d-01,
     2    8.901639685551171d-01,
     3    1.510370844469048d+00,
     4    2.182159307219590d+00,
     5    2.970898115468925d+00/
      data wr5/
     1    1.347575356075990d+01,
     2   -1.617896221170955d+01,
     3    3.933599507145248d+00,
     4   -2.329555498501443d-01,
     5    2.564673612230538d-03/
      data wi5/
     1    3.484098246693387d+01,
     2   -8.352129542827246d+00,
     3   -8.149531822098305d-01,
     4    2.328754348729245d-01,
     5   -5.249948161805266d-03/      
c      
c     soe approximation of the Gaussian with 12 exponentials
c     by symmetry, only 6 are needed
c     maximum absolute error = 3.022080363734858e-10
c      
      data tr6/
     1    2.213645327254045d+00,
     2    2.206422844728666d+00,
     3    2.191500702729535d+00,
     4    2.167727555011547d+00,
     5    2.132617049485016d+00,
     6    2.079609194319604d+00/
      data ti6/
     1    2.696881875222051d-01,
     2    8.137426522157325d-01,
     3    1.372882493475410d+00,
     4    1.961462532263335d+00,
     5    2.604615798736701d+00,
     6    3.365035256316001d+00/
      data wr6/
     1    2.754140321166421d+01,
     2   -3.712394227425998d+01,
     3    1.160506323094613d+01,
     4   -1.037024624392356d+00,
     5    1.444285777933479d-02,
     6    5.759796045271480d-05/
      data wi6/
     1    7.096775937178015d+01,
     2   -1.882232395234807d+01,
     3   -2.747041853363806d+00,
     4    1.238467238588732d+00,
     5   -8.150092347951716d-02,
     6    8.977773646777127d-04/
c      
c     soe approximation of the Gaussian with 14 exponentials
c     by symmetry, only 7 are needed
c     maximum absolute error = 4.396816e-12 occured at x=1.19e-5
c      
      data tr7/
     1    2.384508883129899d0,
     2    2.378748363769416d0,
     3    2.366950815439966d0,
     4    2.348483552895229d0,
     5    2.322129865359463d0,
     6    2.285431589285171d0,
     7    2.232028511136604d0/
      data ti7/
     1    2.503812606417072d-1,
     2    7.543536816315837d-1,
     3    1.268478185247206d0,
     4    1.801538303670637d0,
     5    2.366641441580562d0,
     6    2.987710524123123d0,
     7    3.725542056775757d0/
      data wr7/
     1    5.676243263824121d+01,
     2   -8.315430631762199d+01,
     3    3.096383810675961d+01,
     4   -3.605394150263557d+00,
     5    2.604202645105868d-02,
     6    7.458238052201011d-03,
     7   -7.054162290338943d-05/
      data wi7/
     1    1.459042184345505d+02,
     2   -4.165421433230362d+01,
     3   -7.893835305910772d+00,
     4    4.759586746083758d+00,
     5   -5.382322766556747d-01,
     6    1.690602909533891d-02,
     7   -1.103474906299145d-04/      
c      
c     soe approximation of the Gaussian with 18 exponentials
c     by symmetry, only 9 are needed
c     maximum absolute error = 4.35207e-14 occured at x=6.13e-02
c      

      data tr9/
     1    1.981663215750117d+00,
     2    1.981647628200388d+00,
     3    1.981689205244672d+00,
     4    1.982013944671143d+00,
     5    1.982404603659796d+00,
     6    1.980775824606087d+00,
     7    1.973465871627027d+00,
     8    1.956362674897580d+00,
     9    1.922784087489885d+00/
      data ti9/
     1    1.981670504723005d-01,
     2    5.945399651131158d-01,
     3    9.914045722148965d-01,
     4    1.391493300900841d+00,
     5    1.803543239304497d+00,
     6    2.242874054628137d+00,
     7    2.727161447761292d+00,
     8    3.279508992292855d+00,
     9    3.956014286362498d+00/
      data wr9/
     1    1.543145316345868d+01,
     2   -1.127561215561926d+01,
     3   -6.015030818595092d+00,
     4    2.402423010030009d+00,
     5    6.071597114970148d-01,
     6   -1.475371450228843d-01,
     7   -3.483098287003385d-03,
     8    6.325201288208040d-04,
     9   -5.187590291256931d-06/
      data wi9/
     1    1.543237039314943d+01,
     2    1.127174659689771d+01,
     3   -6.038799787955659d+00,
     4   -2.322703836849708d+00,
     5    7.163573316855557d-01,
     6    8.668250059198450d-02,
     7   -1.639557459051920d-02,
     8    1.991264418638442d-04,
     9    2.414377231363515d-06/
      
      n2=n/2
      if (n .lt. 4 .or. n .gt. 18) then
         write(*,*) 'n must be an integer in [4,14] or 18'
         stop
      endif
      
      if (n .le. 4) then
         do i=1,n2
            ts(i)=tr2(i)+ti2(i)*ima
            ws(i)=wr2(i)+wi2(i)*ima            
         enddo
      elseif (n .eq. 6) then
         do i=1,n2
            ts(i)=tr3(i)+ti3(i)*ima
            ws(i)=wr3(i)+wi3(i)*ima            
         enddo
      elseif (n .eq. 8) then
         do i=1,n2
            ts(i)=tr4(i)+ti4(i)*ima
            ws(i)=wr4(i)+wi4(i)*ima            
         enddo
      elseif (n .eq. 10) then
         do i=1,n2
            ts(i)=tr5(i)+ti5(i)*ima
            ws(i)=wr5(i)+wi5(i)*ima            
         enddo
      elseif (n .eq. 12) then
         do i=1,n2
            ts(i)=tr6(i)+ti6(i)*ima
            ws(i)=wr6(i)+wi6(i)*ima            
         enddo
      elseif (n .eq. 14) then
         do i=1,n2
            ts(i)=tr7(i)+ti7(i)*ima
            ws(i)=wr7(i)+wi7(i)*ima            
         enddo
      elseif (n .eq. 18) then
         do i=1,n2
            ts(i)=tr9(i)+ti9(i)*ima
            ws(i)=wr9(i)+wi9(i)*ima            
         enddo
      endif


      do i=1,n2
         ts(i+n2)=conjg(ts(i))
         ws(i+n2)=conjg(ws(i))
      enddo
      
      do i=1,n
         ts(i)=ts(i)*2
         ws(i)=ws(i)
      enddo

      return
      end subroutine
c      
c      
c
c      
      subroutine get_soe2(n,ws,ts)
c     SOE nodes and weights, n=4-12 uses least squares to improve the accuracy,
c     starting with the initial guesses from Trefethen's code
c
c     returns to the user the SOE quadrature nodes and weights for the 
c     1D Gaussian, i.e.,
c
c     exp(-x^2) = \sum_{i=1}1^n ws(i) exp(-ts(i) |x|)  
c     
c     Note that the nodes and weights appear in conjugate pairs
c     
c
c     input parameters:
c     
c     n - the total number of quadrature nodes and weights
c         n must be an even integer in [4,14]
c      
c     output parameters:
c
c     ts - quadrature nodes 
c     ws - quadrature weights 
c      
c     errmax - maximum absolute error of the SOE approximation
c              n          errmax
c              2          3.516715789022018e-04
c              3          6.510784017060844e-06
c              4          6.853144535190481e-08
c              5          2.005336785515510e-09
c              6          3.729017095110976e-11
c              7          6.112887973586112e-13
c      
c              9          4.35207e-14 
c     
      implicit real *8 (a-h,o-z)
      complex *16 ts(*), ws(*)
      complex *16 ima
      data ima/(0.0d0,1.0d0)/
      
      dimension tr2(2),ti2(2),wr2(2),wi2(2)
      dimension tr3(3),ti3(3),wr3(3),wi3(3)
      dimension tr4(4),ti4(4),wr4(4),wi4(4)
      dimension tr5(5),ti5(5),wr5(5),wi5(5)
      dimension tr6(6),ti6(6),wr6(6),wi6(6)
      dimension tr7(7),ti7(7),wr7(7),wi7(7)

      dimension tr9(9),ti9(9),wr9(9),wi9(9)
c
c     soe approximation of the Gaussian with 4 exponentials
c     by symmetry, only 2 are needed
c     maximum absolute error = 3.516715789022018e-04.
c
      data  tr2/
     1    -2.661798367170939d+00,
     2    -2.532564794761101d+00/
      data ti2/
     1    -8.529355096064692d-01,
     2    -2.699519851786320d+00/
      data wr2/
     1    9.372401876584128d-01,
     2   -4.372616823494260d-01/
      data wi2/
     1    2.192847167149230d+00,
     2   -1.766470350798564d-01/
c      
c     soe approximation of the Gaussian with 6 exponentials
c     by symmetry, only 3 are needed
c     maximum absolute error = 6.510784017060844e-06
c      
      data tr3/
     1    -3.182995116354765d+00,
     2    -3.260753517129289d+00,
     3    -3.217466389941865d+00/
      data ti3/
     1    -7.850177885864261d-01,
     2    -2.344606600656128d+00,
     3    -3.902467582782956d+00/
      data wr3/
     1    1.621465841345878d+00,
     2   -1.267954650617345d+00,
     3    1.464884524546549d-01/
      data wi3/
     1    4.502292586423311d+00,
     2   -9.019856645033295d-01,
     3    2.012930239441965d-02/
c      
c     soe approximation of the Gaussian with 8 exponentials
c     by symmetry, only 4 are needed
c     maximum absolute error = 6.853144535190481e-08
c
      data tr4/
     1    -3.647246911001136d+00,
     2    -3.626225133829404d+00,
     3    -3.603022265862902d+00,
     4    -3.535187295699224d+00/
      data ti4/
     1    -6.529696220312070d-01,
     2    -1.990418453744666d+00,
     3    -3.408768895825578d+00,
     4    -4.990180992026964d+00/
      data wr4/
     1    3.343792210857979d+00,
     2   -3.403876847305482d+00,
     3    5.796032765882220d-01,
     4   -1.951864411560219d-02/
      data wi4/
     1    8.625862610542274d+00,
     2   -1.793605606794936d+00,
     3   -7.242128225115123d-02,
     4    1.126277582434443d-02/
c      
c     soe approximation of the Gaussian with 10 exponentials
c     by symmetry, only 5 are needed
c     maximum absolute error = 2.005336785515510e-09
c
      data tr5/
     1    -4.056700573705696d+00,
     2    -4.037965286896315d+00,
     3    -3.999001176267095d+00,
     4    -3.934063076381675d+00,
     5    -3.821336891615266d+00/
      data ti5/
     1    -5.885148526378086d-01,
     2    -1.780339417934965d+00,
     3    -3.020848791356881d+00,
     4    -4.363368237559863d+00,
     5    -5.944299221990216d+00/
      data wr5/
     1    6.739042108197683d+00,
     2   -8.092480046197947d+00,
     3    1.969157128297097d+00,
     4   -1.169602276245181d-01,
     5    1.241037171898254d-03/
      data wi5/
     1    1.741660958138088d+01,
     2   -4.175484069685802d+00,
     3   -4.067804058685952d-01,
     4    1.165550704739487d-01,
     5   -2.608776038901437d-03/
c      
c     soe approximation of the Gaussian with 12 exponentials
c     by symmetry, only 6 are needed
c     maximum absolute error = 3.729017095110976e-11
c
      data tr6/
     1    -4.427289426763536d+00,
     2    -4.412852687801833d+00,
     3    -4.382977192701978d+00,
     4    -4.335453870889060d+00,
     5    -4.264890006622823d+00,
     6    -4.159216567114453d+00/
      data ti6/
     1    -5.393740774982838d-01,
     2    -1.627485818541268d+00,
     3    -2.745764850165638d+00,
     4    -3.922957815027205d+00,
     5    -5.209121159851328d+00,
     6    -6.730081860478507d+00/
      data wr6/
     1     1.377070839769257d+01,
     2    -1.856201081252789d+01,
     3     5.802530575585394d+00,
     4    -5.184517677117910d-01,
     5     7.194472010677156d-03,
     6     2.913494780561260d-05/
      data wi6/
     1     3.548389585570825d+01,
     2    -9.411123583149262d+00,
     3    -1.373461789441497d+00,
     4     6.191099357882870d-01,
     5    -4.074506010315759d-02,
     6     4.481549603441414d-04/      
c      
c     soe approximation of the Gaussian with 14 exponentials
c     by symmetry, only 7 are needed
c     maximum absolute error = 6.112887973586112e-13
c
      data tr7/
     1    -4.769017777253634d+00,
     2    -4.757496924923473d+00,
     3    -4.733900833344295d+00,
     4    -4.696971190047083d+00,
     5    -4.644248847774636d+00,
     6    -4.570887586163360d+00,
     7    -4.464062381041680d+00/
      data ti7/
     1    -5.007625408668818d-01,
     2    -1.508707160675464d+00,
     3    -2.536957043341105d+00,
     4    -3.603074941435481d+00,
     5    -4.733268195340075d+00,
     6    -5.975411751927892d+00,
     7    -7.451081423130120d+00/
      data wr7/
     1     2.838122801479766d+01,
     2    -4.157715756962441d+01,
     3     1.548192113589467d+01,
     4    -1.802696453672285d+00,
     5     1.301085278245687d-02,
     6     3.729309561777345d-03,
     7    -3.528973992125576d-05/
      data wi7/
     1     7.295210990312378d+01,
     2    -2.082710997019782d+01,
     3    -3.946915218407423d+00,
     4     2.379791870406019d+00,
     5    -2.691218286283854d-01,
     6     8.453146344960496d-03,
     7    -5.518962089533348d-05/      
c      
c     soe approximation of the Gaussian with 18 exponentials
c     by symmetry, only 9 are needed
c     maximum absolute error = 4.35207e-14 occured at x=6.13e-02
c      

      data tr9/
     1    1.981663215750117d+00,
     2    1.981647628200388d+00,
     3    1.981689205244672d+00,
     4    1.982013944671143d+00,
     5    1.982404603659796d+00,
     6    1.980775824606087d+00,
     7    1.973465871627027d+00,
     8    1.956362674897580d+00,
     9    1.922784087489885d+00/
      data ti9/
     1    1.981670504723005d-01,
     2    5.945399651131158d-01,
     3    9.914045722148965d-01,
     4    1.391493300900841d+00,
     5    1.803543239304497d+00,
     6    2.242874054628137d+00,
     7    2.727161447761292d+00,
     8    3.279508992292855d+00,
     9    3.956014286362498d+00/
      data wr9/
     1    1.543145316345868d+01,
     2   -1.127561215561926d+01,
     3   -6.015030818595092d+00,
     4    2.402423010030009d+00,
     5    6.071597114970148d-01,
     6   -1.475371450228843d-01,
     7   -3.483098287003385d-03,
     8    6.325201288208040d-04,
     9   -5.187590291256931d-06/
      data wi9/
     1    1.543237039314943d+01,
     2    1.127174659689771d+01,
     3   -6.038799787955659d+00,
     4   -2.322703836849708d+00,
     5    7.163573316855557d-01,
     6    8.668250059198450d-02,
     7   -1.639557459051920d-02,
     8    1.991264418638442d-04,
     9    2.414377231363515d-06/
      
      n2=n/2
      if (n .lt. 4 .or. n .gt. 18) then
         write(*,*) 'n must be an integer in [4,14] or 18'
         stop
      endif
         
      if (n .eq. 4) then
         do i=1,n2
            ts(i)=tr2(i)+ti2(i)*ima
            ws(i)=wr2(i)+wi2(i)*ima            
         enddo
         errmax = 3.516715789022018d-04
      elseif (n .eq. 6) then
         do i=1,n2
            ts(i)=tr3(i)+ti3(i)*ima
            ws(i)=wr3(i)+wi3(i)*ima            
         enddo
         errmax = 6.510784017060844d-06
      elseif (n .eq. 8) then
         do i=1,n2
            ts(i)=tr4(i)+ti4(i)*ima
            ws(i)=wr4(i)+wi4(i)*ima            
         enddo
         errmax = 6.853144535190481d-08
      elseif (n .eq. 10) then
         do i=1,n2
            ts(i)=tr5(i)+ti5(i)*ima
            ws(i)=wr5(i)+wi5(i)*ima            
         enddo
         errmax = 2.005336785515510d-09
      elseif (n .eq. 12) then
         do i=1,n2
            ts(i)=tr6(i)+ti6(i)*ima
            ws(i)=wr6(i)+wi6(i)*ima            
         enddo
         errmax = 3.729017095110976d-11
      elseif (n .eq. 14) then
         do i=1,n2
            ts(i)=tr7(i)+ti7(i)*ima
            ws(i)=wr7(i)+wi7(i)*ima            
         enddo
         errmax = 6.112887973586112d-13
      elseif (n .eq. 18) then
         do i=1,n2
            ts(i)=tr9(i)+ti9(i)*ima
            ws(i)=wr9(i)+wi9(i)*ima            
         enddo
      endif

      do i=1,n2
         ts(i+n2)=conjg(ts(i))
         ws(i+n2)=conjg(ws(i))
      enddo

      if (n.le.14) then
         do i=1,n
            ts(i)=-ts(i)
            ws(i)=ws(i)*2
         enddo
      elseif (n.eq.18) then
         do i=1,n
            ts(i)=ts(i)*2
            ws(i)=ws(i)
         enddo
      endif
      
      return
      end subroutine
c      
c      
c
c      
c*********************************************************************
C
C shift SOE expansions (mp to loc, mp to mp, loc to loc)
C
C*********************************************************************
      subroutine shiftsoe_translation_matrices(boxsize,nsoe,nmax,
     1           wshift,ws,ts)
C
C     This subroutine computes the complex exponentials in the SOE translations.
C
C     INPUT
C
c     boxsize  = scaled (by 1/sqrt(delta) size of the box at the finest level
C     nsoe     = number of terms in SOE expansion
C     ntr      = number of different translation lengths in the whole scheme 
C     ws,ts    = SOE weights and nodes
C
C     OUTPUT:
C
C     wshift   = table of translation matrices for soe shift used in mp to loc stage
C
      implicit real *8 (a-h,o-z)
      real *8 boxsize
      complex *16, allocatable :: wexp(:,:)
      complex *16 wshift(nsoe*nsoe/2,nmax*nmax)
      complex *16 ws(nsoe),ts(nsoe)
C
      allocate(wexp(nsoe,nmax))
      
      do j2=1,nmax
         do j1=1,nsoe
            wexp(j1,j2) = cdexp(-ts(j1)*j2*boxsize)
         enddo
      enddo

      k=0
      do k1=1,nmax
         do k2=1,nmax
            k=k+1
            j=0
            do j1=1,nsoe/2
               do j2=1,nsoe
                  j=j+1
                  wshift(j,k)=wexp(j1,k1)*wexp(j2,k2)
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
      subroutine g2dshiftsoe_vec(nd,nn,soeall1,
     1              soeall2,wshift,ip)
C
C     This subroutine converts the SOE expansion (soeall1) about
C     the center (CENT1) into an SOE  expansion (soeall2) about 
C     (CENT2) for pp direction.
C
C     INPUT
C
c     nd      = vector length (for vector input)
C     nn      = number of terms in SOE expansion
C     soeall1 = original expansion 
C     wshift  = precomputed shift array
C     ip      = SOE expansion type index
c               0 - pp, 1 - pm, 2 - mp, 3 - mm
c
c
c
c
C     OUTPUT:
C
C     soeall2 = shifted expansion 
C
      implicit none
      integer nd,nn,j,j1,j2,ind,nexp,ip
      complex *16 soeall1(4*nn*nn/2,nd)
      complex *16 soeall2(4*nn*nn/2,nd)
      complex *16 wshift(nn*nn/2)
C
      nexp = nn*nn/2
      do ind=1,nd
         do j=1,nexp
            soeall2(j+ip*nexp,ind) = soeall2(j+ip*nexp,ind)
     1          +soeall1(j+ip*nexp,ind)*wshift(j)
         enddo
      enddo
c
      return
      end
c
C
c
C
      subroutine shiftsoeall_translation_matrices(xmin,nsoe,nmax,
     1           wshift,ws,ts)
C
C     This subroutine precomputes all translation matrices for all
c     SOE translations from child to parent or vice versa.
C
C     INPUT
C
c     xmin     = half of the scaled (by 1/sqrt(delta) size of the box 
c                at the finest level
C     nsoe     = number of terms in SOE expansion
C     nmax     = number of different translation lengths in the whole scheme 
C     ws,ts    = SOE weights and nodes
C
C     OUTPUT:
C
C     wshift   = table of translation matrices for soe shift used in
c                merge mp and split loc stage
C
      implicit real *8 (a-h,o-z)
      real *8 xmin
      complex *16 wshift(4*nsoe*nsoe/2,-nmax:nmax,-nmax:nmax)
      complex *16 ws(nsoe),ts(nsoe)

      complex *16, allocatable:: wwp(:,:),wwm(:,:)
C
      allocate(wwp(nsoe,-nmax:nmax),wwm(nsoe,-nmax:nmax))

      do j2=-nmax,nmax
         do j1=1,nsoe
            wwp(j1,j2)=exp(-ts(j1)*xmin*j2)
            wwm(j1,j2)=exp( ts(j1)*xmin*j2)
         enddo
      enddo
c
      nexp=nsoe*nsoe/2

      do m=-nmax,nmax
      do n=-nmax,nmax
         j = 0
         do j1=1,nsoe/2
         do j2=1,nsoe
            j = j+1
            wshift(j,       m,n) = wwp(j1,m)*wwp(j2,n)
            wshift(j+nexp  ,m,n) = wwp(j1,m)*wwm(j2,n)
            wshift(j+2*nexp,m,n) = wwm(j1,m)*wwp(j2,n)
            wshift(j+3*nexp,m,n) = wwm(j1,m)*wwm(j2,n)
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
      subroutine g2dshiftsoeall_vec(nd,nn,soeall1,
     1              soeall2,wshift)
C
C     This subroutine shifts ALL directional SOE expansions (soeall1) 
C     from the center (CENT1) to SOE expansions (soeall2) about 
C     (CENT2) using the precomputed shift matrix wshift
C
C     INPUT
C
c     nd      = vector length (for vector input)
C     nn      = number of terms in SOE expansion
C     soeall1 = original expansion 
C     wshift  = precomputed translation matrix
c      
C     OUTPUT:
C
C     soeall2 = shifted expansion 
C
      implicit none
      integer nd,nn,j,ind,nexp
      complex *16 wshift(4*nn*nn/2)
      complex *16 soeall1(4*nn*nn/2,nd)
      complex *16 soeall2(4*nn*nn/2,nd)
C
      nexp = 4*nn*nn/2
      do ind=1,nd
      do j=1,nexp
         soeall2(j,ind) = soeall2(j,ind) + soeall1(j,ind)*wshift(j)
      enddo
      enddo
      
      return
      end
c
c
C      
c
C***********************************************************************
      subroutine g2dsoezero_vec(nd,soeall,nsoe)
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
      integer n,nsoe,nd,ii
      complex *16 soeall(4*nsoe*nsoe/2,nd)
c
      do ii=1,nd
      do n=1,4*nsoe*nsoe/2
         soeall(n,ii)=0.0d0
      enddo
      enddo
      return
      end
c
c
C
