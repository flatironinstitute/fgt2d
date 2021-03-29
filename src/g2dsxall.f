c     This file contains the following subroutines
c
c     g2dformsxc_vec : computes the SX expansion about the center CENT
c                     due to sources of strength charge()
c
c     g2dformsxd_vec : computes the SX expansion about the center CENT
c                     due to dipole
c
c     g2dformsxcd_vec : computes the SX expansion about the center CENT
c                     due to charge and dipole
c
C     g2dsxevalp_vec : evaluates the SX expansion 
C                         potential only.
c      
C     g2dsxevalg_vec : evaluates the SX expansion 
C                         potential + gradient
c      
C     g2dsxevalh_vec : evaluates the SX expansion 
C                         potential + gradient + hessian
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
C form SOE/X expansions (charge, dipole, charge+dipole)
C
C*********************************************************************
      subroutine g2dformsxc_vec(nd,delta,sources,ns,charge,center,
     1            npw,wx,tx,nns,ws,ts,spall)
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
C     npw           = number of Fourier plane waves
C     wx,tx         = planewave weights and nodes
C     nns           = number of SOE terms (half of
C                     full quadrature - no complex conjugates)
C     ws,ts         = SOE weights and nodes
C
C     OUTPUT:
C
C     spall    = hybrid expansions in order px,mx,xp,xm
C
      implicit none
      integer npw,nd,ns,nns,i,ii,ind,j1,j2,j,nexp
      real *8 center(2),sources(2,ns),charge(nd,ns)
      real *8 wx(-npw:npw),tx(-npw:npw)
      complex *16 spall(4*(2*npw+1)*nns,nd)
      real *8 x,xp,xq,delta,dsq,y,pi,w
      complex *16 qqq,eye
      complex *16 qqx,qqy,qq1,qq2
      complex *16 ws(nns),ts(nns)
      complex *16 ww1p(100),ww2p(100),ww1m(100),ww2m(100)
      complex *16, allocatable :: wtall(:,:)
      complex *16, allocatable :: ww1c(:,:),ww2c(:,:)
      complex *16 ww1(-100:100)
      complex *16 ww2(-100:100)
C
C     initialize coefficients to zero.
C
      allocate(wtall(4*(2*npw+1)*nns,nd))
c      
      eye = dcmplx(0,1)
c
      dsq = 1.0D0/dsqrt(delta)
      allocate(ww1c(-npw:npw,nd))
      allocate(ww2c(-npw:npw,nd))
c
      nexp = (2*npw+1)*nns
      
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
         qqx = cdexp(-eye*tx(1)*x)
         qqy = cdexp(-eye*tx(1)*y)
         qq1 = 1.0d0
         qq2 = 1.0d0
         ww1(0) = wx(0)
         ww2(0) = wx(0)
         do j1=1,npw
            qq1 = qq1*qqx
            qq2 = qq2*qqy
            ww1(j1) = wx(j1)*qq1
            ww2(j1) = wx(j1)*qq2
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
      subroutine g2dformsxd_vec(nd,delta,sources,ns,rnormal,dipstr,
     1            center,npw,wx,tx,nns,ws,ts,spall)
C
C     This subroutine computes the four hybrid 
C     SOE/plane wave expansions about CENTER
C     due to dipole sources of strength dipstr().
C
C     INPUT
C
c     nd            = vector length (for multiple charges at same locations)
C     delta         = Gaussian variance
C     sources(2,ns) = coordinates of sources
C     ns            = number of sources
C     rnormal       = normals
C     dipstr        = dipole strengths
C     center        = center of the expansion
C     npw           = number of Fourier plane waves
C     wx,tx         = planewave weights and nodes
C     nns           = number of SOE terms (half of
C                     full quadrature - no complex conjugates)
C     ws,ts         = SOE weights and nodes
C
C     OUTPUT:
C
C     spall    = hybrid expansions in order px,mx,xp,xm
C
      implicit none
      integer npw,nd,ns,nns,i,ii,ind,j1,j2,j,nexp
      real *8 center(2),sources(2,ns),dipstr(nd,ns)
      real *8 rnormal(2,ns)
      real *8 wx(-npw:npw),tx(-npw:npw)
      complex *16 spall(4*(2*npw+1)*nns,nd)
      real *8 x,xp,xq,delta,dsq,y,pi,d1,d2
      complex *16 qqq,eye
      complex *16 z0,z1,z2,z3
      complex *16 c0,c1,c2,c3
      complex *16 qqx,qqy,qq1,qq2
      complex *16 ws(nns),ts(nns)
      complex *16 ww1p(100),ww2p(100),ww1m(100),ww2m(100)
      complex *16 ww1px(100),ww2py(100),ww1mx(100),ww2my(100)
      complex *16 ww1(-100:100),ww1x(-100:100)
      complex *16 ww2(-100:100),ww2y(-100:100)
C
      eye = dcmplx(0,1)
c
      dsq = 1.0D0/dsqrt(delta)
c
      do i=1,ns
         x = (sources(1,i) - center(1))*dsq
         y = (sources(2,i) - center(2))*dsq
         do j1=1,nns
            qqq = cdexp(ts(j1)*x)
            ww1p(j1) = ws(j1)*qqq
            ww1px(j1)=ts(j1)*dsq*ww1p(j1)
            ww1m(j1) = ws(j1)/qqq
            ww1mx(j1)=-ts(j1)*dsq*ww1m(j1)
            
            qqq = cdexp(ts(j1)*y)
            ww2p(j1) = ws(j1)*qqq
            ww2py(j1)=ts(j1)*dsq*ww2p(j1)
            ww2m(j1) = ws(j1)/qqq
            ww2my(j1)=-ts(j1)*dsq*ww2m(j1)
         enddo
c
         qqx = cdexp(-eye*tx(1)*x)
         qqy = cdexp(-eye*tx(1)*y)
         qq1 = 1.0d0
         qq2 = 1.0d0
         ww1(0) = wx(0)
         ww2(0) = wx(0)
         ww1x(0)=0
         ww2y(0)=0
         do j1=1,npw
            qq1 = qq1*qqx
            qq2 = qq2*qqy
            ww1(j1) = wx(j1)*qq1
            ww1x(j1)=-eye*tx(j1)*dsq*ww1(j1)
            ww2(j1) = wx(j1)*qq2
            ww2y(j1)=-eye*tx(j1)*dsq*ww2(j1)
            ww1(-j1) = dconjg(ww1(j1))
            ww2(-j1) = dconjg(ww2(j1))
            ww1x(-j1)=dconjg(ww1x(j1))
            ww2y(-j1)=dconjg(ww2y(j1))
         enddo
c
         nexp = (2*npw+1)*nns
         do ind=1,nd
            d1=dipstr(ind,i)*rnormal(1,i)
            d2=dipstr(ind,i)*rnormal(2,i)
            j = 0
            do j1=1,nns
               z0=d1*ww1px(j1)
               z1=d1*ww1mx(j1)
               z2=d1*ww2p(j1)
               z3=d1*ww2m(j1)

               c0=d2*ww1p(j1)
               c1=d2*ww1m(j1)
               c2=d2*ww2py(j1)
               c3=d2*ww2my(j1)
               do j2=-npw,npw
                  j = j+1
                  spall(j,ind) = spall(j,ind) + 
     1                z0*ww2(j2) + c0*ww2y(j2)
                  spall(nexp+j,ind) = spall(nexp+j,ind) +
     1                z1*ww2(j2) + c1*ww2y(j2)
                  spall(2*nexp+j,ind) = spall(2*nexp+j,ind) +
     1                z2*ww1x(j2) + c2*ww1(j2)
                  spall(3*nexp+j,ind) = spall(3*nexp+j,ind) +
     1                z3*ww1x(j2) + c3*ww1(j2)
               enddo
            enddo
         enddo
      enddo

      return
      end
C
C
C
      subroutine g2dformsxcd_vec(nd,delta,sources,ns,charge,rnormal,
     1            dipstr,center,npw,wx,tx,nns,ws,ts,spall)
C
C     This subroutine computes the four hybrid 
C     SOE/plane wave expansions about CENTER
C     due to sources of strength charge()
C     and dipole sources of strength dipstr().
c      
C     INPUT
C
c     nd            = vector length (for multiple charges at same locations)
C     delta         = Gaussian variance
C     sources(2,ns) = coordinates of sources
C     ns            = number of sources
C     charge        = strength of sources
C     rnormal       = normals
C     dipstr        = dipole strengths
C     center        = center of the expansion
C     npw           = number of Fourier plane waves
C     wx,tx         = planewave weights and nodes
C     nns           = number of SOE terms (half of
C                     full quadrature - no complex conjugates)
C     ws,ts         = SOE weights and nodes
C
C     OUTPUT:
C
C     spall         = four SX hybrid expansions are incremented
c                      in the order px,mx,xp,xm
C
      implicit none
      integer npw,nd,ns,nns,i,ii,ind,j1,j2,j,nexp
      real *8 center(2),sources(2,ns),charge(nd,ns),dipstr(nd,ns)
      real *8 rnormal(2,ns)
      real *8 wx(-npw:npw),tx(-npw:npw)
      complex *16 spall(4*(2*npw+1)*nns,nd)
      real *8 x,xp,xq,delta,dsq,y,pi,w,d1,d2
      complex *16 qqq,eye
      complex *16 qqx,qqy,qq1,qq2
      complex *16 z0,z1,z2,z3
      complex *16 c0,c1,c2,c3
      complex *16 ws(nns),ts(nns)
      complex *16 ww1p(100),ww2p(100),ww1m(100),ww2m(100)
      complex *16 ww1px(100),ww2py(100),ww1mx(100),ww2my(100)
      complex *16 ww1(-100:100),ww1x(-100:100)
      complex *16 ww2(-100:100),ww2y(-100:100)
c      
      eye = dcmplx(0,1)
c
      dsq = 1.0D0/dsqrt(delta)
c
      do i=1,ns
         x = (sources(1,i) - center(1))*dsq
         y = (sources(2,i) - center(2))*dsq
         do j1=1,nns
            qqq = cdexp(ts(j1)*x)
            ww1p(j1) = ws(j1)*qqq
            ww1px(j1)=ts(j1)*dsq*ww1p(j1)
            ww1m(j1) = ws(j1)/qqq
            ww1mx(j1)=-ts(j1)*dsq*ww1m(j1)
            
            qqq = cdexp(ts(j1)*y)
            ww2p(j1) = ws(j1)*qqq
            ww2py(j1)=ts(j1)*dsq*ww2p(j1)
            ww2m(j1) = ws(j1)/qqq
            ww2my(j1)=-ts(j1)*dsq*ww2m(j1)
         enddo
c
         qqx = cdexp(-eye*tx(1)*x)
         qqy = cdexp(-eye*tx(1)*y)
         qq1 = 1.0d0
         qq2 = 1.0d0
         ww1(0) = wx(0)
         ww2(0) = wx(0)
         ww1x(0)=0
         ww2y(0)=0
         do j1=1,npw
            qq1 = qq1*qqx
            qq2 = qq2*qqy
            ww1(j1) = wx(j1)*qq1
            ww1x(j1)=-eye*tx(j1)*dsq*ww1(j1)
            ww2(j1) = wx(j1)*qq2
            ww2y(j1)=-eye*tx(j1)*dsq*ww2(j1)
            ww1(-j1) = dconjg(ww1(j1))
            ww2(-j1) = dconjg(ww2(j1))
            ww1x(-j1)=dconjg(ww1x(j1))
            ww2y(-j1)=dconjg(ww2y(j1))
         enddo
c
         nexp = (2*npw+1)*nns
         do ind=1,nd
            d1=dipstr(ind,i)*rnormal(1,i)
            d2=dipstr(ind,i)*rnormal(2,i)
            j = 0
            do j1=1,nns
               z0=d1*ww1px(j1)+charge(ind,i)*ww1p(j1)
               z1=d1*ww1mx(j1)+charge(ind,i)*ww1m(j1)
               z2=d1*ww2p(j1)
               z3=d1*ww2m(j1)

               c0=d2*ww1p(j1)
               c1=d2*ww1m(j1)
               c2=d2*ww2py(j1)+charge(ind,i)*ww2p(j1)
               c3=d2*ww2my(j1)+charge(ind,i)*ww2m(j1)

               do j2=-npw,npw
                  j = j+1
                  spall(j,ind) = spall(j,ind) + 
     1                z0*ww2(j2) + c0*ww2y(j2)
                  spall(nexp+j,ind) = spall(nexp+j,ind) +
     1                z1*ww2(j2) + c1*ww2y(j2)
                  spall(2*nexp+j,ind) = spall(2*nexp+j,ind) +
     1                z2*ww1x(j2) + c2*ww1(j2)
                  spall(3*nexp+j,ind) = spall(3*nexp+j,ind) +
     1                z3*ww1x(j2) + c3*ww1(j2)
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
c*********************************************************************
C
C evaluate SOE/X expansions (potential, pot + grad, pot + grad + hess)
C
C*********************************************************************
      subroutine g2dsxevalp_vec(nd,delta,center,npw,wx,tx,nns,
     1              ws,ts,spall,targ,ntarg,pot)
C
C     This subroutine evaluates the hybrid SOE/plane wave 
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
      real *8 wx(-npw:npw),tx(-npw:npw)
      complex *16 ws(nns),ts(nns)
      
      real *8 x,y,xpoint,pot(nd,ntarg),delta,dsq
      real *8 facx,facy,pi,w,xq
      complex *16 z,eye,cd
      complex *16 qqx,qqy,qq1,qq2
      complex *16 z0,z1,z2,z3
      
      complex *16 spall(4*(2*npw+1)*nns,nd)
      complex *16 ww1p(100),ww2p(100),ww1m(100),ww2m(100)
      complex *16 ww1(-100:100)
      complex *16 ww2(-100:100)
C
      eye = dcmplx(0,1)
      dsq = 1.0D0/dsqrt(delta)
C
      nexp = (2*npw+1)*nns

      do itarg=1,ntarg
         x = (targ(1,itarg) - center(1))*dsq
         y = (targ(2,itarg) - center(2))*dsq

         qqx = cdexp(eye*tx(1)*x)
         qqy = cdexp(eye*tx(1)*y)
         qq1 = 1.0d0
         qq2 = 1.0d0
         ww1(0) = 1
         ww2(0) = 1
         do j1=1,npw
            qq1 = qq1*qqx
            qq2 = qq2*qqy
            ww1(j1) = qq1
            ww2(j1) = qq2
            ww1(-j1) = dconjg(ww1(j1))
            ww2(-j1) = dconjg(ww2(j1))
         enddo
         
         do j1=1,nns
            ww1p(j1) = cdexp(-ts(j1)*x)
            ww2p(j1) = cdexp(-ts(j1)*y)
            ww1m(j1) = 1/ww1p(j1)
            ww2m(j1) = 1/ww2p(j1)
         enddo
c
         do ind=1,nd
            cd = 0
            j = 0
            do j1=1,nns
               z0=0
               z1=0
               z2=0
               z3=0
               do j2=-npw,npw
                  j = j+1
                  z0=z0+spall(j,ind)*ww2(j2)
                  z1=z1+spall(nexp+j,ind)*ww2(j2)
                  z2=z2+spall(2*nexp+j,ind)*ww1(j2)
                  z3=z3+spall(3*nexp+j,ind)*ww1(j2)
               enddo
               cd=cd+z0*ww1p(j1)+z1*ww1m(j1)+z2*ww2p(j1)+z3*ww2m(j1)
            enddo
            pot(ind,itarg) = pot(ind,itarg) + dble(cd)
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
      subroutine g2dsxevalg_vec(nd,delta,center,npw,wx,tx,nns,
     1              ws,ts,sxall,targ,ntarg,pot,grad)
C
C     This subroutine evaluates the hybrid SOE/plane wave 
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
C     wx,tx         = planewave weights and nodes
C     nns           = number of SOE terms (half of
C                     full quadrature - no complex conjugates)
C     sxall         = hybrid SX expansions in order px,mx,xp,xm
C     targ          = target
C     ws,ts         = SOE weights and nodes
C
C     OUTPUT:
C
C     pot           = potential (or vectorized potentials) incremented
C     grad          = gradient (or vectorized gradients) incremented
C
      implicit none
      integer nd,nns,i,ind,j1,j2,j,npw,nexp,itarg,ntarg
      real *8 delta,pot(nd,ntarg),grad(nd,2,ntarg)
      real *8 center(2),targ(2,ntarg)
      real *8 wx(-npw:npw),tx(-npw:npw)
      complex *16 ws(nns),ts(nns)
      
      real *8 x,y,dsq
      complex *16 eye
      complex *16 qqx,qqy,qq1,qq2
      
      complex *16 sxall(4*(2*npw+1)*nns,nd)
      complex *16 ww1p(100),ww2p(100),ww1m(100),ww2m(100)
      complex *16 ww1px(100),ww2py(100),ww1mx(100),ww2my(100)
      complex *16 ww1(-100:100),ww1x(-100:100)
      complex *16 ww2(-100:100),ww2y(-100:100)
      complex *16, allocatable ::  sumz(:),sumzx(:),sumzy(:)
      complex *16, allocatable :: wtall(:),wtxall(:),wtyall(:)
C
      eye = dcmplx(0,1)
      dsq = 1.0D0/dsqrt(delta)
C
      nexp = (2*npw+1)*nns
      allocate(wtall(4*nexp))
      allocate(wtxall(4*nexp))
      allocate(wtyall(4*nexp))

      allocate(sumz(nd))
      allocate(sumzx(nd))
      allocate(sumzy(nd))
      
      do itarg=1,ntarg
         do ind = 1,nd
            sumz(ind) = 0.0d0
            sumzx(ind) = 0.0d0
            sumzy(ind) = 0.0d0
         enddo
         
         x = (targ(1,itarg) - center(1))*dsq
         y = (targ(2,itarg) - center(2))*dsq

         qqx = cdexp(eye*tx(1)*x)
         qqy = cdexp(eye*tx(1)*y)
         qq1 = 1.0d0
         qq2 = 1.0d0
         ww1(0) = 1
         ww2(0) = 1
         ww1x(0)=0
         ww2y(0)=0
         do j1=1,npw
            qq1 = qq1*qqx
            qq2 = qq2*qqy
            
            ww1(j1) = qq1
            ww2(j1) = qq2
            ww1x(j1)= eye*tx(j1)*dsq*ww1(j1)
            ww2y(j1)= eye*tx(j1)*dsq*ww2(j1)
            
            ww1(-j1) = dconjg(ww1(j1))
            ww2(-j1) = dconjg(ww2(j1))
            ww1x(-j1)= dconjg(ww1x(j1))
            ww2y(-j1)= dconjg(ww2y(j1))
         enddo
         
         do j1=1,nns
            ww1p(j1) = cdexp(-ts(j1)*x)
            ww2p(j1) = cdexp(-ts(j1)*y)
            ww1m(j1) = 1/ww1p(j1)
            ww2m(j1) = 1/ww2p(j1)

            ww1px(j1) = -ts(j1)*dsq*ww1p(j1)
            ww2py(j1) = -ts(j1)*dsq*ww2p(j1)
            ww1mx(j1) =  ts(j1)*dsq*ww1m(j1) 
            ww2my(j1) =  ts(j1)*dsq*ww2m(j1)
         enddo
c
         j = 0
         do j1=1,nns
            do j2=-npw,npw
               j = j+1
               wtall(j) = ww1p(j1)*ww2(j2)
               wtall(nexp+j) = ww1m(j1)*ww2(j2)
               wtall(2*nexp+j) = ww2p(j1)*ww1(j2)
               wtall(3*nexp+j) = ww2m(j1)*ww1(j2)

               wtxall(j) = ww1px(j1)*ww2(j2)
               wtxall(nexp+j) = ww1mx(j1)*ww2(j2)
               wtxall(2*nexp+j) = ww2p(j1)*ww1x(j2)
               wtxall(3*nexp+j) = ww2m(j1)*ww1x(j2)

               wtyall(j) = ww1p(j1)*ww2y(j2)
               wtyall(nexp+j) = ww1m(j1)*ww2y(j2)
               wtyall(2*nexp+j) = ww2py(j1)*ww1(j2)
               wtyall(3*nexp+j) = ww2my(j1)*ww1(j2)
            enddo
         enddo

         do ind=1,nd
         do j=1,4*nexp
            sumz(ind) = sumz(ind) +
     1      wtall(j)*sxall(j,ind)
            sumzx(ind) = sumzx(ind) +
     1      wtxall(j)*sxall(j,ind)
            sumzy(ind) = sumzy(ind) +
     1      wtyall(j)*sxall(j,ind)
         enddo
         enddo
c
         do ind = 1,nd
            pot(ind,itarg) = pot(ind,itarg)+dreal(sumz(ind))
            grad(ind,1,itarg) = grad(ind,1,itarg)+dreal(sumzx(ind))
            grad(ind,2,itarg) = grad(ind,2,itarg)+dreal(sumzy(ind))
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
      subroutine g2dsxevalh_vec(nd,delta,center,npw,wx,tx,nns,
     1              ws,ts,sxall,targ,ntarg,pot,grad,hess)
C
C     This subroutine evaluates the hybrid SOE/plane wave 
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
C     nns           = number of SOE terms (half of
C                     full quadrature - no complex conjugates)
C     sxall         = hybrid SX expansions in order px,mx,xp,xm
C     targ          = target
C     ws,ts         = SOE weights and nodes
C
C     OUTPUT:
C     pot           = potential (or vectorized potentials) incremented
C     grad          = gradient (or vectorized gradients) incremented
C     hess          = Hessian (or vectorized Hessians) incremented
C
      implicit none
      integer nd,nns,npw,ntarg
      real *8 delta,center(2),targ(2,ntarg)
      real *8 pot(nd,ntarg),grad(nd,2,ntarg),hess(nd,3,ntarg)
      
      real *8 wx(-npw:npw),tx(-npw:npw)
      complex *16 ws(nns),ts(nns)
      complex *16 sxall(4*(2*npw+1)*nns,nd)

      integer i,ind,j1,j2,j,nexp,itarg
      real *8 x,y,dsq
      complex *16 eye
      complex *16 qqx,qqy,qq1,qq2
      
      complex *16 ww1p(100),ww2p(100),ww1m(100),ww2m(100)
      complex *16 ww1px(100),ww2py(100),ww1mx(100),ww2my(100)
      complex *16 ww1pxx(100),ww2pyy(100),ww1mxx(100),ww2myy(100)
      
      complex *16 ww1(-100:100),ww1x(-100:100),ww1xx(-100:100)
      complex *16 ww2(-100:100),ww2y(-100:100),ww2yy(-100:100)
      
      complex *16, allocatable ::  sumz(:),sumzx(:),sumzy(:)
      complex *16, allocatable ::  sumzxx(:),sumzxy(:),sumzyy(:)
      complex *16, allocatable :: wtall(:),wtxall(:),wtyall(:)
      complex *16, allocatable :: wtxxall(:),wtxyall(:),wtyyall(:)
C
      eye = dcmplx(0,1)
      dsq = 1.0D0/dsqrt(delta)
C
      nexp = (2*npw+1)*nns
      
      allocate(wtall(4*nexp))
      allocate(wtxall(4*nexp))
      allocate(wtyall(4*nexp))
      allocate(wtxxall(4*nexp))
      allocate(wtxyall(4*nexp))
      allocate(wtyyall(4*nexp))

      allocate(sumz(nd))
      allocate(sumzx(nd))
      allocate(sumzy(nd))
      allocate(sumzxx(nd))
      allocate(sumzxy(nd))
      allocate(sumzyy(nd))
      
      do itarg=1,ntarg
         do ind = 1,nd
            sumz(ind) = 0.0d0
            sumzx(ind) = 0.0d0
            sumzy(ind) = 0.0d0
            sumzxx(ind) = 0.0d0
            sumzxy(ind) = 0.0d0
            sumzyy(ind) = 0.0d0
         enddo
         
         x = (targ(1,itarg) - center(1))*dsq
         y = (targ(2,itarg) - center(2))*dsq

         qqx = cdexp(eye*tx(1)*x)
         qqy = cdexp(eye*tx(1)*y)
         qq1 = 1.0d0
         qq2 = 1.0d0
         
         ww1(0) = 1
         ww2(0) = 1
         ww1x(0)=0
         ww2y(0)=0
         ww1xx(0)=0
         ww2yy(0)=0

         do j1=1,npw
            qq1 = qq1*qqx
            qq2 = qq2*qqy
            
            ww1(j1) = qq1
            ww2(j1) = qq2
            ww1x(j1)= eye*tx(j1)*dsq*ww1(j1)
            ww2y(j1)= eye*tx(j1)*dsq*ww2(j1)
            
            ww1(-j1) = dconjg(ww1(j1))
            ww2(-j1) = dconjg(ww2(j1))
            ww1x(-j1)= dconjg(ww1x(j1))
            ww2y(-j1)= dconjg(ww2y(j1))
            
            ww1xx(j1)= eye*tx(j1)*dsq*ww1x(j1)
            ww2yy(j1)= eye*tx(j1)*dsq*ww2y(j1)
            ww1xx(-j1)= dconjg(ww1xx(j1))
            ww2yy(-j1)= dconjg(ww2yy(j1))
         enddo
c
         do j1=1,nns
            ww1p(j1) = cdexp(-ts(j1)*x)
            ww2p(j1) = cdexp(-ts(j1)*y)
            ww1m(j1) = 1/ww1p(j1)
            ww2m(j1) = 1/ww2p(j1)

            ww1px(j1) = -ts(j1)*dsq*ww1p(j1)
            ww2py(j1) = -ts(j1)*dsq*ww2p(j1)
            ww1mx(j1) =  ts(j1)*dsq*ww1m(j1) 
            ww2my(j1) =  ts(j1)*dsq*ww2m(j1)
c
            ww1pxx(j1) = -ts(j1)*dsq*ww1px(j1)
            ww2pyy(j1) = -ts(j1)*dsq*ww2py(j1)
            ww1mxx(j1) =  ts(j1)*dsq*ww1mx(j1) 
            ww2myy(j1) =  ts(j1)*dsq*ww2my(j1)
         enddo
c
         j = 0
         do j1=1,nns
            do j2=-npw,npw
               j = j+1
               wtall(j) = ww1p(j1)*ww2(j2)
               wtall(nexp+j) = ww1m(j1)*ww2(j2)
               wtall(2*nexp+j) = ww2p(j1)*ww1(j2)
               wtall(3*nexp+j) = ww2m(j1)*ww1(j2)

               wtxall(j) = ww1px(j1)*ww2(j2)
               wtxall(nexp+j) = ww1mx(j1)*ww2(j2)
               wtxall(2*nexp+j) = ww2p(j1)*ww1x(j2)
               wtxall(3*nexp+j) = ww2m(j1)*ww1x(j2)

               wtyall(j) = ww1p(j1)*ww2y(j2)
               wtyall(nexp+j) = ww1m(j1)*ww2y(j2)
               wtyall(2*nexp+j) = ww2py(j1)*ww1(j2)
               wtyall(3*nexp+j) = ww2my(j1)*ww1(j2)
c
               wtxxall(j) = ww1pxx(j1)*ww2(j2)
               wtxxall(nexp+j) = ww1mxx(j1)*ww2(j2)
               wtxxall(2*nexp+j) = ww2p(j1)*ww1xx(j2)
               wtxxall(3*nexp+j) = ww2m(j1)*ww1xx(j2)
c
               wtxyall(j) = ww1px(j1)*ww2y(j2)
               wtxyall(nexp+j) = ww1mx(j1)*ww2y(j2)
               wtxyall(2*nexp+j) = ww2py(j1)*ww1x(j2)
               wtxyall(3*nexp+j) = ww2my(j1)*ww1x(j2)
c
               wtyyall(j) = ww1p(j1)*ww2yy(j2)
               wtyyall(nexp+j) = ww1m(j1)*ww2yy(j2)
               wtyyall(2*nexp+j) = ww2pyy(j1)*ww1(j2)
               wtyyall(3*nexp+j) = ww2myy(j1)*ww1(j2)
            enddo
         enddo
c
         do ind=1,nd
         do j=1,4*nexp
            sumz(ind) = sumz(ind) +
     1          wtall(j)*sxall(j,ind)
            sumzx(ind) = sumzx(ind) +
     1          wtxall(j)*sxall(j,ind)
            sumzy(ind) = sumzy(ind) +
     1          wtyall(j)*sxall(j,ind)
            sumzxx(ind) = sumzxx(ind) +
     1          wtxxall(j)*sxall(j,ind)
            sumzxy(ind) = sumzxy(ind) +
     1          wtxyall(j)*sxall(j,ind)
            sumzyy(ind) = sumzyy(ind) +
     1          wtyyall(j)*sxall(j,ind)
         enddo
         enddo
c
         do ind = 1,nd
            pot(ind,itarg) = pot(ind,itarg)+dreal(sumz(ind))
            grad(ind,1,itarg) = grad(ind,1,itarg)+dreal(sumzx(ind))
            grad(ind,2,itarg) = grad(ind,2,itarg)+dreal(sumzy(ind))
            hess(ind,1,itarg) = hess(ind,1,itarg)+dreal(sumzxx(ind))
            hess(ind,2,itarg) = hess(ind,2,itarg)+dreal(sumzxy(ind))
            hess(ind,3,itarg) = hess(ind,3,itarg)+dreal(sumzyy(ind))
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
      subroutine shiftsx_translation_matrices(boxsize,npw,nsoe,nmax,
     1              wshift,ws,ts)
C
C     This subroutine computes the complex exponentials in the SOE/X translations.
C
C     INPUT
C
c     nlevels    = number of levels in the uniform pruned tree
c     boxsize    = scaled (by 1/sqrt(delta) size of the box at the finest level
C     nsoe       = number of terms in SOE expansion
c     npw        = plane wave expansion length
C     nmax       = number of different translation lengths in the whole scheme 
C     ws,ts      = SOE weights and nodes
C
C     OUTPUT:
C
C     wshift     = table of shift arrays for sx translations, which are used in 
C                  mp to loc stage
c      
      implicit real *8 (a-h,o-z)
      real *8 boxsize
      complex *16, allocatable :: wexp(:,:)
      complex *16 wshift((2*npw+1)*nsoe/2,nmax)
      complex *16 ws(nsoe),ts(nsoe)
      complex *16 ztmp
C
      allocate(wexp(nsoe,nmax))
      
      do j1=1,nsoe
         ztmp = cdexp(-ts(j1)*boxsize)
         do j2=1,nmax
            wexp(j1,j2) = ztmp
            ztmp = ztmp*ztmp
         enddo
      enddo

     
      do k=1,nmax
         j=0
         do j1=1,nsoe/2
            do j2=-npw,npw
               j=j+1
               wshift(j,k)=wexp(j1,k)
            enddo
         enddo
      enddo
      
      return
      end
C      
C      
C      
C      
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
      subroutine shiftsxall_translation_matrices(xmin,npw,nns,nmax,
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
      
      complex *16 wshift(4*(2*npw+1)*nns,4,nmax)
      
      complex *16,allocatable:: wwx(:,:)
      complex *16,allocatable:: wwp(:,:),wwm(:,:)

      complex *16 eye,ztmp
C
      eye = dcmplx(0,1)
      
      allocate(wwp(nns,nmax),wwm(nns,nmax))
      
      allocate(wwx(-npw:npw,nmax))

      do j1=-npw,npw
         ztmp = exp(eye*tx(j1)*xmin)
         do j2=1,nmax
            wwx(j1,j2) = ztmp
            ztmp = ztmp*ztmp
         enddo
      enddo

      do j1=1,nns
         ztmp = exp(-ts(j1)*xmin)
         do j2=1,nmax
            wwp(j1,j2) = ztmp
            wwm(j1,j2) = 1.0d0/wwp(j1,j2)
            ztmp = ztmp*ztmp
         enddo
      enddo

      nexp=(2*npw+1)*nns
      
      do n=1,nmax
         j = 0
         do j1=1,nns
         do j2=-npw,npw
            j = j+1
c           pp            
            wshift(j,       1,n) = wwp(j1,n)*wwx(j2,n)
            wshift(j+nexp  ,1,n) = wwm(j1,n)*wwx(j2,n)
            wshift(j+2*nexp,1,n) = wwx(j2,n)*wwp(j1,n)
            wshift(j+3*nexp,1,n) = wwx(j2,n)*wwm(j1,n)
c           pm            
            wshift(j,       2,n) = wwp(j1,n)/wwx(j2,n)
            wshift(j+nexp  ,2,n) = wwm(j1,n)/wwx(j2,n)
            wshift(j+2*nexp,2,n) = wwx(j2,n)/wwp(j1,n)
            wshift(j+3*nexp,2,n) = wwx(j2,n)/wwm(j1,n)
c           mp            
            wshift(j,       3,n) = wwx(j2,n)/wwp(j1,n)
            wshift(j+nexp  ,3,n) = wwx(j2,n)/wwm(j1,n)
            wshift(j+2*nexp,3,n) = wwp(j1,n)/wwx(j2,n)
            wshift(j+3*nexp,3,n) = wwm(j1,n)/wwx(j2,n)
c           mm            
            wshift(j,       4,n) = 1.0d0/(wwx(j2,n)*wwp(j1,n))
            wshift(j+nexp  ,4,n) = 1.0d0/(wwx(j2,n)*wwm(j1,n))
            wshift(j+2*nexp,4,n) = 1.0d0/(wwx(j2,n)*wwp(j1,n))
            wshift(j+3*nexp,4,n) = 1.0d0/(wwx(j2,n)*wwm(j1,n))
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
      subroutine g2dshiftall_vec(nd,nexp,sxall1,sxall2,wshift)
C
C     This subroutine translate all four expansions (sexpall1)
C     about the center (CENT1) to the center (cent2)
C     using a precomputed translation matrix wshift.
C
c     used in merge mp (from child to parent) 
c     or split loc (from parent to child) stages
c
c     applicable to SOE, SOE/X, or SOE/Hermite expansions
c      
C     INPUT
C
c     nd      = vector length (for multiple charges at same locations)
c     nexp    = length of the whole four SOE/X expansions
C     sxall1  = original expansion FFEXP
C     wshift  = precomputed translation matrix
c      
C     OUTPUT:
C
C     sxall2  = shifted expansion 
C
      implicit real *8 (a-h,o-z)

      complex *16 sxall1(nexp,nd)
      complex *16 sxall2(nexp,nd)
      
      complex *16 wshift(nexp)
C
      do ind=1,nd
         do j=1,nexp
            sxall2(j,ind) = sxall2(j,ind)+sxall1(j,ind)*wshift(j)
         enddo
      enddo
c
      return
      end
c
c
c     
c
C
      subroutine g2dshiftsxall0_vec(nd,delta,npw,nns,sxall1,
     1              cent1,sxall2,cent2,ts,tx)
C
C     This subroutine converts the hybrid SOE/X expansions
C     (sxall1) about the center (CENT1) into a hybrid SOE/X
C     expansion (sxall2) about (CENT2).
C
c     used in merge mp or split loc stage
c      
C     INPUT
C
c     nd      = vector length (for multiple charges at same locations)
C     delta   = Gaussian variance
C     pmax    = plane wave integration param
C     npw     = number of terms in plane wave exp
C     nns     = number of SOE terms
C     sxall1  = original expansion FFEXP
C     cent1   = center of pwexp1
C     cent2   = center pwexp2
C     ts      = SOE nodes
C     tx      = X nodes
c      
C     OUTPUT:
C
C     sxall2  = shifted expansion 
C
      implicit real *8 (a-h,o-z)
      complex *16 ts(nns)
      real *8 tx(-npw:npw)
      
      real *8 cent1(2),cent2(2),x,y,xpoint
      complex *16 sxall1(4*(2*npw+1)*nns,nd)
      complex *16 sxall2(4*(2*npw+1)*nns,nd)
      
      complex *16,allocatable:: wshiftpx(:)
      complex *16,allocatable:: wshiftmx(:)
      complex *16,allocatable:: wshiftxp(:)
      complex *16,allocatable:: wshiftxm(:)

      complex *16 ww1p(20)
      complex *16 ww2p(20)
      complex *16 ww2m(20)
      complex *16 ww1m(20)
      complex *16,allocatable:: wwx(:)
      complex *16,allocatable:: wwy(:)

      complex *16 eye
C
      allocate(wwx(-npw:npw))
      allocate(wwy(-npw:npw))
      
      allocate(wshiftpx((2*npw+1)*nns))
      allocate(wshiftmx((2*npw+1)*nns))
      allocate(wshiftxp((2*npw+1)*nns))
      allocate(wshiftxm((2*npw+1)*nns))

      eye = dcmplx(0,1)
      dsq = 1.0D0/dsqrt(delta)
      x = (cent2(1) - cent1(1))*dsq
      y = (cent2(2) - cent1(2))*dsq
      do j1=1,nns
         ww1p(j1) = cdexp(-ts(j1)*x)
         ww1m(j1) = cdexp(ts(j1)*x)
         ww2p(j1) = cdexp(-ts(j1)*y)
         ww2m(j1) = cdexp(ts(j1)*y)
      enddo

      do j2=-npw,npw
         wwx(j2) = exp(eye*tx(j2)*x)
         wwy(j2) = exp(eye*tx(j2)*y)
      enddo
c
      j = 0
      do j1=1,nns
      do j2=-npw,npw
         j = j+1
         wshiftpx(j) = ww1p(j1)*wwy(j2)
         wshiftmx(j) = ww1m(j1)*wwy(j2)
         wshiftxp(j) = wwx(j2)*ww2p(j1)
         wshiftxm(j) = wwx(j2)*ww2m(j1)
      enddo
      enddo
      
      nexp = (2*npw+1)*nns
      do ind=1,nd
         do j=1,nexp
            sxall2(j,ind) = sxall2(j,ind)+
     1                 sxall1(j,ind)*wshiftpx(j)
            sxall2(j+nexp,ind) = sxall2(j+nexp,ind)+
     1                sxall1(j+nexp,ind)*wshiftmx(j)
            sxall2(j+2*nexp,ind) = sxall2(j+2*nexp,ind)+
     1                sxall1(j+2*nexp,ind)*wshiftxp(j)
            sxall2(j+3*nexp,ind) = sxall2(j+3*nexp,ind)+
     1         sxall1(j+3*nexp,ind)*wshiftxm(j)
         enddo
      enddo
c
      return
      end
c
c
c     
c
C
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

