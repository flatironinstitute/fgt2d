c
c
c
c
c
c
c
c
c
c
c
      implicit real *8 (a-h,o-z)
      real *8, allocatable :: sources(:,:),targ(:,:)
      real *8, allocatable :: rnormal(:,:)
      real *8, allocatable :: charges(:,:),dipstr(:,:)
      real *8, allocatable :: pot(:,:),grad(:,:,:),hess(:,:,:)
      real *8, allocatable :: pottarg(:,:),gradtarg(:,:,:),
     1    hesstarg(:,:,:)
      real *8, allocatable :: potex(:,:),gradex(:,:,:),hessex(:,:,:)
      real *8, allocatable :: pottargex(:,:),gradtargex(:,:,:),
     1                             hesstargex(:,:,:)

      complex *16 ima
      data ima/(0.0d0,1.0d0)/

      call prini(6,13)

      done = 1
      pi = atan(done)*4

      nsrc =  1 000 000
      ntarg = nsrc
      nd = 1
      delta = 5.0d-4
      bb = 1.0d0/(2.0d0**6)
      n = 4
      delta = bb*bb/(1.5*1.5)*2.0d0**n
cccc      delta = 0.1*delta
      
      call prin2(' delta = *',delta,1)
      call prinf_long(' nsrc = *',nsrc,1)
c
      allocate(sources(3,nsrc),charges(nd,nsrc),dipstr(nd,nsrc))
      allocate(rnormal(3,nsrc))
      allocate(targ(3,ntarg))
      allocate(pot(nd,nsrc),grad(nd,3,nsrc),hess(nd,6,nsrc))
      allocate(pottarg(nd,ntarg),gradtarg(nd,3,ntarg))
      allocate(hesstarg(nd,6,ntarg))

      rin = 0.3d0
      rwig = 0.15d0
      nwig = 10
      
      do i=1,nsrc
c        nonuniform source distribution
         theta = hkrand(0)*pi
         phi = hkrand(0)*2*pi
         sources(1,i) = (rin+rwig*cos(theta))*sin(theta)*cos(phi)+0.5d0
         sources(2,i) = (rin+rwig*cos(theta))*sin(theta)*sin(phi)+0.5d0
         sources(3,i) = (rin+rwig*cos(theta))*cos(theta)+0.5d0
         
c        uniform source distribution
         sources(1,i) = hkrand(0)
         sources(2,i) = hkrand(0)
         sources(3,i) = hkrand(0)

         
         rnormal(1,i) = hkrand(0)
         rnormal(2,i) = hkrand(0)
         rnormal(3,i) = hkrand(0)
         do ind = 1,nd
            charges(ind,i) = hkrand(0) 
            dipstr(ind,i) = hkrand(0)
         enddo
      enddo

      do ind=1,nd
         dipstr(ind,1) = 0.0d0 
         dipstr(ind,2) = 0.0d0 
         charges(ind,1) = 0.0d0 
         charges(ind,2) = 0.0d0 
         charges(ind,3) = 1.0d0
      enddo
      
      sources(1,1) = 0.0d0
      sources(2,1) = 0.0d0
      sources(3,1) = 0.0d0

      sources(1,2) = 1.0d0
      sources(2,2) = 1.0d0
      sources(3,2) = 1.0d0

      sources(1,3) = 0.05d0
      sources(2,3) = 0.05d0
      sources(3,3) = 0.05d0
 
      do i=1,ntarg
         targ(1,i) = hkrand(0)
         targ(2,i) = hkrand(0)
         targ(3,i) = hkrand(0)
ccc         targ(1,i) = 0.52d0
ccc         targ(2,i) = 0.52d0
      enddo

      nts = min(20,nsrc)
      ntt = min(20,ntarg)

      allocate(potex(nd,nts),gradex(nd,3,nts),hessex(nd,6,nts))
      allocate(pottargex(nd,ntt),gradtargex(nd,3,ntt))
      allocate(hesstargex(nd,6,ntt))

      eps = 1d-6
c
      call dzero(potex,nts*nd)
      call dzero(gradex,3*nts*nd)
      call dzero(hessex,6*nts*nd)
      
      call dzero(pottargex,ntt*nd)
      call dzero(gradtargex,3*ntt*nd)
      call dzero(hesstargex,6*ntt*nd)
      
      ifcharge = 1
      ifdipole = 0
      ifpgh = 1
      ifpghtarg = 0
      call prinf(' ifcharge is *',ifcharge,1)
      call prinf(' ifdipole is *',ifdipole,1)
      call prinf(' ifpgh is *',ifpgh,1)
      call prinf(' ifpghtarg is *',ifpghtarg,1)

      iper = 0
      call fgt3d(nd,delta,eps,nsrc,sources,ifcharge,charges,
     1    ifdipole,rnormal,dipstr,iper,ifpgh,pot,grad,hess,
     2    ntarg,targ,ifpghtarg,pottarg,gradtarg,
     3    hesstarg)

cccc      call prin2('pot=*',pot,nd*nts)
cccc      call prin2('pottarg=*',pottarg,nd*ntt)
cccc      call prin2('gradtarg=*',gradtarg,nd*3*ntt)
cccc      call prin2('hesstarg=*',hesstarg,nd*6*ntt)

      dmax = 1.0d30
      call cpu_time(t1)
      call fgt3dpart_direct_vec(nd,delta,dmax,1,nsrc,1,nts,sources,
     1    ifcharge,charges,ifdipole,rnormal,dipstr,
     2    sources,ifpgh,potex,gradex,hessex)
cccc      call prin2('potex=*',potex,nts*nd)
      call cpu_time(t2)
      print *, 'direct eval time = ', t2-t1
      call fgt3dpart_direct_vec(nd,delta,dmax,1,nsrc,1,ntt,sources,
     1    ifcharge,charges,ifdipole,rnormal,dipstr,
     2    targ,ifpghtarg,pottargex,gradtargex,hesstargex)
cccc      call prin2('pottargex=*',pottargex,ntt*nd)
cccc      call prin2('gradtargex=*',gradtargex,3*ntt*nd)
cccc      call prin2('hesstargex=*',hesstargex,6*ntt(nd)

      if (ifpgh .gt. 0) call derr(potex,pot,nts*nd,errps)
      if (ifpgh .gt. 1) call derr(gradex,grad,3*nts*nd,errgs)
      if (ifpgh .gt. 2) call derr(hessex,hess,6*nts*nd,errhs)
      
      if (ifpghtarg .gt. 0) call derr(pottargex,pottarg,ntt*nd,errpt)
      if (ifpghtarg .gt. 1)
     1    call derr(gradtargex,gradtarg,3*ntt*nd,errgt)
      if (ifpghtarg .gt. 2) 
     1    call derr(hesstargex,hesstarg,6*ntt*nd,errht)
cc      errgs = 0
cc      errgt = 0

cc      errhs = 0
cc      errht = 0
      call errprint(errps,errgs,errhs,errpt,errgt,errht,ifpgh,ifpghtarg)

      stop
      end
c-----------------------------------------------------     
      subroutine dzero(vec,n)
      implicit real *8 (a-h,o-z)
      real *8 vec(*)

      do i=1,n
         vec(i) = 0
      enddo

      return
      end
c------------------------------------
      subroutine derr(vec1,vec2,n,erra)
      implicit real *8 (a-h,o-z)
      real *8 vec1(*),vec2(*)

      ra = 0
      erra = 0
      do i=1,n
         ra = ra + vec1(i)**2
         erra = erra + (vec1(i)-vec2(i))**2
      enddo

      if (sqrt(ra)/n .lt. 1d-10) then
         call prin2('vector norm =*', sqrt(ra)/n,1)
         call prin2('switch to absolute error*',a,0)
         erra = sqrt(erra)/n
      else
         erra = sqrt(erra/ra)
      endif
ccc      erra = sqrt(erra)/n

      return
      end
c----------------------------------
      subroutine errprint(errps,errgs,errhs,errpt,errgt,errht,
     1    ifpgh,ifpghtarg)
      implicit real *8 (a-h,o-z)
 1100 format(3(2x,e11.5))


      write(6,*) 'ifpgh is ', ifpgh
      write(6,*) 'ifpghtarg is ', ifpghtarg
      if (ifpgh .gt. 0) write(6,*) 'error in sources'
      if (ifpgh .gt. 0) call prin2('pot err =*', errps,1)
      if (ifpgh .gt. 1) call prin2('grad err=*', errgs,1)
      if (ifpgh .gt. 2) call prin2('hess err=*', errhs,1)
      write(6,*) 
      if (ifpghtarg .gt. 0) write(6,* ) 'error in targets'
      if (ifpghtarg .gt. 0) call prin2('pot err =*', errpt,1)
      if (ifpghtarg .gt. 1) call prin2('grad err=*', errgt,1)
      if (ifpghtarg .gt. 2) call prin2('hess err=*', errht,1)
      write(6,*)
      write(6,*)'==================='

      return
      end
      
