c
c
c      |       |       |   Box 2: 5 is pp
c      |       |       |   Box 5: 2 is mm
c      |   4   |   5   |   Box 4: 3 is pm
c      |       |       |   Box 3: 4 is mp
c      |---|---|---|---|
c      |       |       |
c      |   2   |   3   |
c      |       |       |
c
c
c
c      |16 |17 |20 |21 |   Box 6: 9,12,13,15,17  pp
c      |---|---|---|---|   Box 7: 12,13 pp; 8,14,16 mp
c      |14 |15 |18 |19 |   Box 8: 15,17 pp; 7,10,11 pm
c      |---|---|---|---|   Box 9: 14,16 mp, 6 mm, 10,11, pm
c      | 8 | 9 |12 |13 |   Box 10: 8,9 mp; 13,19,21 pp
c      |---|---|---|---|   Box 11: 8,9, 12, 18,20 mp
c      | 6 | 7 |10 |11 |   Box 12: 6,7 mm; 11 pm;19,21 pp
c                          Box 13: 6,7,10 mm; 18,20 mp
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
      
      nsrc = 1 000 000
      ntarg = nsrc
      nd = 1
      delta = 5.0d-2
      bb = 1.0d0/(2.0d0**6)
      
c      
c     n = 1,2,3 bad cases
c      
      n = -4
      delta = bb*bb/(1.5*1.5)*2.0d0**n

      
      call prin2(' delta = *',delta,1)
      call prinf_long(' nsrc = *',nsrc,1)
c
      allocate(sources(2,nsrc),charges(nd,nsrc),dipstr(nd,nsrc))
      allocate(rnormal(2,nsrc))
      allocate(targ(2,ntarg))
      allocate(pot(nd,nsrc),grad(nd,2,nsrc),hess(nd,3,nsrc))
      allocate(pottarg(nd,ntarg),gradtarg(nd,2,ntarg))
      allocate(hesstarg(nd,3,ntarg))


      rin = 0.3d0
      rwig = 0.15d0
      nwig = 10
      
      do i=1,nsrc
         
         thet = hkrand(0)*2*pi
c        nonuniform distribution of sources
         sources(1,i) = (rin + rwig*cos(thet))*cos(thet)+0.5d0
         sources(2,i) = (rin + rwig*cos(thet))*sin(thet)+0.5d0

         sources(1,i) = hkrand(0)
         sources(2,i) = hkrand(0)

         rnormal(1,i) = hkrand(0)
         rnormal(2,i) = hkrand(0)

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
      sources(1,2) = 1.0d0
      sources(2,2) = 1.0d0
      sources(1,3) = 0.05d0
      sources(2,3) = 0.05d0

      do i=1,ntarg
         targ(1,i) = hkrand(0)
         targ(2,i) = hkrand(0)
ccc         targ(1,i) = 0.52d0
ccc         targ(2,i) = 0.52d0
      enddo

      nts = min(20,nsrc)
      ntt = min(20,ntarg)
      
cccc      nts = min(1000,nsrc)
cccc      ntt = min(1000,ntarg)
      allocate(potex(nd,nts),gradex(nd,2,nts),hessex(nd,3,nts))
      allocate(pottargex(nd,ntt),gradtargex(nd,2,ntt))
      allocate(hesstargex(nd,3,ntt))

      eps = 0.5d-10
c
      call dzero(potex,nts*nd)
      call dzero(gradex,2*nts*nd)
      call dzero(hessex,3*nts*nd)
      
      call dzero(pottargex,ntt*nd)
      call dzero(gradtargex,2*ntt*nd)
      call dzero(hesstargex,3*ntt*nd)


      ifcharge = 1
      ifdipole = 0
      ifpgh = 1
      ifpghtarg = 0
      call prinf(' ifcharge is *',ifcharge,1)
      call prinf(' ifdipole is *',ifdipole,1)
      call prinf(' ifpgh is *',ifpgh,1)
      call prinf(' ifpghtarg is *',ifpghtarg,1)

      iper = 0
      call fgt2d(nd,delta,eps,nsrc,sources,ifcharge,charges,
     1    ifdipole,rnormal,dipstr,iper,ifpgh,pot,grad,hess,
     2    ntarg,targ,ifpghtarg,pottarg,gradtarg,
     3    hesstarg)
cccc      call prin2('pot=*',pot,nd*nts)
cccc      call prin2('pottarg=*',pottarg,nd*ntt)
cccc      call prin2('gradtarg=*',gradtarg,nd*2*ntt)
cccc      call prin2('hesstarg=*',hesstarg,nd*3*ntt)

      reps = 1.0d-100
      call fgt2dpart_direct_vec(nd,delta,reps,1,nsrc,1,nts,sources,
     1    ifcharge,charges,ifdipole,rnormal,dipstr,
     2    sources,ifpgh,potex,gradex,hessex)
cccc      call prin2('potex=*',potex,nts*nd)

      call fgt2dpart_direct_vec(nd,delta,reps,1,nsrc,1,ntt,sources,
     1    ifcharge,charges,ifdipole,rnormal,dipstr,
     2    targ,ifpghtarg,pottargex,gradtargex,hesstargex)
cccc      call prin2('pottargex=*',pottargex,ntt*nd)
cccc      call prin2('gradtargex=*',gradtargex,2*ntt*nd)
cccc      call prin2('hesstargex=*',hesstargex,3*ntt(nd)

      if (ifpgh .gt. 0) call derr(potex,pot,nts*nd,errps)
      if (ifpgh .gt. 1) call derr(gradex,grad,2*nts*nd,errgs)
      if (ifpgh .gt. 2) call derr(hessex,hess,3*nts*nd,errhs)
      
      if (ifpghtarg .gt. 0) call derr(pottargex,pottarg,ntt*nd,errpt)
      if (ifpghtarg .gt. 1)
     1    call derr(gradtargex,gradtarg,2*ntt*nd,errgt)
      if (ifpghtarg .gt. 2) 
     1    call derr(hesstargex,hesstarg,3*ntt*nd,errht)
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
ccc      

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
      
