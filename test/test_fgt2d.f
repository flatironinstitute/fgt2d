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
      real *8, allocatable :: charges(:),dipstr(:)
      real *8, allocatable :: pot(:),grad(:,:),hess(:,:)
      real *8, allocatable :: pottarg(:),gradtarg(:,:),hesstarg(:,:)
      real *8, allocatable :: potex(:),gradex(:,:),hessex(:,:)
      real *8, allocatable :: pottargex(:),gradtargex(:,:),
     1                             hesstargex(:,:)

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
      n = 3
      delta = bb*bb/(1.5*1.5)*4.0d0**n

      
      call prin2(' delta = *',delta,1)
      call prinf_long(' nsrc = *',nsrc,1)
c
      allocate(sources(2,nsrc),charges(nsrc),dipstr(nsrc))
      allocate(rnormal(2,nsrc))
      allocate(targ(2,ntarg))
      allocate(pot(nsrc),grad(2,nsrc),hess(3,nsrc))
      allocate(pottarg(ntarg),gradtarg(2,ntarg),hesstarg(3,ntarg))

      do i=1,nsrc

         sources(1,i) = hkrand(0)
         sources(2,i) = hkrand(0)
         rnormal(1,i) = hkrand(0)
         rnormal(2,i) = hkrand(0)

         charges(i) = hkrand(0) 

         dipstr(i) = hkrand(0) + ima*hkrand(0)

      enddo
      dipstr(1) = 0.0d0 
      dipstr(2) = 0.0d0 
      charges(1) = 0.0d0 
      charges(2) = 0.0d0 
      charges(3) = 1.0d0 
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

      allocate(potex(nts),gradex(2,nts),hessex(3,nts))
      allocate(pottargex(ntt),gradtargex(2,ntt),hesstargex(3,ntt))

      eps = 0.5d-10
c
      call dzero(potex,nts)
      call dzero(pottargex,ntt)

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

cccc      write(6,*) 'pot is',(pot(i),i=1,nts)
cccc      write(6,*) 'gradx is',(grad(1,i),i=1,nts)
cccc      write(6,*) 'grady is',(grad(2,i),i=1,nts)

      reps = 1.0d-100
      call fgt2dpart_direct_vec(1,delta,reps,1,nsrc,1,nts,sources,
     1            ifcharge,charges,ifdipole,rnormal,dipstr,
     2            sources,ifpgh,potex,gradex,hessex)
cccc      write(6,*) 'potex is',(potex(i),i=1,nts)
cccc      write(6,*) 'gradexx is',(gradex(1,i),i=1,nts)
cccc      write(6,*) 'gradexy is',(gradex(2,i),i=1,nts)

      call fgt2dpart_direct_vec(1,delta,reps,1,nsrc,1,ntt,sources,
     1            ifcharge,charges,ifdipole,rnormal,dipstr,
     2            targ,ifpghtarg,pottargex,gradtargex,hesstargex)
cccc      write(6,*) 'pottarg is',(pottarg(i),i=1,ntt)
cccc      write(6,*) 'pottargex is',(pottargex(i),i=1,ntt)

      if (ifpgh .gt. 0) call derr(potex,pot,nts,errps)
      if (ifpgh .gt. 1) call derr(gradex,grad,2*nts,errgs)
      if (ifpgh .gt. 2) call derr(hessex,hess,3*nts,errhs)
      
      if (ifpghtarg .gt. 0) call derr(pottargex,pottarg,ntt,errpt)
      if (ifpghtarg .gt. 1) call derr(gradtargex,gradtarg,2*ntt,errgt)
      if (ifpghtarg .gt. 2) call derr(hesstargex,hesstarg,3*ntt,errht)

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

      erra = sqrt(erra/ra)
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
      
