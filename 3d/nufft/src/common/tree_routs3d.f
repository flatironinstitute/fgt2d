c
c
c    common routines for generating and processing
c     a level restricted quad tree in 2D
c   
c
c
c
      subroutine tree_refine_boxes(irefinebox,nboxes,
     1  ifirstbox,nbloc,centers,bs,nbctr,nlctr,
     2  ilevel,iparent,nchild,ichild)
      implicit none
      integer nboxes,nbloc,nbctr,nlctr
      real *8 centers(3,nboxes),bs
      integer ilevel(nboxes),iparent(nboxes)
      integer ichild(8,nboxes),nchild(nboxes)
      integer irefinebox(nbloc)
      integer ifirstbox
      integer, allocatable :: isum(:)
      integer ii

      integer i,ibox,nel0,j,l,jbox,nel1,nbl,jj

      allocate(isum(nbloc))
      if(nbloc.gt.0) call cumsum(nbloc,irefinebox,isum)
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,nbl,j,jbox,ii,jj,l)
      do i = 1,nbloc
        ibox = ifirstbox + i-1
        if(irefinebox(i).eq.1) then
          nbl = nbctr + (isum(i)-1)*8
          
          nchild(ibox) = 8
          do j=1,8
            ii = 2
            jj = 2
            if(j.eq.1.or.j.eq.2.or.j.eq.5.or.j.eq.6) ii = 1
            if(j.lt.5) jj = 1
            jbox = nbl+j
            centers(1,jbox) = centers(1,ibox)+(-1)**j*bs/2
            centers(2,jbox) = centers(2,ibox)+(-1)**ii*bs/2 
            centers(3,jbox) = centers(3,ibox)+(-1)**jj*bs/2 
            iparent(jbox) = ibox
            nchild(jbox) = 0
            do l=1,8
              ichild(l,jbox) = -1
            enddo
            ichild(j,ibox) = jbox
            ilevel(jbox) = nlctr 
          enddo
        endif
      enddo
C$OMP END PARALLEL DO      

      if(nbloc.gt.0) nbctr = nbctr + isum(nbloc)*8


      return
      end
c
c
c
c
c
c
       subroutine tree_copy(nb,centers,ilevel,iparent,nchild,ichild,
     1              centers2,ilevel2,iparent2,nchild2,ichild2)

       implicit none
       integer nd,nb,npb
       real *8 centers(3,nb),centers2(3,nb)
       integer ilevel(nb),ilevel2(nb)
       integer iparent(nb),iparent2(nb)
       integer nchild(nb),nchild2(nb)
       integer ichild(8,nb),ichild2(8,nb)

       integer i,j,nel


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
       do i=1,nb
         centers2(1,i) = centers(1,i)
         centers2(2,i) = centers(2,i)
         centers2(3,i) = centers(3,i)
         ilevel2(i) = ilevel(i)
         iparent2(i) = iparent(i)
         nchild2(i) = nchild(i)
         ichild2(1,i) = ichild(1,i)
         ichild2(2,i) = ichild(2,i)
         ichild2(3,i) = ichild(3,i)
         ichild2(4,i) = ichild(4,i)
         ichild2(5,i) = ichild(5,i)
         ichild2(6,i) = ichild(6,i)
         ichild2(7,i) = ichild(7,i)
         ichild2(8,i) = ichild(8,i)
       enddo
C$OMP END PARALLEL DO       
       

       return
       end
c
c
c
c
c

      subroutine computecoll(nlevels,nboxes,laddr,boxsize,
     1                       centers,iparent,nchild,ichild,
     2                       iper,nnbors,nbors)

c     This subroutine computes the colleagues for an adaptive
c     pruned tree. box j is a colleague of box i, if they share a
c     vertex or an edge and the two boxes are at the same
c     level in the tree
c
c     INPUT arguments
c     nlevels     in: integer
c                 Number of levels
c
c     nboxes      in: integer
c                 Total number of boxes
c
c     laddr       in: integer(2,0:nlevels)
c                 indexing array providing access to boxes at
c                 each level. 
c                 the first box on level i is laddr(1,i)
c                 the last box on level i is laddr(2,i)
c
c     boxsize     in: double precision(0:nlevels)
c                 Array of boxsizes
c 
c     centers     in: double precision(3,nboxes)
c                 array of centers of boxes
c   
c     iparent     in: integer(nboxes)
c                 iparent(i) is the box number of the parent of
c                 box i
c
c     nchild      in: integer(nboxes)
c                 nchild(i) is the number of children of box i
c
c     ichild      in: integer(8,nboxes)
c                 ichild(j,i) is the box id of the jth child of
c                 box i
c
c     iper        in: integer
c                 flag for periodic implementations. 
c                 Currently not used. Feature under construction.
c
c----------------------------------------------------------------
c     OUTPUT
c     nnbors      out: integer(nboxes)
c                 nnbors(i) is the number of colleague boxes of
c                 box i
c
c     nbors       out: integer(27,nboxes)
c                 nbors(j,i) is the box id of the jth colleague
c                 box of box i
c---------------------------------------------------------------
      implicit none
      integer nlevels,nboxes,iper
      integer laddr(2,0:nlevels)
      double precision boxsize(0:nlevels)
      double precision centers(3,nboxes)
      integer iparent(nboxes), nchild(nboxes), ichild(8,nboxes)
      integer nnbors(nboxes)
      integer nbors(27,nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,dad
      integer i,j,ifirstbox,ilastbox


c     Setting parameters for level = 0
      nnbors(1) = 1
      nbors(1,1) = 1
      do ilev = 1,nlevels
c        Find the first and the last box at level ilev      
         ifirstbox = laddr(1,ilev)
         ilastbox = laddr(2,ilev)
c        Loop over all boxes to evaluate neighbors, list1 and updating
c        hunglists of targets

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,dad,i,jbox,j,kbox)
         do ibox = ifirstbox,ilastbox
c           Find the parent of the current box         
            dad = iparent(ibox)
c           Loop over the neighbors of the parent box
c           to find out list 1 and list 2
            do i=1,nnbors(dad)
                jbox = nbors(i,dad)
                do j=1,8
c               ichild(j,jbox) is one of the children of the
c               neighbors of the parent of the current
c               box
                   kbox = ichild(j,jbox)
                   if(kbox.gt.0) then
c               Check if kbox is a nearest neighbor or in list 2
                      if((abs(centers(1,kbox)-centers(1,ibox)).le.
     1                   1.05*boxsize(ilev)).and.
     2                   (abs(centers(2,kbox)-centers(2,ibox)).le.
     3                   1.05*boxsize(ilev)).and.
     2                   (abs(centers(3,kbox)-centers(3,ibox)).le.
     3                   1.05*boxsize(ilev))) then
                     
                         nnbors(ibox) = nnbors(ibox)+1
                         nbors(nnbors(ibox),ibox) = kbox
                      endif
                   endif
                enddo
            enddo
c           End of computing colleagues of box i
         enddo
C$OMP END PARALLEL DO         
      enddo

      return
      end
c
c
c
c
c--------------------------------------------------------------------      
      subroutine updateflags(curlev,nboxes,nlevels,laddr,nchild,ichild,
     1                    nnbors,nbors,centers,boxsize,iflag)

c      This subroutine is to check the boxes flagged as flag++
c      and determine which of the boxes need refinement. The flag
c      of the box which need refinement is updated to iflag(box)=1
c      and that of the boxes which do not need refinement is
c      updated to iflag(box) = 0
c
c      INPUT arguments
c      curlev         in: integer
c                     the level for which boxes need to be processed
c
c      nboxes         in: integer
c                     total number of boxes
c
c      nlevels        in: integer
c                     total number of levels
c
c      laddr          in: integer(2,0:nlevels)
c                     boxes from laddr(1,ilev) to laddr(2,ilev)
c                     are at level ilev
c
c      nchild         in: integer(nboxes)
c                     nchild(ibox) is the number of children
c                     of box ibox
c
c      ichild         in: integer(8,nboxes)
c                     ichild(j,ibox) is the box id of the jth
c                     child of box ibox
c
c      nnbors         in: integer(nboxes)
c                     nnbors(ibox) is the number of colleagues
c                     of box ibox
c
c      nbors          in: integer(27,nboxes)
c                     nbors(j,ibox) is the jth colleague of box
c                     ibox
c
c      centers        in: double precision(2,nboxes)
c                     x and y coordinates of the box centers
c
c      boxsize        in: double precision(0:nlevels)
c                     boxsize(i) is the size of the box at level i
c
c      iflag          in/out: integer(nboxes)
c                     iflag(ibox)=3 if it is flag++. iflag(ibox) =1
c                     or 0 at the end of routine depending on
c                     whether box needs to be subdivided or not
c
      implicit none
c     Calling sequence variables
      integer curlev, nboxes, nlevels
      integer laddr(2,0:nlevels),nchild(nboxes),ichild(8,nboxes)
      integer nnbors(nboxes), nbors(27,nboxes)
      integer iflag(nboxes)
      double precision centers(3,nboxes),boxsize(0:nlevels)

c     Temporary variables
      integer i,j,k,l,ibox,jbox,kbox,lbox, ict
      double precision distest,xdis,ydis,zdis

      distest = 1.05d0*(boxsize(curlev) + boxsize(curlev+1))/2.0d0
c     Loop over all boxes at the current level     

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,i,jbox,j,kbox,xdis,ydis)
C$OMP$PRIVATE(zdis,ict)
      do ibox = laddr(1,curlev),laddr(2,curlev)
         if(iflag(ibox).eq.3) then
            iflag(ibox) = 0
c           Loop over colleagues of the current box      
            do i=1,nnbors(ibox)
c              Loop over colleagues of flag++ box        
               jbox = nbors(i,ibox)
              
c              Loop over the children of the colleague box
c              Note we do not need to exclude self from
c              the list of colleagues as a self box which
c              is flag++ does not have any children 
c              and will not enter the next loop
               do j=1,8
                  kbox = ichild(j,jbox)
                  if(kbox.gt.0) then
                     if(nchild(kbox).gt.0) then
                        xdis = centers(1,kbox) - centers(1,ibox)
                        ydis = centers(2,kbox) - centers(2,ibox)
                        zdis = centers(3,kbox) - centers(3,ibox)
                        ict = 0
                        if(abs(xdis).le.distest) ict = ict + 1
                        if(abs(ydis).le.distest) ict = ict + 1
                        if(abs(zdis).le.distest) ict = ict + 1
                        if(ict.eq.3) then
                           iflag(ibox) = 1
                           goto 1111
                        endif
                     endif
                  endif
c                 End of looping over the children of the child
c                 of the colleague box
               enddo
c              End of looping over the children of the colleague box       
            enddo
c           End of looping over colleagues            
 1111       continue        
         endif
c        End of testing if the current box needs to checked for         
      enddo
c     End of looping over boxes at the current level      
C$OMP END PARALLEL DO      

      return
      end
c
c
c
c
c
c

      subroutine tree_refine_boxes_flag(iflag,nboxes,
     1  ifirstbox,nbloc,centers,bs,nbctr,nlctr,
     2  ilevel,iparent,nchild,ichild)
      implicit none
      integer nboxes,nbloc,nbctr,nlctr
      real *8 centers(3,nboxes),bs
      integer ilevel(nboxes),iparent(nboxes)
      integer ichild(8,nboxes),nchild(nboxes)
      integer iflag(nboxes)
      integer ifirstbox
      integer, allocatable :: isum(:),itmp(:)

      integer i,ibox,nel0,j,l,jbox,nel1,nbl
      integer ii,jj

      allocate(isum(nbloc),itmp(nbloc))
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,i)
      do i=1,nbloc
        ibox = ifirstbox+i-1
        itmp(i) = 0
        if(iflag(ibox).gt.0) itmp(i) = 1
      enddo
C$OMP END PARALLEL DO
      if(nbloc.gt.0) call cumsum(nbloc,itmp,isum)
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,nbl,j,jbox,ii,jj,l)
      do i = 1,nbloc
        ibox = ifirstbox + i-1
        if(iflag(ibox).gt.0) then
          nbl = nbctr + (isum(i)-1)*8
          
          nchild(ibox) = 8
          do j=1,8
            ii = 2
            jj = 2
            if(j.eq.1.or.j.eq.2.or.j.eq.5.or.j.eq.6) ii = 1
            if(j.lt.5) jj = 1
            jbox = nbl+j
            centers(1,jbox) = centers(1,ibox)+(-1)**j*bs/2
            centers(2,jbox) = centers(2,ibox)+(-1)**ii*bs/2 
            centers(3,jbox) = centers(3,ibox)+(-1)**jj*bs/2 
            iparent(jbox) = ibox
            nchild(jbox) = 0
            do l=1,8
              ichild(l,jbox) = -1
            enddo
            ichild(j,ibox) = jbox
            ilevel(jbox) = nlctr+1 
            if(iflag(ibox).eq.1) iflag(jbox) = 3
            if(iflag(ibox).eq.2) iflag(jbox) = 0
          enddo
        endif
      enddo
C$OMP END PARALLEL DO      

      if(nbloc.gt.0) nbctr = nbctr + isum(nbloc)*8


      return
      end
c
c
c
c
c
      subroutine computemnlist1(nlevels,nboxes,laddr,boxsize,
     1                   centers,iparent,nchild,
     2                   ichild,isep,nnbors,mnbors,nbors,iper,mnlist1)
c     Compute max nuber of boxes in list1
      implicit none
      integer nlevels,nboxes
      integer iper
      integer laddr(2,0:nlevels)
      double precision boxsize(0:nlevels)
      double precision centers(3,nboxes)
      integer iparent(nboxes),nchild(nboxes),ichild(8,nboxes)
      integer mnbors,isep
      integer nnbors(nboxes),nbors(mnbors,nboxes)
      integer mnlist1
      integer nlist1(nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,i,j,k,l
      integer firstbox,lastbox,dad
      double precision xdis,ydis,zdis,distest


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nboxes
         nlist1(i) = 0
      enddo
C$OMP END PARALLEL DO
      
      nlist1(1) = 1

      do ilev = 1,nlevels
         firstbox = laddr(1,ilev)
         lastbox = laddr(2,ilev)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,dad,i,jbox,j,kbox,distest,xdis,ydis,zdis)
         do ibox = firstbox,lastbox
            dad = iparent(ibox)
c           Compute list1 and list3 of ibox if it is childless
            if(nchild(ibox).eq.0) then
               do i=1,nnbors(ibox)
                  jbox = nbors(i,ibox)

c
cc                     check for list1 at the same level
c
                  if(nchild(jbox).eq.0) then
                     nlist1(ibox) = nlist1(ibox) + 1
                  endif
c
cc                     check for list1 and list3 at one ilev+1
                  if(nchild(jbox).gt.0) then
                     distest = 1.05d0*(boxsize(ilev)+boxsize(ilev+1))/
     1                         2.0d0*isep
                     do j=1,8
                        kbox = ichild(j,jbox)
                        if(kbox.gt.0) then
                           xdis = dabs(centers(1,kbox)-centers(1,ibox))
                           ydis = dabs(centers(2,kbox)-centers(2,ibox))
                           zdis = dabs(centers(3,kbox)-centers(3,ibox))

                           if(xdis.lt.distest.and.ydis.lt.distest.and.
     1                        zdis.lt.distest) then
                              nlist1(ibox) = nlist1(ibox)+1
                           endif
                        endif
                     enddo
                  endif
               enddo
c
cc               compute list1 and list4 for boxes at level ilev-1 
               do i=1,nnbors(dad)
                   jbox = nbors(i,dad)
                   if(nchild(jbox).eq.0) then
                      distest = 1.05d0*(boxsize(ilev)+boxsize(ilev-1))/
     1                         2.0d0*isep
                      xdis = dabs(centers(1,jbox)-centers(1,ibox))
                      ydis = dabs(centers(2,jbox)-centers(2,ibox))
                      zdis = dabs(centers(3,jbox)-centers(3,ibox))
                      if(xdis.lt.distest.and.ydis.lt.distest.and.
     1                  zdis.lt.distest) then
                         nlist1(ibox) = nlist1(ibox)+1
                      endif
                   endif
               enddo
            endif
         enddo
C$OMP END PARALLEL DO         
      enddo

      mnlist1 = 0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) 
C$OMP$REDUCTION(max:mnlist1,mnlist2,mnlist3,mnlist4)      
      do i=1,nboxes
         if(nlist1(i).gt.mnlist1) mnlist1 = nlist1(i)
      enddo
C$OMP END PARALLEL DO      

      return
      end
c---------------------------------------------------------------      
      subroutine compute_modified_list1(nlevels,npwlevel,
     1                   nboxes,laddr,boxsize,
     1                   centers,iparent,nchild,
     2                   ichild,isep,nnbors,mnbors,nbors,iper,nlist1,
     3                   mnlist1,list1)
c     Compute max nuber of boxes in list1
      implicit none
      integer nlevels,npwlevel,nboxes
      integer iper
      integer laddr(2,0:nlevels)
      double precision boxsize(0:nlevels)
      double precision centers(3,nboxes)
      integer iparent(nboxes),nchild(nboxes),ichild(8,nboxes)
      integer mnbors
      integer nnbors(nboxes),nbors(mnbors,nboxes)
      integer mnlist1,isep
      integer nlist1(nboxes)
      integer list1(mnlist1,nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,i,j,k,l
      integer firstbox,lastbox,dad
      double precision xdis,ydis,zdis,distest

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,nboxes
         nlist1(i) = 0
      enddo
C$OMP END PARALLEL DO      
      if(nchild(1).eq.0) then
         nlist1(1) = 1
         list1(1,1) = 1
      else
         nlist1(1) = 0
      endif

      do ilev = 1,nlevels
        firstbox = laddr(1,ilev)
        lastbox = laddr(2,ilev)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,dad,i,jbox,j,kbox,xdis,ydis,zdis,distest)
        do ibox = firstbox,lastbox
          dad = iparent(ibox)
c         Compute list1 of ibox if it is childless
          if(nchild(ibox).eq.0) then
            do i=1,nnbors(ibox)
              jbox = nbors(i,ibox)
c
cc            boxes in list 1 at the same level
c
              if (ilev .eq. npwlevel) then
c               at the cutoff level, colleagues are directly added to list1
                nlist1(ibox) = nlist1(ibox)+1
                list1(nlist1(ibox),ibox) = jbox
              else
                if(nchild(jbox).eq.0) then
                   nlist1(ibox) = nlist1(ibox) + 1
                   list1(nlist1(ibox),ibox) = jbox
                else
c
c                  boxes in list1 at level ilev+1
c
                   distest = 1.05d0*(boxsize(ilev)+boxsize(ilev+1))/
     1                 2.0d0*isep
                   do j=1,8
                      kbox = ichild(j,jbox)
                      if(kbox.gt.0) then
                         xdis = dabs(centers(1,kbox)-centers(1,ibox))
                         ydis = dabs(centers(2,kbox)-centers(2,ibox))
                         zdis = dabs(centers(3,kbox)-centers(3,ibox))

                         if(xdis.lt.distest.and.ydis.lt.distest.and.
     1                       zdis.lt.distest) then
                            nlist1(ibox) = nlist1(ibox)+1
                            list1(nlist1(ibox),ibox) = kbox
                         endif
                      endif
                   enddo

                endif
              endif
            enddo
c     
cc          compute list1 at level ilev-1 
            do i=1,nnbors(dad)
               jbox = nbors(i,dad)
               if(nchild(jbox).eq.0) then
                  distest = 1.05d0*(boxsize(ilev)+boxsize(ilev-1))/
     1                2.0d0*isep
                  xdis = dabs(centers(1,jbox)-centers(1,ibox))
                  ydis = dabs(centers(2,jbox)-centers(2,ibox))
                  zdis = dabs(centers(3,jbox)-centers(3,ibox))
                  if(xdis.lt.distest.and.ydis.lt.distest.and.
     1                zdis.lt.distest) then
                     nlist1(ibox) = nlist1(ibox)+1
                     list1(nlist1(ibox),ibox) = jbox
                  endif
               endif
            enddo
          endif
        enddo
C$OMP END PARALLEL DO         
      enddo

      return
      end
c----------------------------------------------------------------
c      
c      
c      
c      
      subroutine gt3d_computemnlistpw(nlevels,nboxes,itree,ltree,
     1   iptr,centers,
     1   boxsize,iper,mnlistpw)
c
c     determine maximum number of elements in listpw
c
c     NOTE: we use max values
c
      implicit real *8 (a-h,o-z)
      integer *8 ltree, iptr(8)
      integer nlevels,nboxes,itree(ltree)
      real *8 centers(3,nboxes),boxsize(0:nlevels)
      integer mnlistpw

      mnlistpw = 26

      return
      end
c
c
c
c
c
      subroutine gt3d_computelistpw(nlevels,npwlevel,nboxes,
     1    itree,ltree,iptr,centers,nchild,ichild,
     2    boxsize,laddr,mnlistpw,nlistpw,listpw)
c     this subroutine returns lists of pw interaction boxes
c
c     input parameters:
c
c     nlevels = number of levels, level 0 contains boxes of size 6\sqrt{delta}
c     npwlevel   = recursive SOE picture stops at the level npwlevel
c     nboxes = total number of boxes in the tree
c     itree - laddr: tree stuff
c     nchild         in: integer(nboxes)
c                     nchild(ibox) is the number of children
c                     of box ibox
c
c     ichild         in: integer(8,nboxes)
c                     ichild(j,ibox) is the box id of the jth
c                     child of box ibox
c
c
c     mnlistpw = maximum number of a particular type of PW expansions
c      
c     output parameters:
c
c     nlistpw (nboxes) = contains the number of boxes for the PW interaction
c                            for each box
c     listpw (mnlistpw,nboxes) = contains the box ID for the PW interaction
c                            for each box
c      
      implicit real *8 (a-h,o-z)
      integer *8 iptr(8), ltree
      integer nlevels,nboxes
      integer itree(ltree),laddr(2,0:nlevels)
      real *8 centers(3,nboxes),boxsize(0:nlevels)
      integer nchild(nboxes),ichild(8,nboxes)
      integer mnlistpw
      integer nlistpw(nboxes), listpw(mnlistpw,nboxes)

      integer ibox

      do ibox=1,nboxes
         nlistpw(ibox)=0
      enddo

      nlevstart=max(npwlevel,0)
      do ilev=npwlevel,npwlevel
         do ibox = laddr(1,ilev),laddr(2,ilev)
            ncoll = itree(iptr(6)+ibox-1)
            do i=1,ncoll
               jbox = itree(iptr(7) + (ibox-1)*27+i-1)
               if (jbox.ne.ibox .and. nchild(jbox).gt.0) then
                  nlistpw(ibox)=nlistpw(ibox)+1
                  listpw(nlistpw(ibox),ibox) = jbox
               endif
            enddo
         enddo
      enddo

      return
      end
c
c
c
c
c
      
