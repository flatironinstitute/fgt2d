cc Copyright (C) 2020-2021: Leslie Greengard, Shidong Jiang, Manas Rachh 
cc Contact: lgreengard@flatironinstitute.org
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
c    $Date$
c    $Revision$

      subroutine fgt2d(nd,delta,eps,ns,sources,ifcharge,charge,
     1            ifdipole,rnormal,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of FGTs (same source and target locations, 
c                   different charge, dipole strengths)
c   delta         : Gaussian variance 
c   eps           : precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   ifcharge      : flag for including charge interactions
c                   charge interactions included if ifcharge =1
c                   not included otherwise
c   charge(nd,ns) : charge strengths
c   ifdipole      : flag for including dipole interactions
c                   dipole interactions included if ifcharge =1
c                   not included otherwise
c   rnormal(2,ns) : dipole directions
c   dipstr(nd,ns) : dipole strengths
c   iper          : flag for periodic implmentations. Currently unused
c   ifpgh         : flag for computing pot/grad/hess
c                   ifpgh = 1, only potential is computed
c                   ifpgh = 2, potential and gradient are computed
c                   ifpgh = 3, potential, gradient, and hessian 
c                   are computed
c   nt            : number of targets
c   targ(2,nt)    : target locations
c   ifpghtarg     : flag for computing pottarg/gradtarg/hesstarg
c                   ifpghtarg = 1, only potential is computed at targets
c                   ifpghtarg = 2, potential and gradient are 
c                   computed at targets
c                   ifpghtarg = 3, potential, gradient, and hessian are 
c                   computed at targets
c
c   OUTPUT PARAMETERS
c   pot(nd,*)       : potential at the source locations
c   grad(nd,2,*)    : gradients at the source locations
c   hess(nd,3,*)    : hessian at the source locations
c   pottarg(nd,*)   : potential at the target locations
c   gradtarg(nd,2,*): gradient at the target locations
c   hesstarg(nd,3,*): hessian at the target locations
c
      implicit none
c
cc      calling sequence variables
c 
      integer nd
      real *8 eps,delta
      integer ns,nt
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      real *8 sources(2,ns),targ(2,nt)
      real *8 rnormal(2,ns)
      real *8 charge(nd,*),dipstr(nd,*)

      real *8 pot(nd,*),grad(nd,2,*),hess(nd,3,*)
      real *8 pottarg(nd,*),gradtarg(nd,2,*),hesstarg(nd,3,*)

c
cc      Tree variables
c
      integer, allocatable :: itree(:)
      integer iptr(8)
      integer iper,nlmin,nlmax,ifunif
      real *8, allocatable :: tcenters(:,:),boxsize(:)
      integer idivflag,nlevels,nboxes,ndiv,ndiv0
      integer ltree,npwlevel

c
cc     sorted arrays
c
      integer, allocatable :: isrc(:),isrcse(:,:)
      integer, allocatable :: itarg(:),itargse(:,:)
      real *8, allocatable :: sourcesort(:,:)
      real *8, allocatable :: rnormalsort(:,:)
      real *8, allocatable :: targsort(:,:)
      real *8, allocatable :: chargesort(:,:),dipstrsort(:,:)
      real *8, allocatable :: potsort(:,:),gradsort(:,:,:),
     1                             hesssort(:,:,:)
      real *8, allocatable :: pottargsort(:,:),gradtargsort(:,:,:),
     1                              hesstargsort(:,:,:)

c
cc     additional fmm variables

      integer *8 lmptot
      integer ntermmax
      integer, allocatable :: iaddr(:,:)
      integer, allocatable :: ntermsh(:),nlocal(:)
      real *8, allocatable :: rmlexp(:)

c
cc      temporary variables
c
      integer i,ilev,lmptmp,idim
      integer nlocal0,npw,nsoe,nadd,ifprint,ier
      real *8 omp_get_wtime
      real *8 time1,time2,pi,done,pmax,bs0,cen0(2),bsize

      done = 1
      pi = atan(done)*4.0d0
      
      ifprint = 1

      call pts_tree_boxsize0(delta,eps,sources,ns,targ,nt,
     1    npwlevel,bs0,cen0)
      write(6,*) ' npwlevel',npwlevel

C
C     if box is too small as compared with delta, there is no need
C     to do anything - simply form local, then eval it.
c      
      if (npwlevel .lt. -2) then
         call g2dhlterms(2,bs0,delta,eps,nlocal0)
         nlocal0 = nlocal0 + max(ifpgh,ifpghtarg) -1 
         write(6,*) ' nlocal0',nlocal0
      
         call cpu_time(time1)
C$      time1=omp_get_wtime()
         call fgt2d_large_delta(nd,delta,eps,ns,sources,
     1       ifcharge,charge,ifdipole,rnormal,dipstr,iper,
     2       nlocal0,cen0,
     3       ifpgh,pot,grad,hess,nt,targ,
     4       ifpghtarg,pottarg,gradtarg,hesstarg)
         call cpu_time(time2)
C$        time2=omp_get_wtime()
         if( ifprint .eq. 1 ) call prin2('time in fgt main=*',
     1       time2-time1,1)
         return
      endif     
c
c     need to fix the dependence of ndiv on eps
c   
c
      idivflag =0
c     ndiv0 is the maximum number of points per box above the cutoff level
c     it determines the speed of the algorithm when delta goes to zero.
      ndiv0 = 10
c     ndiv is the maximum number of points per box at or below the cutoff level
c     it's determined by numerical experiments on finding the crossover point
c     between direct evaluation and the fast scheme.
      ndiv = 80
c
      ifunif = 0
      iper = 0

c
c     call the tree memory management
c     code to determine number of boxes,
c     number of levels and length of tree
c
      nlmax = npwlevel + 4
      nlmin = 0
      call pts_tree_mem(sources,ns,targ,nt,idivflag,
     1    ndiv,nlmin,nlmax,iper,
     2    ndiv0,npwlevel,bs0,cen0,nlevels,nboxes,ltree) 

c 
      write(6,*) ' nlmax',nlmax
      write(6,*) ' nlevels',nlevels
      write(6,*) ' nboxes',nboxes

      allocate(itree(ltree))
      allocate(boxsize(0:nlevels))
      allocate(tcenters(2,nboxes))
c
c     call the tree code
c
      call pts_tree_build(sources,ns,targ,nt,idivflag,ndiv,
     1    nlmin,iper,nlevels,nboxes,
     2    ndiv0,npwlevel,bs0,cen0,
     3    ltree,itree,iptr,tcenters,boxsize)
cccc      write(6,*) ' boxsize',boxsize(0)
      write(6,*) ' nperbox',ns*4/(3*nboxes)

      allocate(isrc(ns),isrcse(2,nboxes))
      allocate(itarg(nt),itargse(2,nboxes))

      call pts_tree_sort(ns,sources,itree,ltree,nboxes,nlevels,iptr,
     1    tcenters,isrc,isrcse)
cccc      call prinf('isrcse=*',isrcse,20)

      call pts_tree_sort(nt,targ,itree,ltree,nboxes,nlevels,iptr,
     1   tcenters,itarg,itargse)
      allocate(sourcesort(2,ns))
      allocate(targsort(2,nt))


      if(ifcharge.eq.1.and.ifdipole.eq.0) then
        allocate(chargesort(nd,ns),dipstrsort(nd,1))
        allocate(rnormalsort(2,1))
      endif
      if(ifcharge.eq.0.and.ifdipole.eq.1) then
        allocate(chargesort(nd,1),dipstrsort(nd,ns))
        allocate(rnormalsort(2,ns))
      endif
      if(ifcharge.eq.1.and.ifdipole.eq.1) then
        allocate(chargesort(nd,ns),dipstrsort(nd,ns))
        allocate(rnormalsort(2,ns))
      endif

      if(ifpgh.eq.1) then
        allocate(potsort(nd,ns),gradsort(nd,2,1),hesssort(nd,3,1))
      else if(ifpgh.eq.2) then
        allocate(potsort(nd,ns),gradsort(nd,2,ns),hesssort(nd,3,1))
      else if(ifpgh.eq.3) then
        allocate(potsort(nd,ns),gradsort(nd,2,ns),hesssort(nd,3,ns))
      else
        allocate(potsort(nd,1),gradsort(nd,2,1),hesssort(nd,3,1))
      endif
c      
      if(ifpghtarg.eq.1) then
        allocate(pottargsort(nd,nt),gradtargsort(nd,2,1),
     1     hesstargsort(nd,3,1))
      else if(ifpghtarg.eq.2) then
        allocate(pottargsort(nd,nt),gradtargsort(nd,2,nt),
     1      hesstargsort(nd,3,1))
      else if(ifpghtarg.eq.3) then
        allocate(pottargsort(nd,nt),gradtargsort(nd,2,nt),
     1     hesstargsort(nd,3,nt))
      else
        allocate(pottargsort(nd,1),gradtargsort(nd,2,1),
     1     hesstargsort(nd,3,1))
      endif
c
c     initialize potentials,hessians,gradients
c
      if(ifpgh.eq.1) then
        do i=1,ns
          do idim=1,nd
            potsort(idim,i) = 0
          enddo
        enddo
      endif

      if(ifpgh.eq.2) then
        do i=1,ns
          do idim=1,nd
            potsort(idim,i) = 0
            gradsort(idim,1,i) = 0
            gradsort(idim,2,i) = 0
          enddo
        enddo
      endif

      if(ifpgh.eq.3) then
        do i=1,ns
          do idim=1,nd
            potsort(idim,i) = 0
            gradsort(idim,1,i) = 0
            gradsort(idim,2,i) = 0
            hesssort(idim,1,i) = 0
            hesssort(idim,2,i) = 0
            hesssort(idim,3,i) = 0
          enddo
        enddo
      endif


      if(ifpghtarg.eq.1) then
        do i=1,nt
          do idim=1,nd
            pottargsort(idim,i) = 0
          enddo
        enddo
      endif

      if(ifpghtarg.eq.2) then
        do i=1,nt
          do idim=1,nd
            pottargsort(idim,i) = 0
            gradtargsort(idim,1,i) = 0
            gradtargsort(idim,2,i) = 0
          enddo
        enddo
      endif

      if(ifpghtarg.eq.3) then
        do i=1,nt
          do idim=1,nd
            pottargsort(idim,i) = 0
            gradtargsort(idim,1,i) = 0
            gradtargsort(idim,2,i) = 0
            hesstargsort(idim,1,i) = 0
            hesstargsort(idim,2,i) = 0
            hesstargsort(idim,3,i) = 0
          enddo
        enddo
      endif
c
c     compute scaling factor for multipole/local expansions
c     and lengths of multipole and local expansions
c
      allocate(ntermsh(0:nlevels))
      allocate(nlocal(0:nlevels))

      do i=0,nlevels
         nlocal(i)=0
         ntermsh(i)=0
      enddo

      ntermmax = 0
      nadd = ifdipole+max(ifpgh,ifpghtarg)
      do i = max(npwlevel,0),nlevels
         call g2dhlterms(2,boxsize(i),delta,eps,nlocal(i))
         if (ifdipole.eq.1 .or. ifpgh.ge.2 .or. ifpghtarg.ge.2) then
            nlocal(i)=nlocal(i)+nadd
         endif

         ntermsh(i) = nlocal(i)
         if (ntermmax .lt. ntermsh(i)) ntermmax = ntermsh(i)
      enddo
c
c     compute the length of plane wave expansion
      bsize = bs0/(2.0d0**(max(npwlevel,0)))
      call g2dpwterms(bsize,delta,eps,pmax,npw)
      
c     compute the length of SOE expansion
      call g2dsoeterms(eps,nsoe)      
      if (ifdipole.eq.1 .or. ifpgh.eq.3 .or. ifpghtarg.eq.3) then
         if (nsoe .le. 10) then
            nsoe = nsoe+4
         else
            nsoe = 18
         endif
      elseif (ifdipole.eq.1 .or. ifpgh.eq.2 .or. ifpghtarg.eq.2) then
         if (nsoe .le. 12) then
            nsoe = nsoe+2
         else
            nsoe = 18
         endif
      endif


      call prinf(' nlocal =*',nlocal,nlevels+1)
      call prinf(' ntermmax =*',ntermmax,1)
      call prin2(' pmax =*',pmax,1)
      call prinf(' npw =*',npw,1)
      call prinf(' nsoe =*',nsoe,1)
c       
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(3,nboxes).
c
c     iiaddr is pointer to iaddr array, itself contained in workspace.
c
c       ... allocate iaddr and temporary arrays
c
      allocate(iaddr(6,nboxes))
c
c
c     reorder sources
c
      call dreorderf(2,ns,sources,sourcesort,isrc)
      if(ifcharge.eq.1) 
     1    call dreorderf(nd,ns,charge,chargesort,isrc)
      if(ifdipole.eq.1) then
         call dreorderf(nd,ns,dipstr,dipstrsort,isrc)
         call dreorderf(2,ns,rnormal,rnormalsort,isrc)
      endif
c
cc     reorder targets
c
      call dreorderf(2,nt,targ,targsort,itarg)
c
c
c     allocate memory need by multipole, local expansions at all levels
c     
c     irmlexp is pointer for workspace need by various expansions.
c
      call g2dmpalloc(nd,itree,iaddr,
     1    nlevels,npwlevel,lmptot,
     2    ntermsh,nlocal,npw,nsoe)
      if(ifprint .eq. 1) call prinf_long(' lmptot is *',lmptot,1)
c
      allocate(rmlexp(lmptot),stat=ier)
      do i=1,lmptot
         rmlexp(i)=0
      enddo
c
c     call the main FGT routine
c

c     Memory allocation is complete. 
c     Call main FGT routine
c
      call cpu_time(time1)
C$      time1=omp_get_wtime()
      call fgt2dmain(nd,delta,eps,
     $   ns,sourcesort,
     $   ifcharge,chargesort,
     $   ifdipole,rnormalsort,dipstrsort,
     $   nt,targsort,
     $   iaddr,rmlexp,
     $   itree,ltree,iptr,nlevels,npwlevel,
     $   nboxes,iper,boxsize,tcenters,itree(iptr(1)),
     $   isrcse,itargse,ntermsh,nlocal,ntermmax,pmax,npw,nsoe,
     $   ifpgh,potsort,gradsort,hesssort,
     $   ifpghtarg,pottargsort,gradtargsort,
     $   hesstargsort)
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      if( ifprint .eq. 1 ) call prin2('time in fgt main=*',
     1   time2-time1,1)


c
c     resort the output arrays in input order
c
      if(ifpgh.eq.1) then
        call dreorderi(nd,ns,potsort,pot,isrc)
      endif

      if(ifpgh.eq.2) then
        call dreorderi(nd,ns,potsort,pot,isrc)
        call dreorderiv(nd,2,ns,gradsort,grad,isrc)
      endif

      if(ifpgh.eq.3) then
        call dreorderi(nd,ns,potsort,pot,isrc)
        call dreorderiv(nd,2,ns,gradsort,grad,isrc)
        call dreorderiv(nd,3,ns,hesssort,hess,isrc)
      endif

      if(ifpghtarg.eq.1) then
        call dreorderi(nd,nt,pottargsort,pottarg,itarg)
      endif

      if(ifpghtarg.eq.2) then
        call dreorderi(nd,nt,pottargsort,pottarg,itarg)
        call dreorderiv(nd,2,nt,gradtargsort,gradtarg,itarg)
      endif

      if(ifpghtarg.eq.3) then
        call dreorderi(nd,nt,pottargsort,pottarg,itarg)
        call dreorderiv(nd,2,nt,gradtargsort,gradtarg,itarg)
        call dreorderiv(nd,3,nt,hesstargsort,hesstarg,itarg)
      endif

      return
      end
c
c
c
c
c
      subroutine fgt2dmain(nd,delta,eps,
     $     nsource,sourcesort,
     $     ifcharge,chargesort,
     $     ifdipole,rnormalsort,dipstrsort,
     $     ntarget,targetsort,
     $     iaddr,rmlexp,
     $     itree,ltree,iptr,nlevels,npwlevel,
     $     nboxes,iper,boxsize,centers,laddr,
     $     isrcse,itargse,ntermsh,nlocal,ntermmax,pmax,npw,nsoe,
     $     ifpgh,pot,grad,hess,
     $     ifpghtarg,pottarg,gradtarg,hesstarg)
c
c
c   the FGT in R^2: evaluate all pairwise particle
c   interactions
c   and interactions with targets using SOE/X expansions.
c
c
c   \phi(x_i) = \sum_{j\ne i} charge_j e^{-|x_i-x_j|^2/delta)
c   + dipstr_j 2*(x_i - x_j)\cdot rnormal/delta * e^{-|x_i-x_j|^2/delta)
c
c
c   All the source/target/expansion center related quantities
c   are assumed to be tree-sorted
c
c-----------------------------------------------------------------------
c   INPUT PARAMETERS:
c
c   nd:   number of charge densities
c
c   eps:  FGT precision requested
c
c   nsource:     integer:  number of sources
c   sourcesort: real *8 (2,ns):  source locations
c
c   ifcharge:  charge computation flag
c              ifcharge = 1   =>  include charge contribution
c                                     otherwise do not
c   chargesort: complex *16 (nsource): charge strengths
c
c   ifdipole:  dipole computation flag
c              ifdipole = 1   =>  include dipole contribution
c                                     otherwise do not
c   dipstrsort: complex *16 (nsource): dipole strengths
c   ntarget: integer:  number of targets
c   targetsort: real *8 (2,ntarget):  target locations
c   iaddr: (6,nboxes): pointer in rmlexp where multipole
c                      and local expansions for each
c                      box is stored
c                      iaddr(1,ibox) is the
c                      starting index in rmlexp for the 
c                      multipole SOE expansion of ibox
c                      and iaddr(2,ibox) is the
c                      starting index in rmlexp
c                      for the multipole SX expansion of ibox
c                      iaddr(3,ibox) is the
c                      starting index in rmlexp for the 
c                      incoming SOE expansion of ibox
c                      iaddr(4,ibox) is the
c                      starting index in rmlexp
c                      for the incoming SX expansion of ibox
c                      iaddr(5,ibox) is the
c                      starting index in rmlexp
c                      for the local Taylor expansion of ibox
c                      iaddr(6,ibox) is the
c                      starting index in rmlexp
c                      for the Hermite expansion of ibox
c
c
c   Note: there is probably no need to store the Hermite expansion if it's not
c         used after the form stage (step 1 of the algorithm).
c
c   itree    in: integer (ltree)
c             This array contains all the information
c             about the tree
c             Refer to pts_tree2d.f
c
c   ltree    in: integer
c            length of tree
c
c    iptr in: integer(8)
c             iptr is a collection of pointers 
c             which points to where different elements 
c             of the tree are stored in the itree array
c
c     nlevels in: integer
c             number of levels in the tree
c
c     
c     npwlevel in: integer
c             cutoff level at which the X expansion is valid
c
c     
c     nboxes  in: integer
c             number of boxes in the tree
c
c     boxsize in: real*8 (0:nlevels)
c             boxsize(i) is the size of the box from end to end
c             at level i
c     iper    in: integer
c             flag for periodic implementation
c
c     centers in: real *8(2,nboxes)
c                 array containing the centers of all the boxes
c
c     isrcse in: integer(2,nboxes)
c               starting and ending location of sources in ibox
c                in sorted list of sources
c
c     itargse in: integer(2,nboxes)
c               starting and ending location of targets in ibox
c                in sorted list of sources
c
c     ntermsh in:  length of Hermite expansions
c     nlocal  in:  length of Taylor expansions
c     ntermmax in: maximum length of Hermite/Taylor expansions
c     pmax    in:  cutoff limit in the planewave expansion
c     npw     in:  length of planewave expansions
c     nsoe    in:  length of SOE expansions
c
c     ifpgh  in: integer
c             flag for evaluating potential/gradients/hessians 
c             at sources.
c             ifpgh = 1, only potentials will be evaluated
c             ifpgh = 2, potentials/gradients will be evaluated
c             ifpgh = 3, potentials/gradients/hessians will be evaluated
c
c     ifpghtarg  in: integer
c             flag for evaluating potential/gradients/hessians 
c             at targets.
c             ifpghtarg = 1, only potentials will be evaluated
c             ifpghtarg = 2, potentials/gradients will be evaluated
c             ifpghtarg = 3, potentials/gradients/hessians will be evaluated
c
c   OUTPUT
c
c   pot: potential at the source locations
c   grad: gradient at the source locations
c   hess: gradient at the source locations
c  
c   pottarg: potential at the target locations
c   gradtarg: gradient at the target locations
c   hesstarg: gradient at the target locations
c------------------------------------------------------------------

      implicit none

      integer nd
      integer iper
      integer nsource,ntarget
      integer ndiv,nlevels,npwlevel
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      real *8 eps,delta
      real *8 sourcesort(2,nsource)
      real *8 rnormalsort(2,nsource)
      real *8 chargesort(nd,*)
      real *8 dipstrsort(nd,*)
      real *8 targetsort(2,ntarget)
      real *8 pot(nd,*)
      real *8 grad(nd,2,*)
      real *8 hess(nd,3,*)
      real *8 pottarg(nd,*)
      real *8 gradtarg(nd,2,*)
      real *8 hesstarg(nd,3,*)
c
      integer iaddr(6,nboxes)
      real *8 rmlexp(*)
      real *8 pmax
      real *8 timeinfo(10)
      real *8 timelev(0:200)
      real *8 centers(2,*)
c
      integer laddr(2,0:nlevels)
      integer ntermsh(0:nlevels),nlocal(0:nlevels),ntermmax
      integer npw,nsoe,nsoehalf
      integer iptr(8),ltree
      integer itree(ltree)
      integer nboxes
      integer isrcse(2,nboxes),itargse(2,nboxes)
      real *8 boxsize(0:nlevels)
c
c     temp variables
      integer i,j,k,l,idim,ip,isx,ii,j1,j2
      integer ibox,jbox,ilev,npts,nptssrc,nptstarg
      integer nchild,ncoll
c
      integer istart,iend,istarts,iends
      integer isstart,isend,jsstart,jsend
      integer jstart,jend
      integer istarte,iende,istartt,iendt

      integer ifprint

      integer nexp,nshiftmax,nmax,kshift
      integer nsoeall,nsxall
      
      integer nn,jx,jy,nlevstart,nlevend
      real *8 d,time1,time2,omp_get_wtime
      real *8 tt1,tt2,xmin,xmin2,t1,t2
      real *8 pottmp,gradtmp(2),hesstmp(3)

      integer mnlistsoe,mnlistsx

      integer, allocatable :: nlistsoe(:,:),listsoe(:,:,:)
      integer, allocatable :: nlistsx(:,:),listsx(:,:,:)
      integer, allocatable :: ifhung(:),iflocal(:)
      
      real *8, allocatable:: wx(:),tx(:)
      real *8, allocatable :: h2l(:,:),rh2x(:,:),h2ltmp(:,:)
      
      complex *16 ws(100),ts(100)
      complex *16, allocatable :: h2s(:,:),h2x(:,:)
      complex *16, allocatable :: s2l(:,:),x2l(:,:)
      complex *16, allocatable :: wsoeshift(:,:),wsxshift(:,:)
      complex *16, allocatable :: wsoeshiftall(:,:,:),wsxshiftall(:,:,:)
      
      double precision pi
c      
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
c
      ifprint=1
      pi = 4*atan(1.0d0)


      nlevend=nlevels
      if (npwlevel.lt.nlevels) nlevend=npwlevel
      
      allocate(ifhung(nboxes))
      do i=1,nboxes
         ifhung(i)=0
      enddo

      do ilev = 0,nlevend
C
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            nchild = itree(iptr(4)+ibox-1)
            istart = isrcse(1,ibox)
            iend = isrcse(2,ibox)
            npts = iend-istart+1
c           Check if the current box is a nonempty leaf box            
            if(nchild.eq.0.and.npts.gt.0) then
               ifhung(ibox) = 1
            endif
         enddo
      enddo

      do i=1,10
         timeinfo(i)=0
      enddo

      do i=0,nlevels
         timelev(i) = 0
      enddo
      
      if (npwlevel .eq. nlevels) goto 1800




c      
      nsoehalf = nsoe/2

c     nsoe should be an even integer in [4,14] or 18
c     get soe nodes and weights
      if (nsoe .le. 12) then
         call get_soe1(nsoe,ws,ts)
      else
         call get_soe2(nsoe,ws,ts)
      endif
      
c     get planewave nodes and weights
      allocate(wx(-npw:npw))
      allocate(tx(-npw:npw))
      call get_pwnodes(pmax,npw,wx,tx)

c     compute translation matrices for SOE and SOE/X expansions
      xmin  = boxsize(nlevels)/sqrt(delta)
      xmin2 = xmin/2

      nexp = nsoe*nsoe/2
      nshiftmax = 2**(nlevels-max(npwlevel,0)+1)
      
      allocate(wsoeshift(nexp,nshiftmax**2))
      allocate(wsxshift((2*npw+1)*nsoe/2,nshiftmax))

      call shiftsoe_translation_matrices(xmin,nsoe,nshiftmax,
     1    wsoeshift,ws,ts)
      call shiftsx_translation_matrices(xmin,npw,nsoe,nshiftmax,
     1    wsxshift,ws,ts)

      nmax = 2**(nlevels-max(npwlevel,0)+1)

      nsoeall = 4*nsoe*nsoe/2
      nsxall  = 4*(2*npw+1)*nsoe/2
      
      allocate(wsoeshiftall(nsoeall,-nmax:nmax,-nmax:nmax))
      allocate(wsxshiftall(nsxall,-nmax:nmax,-nmax:nmax))
      call shiftsoeall_translation_matrices(xmin2,nsoe,nmax,
     1           wsoeshiftall,ws,ts)
      call shiftsxall_translation_matrices(xmin2,npw,nsoehalf,nmax,
     1           wsxshiftall,ts,tx)
      

c     xmin is used in shiftsoe and shiftsx subroutines to
c     determine the right translation matrices
c      
      xmin  = boxsize(nlevels)
      xmin2 = xmin/2

c
c     1d translation matrices
c
c     Hermite to SOE
      allocate(h2s(nsoe,0:ntermmax))
      call g2dh2smat(ntermmax,nsoe,ws,ts,h2s)

c     SOE to local
      allocate(s2l(nsoe,0:ntermmax))
      call g2ds2lmat(ntermmax,nsoe,ws,ts,s2l)
      
c     Hermite to planewave
      allocate(h2x(-npw:npw,0:ntermmax))
      call g2dh2xmat(ntermmax,npw,wx,tx,h2x)

c     planewave to local
      allocate(x2l(-npw:npw,0:ntermmax))
      call g2dx2lmat(ntermmax,npw,wx,tx,x2l)

c     Hermite to local
      allocate(h2l(0:ntermmax,0:ntermmax))
      call g2dh2lmat(ntermmax,ntermmax,h2l)
      
c
c     check whether we need to evaluate local expansion for boxes at the
c     cutoff level
      allocate(iflocal(nboxes))
      if (npwlevel .ge. 0 .and. npwlevel .le. nlevels) then
         do i=1,nboxes
            iflocal(i)=0
         enddo

         do ilev=npwlevel,npwlevel
            do jbox=laddr(1,ilev),laddr(2,ilev)
               nchild = itree(iptr(4)+jbox-1)
               if (nchild .eq. 0) then
                  ncoll = itree(iptr(6)+jbox-1)
                  do i=1,ncoll
                     ibox = itree(iptr(7) + (jbox-1)*9+i-1)
                     nchild = itree(iptr(4)+ibox-1)
                     if (nchild .gt. 0) then
                        iflocal(jbox)=1
                        exit
                     endif
                  enddo
               endif
            enddo
cccc            call prinf('iflocal=*',iflocal(laddr(1,ilev)),
cccc     1          laddr(2,ilev)-laddr(1,ilev)+1)
         enddo
      endif

      
c
c     compute list info
c
      call gt2d_computemnlists(nlevels,nboxes,itree,ltree,
     1    iptr,centers,
     2    boxsize,iper,mnlistsoe,mnlistsx)

      allocate(nlistsoe(0:3,nboxes),listsoe(mnlistsoe,0:3,nboxes))
      allocate(nlistsx (0:3,nboxes), listsx (mnlistsx,0:3,nboxes))

      call gt2d_computelists_slow(nlevels,npwlevel,nboxes,
     1    itree,ltree,iptr,centers,boxsize,laddr,
     2    mnlistsoe,nlistsoe,listsoe,
     3    mnlistsx, nlistsx, listsx)
      
c
c     ... set all multipole and local expansions to zero
c
cccc      do ilev = 0,nlevels
ccccC$OMP PARALLEL DO DEFAULT (SHARED)
ccccC$OMP$PRIVATE(ibox)
cccc         do ibox = laddr(1,ilev),laddr(2,ilev)
cccc           call g2dsoezero_vec(nd,rmlexp(iaddr(1,ibox)),nsoe)
cccc           call g2dsxzero_vec(nd,rmlexp(iaddr(2,ibox)),
cccc     1          npw,nsoehalf)
cccc           call g2dsoezero_vec(nd,rmlexp(iaddr(3,ibox)),nsoe)
cccc           call g2dsxzero_vec(nd,rmlexp(iaddr(4,ibox)),
cccc     1          npw,nsoehalf)
cccc           call g2dlocalzero_vec(nd,rmlexp(iaddr(5,ibox)),nlocal)
cccc           call g2dhermzero_vec(nd,rmlexp(iaddr(6,ibox)),ntermsh)
cccc
cccc         enddo
ccccC$OMP END PARALLEL DO         
cccc       enddo
c
c

      call prinf('laddr=*',laddr,2*(nlevels+1))

      
      if(ifprint .ge. 1) 
     $   call prinf('=== STEP 1 (form mp) ====*',i,0)
      call cpu_time(time1)
C$        time1=omp_get_wtime()

      
c
c       ... step 1, locate all charges, assign them to boxes, and
c       form multipole expansions


      
      do 1100 ilev = nlevels,max(npwlevel+1,0),-1

         allocate(h2ltmp(0:ntermsh(ilev),0:nlocal(ilev)))
         do j1=0,nlocal(ilev)
            do j2=0,ntermsh(ilev)
               h2ltmp(j2,j1)=h2l(j2,j1)
            enddo
         enddo
C
         if(ifcharge.eq.1.and.ifdipole.eq.0) then
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
            do ibox=laddr(1,ilev),laddr(2,ilev)
               nchild = itree(iptr(4)+ibox-1)
               istart = isrcse(1,ibox)
               iend = isrcse(2,ibox)
               npts = iend-istart+1
c              Check if current box is a leaf box            
               if(nchild.eq.0.and.npts.gt.0) then
c                 form the Hermite expansion
                  call g2dformhc_vec(nd,delta,
     1                sourcesort(1,istart),npts,chargesort(1,istart),
     2                centers(1,ibox),
     3                ntermsh(ilev),rmlexp(iaddr(6,ibox)))
c                 Hermite exp to SOE/X and SOE exp
                  call g2dh2sx_h2s_vec(nd,ntermsh(ilev),
     1                npw,nsoe,h2s,h2x,
     2                rmlexp(iaddr(6,ibox)),
     3                rmlexp(iaddr(2,ibox)),rmlexp(iaddr(1,ibox)))
cccc                call g2dformsc_vec(nd,delta,
cccc     1              sourcesort(1,istart),npts,chargesort(1,istart),
cccc     2              centers(1,ibox),
cccc     3              nsoe,ws,ts,rmlexp(iaddr(1,ibox)))
cccc                call g2dformsxc_vec(nd,delta,
cccc     1              sourcesort(1,istart),npts,chargesort(1,istart),
cccc     2              centers(1,ibox),pmax,npw,
cccc     3              nsoehalf,ws,ts,rmlexp(iaddr(2,ibox)))
                  if (ntermsh(ilev) .gt. 40) then
c                    form local expansion
                     call g2dformlc_vec(nd,delta,
     1                   sourcesort(1,istart),npts,chargesort(1,istart),
     2                   centers(1,ibox),
     3                   nlocal(ilev),rmlexp(iaddr(5,ibox)))
                  else
c                    Hermite exp to local exp
                     call g2dherm2local_vec(nd,nlocal(ilev),
     1                   ntermsh(ilev),h2ltmp,rmlexp(iaddr(6,ibox)),
     2                   rmlexp(iaddr(5,ibox)))
                  endif
               endif
            enddo
C$OMP END PARALLEL DO 
         endif

         if(ifdipole.eq.1.and.ifcharge.eq.0) then
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
            do ibox=laddr(1,ilev),laddr(2,ilev)
               nchild = itree(iptr(4)+ibox-1)
               istart = isrcse(1,ibox)
               iend = isrcse(2,ibox)
               npts = iend-istart+1
c              Check if current box is a leaf box            
               if(nchild.eq.0.and.npts.gt.0) then
                  call g2dformhd_vec(nd,delta,sourcesort(1,istart),
     1                npts,rnormalsort(1,istart),
     2                dipstrsort(1,istart),centers(1,ibox),
     3                ntermsh(ilev),rmlexp(iaddr(6,ibox)))
c                 Hermite exp to SOE/X and SOE exp
                  call g2dh2sx_h2s_vec(nd,ntermsh(ilev),npw,nsoe,
     1                h2s,h2x,rmlexp(iaddr(6,ibox)),
     2                rmlexp(iaddr(2,ibox)),rmlexp(iaddr(1,ibox)))
                  
                  if (ntermsh(ilev) .gt. 40) then
                     call g2dformld_vec(nd,delta,sourcesort(1,istart),
     1                   npts,rnormalsort(1,istart),
     2                   dipstrsort(1,istart),centers(1,ibox),
     3                   nlocal(ilev),rmlexp(iaddr(5,ibox)))
                  else
c                    Hermite expansion to local expansion
                     call g2dherm2local_vec(nd,nlocal(ilev),
     1                   ntermsh(ilev),h2ltmp,rmlexp(iaddr(6,ibox)),
     2                   rmlexp(iaddr(5,ibox)))
                  endif
               endif
            enddo
C$OMP END PARALLEL DO 
         endif

         if(ifdipole.eq.1.and.ifcharge.eq.1) then
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
            do ibox=laddr(1,ilev),laddr(2,ilev)
               nchild = itree(iptr(4)+ibox-1)
               istart = isrcse(1,ibox)
               iend = isrcse(2,ibox)
               npts = iend-istart+1
c              Check if current box is a leaf box            
               if(nchild.eq.0.and.npts.gt.0) then
                  call g2dformhcd_vec(nd,delta,sourcesort(1,istart),
     1                npts,chargesort(1,istart),rnormalsort(1,istart),
     2                dipstrsort(1,istart),centers(1,ibox),
     3                ntermsh(ilev),rmlexp(iaddr(6,ibox)))
c                 Hermite exp to SOE/X and SOE expansions
                  call g2dh2sx_h2s_vec(nd,ntermsh(ilev),npw,nsoe,
     1                h2s,h2x,rmlexp(iaddr(6,ibox)),
     2                rmlexp(iaddr(2,ibox)),rmlexp(iaddr(1,ibox)))
                  if (ntermsh(ilev) .gt. 40) then
                     call g2dformlcd_vec(nd,delta,sourcesort(1,istart),
     1                   npts,chargesort(1,istart),
     2                   rnormalsort(1,istart),dipstrsort(1,istart),
     3                   centers(1,ibox),nlocal(ilev),
     4                   rmlexp(iaddr(5,ibox)))
                  else
c                    Hermite exp to local expansion
                     call g2dherm2local_vec(nd,nlocal(ilev),
     1                   ntermsh(ilev),h2ltmp,
     1                   rmlexp(iaddr(6,ibox)),rmlexp(iaddr(5,ibox)))
                  endif
               endif
            enddo
C     $OMP END PARALLEL DO
         endif
         deallocate(h2ltmp)
c     end of ilev do loop
 1100 continue


      
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(1)=time2-time1


      
      if(ifprint .ge. 1)
     $      call prinf('=== STEP 2 (merge mp) ====*',i,0)
      call cpu_time(time1)
C$    time1=omp_get_wtime()
c

      
      do 1200 ilev=nlevels-1,max(npwlevel,0),-1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,i,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          nchild = itree(iptr(4)+ibox-1)
          do i=1,nchild
            jbox = itree(iptr(5)+4*(ibox-1)+i-1)
            istart = isrcse(1,jbox)
            iend = isrcse(2,jbox)
            npts = iend-istart+1
            if(npts.gt.0) then
              jx= 1.05d0*(centers(1,ibox) - centers(1,jbox))/xmin2
              jy= 1.05d0*(centers(2,ibox) - centers(2,jbox))/xmin2
              call g2dshiftall_vec(nd,nsoeall,rmlexp(iaddr(1,jbox)),
     1             rmlexp(iaddr(1,ibox)),wsoeshiftall(1,jx,jy))
              call g2dshiftall_vec(nd,nsxall,rmlexp(iaddr(2,jbox)),
     1            rmlexp(iaddr(2,ibox)),wsxshiftall(1,jx,jy))
            endif
          enddo
        enddo
C$OMP END PARALLEL DO    
 1200 continue

      
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(2)=time2-time1
c

      
      if(ifprint.ge.1)
     $    call prinf('=== Step 3 (mp to loc) ===*',i,0)
c      ... step 3, convert multipole expansions into local
c       expansions

      call cpu_time(time1)
C$    time1=omp_get_wtime()


      
      do 1300 ilev = max(npwlevel,0),nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,istart,iend,npts,i)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          npts = 0
          if(ifpghtarg.gt.0) then
            istart = itargse(1,ibox)
            iend = itargse(2,ibox)
            npts = npts + iend-istart+1
          endif

          if(ifpgh.gt.0) then
            istart = isrcse(1,ibox)
            iend = isrcse(2,ibox)
            npts = npts + iend-istart+1
          endif
c
          if(npts.gt.0) then
c
c           shift soe expansions
c
            do ip=0,3
            do j=1,nlistsoe(ip,ibox)
              jbox=listsoe(j,ip,ibox)
              jx= 1.05d0*abs(centers(1,jbox) - centers(1,ibox))/xmin
              jy= 1.05d0*abs(centers(2,jbox) - centers(2,ibox))/xmin
              kshift=(jx-1)*nshiftmax+jy
              call g2dshiftsoe_vec(nd,nsoe,rmlexp(iaddr(1,jbox)),
     1            rmlexp(iaddr(3,ibox)),wsoeshift(1,kshift),ip)
            enddo
            enddo
c
c           shift soe/x hybrid expansions
c
            do isx=0,3
            do j=1,nlistsx(isx,ibox)
              jbox=listsx(j,isx,ibox)
              if (isx .lt. 2) then
                 kshift=1.05*abs(centers(1,jbox)-centers(1,ibox))/xmin
              else
                 kshift=1.05*abs(centers(2,jbox)-centers(2,ibox))/xmin
              endif
              call g2dshiftsx_vec(nd,npw,nsoehalf,
     1            rmlexp(iaddr(2,jbox)),rmlexp(iaddr(4,ibox)),
     2            wsxshift(1,kshift),isx)
             enddo
             enddo
          endif
        enddo
C$OMP END PARALLEL DO        
 1300 continue
c

      
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(3) = time2-time1

      call prin2('timeinfo4=*',time2-time1,1)

      
      if(ifprint.ge.1)
     $    call prinf('=== Step 4 (split loc) ===*',i,0)

      call cpu_time(time1)
C$    time1=omp_get_wtime()


      do 1400 ilev = max(npwlevel,0),nlevels-1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,i,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          nchild = itree(iptr(4)+ibox-1)
          npts = 0
          if(ifpghtarg.gt.0) then
            istart = itargse(1,ibox) 
            iend = itargse(2,ibox) 
            npts = npts + iend-istart+1
          endif

          if(ifpgh.gt.0) then
            istart = isrcse(1,ibox) 
            iend = isrcse(2,ibox) 
            npts = npts + iend-istart+1
          endif

          if(npts.gt.0) then
            do i=1,nchild
              jbox = itree(iptr(5)+4*(ibox-1)+i-1)
              jx= 1.05d0*(centers(1,jbox) - centers(1,ibox))/xmin2
              jy= 1.05d0*(centers(2,jbox) - centers(2,ibox))/xmin2
              call g2dshiftall_vec(nd,nsoeall,rmlexp(iaddr(3,ibox)),
     1            rmlexp(iaddr(3,jbox)),wsoeshiftall(1,jx,jy))
              call g2dshiftall_vec(nd,nsxall,rmlexp(iaddr(4,ibox)),
     1            rmlexp(iaddr(4,jbox)),wsxshiftall(1,jx,jy))
            enddo
          endif
        enddo
C$OMP END PARALLEL DO        
 1400 continue

      
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(4) = time2-time1


      
      if(ifprint.ge.1)
     $    call prinf('=== step 5 (eval loc) ===*',i,0)

c     ... step 5, evaluate all local expansions
      call cpu_time(time1)
C$    time1=omp_get_wtime()


      
      do 1500 ilev = max(npwlevel,0),nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istart,iend,i,npts)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            if (ilev .eq. npwlevel .and. iflocal(ibox).eq.0) then
c              do nothing here
            else    
               nchild = itree(iptr(4)+ibox-1)
               if(nchild.eq.0) then
                  istartt = itargse(1,ibox) 
                  iendt = itargse(2,ibox)
                  nptstarg = iendt-istartt + 1

                  istarts = isrcse(1,ibox)
                  iends = isrcse(2,ibox)
                  nptssrc = iends-istarts+1

                  if (nptssrc + nptstarg .gt. 0) then
c                    convert SOE expansions to local expansions
                     call g2dsoe2local_vec(nd,nlocal(ilev),nsoe,s2l,
     1                   rmlexp(iaddr(3,ibox)),rmlexp(iaddr(5,ibox)))
c                    convert SOE/X expansions to local expansions
                     call g2dsx2local_vec(nd,nlocal(ilev),npw,nsoe,
     1                   x2l,s2l,rmlexp(iaddr(4,ibox)),
     2                   rmlexp(iaddr(5,ibox)))

c                    evaluate local expansion at targets
                     if(ifpghtarg.eq.1) then
                        call g2dlevalp_vec(nd,delta,centers(1,ibox),
     1                      nlocal(ilev),rmlexp(iaddr(5,ibox)),
     2                      targetsort(1,istartt),
     3                      nptstarg,pottarg(1,istartt))
                     endif
                     if(ifpghtarg.eq.2) then
                        call g2dlevalg_vec(nd,delta,centers(1,ibox),
     1                      nlocal(ilev),rmlexp(iaddr(5,ibox)),
     2                      targetsort(1,istartt),nptstarg,
     3                      pottarg(1,istartt),gradtarg(1,1,istartt))
                     endif
                     if(ifpghtarg.eq.3) then
                        call g2dlevalh_vec(nd,delta,centers(1,ibox),
     1                      nlocal(ilev),rmlexp(iaddr(5,ibox)),
     2                      targetsort(1,istartt),nptstarg,
     3                      pottarg(1,istartt),gradtarg(1,1,istartt),
     4                      hesstarg(1,1,istartt))
                     endif
c     
c                    evaluate local expansion at sources
                     if(ifpgh.eq.1) then
                        call g2dlevalp_vec(nd,delta,centers(1,ibox),
     1                      nlocal(ilev),rmlexp(iaddr(5,ibox)),
     2                      sourcesort(1,istarts),nptssrc,
     3                      pot(1,istarts))
                     endif
                     if(ifpgh.eq.2) then
                        call g2dlevalg_vec(nd,delta,centers(1,ibox),
     1                      nlocal(ilev),rmlexp(iaddr(5,ibox)),
     2                      sourcesort(1,istarts),nptssrc,
     3                      pot(1,istarts),grad(1,1,istarts))
                     endif
                     if(ifpgh.eq.3) then
                        call g2dlevalh_vec(nd,delta,centers(1,ibox),
     1                      nlocal(ilev),rmlexp(iaddr(5,ibox)),
     2                      sourcesort(1,istarts),nptssrc,
     3                      pot(1,istarts),grad(1,1,istarts),
     4                      hess(1,1,istarts))
                     endif
                  endif
               endif
            endif
         enddo
C$OMP END PARALLEL DO        
 1500 continue

      
      call cpu_time(time2)
C$    time2 = omp_get_wtime()      
      timeinfo(5) = time2 - time1
      
      
 1800 continue

      if(ifprint .ge. 1)
     $     call prinf('=== STEP 6 (direct) =====*',i,0)
c
cc
      call cpu_time(time1)
C$    time1=omp_get_wtime()  
      do 2000 ilev = 0,nlevend
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,istartt,iendt,i,jstart,jend,istarte,iende)
C$OMP$PRIVATE(istarts,iends,nptssrc,nptstarg)
C$OMP$SCHEDULE(DYNAMIC)  
         do jbox = laddr(1,ilev),laddr(2,ilev)
            if (ifhung(jbox) .eq. 1) then
               
               jstart = isrcse(1,jbox)
               jend = isrcse(2,jbox)
               

               ncoll = itree(iptr(6)+jbox-1)
               do i=1,ncoll
                  ibox = itree(iptr(7) + (jbox-1)*9+i-1)

                  istarts = isrcse(1,ibox)
                  iends = isrcse(2,ibox)
                  nptssrc = iends-istarts + 1
                  
                  istartt = itargse(1,ibox)
                  iendt = itargse(2,ibox)
                  nptstarg = iendt-istartt + 1

                  if (nptstarg .gt. 0) then
                     call fgt2dpart_direct_vec(nd,delta,eps/90,
     1                   jstart,jend,istartt,iendt,sourcesort,
     2                   ifcharge,chargesort,
     3                   ifdipole,rnormalsort,dipstrsort,targetsort,
     4                   ifpghtarg,pottarg,gradtarg,hesstarg)
                  endif
                  if (nptssrc .gt. 0) then
                     call fgt2dpart_direct_vec(nd,delta,eps/90,
     1                   jstart,jend,istarts,iends,sourcesort,
     2                   ifcharge,chargesort,
     3                   ifdipole,rnormalsort,dipstrsort,sourcesort,
     4                   ifpgh,pot,grad,hess)
                  endif
               enddo
            endif
         enddo
C$OMP END PARALLEL DO         
 2000 continue

      call cpu_time(time2)
C$    time2=omp_get_wtime()  
      timeinfo(6) = time2-time1
      if(ifprint.ge.1) call prin2('timeinfo=*',timeinfo,6)
      d = 0
      do i = 1,6
         d = d + timeinfo(i)
      enddo

      if(ifprint.ge.1) call prin2('sum(timeinfo)=*',d,1)
cccc      if(ifprint.ge.1) call prin2('timlev=*',timelev,nlevels+1)

      return
      end
c
c
c
c
c------------------------------------------------------------------     
      subroutine fgt2dpart_direct_vec(nd,delta,eps,istart,iend,
     $    jstart,jend,source,ifcharge,charge,
     2    ifdipole,rnormal,dipstr,
     $    targ,ifpgh,pot,grad,hess)
c--------------------------------------------------------------------
c     This subroutine adds the contribution due to sources
c     istart to iend in the source array to the fields at targets
c     jstart to jend in the target array.
c
c     INPUT arguments
c-------------------------------------------------------------------
c     nd           in: integer
c                  number of charge densities
c
c     delta        in: Gaussian variance
c
c     istart       in:Integer
c                  Starting index in source array whose expansions
c                  we wish to add
c
c     iend         in:Integer
c                  Last index in source array whose expansions
c                  we wish to add
c
c     jstart       in: Integer
c                  First index in target array at which we
c                  wish to update the potential and gradients
c 
c     jend         in:Integer
c                  Last index in target array at which we wish
c                  to update the potential and gradients
c
c     source       in: real *8(2,ns)
c                  Source locations
c
c     ifcharge     in: Integer
c                  flag for including expansions due to charges
c                  The expansion due to charges will be included
c                  if ifcharge == 1
c
c     charge       in: complex *16
c                  Charge at the source locations
c
c     ifdipole     in: Integer
c                 flag for including expansions due to dipoles
c                 The expansion due to dipoles will be included
c                 if ifdipole == 1
c
c     dipstr        in: complex *16(ns)
c                 dipole strengths at the source locations
c
c     targ        in: real *8(2,nt)
c                 target locations
c
c     ifpgh        in: Integer
c                  Flag for computing the potential/gradient/hessian.
c                  ifpgh = 1, only potential is computed
c                  ifpgh = 2, potential/gradient are computed
c                  ifpgh = 3, potential/gradient/hessian are computed
c
c------------------------------------------------------------
c     OUTPUT
c
c     pot          potential incremented at targets
c     grad         gradients incremented at targets
c     hess         Hessians  incremented at targets
c-------------------------------------------------------               
        implicit none
c
        integer istart,iend,jstart,jend,ns,j,i
        integer ifcharge,ifdipole
        integer nd
        integer ifpgh
c
        real *8 source(2,*)
        real *8 rnormal(2,*)
        real *8 charge(nd,*),dipstr(nd,*)
        real *8 targ(2,*),delta,eps
        real *8 pot(nd,*)
        real *8 grad(nd,2,*)
        real *8 hess(nd,3,*)
c
        ns = iend - istart + 1
        if(ifcharge.eq.1.and.ifdipole.eq.0) then
          if(ifpgh.eq.1) then
             do j=jstart,jend
               call g2d_directcp_vec(nd,delta,eps,source(1,istart),ns,
     1            charge(1,istart),targ(1,j),pot(1,j))
             enddo
          endif

          if(ifpgh.eq.2) then
             do j=jstart,jend
               call g2d_directcg_vec(nd,delta,eps,source(1,istart),ns,
     1            charge(1,istart),targ(1,j),pot(1,j),grad(1,1,j))
             enddo
          endif
          if(ifpgh.eq.3) then
             do j=jstart,jend
               call g2d_directch_vec(nd,delta,eps,source(1,istart),ns,
     1            charge(1,istart),targ(1,j),pot(1,j),grad(1,1,j),
     2            hess(1,1,j))
             enddo
          endif
        endif

        if(ifcharge.eq.0.and.ifdipole.eq.1) then
          if(ifpgh.eq.1) then
             do j=jstart,jend
               call g2d_directdp_vec(nd,delta,eps,source(1,istart),ns,
     1          rnormal(1,istart),dipstr(1,istart),targ(1,j),pot(1,j))
             enddo
          endif

          if(ifpgh.eq.2) then
             do j=jstart,jend
               call g2d_directdg_vec(nd,delta,eps,source(1,istart),ns,
     1            rnormal(1,istart),dipstr(1,istart),targ(1,j),
     2            pot(1,j),grad(1,1,j))
             enddo
          endif
          if(ifpgh.eq.3) then
             do j=jstart,jend
               call g2d_directdh_vec(nd,delta,eps,source(1,istart),ns,
     1            rnormal(1,istart),dipstr(1,istart),targ(1,j),
     2            pot(1,j),grad(1,1,j),hess(1,1,j))
             enddo
          endif
        endif

        if(ifcharge.eq.1.and.ifdipole.eq.1) then
          if(ifpgh.eq.1) then
             do j=jstart,jend
               call g2d_directcdp_vec(nd,delta,eps,source(1,istart),ns,
     1            charge(1,istart),rnormal(1,istart),dipstr(1,istart),
     2            targ(1,j),pot(1,j))
             enddo
          endif

          if(ifpgh.eq.2) then
             do j=jstart,jend
               call g2d_directcdg_vec(nd,delta,eps,source(1,istart),ns,
     1            charge(1,istart),rnormal(1,istart),dipstr(1,istart),
     2            targ(1,j),pot(1,j),grad(1,1,j))
             enddo
          endif
          if(ifpgh.eq.3) then
             do j=jstart,jend
               call g2d_directcdh_vec(nd,delta,eps,source(1,istart),ns,
     1            charge(1,istart),rnormal(1,istart),dipstr(1,istart),
     2            targ(1,j),pot(1,j),grad(1,1,j),hess(1,1,j))
             enddo
          endif
        endif


c
        return
        end
c
c
c
c
c------------------------------------------------------------------    
      subroutine g2dmpalloc(nd,laddr,iaddr,
     1    nlevels,npwlevel,lmptot,
     2    ntermsh,nlocal,npw,nsoe)
c     This subroutine determines the size of the array
c     to be allocated for multipole/local expansions
c
c     Input arguments
c     nd          in: integer
c                 number of expansions
c
c     laddr       in: Integer(2,0:nlevels)
c                 indexing array providing access to boxes at each
c                 level
c
c     nlevels     in: Integer
c                 Total numner of levels
c     
c     npwlevel    in: Integer
c                 cutoff level where the plane wave expansion is
c                 valid at or below ilev = npwlevel
c
c     ntermsh     in: Integer(0:nlevels)
c                 Number of terms requried in the Hermite expansions at each
c                 level
c     nlocal      in: Integer(0:nlevels)
c                 Number of terms requried in the local expansions at each
c                 level
c     npw         in: Integer
c                 Number of terms in the plane wave expansion
c     nsoe        in: Integer
c                 Number of terms in the SOE expansion
c
c
c
c
c
c------------------------------------------------------------------
c     Output arguments
c     iaddr: (6,nboxes): pointer in rmlexp where multipole
c                      and local expansions for each
c                      box is stored
c                      iaddr(1,ibox) is the
c                      starting index in rmlexp for the 
c                      multipole SOE expansion of ibox
c                      and iaddr(2,ibox) is the
c                      starting index in rmlexp
c                      for the multipole SX expansion of ibox
c                      iaddr(3,ibox) is the
c                      starting index in rmlexp for the 
c                      incoming SOE expansion of ibox
c                      iaddr(4,ibox) is the
c                      starting index in rmlexp
c                      for the incoming SX expansion of ibox
c                      iaddr(5,ibox) is the
c                      starting index in rmlexp
c                      for the local Taylor expansion of ibox
c                      iaddr(6,ibox) is the
c                      starting index in rmlexp
c                      for the Hermite expansion of ibox
c     lmptot      out: Integer
c                 Total length of expansions array required
c------------------------------------------------------------------

      implicit none
      integer nlevels,npwlevel,npw,nsoe,nd
      integer ntermsh(0:nlevels),nlocal(0:nlevels)
      integer iaddr(6,*), laddr(2,0:nlevels)
      integer *8 lmptot
      integer ibox,i,istart,nn,itmp,nlevstart
c
      istart = 1
      if (npwlevel .eq. nlevels) return
      
      nlevstart = 0
      if (npwlevel .ge. 0) nlevstart = npwlevel

      do i = nlevstart,nlevels

         nn = (4*nsoe*nsoe/2)*2*nd
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the multipole SOE expansion         
c
           itmp = ibox - laddr(1,i)
           iaddr(1,ibox) = istart + itmp*nn
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
      enddo
c
      do i = nlevstart,nlevels

         nn = (4*(2*npw+1)*nsoe/2)*2*nd
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the multipole SX expansion         
c
           itmp = ibox - laddr(1,i)
           iaddr(2,ibox) = istart + itmp*nn
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
      enddo


      
      
c     
c     Allocate memory for the local SOE and SOE/X expansions
c
      do i=nlevstart,nlevels
         nn = (4*nsoe*nsoe/2)*2*nd
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
             itmp = ibox - laddr(1,i)
             iaddr(3,ibox) = istart + itmp*nn 
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
      enddo
c
       do i=nlevstart,nlevels
         nn = (4*(2*npw+1)*nsoe/2)*2*nd
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
             itmp = ibox - laddr(1,i)
             iaddr(4,ibox) = istart + itmp*nn 
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
      enddo




      
c
c     Allocate memory for the local Taylor expansions
c
      do i=nlevstart,nlevels
         nn = (nlocal(i)+1)*(nlocal(i)+1)*nd
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
             itmp = ibox - laddr(1,i)
             iaddr(5,ibox) = istart + itmp*nn 
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
      enddo



      
c
c     Allocate memory for the Hermite expansions
c
      nlevstart = 0
      if (npwlevel .ge. 0) nlevstart = npwlevel+1

      do i=nlevstart,nlevels
         nn = (ntermsh(i)+1)*(ntermsh(i)+1)*nd
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
             itmp = ibox - laddr(1,i)
             iaddr(6,ibox) = istart + itmp*nn 
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
      enddo
      
      lmptot = istart

      return
      end
c----------------------------------------------------------------     
c
c
c
c
      subroutine fgt2d_large_delta(nd,delta,eps,ns,sources,
     1    ifcharge,charge,ifdipole,rnormal,dipstr,iper,
     2    nlocal,center,
     3    ifpgh,pot,grad,hess,nt,targ,
     4    ifpghtarg,pottarg,gradtarg,hesstarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of FGTs (same source and target locations, 
c                   different charge, dipole strengths)
c   delta         : Gaussian variance 
c   eps           : precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   ifcharge      : flag for including charge interactions
c                   charge interactions included if ifcharge =1
c                   not included otherwise
c   charge(nd,ns) : charge strengths
c   ifdipole      : flag for including dipole interactions
c                   dipole interactions included if ifcharge =1
c                   not included otherwise
c   rnormal(2,ns) : dipole directions
c   dipstr(nd,ns) : dipole strengths
c   iper          : flag for periodic implmentations. Currently unused
c   nlocal        : number of terms in the local expansion
c   center(2)     : the center of the bounding box
c   ifpgh         : flag for computing pot/grad/hess
c                   ifpgh = 1, only potential is computed
c                   ifpgh = 2, potential and gradient are computed
c                   ifpgh = 3, potential, gradient, and hessian 
c                   are computed
c   nt            : number of targets
c   targ(2,nt)    : target locations
c   ifpghtarg     : flag for computing pottarg/gradtarg/hesstarg
c                   ifpghtarg = 1, only potential is computed at targets
c                   ifpghtarg = 2, potential and gradient are 
c                   computed at targets
c                   ifpghtarg = 3, potential, gradient, and hessian are 
c                   computed at targets
c
c   OUTPUT PARAMETERS
c   pot(nd,*)       : potential at the source locations
c   grad(nd,2,*)    : gradients at the source locations
c   hess(nd,3,*)    : hessian at the source locations
c   pottarg(nd,*)   : potential at the target locations
c   gradtarg(nd,2,*): gradient at the target locations
c   hesstarg(nd,3,*): hessian at the target locations
c
      implicit none
c
cc      calling sequence variables
c 
      integer nd,ns,nt
      integer iper,ifcharge,ifdipole
      integer ifpgh,ifpghtarg,nlocal

      real *8 omp_get_wtime
      real *8 eps,delta,bgf,bsize
      real *8 sources(2,ns),targ(2,nt)
      real *8 rnormal(2,ns)
      real *8 charge(nd,*),dipstr(nd,*)
      real *8 center(2)

      real *8 pot(nd,*),grad(nd,2,*),hess(nd,3,*)
      real *8 pottarg(nd,*),gradtarg(nd,2,*),hesstarg(nd,3,*)


      real *8, allocatable :: local(:,:,:)

      allocate(local(0:nlocal,0:nlocal,nd))
c
c     form local
c
      if (ifcharge.eq.1 .and. ifdipole.eq.0) then
         call g2dformlc_vec(nd,delta,sources,ns,charge,center,
     1       nlocal,local)
      elseif (ifcharge.eq.0 .and. ifdipole.eq.1) then
         call  g2dformld_vec(nd,delta,sources,ns,rnormal,dipstr,
     1       center,nlocal,local)
      elseif (ifcharge.eq.1 .and. ifdipole.eq.1) then
         call g2dformlcd_vec(nd,delta,sources,ns,charge,rnormal,
     1       dipstr,center,nlocal,local)
      endif

c
c     local eval
c
c     targets
c      
      if(ifpghtarg.eq.1) then
         call g2dlevalp_vec(nd,delta,center,nlocal,local,
     1       targ,nt,pottarg)
      endif
      if(ifpghtarg.eq.2) then
         call g2dlevalg_vec(nd,delta,center,nlocal,local,
     1       targ,nt,pottarg,gradtarg)
      endif
      if(ifpghtarg.eq.3) then
         call g2dlevalh_vec(nd,delta,center,nlocal,local,
     1       targ,nt,pottarg,gradtarg,hesstarg)
      endif
c
c     sources
c      
      if(ifpgh.eq.1) then
         call g2dlevalp_vec(nd,delta,center,nlocal,local,
     1       sources,ns,pot)
      endif
      if(ifpgh.eq.2) then
         call g2dlevalg_vec(nd,delta,center,nlocal,local,
     1       sources,ns,pot,grad)
      endif
      if(ifpgh.eq.3) then
         call g2dlevalh_vec(nd,delta,center,nlocal,local,
     1       sources,ns,pot,grad,hess)
      endif
      
      return
      end
c
c
c
c
