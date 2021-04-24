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

      subroutine fgt3d(nd,delta,eps,ns,sources,ifcharge,charge,
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
c   sources(3,ns) : source locations
c   ifcharge      : flag for including charge interactions
c                   charge interactions included if ifcharge =1
c                   not included otherwise
c   charge(nd,ns) : charge strengths
c   ifdipole      : flag for including dipole interactions
c                   dipole interactions included if ifcharge =1
c                   not included otherwise
c   rnormal(3,ns) : dipole directions
c   dipstr(nd,ns) : dipole strengths
c   iper          : flag for periodic implmentations. Currently unused
c   ifpgh         : flag for computing pot/grad/hess
c                   ifpgh = 1, only potential is computed
c                   ifpgh = 2, potential and gradient are computed
c                   ifpgh = 3, potential, gradient, and hessian 
c                   are computed
c   nt            : number of targets
c   targ(3,nt)    : target locations
c   ifpghtarg     : flag for computing pottarg/gradtarg/hesstarg
c                   ifpghtarg = 1, only potential is computed at targets
c                   ifpghtarg = 2, potential and gradient are 
c                   computed at targets
c                   ifpghtarg = 3, potential, gradient, and hessian are 
c                   computed at targets
c
c   OUTPUT PARAMETERS
c   pot(nd,*)       : potential at the source locations
c   grad(nd,3,*)    : gradients at the source locations
c   hess(nd,6,*)    : hessian at the source locations
c   pottarg(nd,*)   : potential at the target locations
c   gradtarg(nd,3,*): gradient at the target locations
c   hesstarg(nd,6,*): hessian at the target locations
c
c   Note: hessians are in the order xx, yy, zz, xy, xz, yz
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
      real *8 sources(3,ns),targ(3,nt)
      real *8 rnormal(3,ns)
      real *8 charge(nd,*),dipstr(nd,*)

      real *8 pot(nd,*),grad(nd,3,*),hess(nd,6,*)
      real *8 pottarg(nd,*),gradtarg(nd,3,*),hesstarg(nd,6,*)

c
cc      Tree variables
c
      integer, allocatable :: itree(:)
      integer *8 iptr(8)
      integer iper,nlmin,nlmax,ifunif
      real *8, allocatable :: tcenters(:,:),boxsize(:)
      integer idivflag,nlevels,nboxes,ndiv,ndiv0,npwlevel
      integer *8 ltree

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
      integer i,ilev,lmptmp,idim,ndim
      integer nlocal0,npw,nadd,ifprint,ier,nlevstart
      real *8 omp_get_wtime
      real *8 time1,time2,pi,done,pmax,bs0,cen0(3),bsize,pweps

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
      ndim = 3
      if (npwlevel .lt. -2) then
         call g3dhlterms(ndim,bs0,delta,eps,nlocal0)
         nlocal0 = nlocal0 + max(ifpgh,ifpghtarg) -1 
         write(6,*) ' nlocal0',nlocal0
      
         call cpu_time(time1)
C$      time1=omp_get_wtime()
         call fgt3d_large_delta(nd,delta,eps,ns,sources,
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
      ndiv0 = 20
c     ndiv is the maximum number of points per box at or below the cutoff level
c     it's determined by numerical experiments on finding the crossover point
c     between direct evaluation and the fast scheme.
      ndiv = 100
c
      ifunif = 0
      iper = 0

c
c     call the tree memory management
c     code to determine number of boxes,
c     number of levels and length of tree
c
c     1. determine the maximum level - it seems that the best strategy for point FGT
c     is to stop refinement as long as the Hermite expansion length is less
c     than 20 for high precision calculation
      nlmax = npwlevel + 4
      nlevstart=max(npwlevel,0)
      
      bsize = bs0/2**nlevstart
      do i = nlevstart,100
         call g3dhlterms(ndim,bsize,delta,eps,nlocal0)
         if (nlocal0 .le. 13) then
            nlmax = i
            exit
         endif
         bsize = bsize/2
      enddo

cccc      npwlevel=npwlevel-1
      if (npwlevel .ge. 0) nlmax=npwlevel+1

      write(6,*) ' nlmax',nlmax

      
      nlmin = 0
      call pts_tree_mem(sources,ns,targ,nt,idivflag,
     1    ndiv,nlmin,nlmax,ifunif,iper,
     2    ndiv0,npwlevel,bs0,cen0,
     3    nlevels,nboxes,ltree) 

c 
      write(6,*) ' nlevels',nlevels
      write(6,*) ' nboxes',nboxes

      allocate(itree(ltree))
      allocate(boxsize(0:nlevels))
      allocate(tcenters(3,nboxes))
c
c     call the tree code
c
      call pts_tree_build(sources,ns,targ,nt,idivflag,ndiv,
     1    nlmin,nlmax,ifunif,iper,nlevels,nboxes,
     2    ndiv0,npwlevel,bs0,cen0,
     3    ltree,itree,iptr,tcenters,boxsize)
cccc      write(6,*) ' boxsize',boxsize(0)
cccc      write(6,*) ' cutoff length', boxsize(0)/2**npwlevel/sqrt(delta)
cccc      write(6,*) ' nperbox',ns*4/(3*nboxes)

      allocate(isrc(ns),isrcse(2,nboxes))
      allocate(itarg(nt),itargse(2,nboxes))

      call pts_tree_sort(ns,sources,itree,ltree,nboxes,nlevels,iptr,
     1    tcenters,isrc,isrcse)
cccc      call prinf('isrcse=*',isrcse,20)

      call pts_tree_sort(nt,targ,itree,ltree,nboxes,nlevels,iptr,
     1   tcenters,itarg,itargse)
      allocate(sourcesort(ndim,ns))
      allocate(targsort(ndim,nt))


      if(ifcharge.eq.1.and.ifdipole.eq.0) then
        allocate(chargesort(nd,ns),dipstrsort(nd,1))
        allocate(rnormalsort(ndim,1))
      endif
      if(ifcharge.eq.0.and.ifdipole.eq.1) then
        allocate(chargesort(nd,1),dipstrsort(nd,ns))
        allocate(rnormalsort(ndim,ns))
      endif
      if(ifcharge.eq.1.and.ifdipole.eq.1) then
        allocate(chargesort(nd,ns),dipstrsort(nd,ns))
        allocate(rnormalsort(ndim,ns))
      endif

      if(ifpgh.eq.1) then
        allocate(potsort(nd,ns),gradsort(nd,3,1),hesssort(nd,6,1))
      else if(ifpgh.eq.2) then
        allocate(potsort(nd,ns),gradsort(nd,3,ns),hesssort(nd,6,1))
      else if(ifpgh.eq.3) then
        allocate(potsort(nd,ns),gradsort(nd,3,ns),hesssort(nd,6,ns))
      else
        allocate(potsort(nd,1),gradsort(nd,3,1),hesssort(nd,6,1))
      endif
c      
      if(ifpghtarg.eq.1) then
        allocate(pottargsort(nd,nt),gradtargsort(nd,3,1),
     1     hesstargsort(nd,6,1))
      else if(ifpghtarg.eq.2) then
        allocate(pottargsort(nd,nt),gradtargsort(nd,3,nt),
     1      hesstargsort(nd,6,1))
      else if(ifpghtarg.eq.3) then
        allocate(pottargsort(nd,nt),gradtargsort(nd,3,nt),
     1     hesstargsort(nd,6,nt))
      else
        allocate(pottargsort(nd,1),gradtargsort(nd,3,1),
     1     hesstargsort(nd,6,1))
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
            gradsort(idim,3,i) = 0
          enddo
        enddo
      endif

      if(ifpgh.eq.3) then
        do i=1,ns
          do idim=1,nd
            potsort(idim,i) = 0
            gradsort(idim,1,i) = 0
            gradsort(idim,2,i) = 0
            gradsort(idim,3,i) = 0
            hesssort(idim,1,i) = 0
            hesssort(idim,2,i) = 0
            hesssort(idim,3,i) = 0
            hesssort(idim,4,i) = 0
            hesssort(idim,5,i) = 0
            hesssort(idim,6,i) = 0
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
            gradtargsort(idim,3,i) = 0
          enddo
        enddo
      endif

      if(ifpghtarg.eq.3) then
        do i=1,nt
          do idim=1,nd
            pottargsort(idim,i) = 0
            gradtargsort(idim,1,i) = 0
            gradtargsort(idim,2,i) = 0
            gradtargsort(idim,3,i) = 0
            hesstargsort(idim,1,i) = 0
            hesstargsort(idim,2,i) = 0
            hesstargsort(idim,3,i) = 0
            hesstargsort(idim,4,i) = 0
            hesstargsort(idim,5,i) = 0
            hesstargsort(idim,6,i) = 0
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
         call g3dhlterms(ndim,boxsize(i),delta,eps,nlocal(i))
         if (ifdipole.eq.1 .or. ifpgh.ge.2 .or. ifpghtarg.ge.2) then
            nlocal(i)=nlocal(i)+nadd
cccc            nlocal(i)=nlocal(i)*1.1
         endif

         ntermsh(i) = nlocal(i)
         if (ntermmax .lt. ntermsh(i)) ntermmax = ntermsh(i)
      enddo
c
c     compute the length of plane wave expansion
      
      bsize = 2*bs0/(2.0d0**(max(npwlevel,0)))
      pweps = eps
      if (nadd .gt. 2) pweps=pweps/10
      call g3dpwterms(bsize,delta,pweps,pmax,npw)
      
      call prinf(' nlocal =*',nlocal,nlevels+1)
cccc      call prinf(' ntermmax =*',ntermmax,1)
cccc      call prin2(' pmax =*',pmax,1)
      call prinf(' npw =*',npw,1)
c       
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(3,nboxes).
c
c     iiaddr is pointer to iaddr array, itself contained in workspace.
c
c       ... allocate iaddr and temporary arrays
c
      allocate(iaddr(2,nboxes))
c
c
c     reorder sources
c
      call dreorderf(3,ns,sources,sourcesort,isrc)
      if(ifcharge.eq.1) 
     1    call dreorderf(nd,ns,charge,chargesort,isrc)
      if(ifdipole.eq.1) then
         call dreorderf(nd,ns,dipstr,dipstrsort,isrc)
         call dreorderf(3,ns,rnormal,rnormalsort,isrc)
      endif
c
cc     reorder targets
c
      call dreorderf(3,nt,targ,targsort,itarg)
c
c
c     allocate memory need by multipole, local expansions at all levels
c     
c     irmlexp is pointer for workspace need by various expansions.
c
      call g3dmpalloc(nd,itree,iaddr,
     1    nlevels,npwlevel,lmptot,
     2    ntermsh,nlocal,npw)
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
      call fgt3dmain(nd,delta,eps,
     $   ns,sourcesort,
     $   ifcharge,chargesort,
     $   ifdipole,rnormalsort,dipstrsort,
     $   nt,targsort,
     $   iaddr,rmlexp,
     $   itree,ltree,iptr,nlevels,npwlevel,ndiv,
     $   nboxes,iper,boxsize,tcenters,itree(iptr(1)),
     $   isrcse,itargse,ntermsh,nlocal,ntermmax,pmax,npw,
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
        call dreorderiv(nd,3,ns,gradsort,grad,isrc)
      endif

      if(ifpgh.eq.3) then
        call dreorderi(nd,ns,potsort,pot,isrc)
        call dreorderiv(nd,3,ns,gradsort,grad,isrc)
        call dreorderiv(nd,6,ns,hesssort,hess,isrc)
      endif

      if(ifpghtarg.eq.1) then
        call dreorderi(nd,nt,pottargsort,pottarg,itarg)
      endif

      if(ifpghtarg.eq.2) then
        call dreorderi(nd,nt,pottargsort,pottarg,itarg)
        call dreorderiv(nd,3,nt,gradtargsort,gradtarg,itarg)
      endif

      if(ifpghtarg.eq.3) then
        call dreorderi(nd,nt,pottargsort,pottarg,itarg)
        call dreorderiv(nd,3,nt,gradtargsort,gradtarg,itarg)
        call dreorderiv(nd,6,nt,hesstargsort,hesstarg,itarg)
      endif

      return
      end
c
c
c
c
c
      subroutine fgt3dmain(nd,delta,eps,
     $     nsource,sourcesort,
     $     ifcharge,chargesort,
     $     ifdipole,rnormalsort,dipstrsort,
     $     ntarget,targetsort,
     $     iaddr,rmlexp,
     $     itree,ltree,iptr,nlevels,npwlevel,ndiv,
     $     nboxes,iper,boxsize,centers,laddr,
     $     isrcse,itargse,ntermsh,nlocal,ntermmax,pmax,npw,
     $     ifpgh,pot,grad,hess,
     $     ifpghtarg,pottarg,gradtarg,hesstarg)
c
c
c   the FGT in R^3: evaluate all pairwise particle
c   interactions
c   and interactions with targets using PW expansions.
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
c   sourcesort: real *8 (3,ns):  source locations
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
c   targetsort: real *8 (3,ntarget):  target locations
c   iaddr: (2,nboxes): pointer in rmlexp where multipole
c                      and local expansions for each
c                      box is stored
c                      iaddr(1,ibox) is the
c                      starting index in rmlexp for the 
c                      multipole PW expansion of ibox
c                      and iaddr(2,ibox) is the
c                      starting index in rmlexp
c                      for the local PW expansion of ibox
c
c   itree    in: integer (ltree)
c             This array contains all the information
c             about the tree
c             Refer to pts_tree2d.f
c
c   ltree    in: integer *8
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
c     centers in: real *8(3,nboxes)
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
      real *8 sourcesort(3,nsource)
      real *8 rnormalsort(3,nsource)
      real *8 chargesort(nd,*)
      real *8 dipstrsort(nd,*)
      real *8 targetsort(3,ntarget)
      real *8 pot(nd,*)
      real *8 grad(nd,3,*)
      real *8 hess(nd,6,*)
      real *8 pottarg(nd,*)
      real *8 gradtarg(nd,3,*)
      real *8 hesstarg(nd,6,*)
c
      integer iaddr(2,nboxes)
      real *8 rmlexp(*)
      real *8 pmax
      real *8 timeinfo(10)
      real *8 timelev(0:200)
      real *8 centers(3,*)
c
      integer laddr(2,0:nlevels)
      integer ntermsh(0:nlevels),nlocal(0:nlevels),ntermmax
      integer npw
      integer *8 iptr(8)
      integer *8 ltree
      integer itree(ltree)
      integer nboxes
      integer isrcse(2,nboxes),itargse(2,nboxes)
      real *8 boxsize(0:nlevels)
c
c     temp variables
      integer i,j,k,l,idim,ip,isx,ii,j1,j2
      integer ibox,jbox,ilev,klev,npts,nptssrc,nptstarg,nptsj
      integer nchild,ncoll,nb
      integer isep, mnbors

      integer mnlist1,mnlist2,mnlist3,mnlist4,mnlistpw
      integer nlist1
      integer, allocatable :: list1(:,:)
      integer, allocatable :: nlist1s(:)

      integer, allocatable :: nlistpw(:), listpw(:,:)
      
c
      integer istart,iend,istarts,iends
      integer isstart,isend,jsstart,jsend
      integer jstart,jend
      integer istarte,iende,istartt,iendt

      integer ifprint

      integer nexp,nmax,ncutlevbox
      integer ncdir,ncddir,nddir,npdir,ngdir,nhdir
      
      integer nn,jx,jy,jz,nlevstart,nlevend

      real *8 d,time1,time2,omp_get_wtime
      real *8 dx,dy,dz
      real *8 tt1,tt2,xmin,xmin2,t1,t2,dt,dtt
      real *8 pottmp,gradtmp(2),hesstmp(3)
      real *8 dmax
      
      integer, allocatable :: ifhung(:),ifpwexp(:)
      integer ndirect
      
      real *8, allocatable :: ws(:),ts(:)
      real *8, allocatable :: wnufft(:)
      complex *16, allocatable :: wpwshift(:,:,:,:)
      
      double precision pi
      complex *16 eye
C
      eye = dcmplx(0,1)
      
c      
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
c
      ifprint=1
      pi = 4*atan(1.0d0)

      ncutlevbox=laddr(2,npwlevel)-laddr(1,npwlevel)+1
      write(6,*) ' n per box on the cutoff level',nsource/ncutlevbox      
      

      nlevend=nlevels
      if (npwlevel.le.nlevels) nlevend=npwlevel
      
      allocate(ifhung(nboxes))
      do i=1,nboxes
         ifhung(i)=0
      enddo

      ndirect = 0
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
c     Check if the current box is a nonempty leaf box            
            if(nchild.eq.0.and.npts.gt.0.and.npts.le.ndiv) then
               ifhung(ibox) = 1
               ndirect = ndirect + 1
            endif
         enddo
      enddo
      call prinf('number of direct evaluation source boxes=*',
     1    ndirect,1)
      
      do i=1,10
         timeinfo(i)=0
      enddo

      do i=0,nlevels
         timelev(i) = 0
      enddo

c     compute list info
c
      mnbors = 27
      isep = 1
      call computemnlist1(nlevels,nboxes,itree(iptr(1)),boxsize,
     1  centers,itree(iptr(3)),itree(iptr(4)),
     2  itree(iptr(5)),isep,itree(iptr(6)),mnbors,
     2  itree(iptr(7)),iper,mnlist1)
      
      allocate(list1(mnlist1,nboxes),nlist1s(nboxes))

c     modified list1 for direct evaluation
      call compute_modified_list1(nlevels,npwlevel,
     1  nboxes,itree(iptr(1)),boxsize,
     1  centers,itree(iptr(3)),itree(iptr(4)),
     2  itree(iptr(5)),isep,itree(iptr(6)),mnbors,
     3  itree(iptr(7)),iper,nlist1s,mnlist1,list1)

c     direct evaluation if the cutoff level is >= the maximum level 
      if (npwlevel .ge. nlevels) goto 1800
c      
c     get planewave nodes and weights
      allocate(ws(-npw/2:npw/2-1))
      allocate(ts(-npw/2:npw/2-1))
      call get_pwnodes(pmax,npw,ws,ts)

c     compute translation matrices for PW expansions
      xmin  = boxsize(npwlevel)/sqrt(delta)

      nmax = 1
      nexp = npw*npw*npw/2
      allocate(wpwshift(nexp,-nmax:nmax,-nmax:nmax,-nmax:nmax))
      call pw_translation_matrices(xmin,npw,ts,nmax,wpwshift)

      allocate(wnufft(nexp))
      call nufft_weights(npw,ws,ts,nexp,wnufft)

c     xmin is used in shiftpw subroutines to
c     determine the right translation matrices
c      
      xmin  = boxsize(npwlevel)
c
c     compute list info
c
      call gt3d_computemnlistpw(nlevels,nboxes,itree,ltree,
     1    iptr,centers,
     2    boxsize,iper,mnlistpw)
      allocate(nlistpw(nboxes),listpw(mnlistpw,nboxes))
c     listpw contains source boxes in the pw interaction
      call gt3d_computelistpw(nlevels,npwlevel,nboxes,
     1    itree,ltree,iptr,centers,
     2    itree(iptr(4)),itree(iptr(5)),
     1    boxsize,laddr,
     2    mnlistpw,nlistpw,listpw)
c
c     check whether we need to evaluate local PW expansion for boxes 
c     at the cutoff level
      allocate(ifpwexp(nboxes))
      if (npwlevel .ge. 0 .and. npwlevel .le. nlevels) then
         do i=1,nboxes
            ifpwexp(i)=0
         enddo

         do ilev=npwlevel,npwlevel
            do ibox=laddr(1,ilev),laddr(2,ilev)
               istart = isrcse(1,ibox)
               iend = isrcse(2,ibox)
               npts = iend-istart+1
               if (nlistpw(ibox).gt.0 .or. npts.gt.ndiv) then
                  ifpwexp(ibox)=1
               endif
            enddo
         enddo
      endif
c
c     ... set all multipole and local expansions to zero
c
cccc      do ilev = 0,nlevels
ccccC$OMP PARALLEL DO DEFAULT (SHARED)
ccccC$OMP$PRIVATE(ibox)
cccc         do ibox = laddr(1,ilev),laddr(2,ilev)
cccc           call g3dsoezero_vec(nd,rmlexp(iaddr(1,ibox)),npw)
cccc           call g3dsxzero_vec(nd,rmlexp(iaddr(2,ibox)),
cccc     1          npw,npwhalf)
cccc           call g3dsoezero_vec(nd,rmlexp(iaddr(3,ibox)),npw)
cccc           call g3dsxzero_vec(nd,rmlexp(iaddr(4,ibox)),
cccc     1          npw,npwhalf)
cccc           call g3dlocalzero_vec(nd,rmlexp(iaddr(5,ibox)),nlocal)
cccc           call g3dhermzero_vec(nd,rmlexp(iaddr(6,ibox)),ntermsh)
cccc
cccc         enddo
ccccC$OMP END PARALLEL DO         
cccc       enddo
c
c
      ncdir = 550
      nddir = 750
      ncddir = 1000

      npdir = 450
      ngdir = 900
      nhdir = 1200

      call prinf('laddr=*',laddr,2*(nlevels+1))
      
      if(ifprint .ge. 1) 
     $   call prinf('=== STEP 1 (form mp) ====*',i,0)
      call cpu_time(time1)
C$        time1=omp_get_wtime()
c
c       ... step 1, form multipole pw expansions at the cutoff level
c       
      do 1100 ilev = npwlevel,npwlevel
         nb=0
         dt=0
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
c              Check if current box needs to form pw exp         
               if(npts.gt.ndiv) then
                  nb=nb+1
                  call cpu_time(t1)
c                 form the pw expansion
                  if (npts.gt.ncdir) then
                     call g3dformpwc_fast_vec(nd,delta,eps,
     1                   sourcesort(1,istart),npts,
     2                   chargesort(1,istart),centers(1,ibox),
     3                   npw,ws,ts,wnufft,rmlexp(iaddr(1,ibox)))
                  else
                     call g3dformpwc_vec(nd,delta,
     1                   sourcesort(1,istart),npts,
     2                   chargesort(1,istart),centers(1,ibox),
     3                   npw,ws,ts,rmlexp(iaddr(1,ibox)))
                  endif
                  call cpu_time(t2)
                  dt=dt+t2-t1
c                 copy the multipole PW exp into local PW exp
c                 for self interaction 
                  call g3dcopypwexp_vec(nd,nexp,rmlexp(iaddr(1,ibox)),
     1                rmlexp(iaddr(2,ibox)))
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
c              Check if current box needs to form pw exp           
               if(npts.gt.ndiv) then
c                 form the pw expansion
                  if (npts.gt.nddir) then
                     call g3dformpwd_fast_vec(nd,delta,eps,
     1                   sourcesort(1,istart),npts,
     2                   rnormalsort(1,istart),dipstrsort(1,istart),
     3                   centers(1,ibox),
     4                   npw,ws,ts,wnufft,rmlexp(iaddr(1,ibox)))
                  else
                     call g3dformpwd_vec(nd,delta,
     1                   sourcesort(1,istart),npts,
     2                   rnormalsort(1,istart),dipstrsort(1,istart),
     3                   centers(1,ibox),
     4                   npw,ws,ts,rmlexp(iaddr(1,ibox)))
                  endif
c                 copy the multipole PW exp into local PW exp
                  call g3dcopypwexp_vec(nd,nexp,rmlexp(iaddr(1,ibox)),
     1                rmlexp(iaddr(2,ibox)))
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
c              Check if current box needs to form pw exp          
               if(npts.gt.ndiv) then
                  if (npts.gt.ncddir) then
                     call g3dformpwcd_fast_vec(nd,delta,eps,
     1                   sourcesort(1,istart),npts,chargesort(1,istart),
     2                   rnormalsort(1,istart),dipstrsort(1,istart),
     3                   centers(1,ibox),
     4                   npw,ws,ts,wnufft,rmlexp(iaddr(1,ibox)))
                  else
                     call g3dformpwcd_vec(nd,delta,
     1                   sourcesort(1,istart),npts,chargesort(1,istart),
     2                   rnormalsort(1,istart),dipstrsort(1,istart),
     3                   centers(1,ibox),
     4                   npw,ws,ts,rmlexp(iaddr(1,ibox)))
                  endif                     
c                 copy the multipole PW exp into local PW exp
                  call g3dcopypwexp_vec(nd,nexp,rmlexp(iaddr(1,ibox)),
     1                rmlexp(iaddr(2,ibox)))
               endif
            enddo
C     $OMP END PARALLEL DO
         endif
 111     format ('ilev=', i1,4x, 'nb=',i6, 4x,'formpw=', f6.2)
         write(6,111) ilev,nb,dt
c     end of ilev do loop
 1100 continue


      
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(1)=time2-time1


      
      if(ifprint.ge.1)
     $    call prinf('=== Step 2 (mp to loc) ===*',i,0)
c      ... step 2, convert multipole pw expansions into local
c       pw expansions

      call cpu_time(time1)
C$    time1=omp_get_wtime()
      
      do 1300 ilev = npwlevel,npwlevel
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
c           shift PW expansions
c
            do j=1,nlistpw(ibox)
              jbox=listpw(j,ibox)
              jstart = isrcse(1,jbox)
              jend = isrcse(2,jbox)
              nptsj = jend-jstart+1
              if (nptsj .gt. ndiv) then
                jx= nint((centers(1,ibox) - centers(1,jbox))/xmin)
                jy= nint((centers(2,ibox) - centers(2,jbox))/xmin)
                jz= nint((centers(3,ibox) - centers(3,jbox))/xmin)
              
                call g3dshiftpw_vec(nd,nexp,rmlexp(iaddr(1,jbox)),
     1              rmlexp(iaddr(2,ibox)),wpwshift(1,jx,jy,jz))
              endif
            enddo
          endif
        enddo
C$OMP END PARALLEL DO        
 1300 continue
c

      
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(2) = time2-time1

cccc      call prin2('timeinfo2=*',time2-time1,1)

      if(ifprint.ge.1)
     $    call prinf('=== step 3 (eval loc) ===*',i,0)

c     ... step 5, evaluate all local pw expansions
      call cpu_time(time1)
C$    time1=omp_get_wtime()
      
      do 1500 ilev = npwlevel,npwlevel
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istart,iend,i,npts)
C$OMP$SCHEDULE(DYNAMIC)
         call cpu_time(t1)
         nb=0
         do ibox = laddr(1,ilev),laddr(2,ilev)
            if (ifpwexp(ibox).eq.1) then
               istartt = itargse(1,ibox) 
               iendt = itargse(2,ibox)
               nptstarg = iendt-istartt + 1
               
               istarts = isrcse(1,ibox)
               iends = isrcse(2,ibox)
               nptssrc = iends-istarts+1

               nb=nb+1
c     evaluate local expansion at targets
               if(nptstarg.gt.0) then
                  if (ifpghtarg.eq.1.and.nptstarg.gt.npdir) then
                     call g3dpwevalp_fast_vec(nd,delta,eps,
     1                   centers(1,ibox),npw,ws,ts,
     2                   rmlexp(iaddr(2,ibox)),
     3                   targetsort(1,istartt),
     4                   nptstarg,pottarg(1,istartt))
                  elseif (ifpghtarg.eq.1.and.nptstarg.le.npdir) then
                     call g3dpwevalp_vec(nd,delta,
     1                   centers(1,ibox),npw,ws,ts,
     2                   rmlexp(iaddr(2,ibox)),
     3                   targetsort(1,istartt),
     4                   nptstarg,pottarg(1,istartt))
                  endif
                  if (ifpghtarg.eq.2.and.nptstarg.gt.ngdir) then
                     call g3dpwevalg_fast_vec(nd,delta,eps,
     1                   centers(1,ibox),npw,ws,ts,
     2                   rmlexp(iaddr(2,ibox)),
     3                   targetsort(1,istartt),
     4                   nptstarg,pottarg(1,istartt),
     5                   gradtarg(1,1,istartt))
                  elseif (ifpghtarg.eq.2.and.nptstarg.le.ngdir) then
                     call g3dpwevalg_vec(nd,delta,
     1                   centers(1,ibox),npw,ws,ts,
     2                   rmlexp(iaddr(2,ibox)),
     3                   targetsort(1,istartt),
     4                   nptstarg,pottarg(1,istartt),
     5                   gradtarg(1,1,istartt))
                  endif
                  if (ifpghtarg.eq.3.and.nptstarg.gt.nhdir) then
                     call g3dpwevalh_fast_vec(nd,delta,eps,
     1                   centers(1,ibox),npw,ws,ts,
     2                   rmlexp(iaddr(2,ibox)),
     3                   targetsort(1,istartt),
     4                   nptstarg,pottarg(1,istartt),
     5                   gradtarg(1,1,istartt),
     6                   hesstarg(1,1,istartt))
                  elseif (ifpghtarg.eq.3.and.nptstarg.le.nhdir) then
                     call g3dpwevalh_vec(nd,delta,
     1                   centers(1,ibox),npw,ws,ts,
     2                   rmlexp(iaddr(2,ibox)),
     3                   targetsort(1,istartt),
     4                   nptstarg,pottarg(1,istartt),
     5                   gradtarg(1,1,istartt),
     6                   hesstarg(1,1,istartt))
                  endif
               endif
c     
c     evaluate local expansion at sources
               if (nptssrc.gt.0) then
                  if (ifpgh.eq.1.) then
                     if(nptssrc.gt.npdir) then
                     call g3dpwevalp_fast_vec(nd,delta,eps,
     1                   centers(1,ibox),npw,ws,ts,
     2                   rmlexp(iaddr(2,ibox)),
     3                   sourcesort(1,istarts),nptssrc,
     4                   pot(1,istarts))
                     else
                     call g3dpwevalp_vec(nd,delta,
     1                   centers(1,ibox),npw,ws,ts,
     2                   rmlexp(iaddr(2,ibox)),
     3                   sourcesort(1,istarts),nptssrc,
     4                      pot(1,istarts))
                     endif
                  endif
                  if (ifpgh.eq.2.and.nptssrc.gt.ngdir) then
                     call g3dpwevalg_fast_vec(nd,delta,eps,
     1                   centers(1,ibox),npw,ws,ts,
     2                   rmlexp(iaddr(2,ibox)),
     3                   sourcesort(1,istarts),nptssrc,
     4                   pot(1,istarts),grad(1,1,istarts))
                  elseif (ifpgh.eq.2.and.nptssrc.le.ngdir) then
                     call g3dpwevalg_vec(nd,delta,
     1                   centers(1,ibox),npw,ws,ts,
     2                   rmlexp(iaddr(2,ibox)),
     3                   sourcesort(1,istarts),nptssrc,
     4                   pot(1,istarts),grad(1,1,istarts))
                  endif
                  if (ifpgh.eq.3.and.nptssrc.gt.nhdir) then
                     call g3dpwevalh_fast_vec(nd,delta,eps,
     1                   centers(1,ibox),npw,ws,ts,
     2                   rmlexp(iaddr(2,ibox)),
     3                   sourcesort(1,istarts),nptssrc,
     4                   pot(1,istarts),grad(1,1,istarts),
     5                   hess(1,1,istarts))
                  elseif (ifpgh.eq.3.and.nptssrc.le.nhdir) then
                     call g3dpwevalh_vec(nd,delta,
     1                   centers(1,ibox),npw,ws,ts,
     2                   rmlexp(iaddr(2,ibox)),
     3                   sourcesort(1,istarts),nptssrc,
     4                   pot(1,istarts),grad(1,1,istarts),
     5                   hess(1,1,istarts))
                  endif
               endif
            endif
         enddo
         call cpu_time(t2)
 222     format ('ilev=', i1,4x, 'nb=',i6, 4x,'pweval=', f6.2)
         write(6,222) ilev,nb,t2-t1
C$OMP END PARALLEL DO        
 1500 continue

      
      call cpu_time(time2)
C$    time2 = omp_get_wtime()      
      timeinfo(3) = time2 - time1
      
      
 1800 continue

      if(ifprint .ge. 1)
     $     call prinf('=== STEP 4 (direct) =====*',i,0)
c
cc
      call cpu_time(time1)
      dmax = log(1.0d0/eps)*delta
      
C$    time1=omp_get_wtime()  
      do 2000 ilev = 0,nlevend
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,istartt,iendt,i,jstart,jend,istarte,iende)
C$OMP$PRIVATE(istarts,iends,nptssrc,nptstarg)
C$OMP$SCHEDULE(DYNAMIC)  
         do jbox = laddr(1,ilev),laddr(2,ilev)
c        jbox is the source box            
            if (ifhung(jbox) .eq. 1) then
               
               jstart = isrcse(1,jbox)
               jend = isrcse(2,jbox)
               
               nlist1 = nlist1s(jbox)
               do i=1,nlist1
cccc              ibox is the target box
                  ibox = list1(i,jbox)

                  istarts = isrcse(1,ibox)
                  iends = isrcse(2,ibox)
                  nptssrc = iends-istarts + 1
                  
                  istartt = itargse(1,ibox)
                  iendt = itargse(2,ibox)
                  nptstarg = iendt-istartt + 1

                  if (nptstarg .gt. 0) then
                     call fgt3dpart_direct_vec(nd,delta,dmax,
     1                   jstart,jend,istartt,iendt,sourcesort,
     2                   ifcharge,chargesort,
     3                   ifdipole,rnormalsort,dipstrsort,targetsort,
     4                   ifpghtarg,pottarg,gradtarg,hesstarg)
                  endif
                  if (nptssrc .gt. 0) then
                     call fgt3dpart_direct_vec(nd,delta,dmax,
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
      timeinfo(4) = time2-time1
      if(ifprint.ge.1) call prin2('timeinfo=*',timeinfo,4)
      d = 0
      do i = 1,4
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
      subroutine fgt3dpart_direct_vec(nd,delta,dmax,istart,iend,
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
c     source       in: real *8(3,ns)
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
c     targ        in: real *8(3,nt)
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
        integer istart,iend,jstart,jend,ns,j,i,ntarg
        integer ifcharge,ifdipole
        integer nd
        integer ifpgh
c
        real *8 source(3,*)
        real *8 rnormal(3,*)
        real *8 charge(nd,*),dipstr(nd,*)
        real *8 targ(3,*),delta,eps,dmax
        real *8 pot(nd,*)
        real *8 grad(nd,3,*)
        real *8 hess(nd,6,*)
c
        
        ns = iend - istart + 1
        ntarg = jend-jstart+1
        
        if(ifcharge.eq.1.and.ifdipole.eq.0) then
          if(ifpgh.eq.1) then
             call g3d_directcp_vec(nd,delta,dmax,source(1,istart),ns,
     1           charge(1,istart),targ(1,jstart),ntarg,pot(1,jstart))
          endif

          if(ifpgh.eq.2) then
             call g3d_directcg_vec(nd,delta,dmax,source(1,istart),ns,
     1           charge(1,istart),targ(1,jstart),ntarg,pot(1,jstart),
     2           grad(1,1,jstart))
          endif
          if(ifpgh.eq.3) then 
             call g3d_directch_vec(nd,delta,dmax,source(1,istart),ns,
     1           charge(1,istart),targ(1,jstart),ntarg,pot(1,jstart),
     2           grad(1,1,jstart),hess(1,1,jstart))
          endif
        endif

        if(ifcharge.eq.0.and.ifdipole.eq.1) then
          if(ifpgh.eq.1) then
             call g3d_directdp_vec(nd,delta,dmax,source(1,istart),ns,
     1           rnormal(1,istart),dipstr(1,istart),targ(1,jstart),
     2           ntarg,pot(1,jstart))
          endif

          if(ifpgh.eq.2) then
             call g3d_directdg_vec(nd,delta,dmax,source(1,istart),ns,
     1           rnormal(1,istart),dipstr(1,istart),targ(1,jstart),
     2           ntarg,pot(1,jstart),grad(1,1,jstart))
          endif
          if(ifpgh.eq.3) then
             call g3d_directdh_vec(nd,delta,dmax,source(1,istart),ns,
     1           rnormal(1,istart),dipstr(1,istart),targ(1,jstart),
     2           ntarg,pot(1,jstart),grad(1,1,jstart),hess(1,1,jstart))
          endif
        endif 

        if(ifcharge.eq.1.and.ifdipole.eq.1) then
          if(ifpgh.eq.1) then
             call g3d_directcdp_vec(nd,delta,dmax,source(1,istart),ns,
     1           charge(1,istart),rnormal(1,istart),dipstr(1,istart),
     2           targ(1,jstart),ntarg,pot(1,jstart))
          endif

          if(ifpgh.eq.2) then
             call g3d_directcdg_vec(nd,delta,dmax,source(1,istart),ns,
     1           charge(1,istart),rnormal(1,istart),dipstr(1,istart),
     2           targ(1,jstart),ntarg,pot(1,jstart),grad(1,1,jstart))
          endif
          if(ifpgh.eq.3) then
             call g3d_directcdh_vec(nd,delta,dmax,source(1,istart),ns,
     1           charge(1,istart),rnormal(1,istart),dipstr(1,istart),
     2           targ(1,jstart),ntarg,pot(1,jstart),grad(1,1,jstart),
     3           hess(1,1,jstart))
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
      subroutine g3dmpalloc(nd,laddr,iaddr,
     1    nlevels,npwlevel,lmptot,
     2    ntermsh,nlocal,npw)
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
c
c
c------------------------------------------------------------------
c     Output arguments
c     iaddr: (2,nboxes): pointer in rmlexp where multipole
c                      and local expansions for each
c                      box is stored
c                      iaddr(1,ibox) is the
c                      starting index in rmlexp for the 
c                      multipole PW expansion of ibox
c                      and iaddr(2,ibox) is the
c                      starting index in rmlexp
c                      for the local PW expansion of ibox
c     lmptot      out: Integer
c                 Total length of expansions array required
c------------------------------------------------------------------

      implicit none
      integer nlevels,npwlevel,npw,nd
      integer ntermsh(0:nlevels),nlocal(0:nlevels)
      integer iaddr(2,*), laddr(2,0:nlevels)
      integer *8 lmptot
      integer ibox,i,istart,nn,itmp,nlevstart
c
      istart = 1
      if (npwlevel .eq. nlevels) then
         lmptot=0
         return
      endif
      
      nlevstart = 0
      if (npwlevel .ge. 0) nlevstart = npwlevel

      do i = nlevstart,nlevstart
cccc         print *, 'nboxes at npwlevel=',laddr(2,i)-laddr(1,i)+1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the multipole PW expansion         
c
            nn = (npw*npw*npw/2)*2*nd
           itmp = ibox - laddr(1,i)
           iaddr(1,ibox) = istart + itmp*nn
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
      enddo
c
      do i = nlevstart,nlevstart

         nn = (npw*npw*npw/2)*2*nd
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the local PW expansion         
c
           itmp = ibox - laddr(1,i)
           iaddr(2,ibox) = istart + itmp*nn
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
      subroutine fgt3d_large_delta(nd,delta,eps,ns,sources,
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
c   sources(3,ns) : source locations
c   ifcharge      : flag for including charge interactions
c                   charge interactions included if ifcharge =1
c                   not included otherwise
c   charge(nd,ns) : charge strengths
c   ifdipole      : flag for including dipole interactions
c                   dipole interactions included if ifcharge =1
c                   not included otherwise
c   rnormal(3,ns) : dipole directions
c   dipstr(nd,ns) : dipole strengths
c   iper          : flag for periodic implmentations. Currently unused
c   nlocal        : number of terms in the local expansion
c   center(3)     : the center of the bounding box
c   ifpgh         : flag for computing pot/grad/hess
c                   ifpgh = 1, only potential is computed
c                   ifpgh = 2, potential and gradient are computed
c                   ifpgh = 3, potential, gradient, and hessian 
c                   are computed
c   nt            : number of targets
c   targ(3,nt)    : target locations
c   ifpghtarg     : flag for computing pottarg/gradtarg/hesstarg
c                   ifpghtarg = 1, only potential is computed at targets
c                   ifpghtarg = 2, potential and gradient are 
c                   computed at targets
c                   ifpghtarg = 3, potential, gradient, and hessian are 
c                   computed at targets
c
c   OUTPUT PARAMETERS
c   pot(nd,*)       : potential at the source locations
c   grad(nd,3,*)    : gradients at the source locations
c   hess(nd,6,*)    : hessian at the source locations
c   pottarg(nd,*)   : potential at the target locations
c   gradtarg(nd,3,*): gradient at the target locations
c   hesstarg(nd,6,*): hessian at the target locations
c
      implicit none
c
cc      calling sequence variables
c 
      integer nd,ns,nt,i,j,k,ind
      integer iper,ifcharge,ifdipole
      integer ifpgh,ifpghtarg,nlocal

      real *8 omp_get_wtime
      real *8 eps,delta,bgf,bsize
      real *8 sources(3,ns),targ(3,nt)
      real *8 rnormal(3,ns)
      real *8 charge(nd,*),dipstr(nd,*)
      real *8 center(3)

      real *8 pot(nd,*),grad(nd,3,*),hess(nd,6,*)
      real *8 pottarg(nd,*),gradtarg(nd,3,*),hesstarg(nd,6,*)


      real *8, allocatable :: local(:,:,:,:)

      allocate(local(0:nlocal,0:nlocal,0:nlocal,nd))
      do ind=1,nd
         do i=0,nlocal
            do j=0,nlocal
               do k=0,nlocal
                  local(k,j,i,ind)=0
               enddo
            enddo
         enddo
      enddo
      
c
c     form local
c
      if (ifcharge.eq.1 .and. ifdipole.eq.0) then
         call g3dformlc_vec(nd,delta,sources,ns,charge,center,
     1       nlocal,local)
      elseif (ifcharge.eq.0 .and. ifdipole.eq.1) then
         call  g3dformld_vec(nd,delta,sources,ns,rnormal,dipstr,
     1       center,nlocal,local)
      elseif (ifcharge.eq.1 .and. ifdipole.eq.1) then
         call g3dformlcd_vec(nd,delta,sources,ns,charge,rnormal,
     1       dipstr,center,nlocal,local)
      endif

c
c     local eval
c
c     targets
c      
      if(ifpghtarg.eq.1) then
         call g3dlevalp_vec(nd,delta,center,nlocal,local,
     1       targ,nt,pottarg)
      endif
      if(ifpghtarg.eq.2) then
         call g3dlevalg_vec(nd,delta,center,nlocal,local,
     1       targ,nt,pottarg,gradtarg)
      endif
      if(ifpghtarg.eq.3) then
         call g3dlevalh_vec(nd,delta,center,nlocal,local,
     1       targ,nt,pottarg,gradtarg,hesstarg)
      endif
c
c     sources
c      
      if(ifpgh.eq.1) then
         call g3dlevalp_vec(nd,delta,center,nlocal,local,
     1       sources,ns,pot)
      endif
      if(ifpgh.eq.2) then
         call g3dlevalg_vec(nd,delta,center,nlocal,local,
     1       sources,ns,pot,grad)
      endif
      if(ifpgh.eq.3) then
         call g3dlevalh_vec(nd,delta,center,nlocal,local,
     1       sources,ns,pot,grad,hess)
      endif
      
      return
      end
c
c
c
c
