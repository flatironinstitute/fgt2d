cc Copyright (C) 2018-2019: Leslie Greengard, Manas Rachh, Shidong Jiang
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
      real *8 omp_get_wtime
      real *8 eps,delta,bgf,bsize
      integer ns,nt,ntermsh,nlocal,npw,nsoe
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
      integer iper,nlmin,nlmax,ifunif,ngpl,nlgf,ngrid
      real *8, allocatable :: tcenters(:,:),boxsize(:)
      integer idivflag,nlevels,nboxes,ndiv,npwlevel
      integer ltree

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

      integer lmptot
      integer, allocatable :: iaddr(:,:)
      real *8, allocatable :: rmlexp(:)

c
cc      temporary variables
c
      integer i,ilev,lmptmp,idim
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg,ifprint,ier
      real *8 time1,time2,pi,done

      done = 1
      pi = atan(done)*4.0d0

c
c     Need to fix ndiv in Laplace FMM
c   
c
      idivflag =0
      ndiv = 4
c
      ltree = 0
      nboxes = 0
      nlevels = 0
c
      ifunif = 1
      iper = 0

      ifprint = 1
c
c     call the tree memory management
c     code to determine number of boxes,
c     number of levels and length of tree
c
c    
c     ifunif = 1
c
c     grandparent level = finest level with bsize > bgf
c
      bgf = 3*sqrt(delta)
      if (eps.lt. 1.0d-3) bgf = 4*sqrt(delta)
      if (eps.lt. 1.0d-6) bgf = 5*sqrt(delta)
      if (eps.lt. 1.0d-10) bgf = sqrt(28*delta)
      write(6,*) ' bgf is ',bgf
c
c     divide the unit square into ngrid x ngrid boxes at level 0
c     bgf should be a parameter of pts_tree_mem so that the outer
c     box size is an integer of bgf to avoid possible loss of accuracy
c
      ngrid = ceiling(1/bgf)
c     nlmax is the number of levels, either 2 or 3
      nlmax = 1
      if (ns/(ngrid*ngrid*16) .gt. 120) nlmax=3
c      
      call pts_tree_mem(sources,ns,targ,nt,idivflag,ndiv,ngrid,nlmax,
     1  ifunif,iper,nlevels,nboxes,ltree)

c     setting npwlevel = 0 goes back to Leslie's picture
      npwlevel = nlevels
      if (nlevels .eq. 3) npwlevel=2
c 
      write(6,*) ' ngrid',ngrid
      write(6,*) ' nlmax',nlmax
      write(6,*) ' nlevels',nlevels
      write(6,*) ' npwlevel',npwlevel
      write(6,*) ' nboxes',nboxes
      write(6,*) ' nperbox',ns/(ngrid*ngrid*(2**(2*nlevels)))

      allocate(itree(ltree))
      allocate(boxsize(0:nlevels))
      allocate(tcenters(2,nboxes))

c
c     call the tree code
c
      call pts_tree_build(sources,ns,targ,nt,idivflag,ndiv,ngrid,nlmax,
     1  ifunif,iper,nlevels,nboxes,ltree,itree,iptr,tcenters,boxsize)
      write(6,*) ' boxsize',boxsize(0)

      allocate(isrc(ns),isrcse(2,nboxes))
      allocate(itarg(nt),itargse(2,nboxes))

      call pts_tree_sort(ns,sources,itree,ltree,nboxes,nlevels,iptr,
     1   tcenters,ngrid,boxsize,isrc,isrcse)
ccc      write(6,*)' isrcse',isrcse

      call pts_tree_sort(nt,targ,itree,ltree,nboxes,nlevels,iptr,
     1   tcenters,ngrid,boxsize,itarg,itargse)
ccc      write(6,*)' itargse',itargse
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
      if (nlevels .eq. 3) then
         if (ifdipole.eq.1 .or. ifpgh.ge.2 .or. ifpghtarg.ge.2) then
            ntermsh = 16
            nlocal = 16
         else
            ntermsh = 14
            nlocal = 14
         endif
      elseif (nlevels .eq. 2) then
         if (ifdipole.eq.1 .or. ifpgh.ge.2 .or. ifpghtarg.ge.2) then
            ntermsh = 26
            nlocal = 26
         else
            ntermsh = 22
            nlocal = 22
         endif
      elseif (nlevels .eq. 1) then
         if (ifdipole.eq.1 .or. ifpgh.ge.2 .or. ifpghtarg.ge.2) then
            ntermsh = 38
            nlocal = 38
         else
            ntermsh = 34
            nlocal = 34
         endif
      elseif (nlevels .eq. 0) then
         if (ifdipole.eq.1 .or. ifpgh.ge.2 .or. ifpghtarg.ge.2) then
            ntermsh = 76
            nlocal = 76
         else
            ntermsh = 70
            nlocal = 70
         endif
      endif

      if (npwlevel .eq. 2) then
         npw = 11
      elseif (npwlevel .eq. 1) then
         npw = 14
      elseif (npwlevel .eq. 0) then
         npw = 20
      endif

c     use 14 for better accuracy
c     use 18 for even better accuracy in hessian computation
      if (ifdipole.eq.1 .or. ifpgh.ge.2 .or. ifpghtarg.ge.2) then
         nsoe = 14
cccc         nsoe = 18
      else
         nsoe = 12
      endif
      
      call prinf(' ntermsh =*',ntermsh,1)
      call prinf(' nlocal =*',nlocal,1)
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
      call g2dmpalloc(nd,itree,iaddr,nlevels,lmptot,
     1    ntermsh,nlocal,npw,nsoe)
      if(ifprint .eq. 1) call prinf_long(' lmptot is *',lmptot,1)
c
      allocate(rmlexp(lmptot),stat=ier)
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
     $   itree,ltree,iptr,ndiv,nlevels,npwlevel,
     $   nboxes,iper,boxsize,tcenters,itree(iptr(1)),
     $   isrcse,itargse,ntermsh,nlocal,npw,nsoe,
     $   ifpgh,potsort,gradsort,hesssort,
     $   ifpghtarg,pottargsort,gradtargsort,
     $   hesstargsort)
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      if( ifprint .eq. 1 ) call prin2('time in fmm main=*',
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
     $     itree,ltree,iptr,ndiv,nlevels,npwlevel,
     $     nboxes,iper,boxsize,centers,laddr,
     $     isrcse,itargse,ntermsh,nlocal,npw,nsoe,
     $     ifpgh,pot,grad,hess,
     $     ifpghtarg,pottarg,gradtarg,hesstarg)
c
c
c   the FGT in R^2: evaluate all pairwise particle
c   interactions
c   and interactions with targets using Hermite/SOE expansions.
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
c                      SOE expansion of ibox
c                      and iaddr(2,ibox) is the
c                      starting index in rmlexp
c                      for the SX expansion of ibox
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
c     ndiv    in: integer
c             Max number of points per box at the level where only the interactions
c             between its neighbors need to be considered
c
c
c     nlevels in: integer
c             number of levels in the tree
c
c     
c     npwlevel   in: integer
c             level at which the X expansion is valid
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

      integer nd,ii
      integer iper,ncoll
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
      complex *16 ws(100),ts(100)
      real *8,allocatable:: wx(:),tx(:)
c
      integer iaddr(6,nboxes),lmptmp
      real *8 rmlexp(*)
      real *8 timeinfo(10)
      real *8 timelev(0:200)
      real *8 centers(2,*)
c
      integer laddr(2,0:nlevels)
      integer ntermsh,nlocal,npw,nsoe,nsoehalf
      integer iptr(8),ltree
      integer itree(ltree)
      integer nboxes
      integer isrcse(2,nboxes),itargse(2,nboxes)
      real *8 boxsize(0:nlevels)
c
c     temp variables
      integer i,j,k,l,idim,ip,isx
      integer ibox,jbox,ilev,npts,nptssrc,nptstarg
      integer nchild
c
      integer istart,iend,istarts,iends
      integer isstart,isend,jsstart,jsend
      integer jstart,jend
      integer istarte,iende,istartt,iendt

      integer ifprint

      integer nexp,nshiftmax,nmax,kshift
      integer nsoeexp,nsxexp
      
      integer ifhesstarg,nn,jx,jy
      real *8 d,time1,time2,omp_get_wtime
      real *8 tt1,tt2,xmin,xmin2,t1,t2,pmax
      real *8 pottmp,gradtmp(2),hesstmp(3)

      integer mnlistsoe,mnlistsx
      integer, allocatable :: nlistsoe(:,:),listsoe(:,:,:)
      integer, allocatable :: nlistsx(:,:),listsx(:,:,:)
      integer, allocatable :: ifhung(:)
      
      real *8, allocatable :: h2l(:,:),rh2x(:,:)
      
      complex *16, allocatable :: h2s(:,:),h2x(:,:)
      complex *16, allocatable :: s2l(:,:),x2l(:,:)

      complex *16, allocatable :: wsoeshift(:,:),wsxshift(:,:)
      complex *16, allocatable :: wsoeshiftall(:,:,:),wsxshiftall(:,:,:)
      
      double precision pi
      
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
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
      call get_pwnodes(npw,wx,tx)
      pmax = 10.6d0

c     compute translation matrices for SOE and SOE/X expansions
      xmin = boxsize(nlevels)/sqrt(delta)

      nexp = nsoe*nsoe/2
      nshiftmax = 2**nlevels
      
      allocate(wsoeshift(nexp,nshiftmax**2))
      allocate(wsxshift((2*npw+1)*nsoe/2,nshiftmax))

      call shiftsoe_translation_matrices(xmin,nsoe,nshiftmax,
     1    wsoeshift,ws,ts)
      call shiftsx_translation_matrices(xmin,npw,nsoe,nshiftmax,
     1    wsxshift,ws,ts)

      xmin2 = xmin/2
      nmax = 2**nlevels

      nsoeexp = 4*nsoe*nsoe/2
      nsxexp  = 4*(2*npw+1)*nsoe/2
      
      allocate(wsoeshiftall(nsoeexp,-nmax:nmax,-nmax:nmax))
      allocate(wsxshiftall(nsxexp,-nmax:nmax,-nmax:nmax))
      call shiftsoeall_translation_matrices(xmin2,nsoe,nmax,
     1           wsoeshiftall,ws,ts)
      call shiftsxall_translation_matrices(xmin2,npw,nsoehalf,nmax,
     1           wsxshiftall,ts,tx)
      

c     xmin is used in shiftsoe and shiftsx subroutines to
c     determine the right translation matrices
c      
      xmin=boxsize(nlevels)
      xmin2=xmin/2

c
c     1d translation matrices
c
c     Hermite to SOE
      allocate(h2s(0:ntermsh,nsoe))
      call g2dh2smat(ntermsh,nsoe,ws,ts,h2s)

c     SOE to local
      allocate(s2l(nsoe,0:nlocal))
      call g2ds2lmat(nlocal,nsoe,ws,ts,s2l)
      
      
c     Hermite to planewave
      allocate(h2x(0:ntermsh,-npw:npw))
      call g2dh2xmat(ntermsh,npw,wx,tx,h2x)
c     real Hermite to planewave translation matrix      
      allocate(rh2x(0:ntermsh,-npw:npw))
      call g2drh2xmat(ntermsh,npw,wx,tx,rh2x)

c     planewave to local
      allocate(x2l(-npw:npw,0:nlocal))
      call g2dx2lmat(nlocal,npw,wx,tx,x2l)

c     Hermite to local
      allocate(h2l(0:ntermsh,0:nlocal))
      call g2dh2lmat(nlocal,ntermsh,h2l)
      
c
      ifprint=1
      pi = 4*atan(1.0d0)

c
      allocate(ifhung(nboxes))
      do i=1,nboxes
         ifhung(i)=0
      enddo

      do ilev = 0,0
C
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            nchild = itree(iptr(4)+ibox-1)
            istart = isrcse(1,ibox)
            iend = isrcse(2,ibox)
            npts = iend-istart+1
cccc  if (nchild .eq. 0) call prinf('npts=*',npts,1)
c     Check if current box is a leaf box            
            if(nchild.eq.0.and.npts.le.ndiv) then
               ifhung(ibox) = 1
            endif
         enddo
      enddo

cccc      call prinf('ifhung=*',ifhung,nboxes)
      
      do i=0,nlevels
         timelev(i) = 0
      enddo
c
c     compute list info
c
      t1=second()
      call gt2d_computemnlists(nlevels,nboxes,itree,ltree,
     1    iptr,centers,
     2    boxsize,iper,mnlistsoe,mnlistsx)

      allocate(nlistsoe(0:3,nboxes),listsoe(mnlistsoe,0:3,nboxes))
      allocate(nlistsx (0:3,nboxes), listsx (mnlistsx,0:3,nboxes))

      call gt2d_computelists_slow(nlevels,npwlevel,nboxes,
     1    itree,ltree,iptr,centers,boxsize,laddr,
     2    mnlistsoe,nlistsoe,listsoe,
     3    mnlistsx, nlistsx, listsx)
      
      t2=second()
      call prin2('time on computelists*',t2-t1,1)
      
      do i=1,10
         timeinfo(i)=0
      enddo
c
c     ... set all multipole and local expansions to zero
c
      do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox)
         do ibox = laddr(1,ilev),laddr(2,ilev)
           call g2dsoezero_vec(nd,rmlexp(iaddr(1,ibox)),nsoe)
           call g2dsxzero_vec(nd,rmlexp(iaddr(2,ibox)),
     1          npw,nsoehalf)
           call g2dsoezero_vec(nd,rmlexp(iaddr(3,ibox)),nsoe)
           call g2dsxzero_vec(nd,rmlexp(iaddr(4,ibox)),
     1          npw,nsoehalf)
           call g2dlocalzero_vec(nd,rmlexp(iaddr(5,ibox)),nlocal)
           call g2dhermzero_vec(nd,rmlexp(iaddr(6,ibox)),ntermsh)

         enddo
C$OMP END PARALLEL DO         
       enddo
c
c
      if(ifprint .ge. 1) 
     $   call prinf('=== STEP 1 (form mp) ====*',i,0)
        call cpu_time(time1)
C$        time1=omp_get_wtime()
c
c       ... step 1, locate all charges, assign them to boxes, and
c       form multipole expansions

      do ilev = nlevels,nlevels
C
        if(ifcharge.eq.1.and.ifdipole.eq.0) then
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
          do ibox=laddr(1,ilev),laddr(2,ilev)
             if (ifhung(ibox) .eq. 0) then 
             nchild = itree(iptr(4)+ibox-1)
             istart = isrcse(1,ibox)
             iend = isrcse(2,ibox)
             npts = iend-istart+1
cccc             if (nchild .eq. 0) call prinf('npts=*',npts,1)
c            Check if current box is a leaf box            
             if(nchild.eq.0.and.npts.gt.0) then
c
c               form the Hermite expansion
c                
                call g2dformhc_vec(nd,delta,
     1             sourcesort(1,istart),npts,chargesort(1,istart),
     2             centers(1,ibox),
     3             ntermsh,rmlexp(iaddr(6,ibox)))
c
c convert the Hermite expansion to SOE/X expansions and SOE expansions
c                
                call g2dh2sx_h2s_vec(nd,ntermsh,npw,nsoe,h2s,h2x,
c                call g2dh2sx_h2s_real_vec(nd,ntermsh,npw,nsoe,h2s,rh2x,
     1              rmlexp(iaddr(6,ibox)),
     2              rmlexp(iaddr(2,ibox)),rmlexp(iaddr(1,ibox)))
cccc                call g2dformsc_vec(nd,delta,
cccc     1              sourcesort(1,istart),npts,chargesort(1,istart),
cccc     2              centers(1,ibox),
cccc     3              nsoe,ws,ts,rmlexp(iaddr(1,ibox)))
cccc                call g2dformsxc_vec(nd,delta,
cccc     1              sourcesort(1,istart),npts,chargesort(1,istart),
cccc     2              centers(1,ibox),pmax,npw,
cccc     3              nsoehalf,ws,ts,rmlexp(iaddr(2,ibox)))
                if (ntermsh .gt. 40) then
c
c               form local Taylor expansion
c                
                   call g2dformlc_vec(nd,delta,
     1                 sourcesort(1,istart),npts,chargesort(1,istart),
     2                 centers(1,ibox),
     3                 nlocal,rmlexp(iaddr(5,ibox)))
                else
c     
c convert the Hermite expansion to local Taylor expansion
c
                   call g2dherm2local_vec(nd,nlocal,ntermsh,h2l,
     1                 rmlexp(iaddr(6,ibox)),rmlexp(iaddr(5,ibox)))
                endif
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
             if (ifhung(ibox) .eq. 0) then
             nchild = itree(iptr(4)+ibox-1)
             istart = isrcse(1,ibox)
             iend = isrcse(2,ibox)
             npts = iend-istart+1
c              Check if current box is a leaf box            
             if(nchild.eq.0.and.npts.gt.0) then
                call g2dformhd_vec(nd,delta,
     1             sourcesort(1,istart),npts,rnormalsort(1,istart),
     2             dipstrsort(1,istart),centers(1,ibox),
     3             ntermsh,rmlexp(iaddr(6,ibox)))
c
c convert the Hermite expansion to SOE/X expansions and SOE expansions
c                
                call g2dh2sx_h2s_vec(nd,ntermsh,npw,nsoe,h2s,h2x,
cc                
c                          
cc              uncomment to see the performance of real h2x  
c                call g2dh2sx_h2s_real_vec(nd,ntermsh,npw,nsoe,h2s,rh2x,
     1              rmlexp(iaddr(6,ibox)),
     2              rmlexp(iaddr(2,ibox)),rmlexp(iaddr(1,ibox)))
                if (ntermsh .gt. 40) then
                   call g2dformld_vec(nd,delta,
     1                 sourcesort(1,istart),npts,rnormalsort(1,istart),
     2                 dipstrsort(1,istart),centers(1,ibox),
     3                 nlocal,rmlexp(iaddr(5,ibox)))
                else
                
c
c convert the Hermite expansion to local Taylor expansion
c                
                   call g2dherm2local_vec(nd,nlocal,ntermsh,h2l,
     1                 rmlexp(iaddr(6,ibox)),rmlexp(iaddr(5,ibox)))
                endif
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
             if (ifhung(ibox) .eq. 0) then
             nchild = itree(iptr(4)+ibox-1)
             istart = isrcse(1,ibox)
             iend = isrcse(2,ibox)
             npts = iend-istart+1
c             Check if current box is a leaf box            
             if(nchild.eq.0.and.npts.gt.0) then
                call g2dformhcd_vec(nd,delta,
     1             sourcesort(1,istart),npts,chargesort(1,istart),
     2             rnormalsort(1,istart),dipstrsort(1,istart),
     3             centers(1,ibox),ntermsh,rmlexp(iaddr(6,ibox)))
c
c convert the Hermite expansion to SOE/X expansions and SOE expansions
c                
                call g2dh2sx_h2s_vec(nd,ntermsh,npw,nsoe,h2s,h2x,
c                call g2dh2sx_h2s_real_vec(nd,ntermsh,npw,nsoe,h2s,rh2x,
     1              rmlexp(iaddr(6,ibox)),
     2              rmlexp(iaddr(2,ibox)),rmlexp(iaddr(1,ibox)))
                if (ntermsh .gt. 40) then
                   call g2dformlcd_vec(nd,delta,
     1                 sourcesort(1,istart),npts,chargesort(1,istart),
     2                 rnormalsort(1,istart),dipstrsort(1,istart),
     3                 centers(1,ibox),nlocal,rmlexp(iaddr(5,ibox)))
                else
c     
c convert the Hermite expansion to local Taylor expansion
c
                   call g2dherm2local_vec(nd,nlocal,ntermsh,h2l,
     1                 rmlexp(iaddr(6,ibox)),rmlexp(iaddr(5,ibox)))
                endif
             endif
             endif
          enddo
C$OMP END PARALLEL DO 
        endif
      enddo

      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(1)=time2-time1

      if(ifprint .ge. 1)
     $      call prinf('=== STEP 2 (merge mp) ====*',i,0)
      call cpu_time(time1)
C$    time1=omp_get_wtime()
c
      do ilev=nlevels-1,0,-1
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
              call g2dshiftall_vec(nd,nsoeexp,rmlexp(iaddr(1,jbox)),
     1             rmlexp(iaddr(1,ibox)),wsoeshiftall(1,jx,jy))
              if (ilev .ge. npwlevel) then
                 call g2dshiftall_vec(nd,nsxexp,rmlexp(iaddr(2,jbox)),
     1               rmlexp(iaddr(2,ibox)),wsxshiftall(1,jx,jy))
              endif
            endif
          enddo
        enddo
C$OMP END PARALLEL DO    
      enddo
111   continue
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

      do ilev = 0,nlevels
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
            if (ilev .ge. npwlevel) then
              do isx=0,3
              do j=1,nlistsx(isx,ibox)
                jbox=listsx(j,isx,ibox)
                if (isx .lt. 2) then
                  kshift=1.05*abs(centers(1,jbox)-centers(1,ibox))/xmin
                else
                  kshift=1.05*abs(centers(2,jbox)-centers(2,ibox))/xmin
                endif
                call g2dshiftsx_vec(nd,npw,nsoehalf,
     1              rmlexp(iaddr(2,jbox)),rmlexp(iaddr(4,ibox)),
     2              wsxshift(1,kshift),isx)
              enddo
              enddo
            endif
          endif
        enddo
C$OMP END PARALLEL DO        
      enddo
c
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(3) = time2-time1

      call prin2('timeinfo4=*',time2-time1,1)

      if(ifprint.ge.1)
     $    call prinf('=== Step 4 (split loc) ===*',i,0)

      call cpu_time(time1)
C$    time1=omp_get_wtime()
      do ilev = 0,nlevels-1
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
              call g2dshiftall_vec(nd,nsoeexp,rmlexp(iaddr(3,ibox)),
     1             rmlexp(iaddr(3,jbox)),wsoeshiftall(1,jx,jy))
              if (ilev .ge. npwlevel) then
                 call g2dshiftall_vec(nd,nsxexp,rmlexp(iaddr(4,ibox)),
     1               rmlexp(iaddr(4,jbox)),wsxshiftall(1,jx,jy))
              endif
            enddo
          endif
        enddo
C$OMP END PARALLEL DO        
      enddo
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(4) = time2-time1

      if(ifprint.ge.1)
     $    call prinf('=== step 5 (eval loc) ===*',i,0)

c     ... step 5, evaluate all local expansions
      call cpu_time(time1)
C$    time1=omp_get_wtime()
      do ilev = nlevels,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istart,iend,i,npts)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          nchild = itree(iptr(4)+ibox-1)
          if(nchild.eq.0) then
c
c     evaluate local expansion at targets
c                
            istart = itargse(1,ibox) 
            iend = itargse(2,ibox)
            npts = iend-istart + 1

c     
c     convert SOE expansions to local expansions
c
            call g2dsoe2local_vec(nd,nlocal,nsoe,s2l,
     1          rmlexp(iaddr(3,ibox)),rmlexp(iaddr(5,ibox)))
c            
c     convert SOE/X expansions to local expansions
c
            call g2dsx2local_vec(nd,nlocal,npw,nsoe,x2l,s2l,
     1          rmlexp(iaddr(4,ibox)),rmlexp(iaddr(5,ibox)))

c
c     use Taylor expansion
c
            if(ifpghtarg.eq.1) then
               call g2dlevalp_vec(nd,delta,centers(1,ibox),
     1             nlocal,rmlexp(iaddr(5,ibox)),targetsort(1,istart),
     2             npts,pottarg(1,istart))
cccc               call g2dsevalp_vec(nd,delta,centers(1,ibox),nsoe,ws,ts,
cccc     1             rmlexp(iaddr(3,ibox)),targetsort(1,istart),
cccc     2             npts,pottarg(1,istart))
cccc               
cccc               call g2dsxevalp_vec(nd,delta,centers(1,ibox),pmax,npw,
cccc     1             nsoehalf,rmlexp(iaddr(4,ibox)),targetsort(1,istart),
cccc     2             npts,ws,ts,pottarg(1,istart))
               
            endif
            if(ifpghtarg.eq.2) then
               call g2dlevalg_vec(nd,delta,centers(1,ibox),
     1             nlocal,rmlexp(iaddr(5,ibox)),targetsort(1,istart),
     2             npts,pottarg(1,istart),gradtarg(1,1,istart))
            endif
            if(ifpghtarg.eq.3) then
               call g2dlevalh_vec(nd,delta,centers(1,ibox),
     1             nlocal,rmlexp(iaddr(5,ibox)),targetsort(1,istart),
     2             npts,pottarg(1,istart),gradtarg(1,1,istart),
     3             hesstarg(1,1,istart))
            endif
c     
c     evaluate local expansion at sources
c
            istart = isrcse(1,ibox)
            iend = isrcse(2,ibox)
            npts = iend-istart+1
cccc            call prinf('source points=*',npts,1)
            
c
c             use Taylor expansion
c
            if(ifpgh.eq.1) then
               call g2dlevalp_vec(nd,delta,centers(1,ibox),
     1             nlocal,rmlexp(iaddr(5,ibox)),sourcesort(1,istart),
     2             npts,pot(1,istart))
cccc               call g2dsevalp_vec(nd,delta,centers(1,ibox),nsoe,ws,ts,
cccc     1             rmlexp(iaddr(3,ibox)),sourcesort(1,istart),
cccc     2             npts,pot(1,istart))
cccc               call g2dsxevalp_vec(nd,delta,centers(1,ibox),pmax,npw,
cccc     1             nsoehalf,rmlexp(iaddr(4,ibox)),sourcesort(1,istart),
cccc     2             npts,ws,ts,pot(1,istart))
            endif
            if(ifpgh.eq.2) then
               call g2dlevalg_vec(nd,delta,centers(1,ibox),
     1             nlocal,rmlexp(iaddr(5,ibox)),sourcesort(1,istart),
     2             npts,pot(1,istart),grad(1,1,istart))
            endif
            if(ifpgh.eq.3) then
               call g2dlevalh_vec(nd,delta,centers(1,ibox),
     1             nlocal,rmlexp(iaddr(5,ibox)),sourcesort(1,istart),
     2             npts,pot(1,istart),grad(1,1,istart),
     3             hess(1,1,istart))
            endif
          endif
        enddo
C$OMP END PARALLEL DO        
      enddo

      call cpu_time(time2)
C$    time2 = omp_get_wtime()      
      timeinfo(5) = time2 - time1

      if(ifprint .ge. 1)
     $     call prinf('=== STEP 6 (direct) =====*',i,0)

c
cc
      call cpu_time(time1)
C$    time1=omp_get_wtime()  
      do ilev = 0,0
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
cccc               ibox=jbox
                  istarts = isrcse(1,ibox)
                  iends = isrcse(2,ibox)
                  nptssrc = iends-istarts + 1
                  
                  istartt = itargse(1,ibox)
                  iendt = itargse(2,ibox)
                  nptstarg = iendt-istartt + 1
ccccc
ccccc     use Hermite expansion for self interactions
ccccc     
cccc            call g2dhermevalp_vec(nd,delta,centers(1,ibox),
cccc     1          ntermsh,rmlexp(iaddr(6,ibox)),targetsort(1,istartt),
cccc     2          nptstarg,pottarg(1,istartt))
cccc
cccc            call g2dhermevalp_vec(nd,delta,centers(1,ibox),
cccc     1          ntermsh,rmlexp(iaddr(6,ibox)),sourcesort(1,istarts),
cccc     2          nptssrc,pot(1,istarts)) 
c
            
                  call fgt2dpart_direct_vec(nd,delta,eps/90,jstart,jend,
     1                istartt,iendt,sourcesort,ifcharge,chargesort,
     2                ifdipole,rnormalsort,dipstrsort,targetsort,
     3                ifpghtarg,pottarg,gradtarg,hesstarg)
         
                  call fgt2dpart_direct_vec(nd,delta,eps/90,jstart,jend,
     1                istarts,iends,sourcesort,ifcharge,chargesort,
     2                ifdipole,rnormalsort,dipstrsort,sourcesort,
     3                ifpgh,pot,grad,hess)
               enddo
            endif
         enddo
C$OMP END PARALLEL DO         
      enddo

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
c------------------------------------------------------------------    
      subroutine g2dmpalloc(nd,laddr,iaddr,nlevels,lmptot,
     1                          ntermsh,nlocal,npw,nsoe)
c     This subroutine determines the size of the array
c     to be allocated for the multipole expansions
c     iaddr(1,i) points to the starting location of the multipole
c     expansion of box i and iaddr(2,i) points to the local
c     expansion of box i
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
c     nterms      in: Integer(0:nlevels)
c                 Number of terms requried in expansions at each
c                 level
c
c------------------------------------------------------------------
c     Output arguments
c     iaddr       out: Integer(2,nboxes)
c                 Points the multipole and local expansions in box i
c 
c     lmptot      out: Integer
c                 Total length of expansions array required
c------------------------------------------------------------------

      implicit none
      integer nlevels,ntermsh,nlocal,npw,nsoe,nd,nsig,nt1,nt2,next235
      integer iaddr(6,*), lmptot, laddr(2,0:nlevels)
      integer ibox,i,iptr,istart,nn,itmp
      real *8 ddn
c
      istart = 1
      do i = 0,nlevels

         nn = (4*nsoe*nsoe/2)*2*nd
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the multipole expansion         
c
           itmp = ibox - laddr(1,i)
           iaddr(1,ibox) = istart + itmp*nn
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
      enddo
c
      do i = 0,nlevels

         nn = (4*(2*npw+1)*nsoe/2)*2*nd
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the multipole expansion         
c
           itmp = ibox - laddr(1,i)
           iaddr(2,ibox) = istart + itmp*nn
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
      enddo
c
c     Allocate memory for the local SOE and SOE/Hermite expansions
c
       do i=0,nlevels
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
       do i=0,nlevels
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
       do i=0,nlevels
         nn = (nlocal+1)*(nlocal+1)*nd
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
c     Allocate memory for the local Taylor expansions
c
       do i=0,nlevels
         nn = (ntermsh+1)*(ntermsh+1)*nd
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
