*----------------------------------------------------------------------*
      logical function next_msgamdist(first,
     &     ms_a,ms_c,igam_a,igam_c,iocc,nsym,
     &     ms_dist,igam_dist)
*----------------------------------------------------------------------*
*     increment an entire set of Ms and IRREP distributions
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'multd2h.h'
      
      logical, intent(in) ::
     &     first
      integer, intent(in) ::
     &     ms_a, ms_c, igam_a, igam_c,
     &     iocc(ngastp,2), nsym
      integer, intent(inout) ::
     &     ms_dist(ngastp,2), igam_dist(ngastp,2)
      logical, external ::
     &     first_msdist4, first_gamdist4,
     &     next_msdist4, next_gamdist4

c dbg
c      print *,'First',first
c dbg

      if (first) then
        ! initialize
        next_msgamdist = first_msdist4(ms_dist(1,1),ms_c,iocc(1,1))
        next_msgamdist = next_msgamdist.and.
     &      first_msdist4(ms_dist(1,2),ms_a,iocc(1,2))
        next_msgamdist = next_msgamdist.and.        
     &      first_gamdist4(igam_dist(1,1),igam_c,nsym,iocc(1,1))
        next_msgamdist = next_msgamdist.and.        
     &      first_gamdist4(igam_dist(1,2),igam_a,nsym,iocc(1,2))
      else          
        ! innermost index: increment IRREP distribution of C
        if (next_gamdist4(igam_dist(1,1),igam_c,nsym,iocc(1,1))) then
          next_msgamdist = .true.

        ! increment Ms distribution of C
        else if (next_msdist4(ms_dist(1,1),ms_c,iocc(1,1))) then
          ! it's the same call to first_gamdist as above, so
          ! actually the result *must* be true (not checked)
          next_msgamdist =
     &         first_gamdist4(igam_dist(1,1),igam_c,nsym,iocc(1,1))
        ! increment IRREP distribution of A
        else if (next_gamdist4(igam_dist(1,2),igam_a,nsym,iocc(1,2)))
     &         then
          next_msgamdist = first_msdist4(ms_dist(1,1),ms_c,iocc(1,1))
          next_msgamdist = next_msgamdist.and.
     &         first_gamdist4(igam_dist(1,1),igam_c,nsym,iocc(1,1))

        ! increment Ms distribution of A
        else if (next_msdist4(ms_dist(1,2),ms_a,iocc(1,2))) then
          next_msgamdist =
     &         first_gamdist4(igam_dist(1,2),igam_a,nsym,iocc(1,2))
          next_msgamdist = next_msgamdist.and.
     &         first_msdist4(ms_dist(1,1),ms_c,iocc(1,1))
          next_msgamdist = next_msgamdist.and.
     &         first_gamdist4(igam_dist(1,1),igam_c,nsym,iocc(1,1))

        else
          next_msgamdist = .false.
        end if

      end if

      return
      end
*----------------------------------------------------------------------*
*     a few aux-routines follow:
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      logical function next_msdist4(msdis,mst,ihpvdis)
*----------------------------------------------------------------------*
*     function to generate MS distributions (in the 
*     sequence (X,P,H,V), i.e. X running fastest; note that msdis is
*     stored as (H,P,V,X))
*     fixed to 4 space types 
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     ihpvdis(4), mst
      integer, intent(inout) ::
     &     msdis(4)

      next_msdist4 = .true.
      do
        if (msdis(4).gt.-ihpvdis(4)) then
          ! decrement Ms of X:
          msdis(4) = msdis(4)-2
        else if (msdis(2).gt.-ihpvdis(2)) then
          ! decrement Ms of P: 
          msdis(2) = msdis(2)-2
          ! start X with maximum Ms
          msdis(4) = ihpvdis(4)
        else if (msdis(3).gt.-ihpvdis(3)) then
          ! else enter a new Ms block given by V
          msdis(3) = msdis(3)-2
          ! X, P start with maximum Ms
          msdis(2) = ihpvdis(2)
          msdis(4) = ihpvdis(4)
        else
          ! last Ms block of V over --> ende gelaende
          next_msdist4 = .false.
          exit
        end if
        ! Ms value of H block is dictated by overall Ms
        msdis(1) = mst-msdis(2)-msdis(3)-msdis(4)
        ! if this is a possible value, we are done ...
        if (abs(msdis(1)).le.ihpvdis(1)) exit
      end do

      return
      end
*----------------------------------------------------------------------*
      logical function first_msdist4(msdis,mst,ihpvdis)
*----------------------------------------------------------------------*
*     functions to generate first MS distribution 
*     fixed to 4 space types 
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     ihpvdis(4), mst
      integer, intent(inout) ::
     &     msdis(4)

c dbg
c      print *,'First msgam'
c dbg

      first_msdist4 = .true.
      ! start with max value for X, P and V
      msdis(2) = ihpvdis(2)
      msdis(3) = ihpvdis(3)
      msdis(4) = ihpvdis(4)
      msdis(1) = mst-msdis(2)-msdis(3)-msdis(4)
c dbg
c      print *,'msdis(4),ihpvdis(4)',msdis(4),ihpvdis(4)
c      print *,'msdis(3),ihpvdis(3)',msdis(3),ihpvdis(3)
c      print *,'msdis(2),ihpvdis(2)',msdis(2),ihpvdis(2)
c      print *,'msdis(1),ihpvdis(1)',msdis(1),ihpvdis(1)
c dbg
      ! if this is an acceptable value for msdis --> exit
      if (abs(msdis(1)).le.ihpvdis(1)) return
      
      ! else, run through possible combinations:
      do
        if (msdis(4).gt.-ihpvdis(4)) then
          ! decrement Ms of X: 
          msdis(4) = msdis(4)-2
        else if (msdis(2).gt.-ihpvdis(2)) then
          ! decrement Ms of P: 
          msdis(2) = msdis(2)-2
          ! max. value for X:
          msdis(4) = ihpvdis(4)
        else if (msdis(3).gt.-ihpvdis(3)) then
          ! else enter a new Ms block given by V
          msdis(3) = msdis(3)-2
          ! X, P start with maximum Ms
          msdis(2) = ihpvdis(2)
          msdis(4) = ihpvdis(4)
        else
          ! last Ms block of V over --> ende gelaende
          first_msdist4 = .false.
          exit
        end if
        ! Ms value of H block is dictated by overall Ms
        msdis(1) = mst-msdis(2)-msdis(3)-msdis(4)
        ! if this is a possible value, we are done ...
        if (abs(msdis(1)).le.ihpvdis(1)) exit
      end do

      return
      end

*----------------------------------------------------------------------*
      logical function next_gamdist4(igamdis,igamt,nsym,ihpvdis)
*----------------------------------------------------------------------*
*     functions to generate IRREP distributions (in the 
*     sequence (X,P,H,V), i.e. X running fastest; NOTE that ihpvdis
*     is stored in the sequence (H,P,V,X) )
*     fixed to 4 space types 
*----------------------------------------------------------------------*
      implicit none

      include 'multd2h.h'

      integer, intent(in) ::
     &     igamt, nsym, ihpvdis(4)
      integer, intent(inout) ::
     &     igamdis(4)

      integer ::
     &     igam234

      next_gamdist4 = .true.
      do
        if (igamdis(4).lt.nsym.and.ihpvdis(4).gt.0) then
          ! increment IRREP for X
          igamdis(4) = igamdis(4)+1
          igam234 = multd2h(igamdis(2),igamdis(3))
          igam234 = multd2h(igam234,igamdis(4))
        else if (igamdis(2).lt.nsym.and.ihpvdis(2).gt.0) then
          ! increment IRREP for P
          igamdis(2) = igamdis(2)+1
          ! reset X to the lowest value
          igamdis(4) = 1
          igam234 = multd2h(igamdis(2),igamdis(3))
        else if (igamdis(3).lt.nsym.and.ihpvdis(3).gt.0) then
          ! else, we increment IRREP for V ...
          igamdis(3) = igamdis(3)+1
          ! and reset X and P to the lowest value
          igamdis(2) = 1
          igamdis(4) = 1
          igam234 = igamdis(3)
        else
          next_gamdist4 = .false.
          exit
        end if
        ! IRREP of H is dictated by overall symmetry
        igamdis(1) = multd2h(igamt,igam234)
        ! allow only if either totally symmetric or occupied:
        if (igamdis(1).eq.1.or.ihpvdis(1).gt.0) exit
      end do

      return
      end

*----------------------------------------------------------------------*
      logical function first_gamdist4(igamdis,igamt,nsym,ihpvdis)
*----------------------------------------------------------------------*
*     functions to generate first IRREP distributions 
*     fixed to 4 space types 
*----------------------------------------------------------------------*
      implicit none

      include 'multd2h.h'

      integer, intent(in) ::
     &     igamt, nsym, ihpvdis(4)
      integer, intent(inout) ::
     &     igamdis(4)

      integer ::
     &     igam234

      first_gamdist4 = .true.
      ! start out with IRREP 1 for P and V
      ! so IRREP for H must be igamt
      igamdis(1) = igamt 
      igamdis(2) = 1
      igamdis(3) = 1
      igamdis(4) = 1
      ! acceptable?
      if (igamdis(1).eq.1.or.ihpvdis(1).gt.0) return

      ! else, we have to increment until we find something
      do
        if (igamdis(4).lt.nsym.and.ihpvdis(4).gt.0) then
          ! increment IRREP for X
          igamdis(4) = igamdis(4)+1
          igam234 = multd2h(igamdis(2),igamdis(3))
          igam234 = multd2h(igam234,igamdis(4))
        else if (igamdis(2).lt.nsym.and.ihpvdis(2).gt.0) then
          ! increment IRREP for P
          igamdis(2) = igamdis(2)+1
          ! reset X
          igamdis(4) = 1
          igam234 = multd2h(igamdis(2),igamdis(3))
        else if (igamdis(3).lt.nsym.and.ihpvdis(3).gt.0) then
          ! else, we increment IRREP for V ...
          igamdis(3) = igamdis(3)+1
          ! and reset X and P to the lowest value
          igamdis(2) = 1
          igamdis(4) = 1
          igam234 = igamdis(3)
        else
          first_gamdist4 = .false.
          exit
        end if
        ! IRREP of H is dictated by overall symmetry
        igamdis(1) = multd2h(igamt,igam234)
        ! allow only if either totally symmetric or occupied:
        if (igamdis(1).eq.1.or.ihpvdis(1).gt.0) exit
      end do

      return
      end
*----------------------------------------------------------------------*
*     old 3-index routines follow:
*----------------------------------------------------------------------*
      logical function next_msdist3(msdis,mst,ihpvdis)
*----------------------------------------------------------------------*
*     function to generate MS distributions (in the 
*     sequence (P,H,V), i.e. P running fastest; note that msdis is
*     stored as (H,P,V))
*     fixed to 3 space types 
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     ihpvdis(3), mst
      integer, intent(inout) ::
     &     msdis(3)

      next_msdist3 = .true.
      do
        if (msdis(2).gt.-ihpvdis(2)) then
          ! decrement Ms of P: 
          msdis(2) = msdis(2)-2
        else if (msdis(3).gt.-ihpvdis(3)) then
          ! else enter a new Ms block given by V
          msdis(3) = msdis(3)-2
          ! P starts with maximum Ms
          msdis(2) = ihpvdis(2)
        else
          ! last Ms block of V over --> ende gelaende
          next_msdist3 = .false.
          exit
        end if
        ! Ms value of H block is dictated by overall Ms
        msdis(1) = mst-msdis(2)-msdis(3)
        ! if this is a possible value, we are done ...
        if (abs(msdis(1)).le.ihpvdis(1)) exit
      end do

      return
      end
*----------------------------------------------------------------------*
      logical function first_msdist3(msdis,mst,ihpvdis)
*----------------------------------------------------------------------*
*     functions to generate first MS distribution 
*     fixed to 3 space types 
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     ihpvdis(3), mst
      integer, intent(inout) ::
     &     msdis(3)

      first_msdist3 = .true.
      ! start with max value for P and V
      msdis(2) = ihpvdis(2)
      msdis(3) = ihpvdis(3)
      msdis(1) = mst-msdis(2)-msdis(3)
      ! if this is an acceptable value for msdis --> exit
      if (abs(msdis(1)).le.ihpvdis(1)) return
      
      ! else, run through possible combinations:
      do
        if (msdis(2).gt.-ihpvdis(2)) then
          ! decrement Ms of P: 
          msdis(2) = msdis(2)-2
        else if (msdis(3).gt.ihpvdis(3)) then
          ! else enter a new Ms block given by V
          msdis(3) = msdis(3)-2
          ! P starts with maximum Ms
          msdis(2) = ihpvdis(2)
        else
          ! last Ms block of V over --> ende gelaende
          first_msdist3 = .false.
          exit
        end if
        ! Ms value of H block is dictated by overall Ms
        msdis(1) = mst-msdis(2)-msdis(3)
        ! if this is a possible value, we are done ...
        if (abs(msdis(1)).le.ihpvdis(1)) exit
      end do

      return
      end

*----------------------------------------------------------------------*
      logical function next_gamdist3(igamdis,igamt,nsym,ihpvdis)
*----------------------------------------------------------------------*
*     functions to generate IRREP distributions (in the 
*     sequence (P,H,V), i.e. P running fastest; NOTE that ihpvdis
*     is stored in the sequence (H,P,V) )
*     fixed to 3 space types 
*----------------------------------------------------------------------*
      implicit none

      include 'multd2h.h'

      integer, intent(in) ::
     &     igamt, nsym, ihpvdis(3)
      integer, intent(inout) ::
     &     igamdis(3)

      integer ::
     &     igam23

      next_gamdist3 = .true.
      do
        if (igamdis(2).lt.nsym.and.ihpvdis(2).gt.0) then
          ! increment IRREP for P
          igamdis(2) = igamdis(2)+1
          igam23 = multd2h(igamdis(2),igamdis(3))
        else if (igamdis(3).lt.nsym.and.ihpvdis(3).gt.0) then
          ! else, we increment IRREP for V ...
          igamdis(3) = igamdis(3)+1
          ! and reset P to the lowest value
          igamdis(2) = 1
          igam23 = igamdis(3)
        else
          next_gamdist3 = .false.
          exit
        end if
        ! IRREP of H is dictated by overall symmetry
        igamdis(1) = multd2h(igamt,igam23)
        ! allow only if either totally symmetric or occupied:
        if (igamdis(1).eq.1.or.ihpvdis(1).gt.0) exit
      end do

      return
      end

*----------------------------------------------------------------------*
      logical function first_gamdist3(igamdis,igamt,nsym,ihpvdis)
*----------------------------------------------------------------------*
*     functions to generate first IRREP distributions 
*     fixed to 3 space types 
*----------------------------------------------------------------------*
      implicit none

      include 'multd2h.h'

      integer, intent(in) ::
     &     igamt, nsym, ihpvdis(3)
      integer, intent(inout) ::
     &     igamdis(3)

      integer ::
     &     igam23

      first_gamdist3 = .true.
      ! start out with IRREP 1 for P and V
      ! so IRREP for H must be igamt
      igamdis(1) = igamt 
      igamdis(2) = 1
      igamdis(3) = 1
      ! acceptable?
      if (igamdis(1).eq.1.or.ihpvdis(1).gt.0) return

      ! else, we have to increment until we find something
      do
        if (igamdis(2).lt.nsym.and.ihpvdis(2).gt.0) then
          ! increment IRREP for P
          igamdis(2) = igamdis(2)+1
          igam23 = multd2h(igamdis(2),igamdis(3))
        else if (igamdis(3).lt.nsym.and.ihpvdis(3).gt.0) then
          ! else, we increment IRREP for V ...
          igamdis(3) = igamdis(3)+1
          ! and reset P to the lowest value
          igamdis(2) = 1
          igam23 = igamdis(3)
        else
          first_gamdist3 = .false.
          exit
        end if
        ! IRREP of H is dictated by overall symmetry
        igamdis(1) = multd2h(igamt,igam23)
        ! allow only if either totally symmetric or occupied:
        if (igamdis(1).eq.1.or.ihpvdis(1).gt.0) exit
      end do

      return
      end
