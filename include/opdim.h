*----------------------------------------------------------------------*
*     some global definitions regarding operator (occup.) dimensions
*----------------------------------------------------------------------*
      ! constants
      integer, parameter ::
     &     ihole = 1,        ! hole space is always 1
     &     ipart = 2,        ! particle space is always 2
     &     ngastp = 4,       ! number of GAS types (2,3,4)
     &     pack_base = 256   ! make 2*ngastp fit into integer(8)
c     &     pack_base = 10    ! better choice for debugging

      ! global variable
      integer ::
     &     ivale,          ! index of valence space (usually 3)
     &     iextr           ! index of external space (3 or 4)
      common /opdim/
     &     ivale,iextr

