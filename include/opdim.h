*----------------------------------------------------------------------*
*     some global definitions regarding operator (occup.) dimensions
*----------------------------------------------------------------------*
      ! constants
      integer, parameter ::
     &     ihole = 1,        ! hole space is always 1
     &     ipart = 2,        ! particle space is always 2
     &     ivale = 3,        ! valence space is always 3
     &     iextr = 4,        ! external space is always 4
     &     ngastp = 4,       ! number of GAS types (2,3,4)
c     &     pack_base = 256   ! make 2*ngastp fit into integer(8)
c     &     pack_base = 10    ! better choice for debugging
     &     pack_base = 10

