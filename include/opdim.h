*----------------------------------------------------------------------*
*     some global definitions regarding operator (occup.) dimensions
*----------------------------------------------------------------------*
      ! constants
      integer, parameter ::
     &     ihole = 1,        ! hole space is always 1
     &     ipart = 2         ! particle space is always 2

      ! global variable
      integer ::
     &     ngastp,         ! number of GAS types (2,3,4)
     &     ivale,          ! index of valence space (usually 3)
     &     iextr           ! index of external space (3 or 4)
      common /opdim/
     &     ngastp,ivale,iextr

