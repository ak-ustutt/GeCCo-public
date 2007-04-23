*----------------------------------------------------------------------*
      integer function mem_alloc_real(xarr,nalloc,name)
*----------------------------------------------------------------------*
*     allocate a real(8) array of length nalloc in the current section
*     xarr is a pointer to that array
*     make sure to include the interface file ifc_memman.h
*----------------------------------------------------------------------*
      use memman
      implicit none

      real(8), pointer ::
     &     xarr(:)
      integer, intent(in) ::
     &     nalloc
      character, intent(in) ::
     &     name*(*)

      mem_alloc_real = memman_alloc(mtyp_rl8,nalloc,name,xpnt=xarr)

      return
      end
