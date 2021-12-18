*----------------------------------------------------------------------*
      integer function mem_alloc_int(iarr,nalloc,name)
*----------------------------------------------------------------------*
*     allocate an integer array of length nalloc in the current section
*     xarr is a pointer to that array
*     make sure to include the interface file ifc_memman.h
*----------------------------------------------------------------------*
      use memman
      implicit none

      integer, pointer ::
     &     iarr(:)
      integer, intent(in) ::
     &     nalloc
      character, intent(in) ::
     &     name*(*)

      mem_alloc_int = memman_alloc(mtyp_int,nalloc,name,ipnt=iarr)

      return
      end
