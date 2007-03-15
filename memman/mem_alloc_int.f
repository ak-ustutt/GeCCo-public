      integer function mem_alloc_int(iarr,nalloc,name)
      use memman
      implicit none

      integer, pointer, intent(out) ::
     &     iarr(:)
      integer, intent(in) ::
     &     nalloc
      character, intent(in) ::
     &     name*(*)

      mem_alloc_int = memman_alloc(mtyp_int,nalloc,name,ipnt=iarr)

      return
      end
