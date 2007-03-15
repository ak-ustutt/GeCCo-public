      integer function mem_alloc_real(xarr,nalloc,name)
      use memman
      implicit none

      real(8), pointer, intent(out) ::
     &     xarr(:)
      integer, intent(in) ::
     &     nalloc
      character, intent(in) ::
     &     name*(*)

      mem_alloc_real = memman_alloc(mtyp_rl8,nalloc,name,xpnt=xarr)

      return
      end
