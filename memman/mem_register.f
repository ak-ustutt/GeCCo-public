      integer function mem_register(nalloc,name)
      use memman
      implicit none

      integer, intent(in) ::
     &     nalloc
      character, intent(in) ::
     &     name*(*)

      mem_register = memman_alloc(mtyp_reg,nalloc,name)

      return
      end
