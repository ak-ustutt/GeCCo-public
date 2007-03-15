      interface
      integer function mem_alloc_int(iarr,nalloc,name)
c      use memman
      implicit none
      integer, pointer, intent(out) ::
     &     iarr(:)
      integer, intent(in) ::
     &     nalloc
      character, intent(in) ::
     &     name*(*)
      end function

      integer function mem_alloc_real(xarr,nalloc,name)
c      use memman
      implicit none
      real(8), pointer, intent(out) ::
     &     xarr(:)
      integer, intent(in) ::
     &     nalloc
      character, intent(in) ::
     &     name*(*)
      end function

      subroutine mem_clean()
c      use memman
      implicit none
      end subroutine

      integer function mem_dealloc(name)
c      use memman
      implicit none      
      character, intent(in), optional ::
     &     name*(*)
      end function

      integer function mem_flushmark(name)
c      use memman
      implicit none      
      character, intent(in), optional ::
     &     name*(*)
      end function

      subroutine mem_init(mem_free_init)
c      use memman
      implicit none
      integer, intent(in) ::
     &     mem_free_init
      end subroutine

      subroutine mem_map(check)
c      use memman
      implicit none
      include 'stdunit.h'
      logical, intent(in) ::
     &     check
      end subroutine

      integer function mem_register(nalloc,name)
c      use memman
      implicit none
      integer, intent(in) ::
     &     nalloc
      character, intent(in) ::
     &     name*(*)
      end function

      integer function mem_setmark(name)
c      use memman
      implicit none      
      character, intent(in) ::
     &     name*(*)
      end function

      end interface
