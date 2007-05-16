      interface
      integer function mem_alloc_int(iarr,nalloc,name)
      implicit none
      integer, pointer ::
     &     iarr(:)
      integer, intent(in) ::
     &     nalloc
      character, intent(in) ::
     &     name*(*)
      end function

      integer function mem_alloc_real(xarr,nalloc,name)
      implicit none
      real(8), pointer ::
     &     xarr(:)
      integer, intent(in) ::
     &     nalloc
      character, intent(in) ::
     &     name*(*)
      end function

      subroutine mem_clean()
      implicit none
      end subroutine

      integer function mem_dealloc(name)
      implicit none      
      character, intent(in), optional ::
     &     name*(*)
      end function

      integer function mem_flushmark(name)
      implicit none      
      character, intent(in), optional ::
     &     name*(*)
      end function

      subroutine mem_init(mem_free_init)
      implicit none
      integer, intent(in) ::
     &     mem_free_init
      end subroutine

      subroutine mem_map(check)
      implicit none
      logical, intent(in) ::
     &     check
      end subroutine

      subroutine mem_check(label)
      implicit none
      character, intent(in) ::
     &     label*(*)
      end subroutine

      integer function mem_register(nalloc,name)
      implicit none
      integer, intent(in) ::
     &     nalloc
      character, intent(in) ::
     &     name*(*)
      end function

      integer function mem_setmark(name)
      implicit none      
      character, intent(in) ::
     &     name*(*)
      end function

      integer function mem_gotomark(name)
      implicit none
      character, intent(in) ::
     &     name*(*)
      end function

      integer function mem_gotolastmark()
      implicit none
      end function

      subroutine mem_popmark()
      end subroutine

      subroutine mem_pushmark()
      end subroutine

      end interface
