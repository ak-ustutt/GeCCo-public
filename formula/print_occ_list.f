      subroutine print_occ_list(lulog,olist)

      implicit none

      include 'opdim.h'
      include 'def_occ_list.h'

      type(occ_list), intent(in) ::
     &     olist
      integer, intent(in) ::
     &     lulog

      integer ::
     &     idx

      write(lulog,*) ' vtx1 vtx2 start  end'
      write(lulog,*) '----------------------'
      if (olist%n_vtx_inf.eq.0) then
        write(lulog,*) '  empty'
        return
      end if
      write(lulog,'(x,i4,x,i4,2x,i4,x,i4)')
     &     olist%vtx_inf(1:4,1:olist%n_vtx_inf)
        
      do idx = 1, olist%n_occ
        call wrt_occ(lulog,olist%occ(:,:,idx))
      end do

      return
      end
