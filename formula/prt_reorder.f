*----------------------------------------------------------------------*
      subroutine prt_reorder(luout,reo)
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      integer, intent(in) ::
     &     luout
      type(reorder) ::
     &     reo

      integer ::
     &     ij

      if (.not.reo%tra_in) write(luout,'(x,a,i4)')
     &     trim(reo%label_in),reo%iblk_in
      if (     reo%tra_in) write(luout,'(x,a,i4)')
     &     trim(reo%label_in)//'^+',reo%iblk_in

      do ij = 1, reo%nj_in
        call wrt_occ_rstr(luout,ij,
     &       reo%occ_opin(1:,1:,ij),
     &       reo%rst_opin(1:,1:,1:,1:,1:,ij),reo%ngas,reo%nspin)
      end do

      write(luout,*) 'reorder: sign = ',reo%sign
      write(luout,'(3x,8(i3,"->",i3))') reo%from_to(1:2,1:reo%nreo)
      call wrt_occ_n(luout,reo%occ_shift,reo%nreo)

      if (.not.reo%tra_out) write(luout,'(x,"==> ",a,i4)')
     &     trim(reo%label_out),reo%iblk_out
      if (     reo%tra_out) write(luout,'(x,"==> ",a,i4)')
     &     trim(reo%label_out)//'^+',reo%iblk_out

      do ij = 1, reo%nj_out
        call wrt_occ_rstr(luout,ij,
     &       reo%occ_opout(1:,1:,ij),
     &       reo%rst_opout(1:,1:,1:,1:,1:,ij),reo%ngas,reo%nspin)
      end do      

      return
      end
