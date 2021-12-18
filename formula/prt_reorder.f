*----------------------------------------------------------------------*
      subroutine prt_reorder(lulog,reo)
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      integer, intent(in) ::
     &     lulog
      type(reorder) ::
     &     reo

      integer ::
     &     ij

      if (.not.reo%tra_in) write(lulog,'(x,a,i4)')
     &     trim(reo%label_in),reo%iblk_in
      if (     reo%tra_in) write(lulog,'(x,a,i4)')
     &     trim(reo%label_in)//'^+',reo%iblk_in

      do ij = 1, reo%nj_in
        call wrt_occ_rstr(lulog,ij,
     &       reo%occ_opin(1:,1:,ij),
     &       reo%rst_opin(1:,1:,1:,1:,1:,ij),reo%ngas,reo%nspin)
      end do

      write(lulog,*) 'reorder: sign = ',reo%sign
      write(lulog,'(3x,8(i3,"->",i3))') reo%from_to(1:2,1:reo%nreo)
      call wrt_occ_n(lulog,reo%occ_shift,reo%nreo)
      if (reo%nreo_i0.gt.0) write(lulog,'(x,a,3x,8(i3,"->",i3))')
     &       'merge: ',reo%from_to(1:2,reo%nreo+1:reo%nreo+reo%nreo_i0)

      if (.not.reo%tra_out) write(lulog,'(x,"==> ",a,i4)')
     &     trim(reo%label_out),reo%iblk_out
      if (     reo%tra_out) write(lulog,'(x,"==> ",a,i4)')
     &     trim(reo%label_out)//'^+',reo%iblk_out

      do ij = 1, reo%nj_out
        call wrt_occ_rstr(lulog,ij,
     &       reo%occ_opout(1:,1:,ij),
     &       reo%rst_opout(1:,1:,1:,1:,1:,ij),reo%ngas,reo%nspin)
      end do      

      return
      end
