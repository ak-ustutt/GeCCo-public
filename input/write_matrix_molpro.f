      subroutine write_matrix_molpro(xmat,fout,mat_name,mat_type,
     &                               new,isym,nblock,nsym)
*
*     store a matrix in molpro's matrix format
*
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'multd2h.h'

      type(filinf), intent(inout) :: fout
      character(len=*), intent(in) :: mat_name, mat_type
      real(8), intent(in) :: xmat(*)
      logical, intent(in) :: new
      integer, intent(in) :: isym, nsym, nblock(nsym)

      integer :: lout, idxoff, csym, rsym, ncol, nrow, icol

      call file_open(fout)
      lout = fout%unit

      if (.not.new) then
         do
             read(lout,*,end=12)
         end do
  12     continue
      end if

      write(lout,"(a)") "BEGIN_DATA,"
      write(lout,"(a,a,4x,a,a,i1)") "# MATRIX ",trim(mat_name),
     &                           trim(mat_type)," SYMMETRY=",isym 

      idxoff = 0
      do csym = 1, nsym
         ncol = nblock(csym)
         rsym = multd2h(csym,isym)
         nrow = nblock(rsym)
         do icol = 1, ncol
            write(lout,"(5(e19.11,','))") xmat(idxoff+1:idxoff+nrow)
            idxoff = idxoff + nrow
         end do
      end do

      write(lout,"(a)") "END_DATA,"

      call file_close_keep(fout)

      return

      end
