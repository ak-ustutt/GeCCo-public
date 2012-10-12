*----------------------------------------------------------------------*
      subroutine rw_opdef_kernel(irw,lu,op,incore,
     &                           lab_p1,lab_p2,tra,tra1,tra2)
*----------------------------------------------------------------------*
*     read/write contraction from/to lu 
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_operator.h'
      include 'ifc_memman.h'
      include 'par_globalmarks.h'

      integer, intent(in) ::
     &     lu, irw
      type(operator), intent(inout) ::
     &     op
      character(len=*), intent(inout) ::
     &     lab_p1, lab_p2
      logical, intent(inout) ::
     &     tra, tra1, tra2
      integer, intent(inout) ::
     &     incore

      integer ::
     &     nblk, nj, ngas, nspin, lenlab, lenp1lab, lenp2lab, ifree
      

      if (irw.gt.0) then

        ! dimension record
        read(lu,end=100) nblk, nj, ngas, nspin

        ! allocate in operator section:
        call mem_pushmark()
        ifree = mem_gotomark(operator_def)
        allocate(op%ihpvca_occ(ngastp,2,nblk*nj),
     &           op%ica_occ(2,nblk),
     &           op%igasca_restr(2,ngas,2,2,nspin,nblk*nj),
     &           op%formal_blk(nblk))

        ! init_operator_0 has been called before ...
        op%ngas = ngas
        op%nspin = nspin
        op%n_occ_cls = nblk
        op%njoined = nj
        
        read(lu,end=100) lenlab,op%name(1:lenlab),tra,
     &       lenp1lab,lab_p1(1:lenp1lab),tra1,
     &       lenp2lab,lab_p2(1:lenp2lab),tra2,
     &       op%ihpvca_occ,op%igasca_restr,incore

        ifree = mem_register(2*ngastp*nblk*nj
     &                      +2*nblk
     &                      +8*ngas*nblk*nj*nspin
     &                      +nblk,
     &       trim(op%name))

        call occ2caocc(op%ica_occ,op%ihpvca_occ,nj,nblk)

        op%formal_blk = .false.

        call mem_popmark()

      else

        write(lu,err=200) op%n_occ_cls, op%njoined, op%ngas, op%nspin

        lenlab = len_trim(op%name)
        lenp1lab = len_trim(lab_p1)
        lenp2lab = len_trim(lab_p2)
        write(lu,err=200) lenlab,op%name(1:lenlab),tra,
     &       lenp1lab,lab_p1(1:lenp1lab),tra1,
     &       lenp2lab,lab_p2(1:lenp2lab),tra2,
     &       op%ihpvca_occ,op%igasca_restr,incore

      end if

      return

      ! error handling
      ! unexpected EOF on read
 100  call quit(1,'rw_opdef_kernel',
     &     'unexpected EOF while reading formula file')

      ! write errors:
 200  call quit(1,'rw_opdef_kernel','error writing contraction')

      end
