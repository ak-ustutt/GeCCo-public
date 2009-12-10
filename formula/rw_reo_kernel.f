*----------------------------------------------------------------------*
      subroutine rw_reo_kernel(irw,lu,reo)
*----------------------------------------------------------------------*
*     read/write contraction from/to lu 
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'

      integer, intent(in) ::
     &     lu, irw
      type(reorder), intent(inout) ::
     &     reo

      integer ::
     &     nreo, nj_in, nj_out,
     &     ngas, nspin, 
     &     lenlab1, lenlab2, len_m1, len_m1i, len_m2, len_m2i

      integer, external ::
     &     len_merge_map

      if (irw.gt.0) then

        ! dimension record
        read(lu,end=100) nj_out, nj_in, nreo, ngas, nspin,
     &       len_m1, len_m1i, len_m2, len_m2i

        allocate(reo%occ_opin(ngastp,2,nj_in))
        allocate(reo%rst_opin(2,ngas,2,2,nspin,nj_in))
        allocate(reo%occ_opout(ngastp,2,nj_out))
        allocate(reo%rst_opout(2,ngas,2,2,nspin,nj_out))
        allocate(reo%from_to(2,nreo))
        allocate(reo%occ_shift(ngastp,2,nreo))
        allocate(reo%occ_op0(ngastp,2,nj_in))
        allocate(reo%merge_stp1(len_m1))
        allocate(reo%merge_stp1inv(len_m1i))
        allocate(reo%merge_stp2(len_m2))
        allocate(reo%merge_stp2inv(len_m2i))

        reo%nj_out = nj_out
        reo%nj_in  = nj_in
        reo%nreo   = nreo
        reo%ngas = ngas
        reo%nspin = nspin
        reo%label_in (1:maxlen_bc_label) = ' '
        reo%label_out(1:maxlen_bc_label) = ' '
        
        read(lu,end=100) 
     &       lenlab1,reo%label_in (1:lenlab1),
     &       lenlab2,reo%label_out(1:lenlab2),
     &       reo%tra_in,reo%tra_out,
     &       reo%sign,
     &       reo%iblk_in,reo%iblk_out,
     &       reo%occ_opin, reo%rst_opin,
     &       reo%occ_opout,reo%rst_opout,
     &       reo%occ_op0, reo%occ_shift, reo%from_to,
     &       reo%merge_stp1,reo%merge_stp1inv,
     &       reo%merge_stp2,reo%merge_stp2inv

c dbg
c        print *,'merge maps read:'
c        print *,' 1:  ',reo%merge_stp1
c        print *,' 1i: ',reo%merge_stp1inv
c        print *,' 2:  ',reo%merge_stp2
c        print *,' 2i: ',reo%merge_stp2inv
c dbg
      else

        ngas        = reo%ngas
        nspin       = reo%nspin
        nj_in       = reo%nj_in
        nj_out      = reo%nj_out
        nreo        = reo%nreo

        len_m1  = len_merge_map(reo%merge_stp1,nj_in)
        len_m1i = len_merge_map(reo%merge_stp1inv,nj_in)
        len_m2  = len_merge_map(reo%merge_stp2,nj_out)
        len_m2i = len_merge_map(reo%merge_stp2inv,nj_out)

        ! dimension record
        write(lu,err=200) nj_out, nj_in, nreo, ngas, nspin,
     &       len_m1, len_m1i, len_m2, len_m2i

        lenlab1 = len_trim(reo%label_in)
        lenlab2 = len_trim(reo%label_out)
        write(lu,err=200) 
     &       lenlab1,reo%label_in (1:lenlab1),
     &       lenlab2,reo%label_out(1:lenlab2),
     &       reo%tra_in,reo%tra_out,
     &       reo%sign,
     &       reo%iblk_in,reo%iblk_out,
     &       reo%occ_opin, reo%rst_opin,
     &       reo%occ_opout,reo%rst_opout,
     &       reo%occ_op0, reo%occ_shift, reo%from_to,
     &       reo%merge_stp1,reo%merge_stp1inv,
     &       reo%merge_stp2,reo%merge_stp2inv
        
c dbg
c        print *,'merge maps written:'
c        print *,' 1:  ',reo%merge_stp1
c        print *,' 1i: ',reo%merge_stp1inv
c        print *,' 2:  ',reo%merge_stp2
c        print *,' 2i: ',reo%merge_stp2inv
c dbg

      end if

      return

      ! error handling
      ! unexpected EOF on read
 100  call quit(1,'rw_reo_kernel',
     &     'unexpected EOF while reading formula file')

      ! write errors:
 200  call quit(1,'rw_reo_kernel','error writing contraction')

      end
