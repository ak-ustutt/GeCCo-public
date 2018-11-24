*----------------------------------------------------------------------*
      subroutine rw_bcontr_kernel(irw,lu,bcontr)
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
      type(binary_contr), intent(inout) ::
     &     bcontr

      integer ::
     &     n_operands, nj_cnt, nj_res, nj_op1, nj_op2,
     &     ngas, nspin, len_m1, len_m2, len_m12, len_m21,
     &     lenlab1, lenlab2, lenlab3

      integer, external ::
     &     len_merge_map


      if (irw.gt.0) then

        ! dimension record
        read(lu,end=100) n_operands, nj_cnt, nj_res, nj_op1, nj_op2,
     &       ngas, nspin, len_m1, len_m2, len_m12, len_m21

        allocate(bcontr%occ_res(ngastp,2,nj_res))
        allocate(bcontr%rst_res(2,ngas,2,2,nspin,nj_res))
        allocate(bcontr%occ_op1(ngastp,2,nj_op1))
        allocate(bcontr%rst_op1(2,ngas,2,2,nspin,nj_op1))
        allocate(bcontr%occ_ex1(ngastp,2,nj_op1))
        allocate(bcontr%rst_ex1(2,ngas,2,2,nspin,nj_op1))
        if (n_operands.eq.2) then
          allocate(bcontr%occ_op2(ngastp,2,nj_op2))
          allocate(bcontr%rst_op2(2,ngas,2,2,nspin,nj_op2))
          allocate(bcontr%occ_ex2(ngastp,2,nj_op2))
          allocate(bcontr%rst_ex2(2,ngas,2,2,nspin,nj_op2))
        end if
        if (nj_cnt.gt.0) then
          allocate(bcontr%occ_cnt(ngastp,2,nj_cnt))
          allocate(bcontr%rst_cnt(2,ngas,2,2,nspin,nj_cnt))
          allocate(bcontr%merge_op1(len_m1))
          if (n_operands.eq.2) then
            allocate(bcontr%merge_op2(len_m2))
          end if
          allocate(bcontr%merge_op1op2(len_m12))
          allocate(bcontr%merge_op2op1(len_m21))
        end if

        bcontr%n_operands = n_operands
        bcontr%n_cnt = nj_cnt
        bcontr%ngas = ngas
        bcontr%nspin = nspin
        bcontr%label_res(1:maxlen_bc_label) = ' '
        bcontr%label_op1(1:maxlen_bc_label) = ' '
        bcontr%label_op2(1:maxlen_bc_label) = ' '
        bcontr%nj_res = nj_res
        bcontr%nj_op1 = nj_op1
        bcontr%nj_op2 = nj_op2
        
        read(lu,end=100) bcontr%fact,bcontr%fact_itf,
     &       lenlab1,bcontr%label_res(1:lenlab1),
     &       lenlab2,bcontr%label_op1(1:lenlab2),
     &       lenlab3,bcontr%label_op2(1:lenlab3),
     &       bcontr%iblk_res,bcontr%iblk_op1,bcontr%iblk_op2,
     &       bcontr%tra_res,bcontr%tra_op1,bcontr%tra_op2,
     &       bcontr%occ_res,bcontr%rst_res,
     &       bcontr%occ_op1,bcontr%rst_op1

        ! extra records
        if (n_operands.eq.2)
     &       read(lu,end=100) bcontr%occ_op2,bcontr%rst_op2
        if (nj_cnt.gt.0)
     &       read(lu,end=100)
     &       bcontr%occ_ex1,bcontr%occ_cnt,
     &       bcontr%rst_ex1,bcontr%rst_cnt,
     &       bcontr%merge_op1(1:len_m1),
     &       bcontr%merge_op1op2(1:len_m12),
     &       bcontr%merge_op2op1(1:len_m21)
        if (nj_cnt.gt.0.and.n_operands.eq.2)
     &       read(lu,end=100)
     &       bcontr%occ_ex2,bcontr%rst_ex2,bcontr%merge_op2(1:len_m2)

      else

        n_operands  = bcontr%n_operands
        nj_cnt      = bcontr%n_cnt
        ngas        = bcontr%ngas
        nspin       = bcontr%nspin
        nj_res      = bcontr%nj_res
        nj_op1      = bcontr%nj_op1
        nj_op2      = bcontr%nj_op2

        if (nj_cnt.gt.0) then
          len_m1 = len_merge_map(bcontr%merge_op1,nj_op1)
          if (n_operands.eq.2) then
            len_m2 = len_merge_map(bcontr%merge_op2,nj_op2)
          else
            len_m2 = 0
          end if
          len_m12 = len_merge_map(bcontr%merge_op1op2,nj_res)
          len_m21 = len_merge_map(bcontr%merge_op2op1,nj_res)
        else
          len_m1 = 0
          len_m2 = 0
          len_m12 = 0
          len_m21 = 0
        end if

        ! dimension record
        write(lu,err=200) n_operands, nj_cnt, nj_res, nj_op1, nj_op2,
     &       ngas, nspin, len_m1, len_m2, len_m12, len_m21

        lenlab1 = len_trim(bcontr%label_res)
        lenlab2 = len_trim(bcontr%label_op1)
        lenlab3 = len_trim(bcontr%label_op2)
        write(lu,err=200) bcontr%fact,bcontr%fact_itf,
     &       lenlab1,bcontr%label_res(1:lenlab1),
     &       lenlab2,bcontr%label_op1(1:lenlab2),
     &       lenlab3,bcontr%label_op2(1:lenlab3),
     &       bcontr%iblk_res,bcontr%iblk_op1,bcontr%iblk_op2,
     &       bcontr%tra_res,bcontr%tra_op1,bcontr%tra_op2,
     &       bcontr%occ_res,bcontr%rst_res,
     &       bcontr%occ_op1,bcontr%rst_op1
        
        ! extra records
        if (n_operands.eq.2)
     &       write(lu,err=200) bcontr%occ_op2,bcontr%rst_op2
        if (nj_cnt.gt.0)
     &       write(lu,err=200)
     &       bcontr%occ_ex1,bcontr%occ_cnt,
     &       bcontr%rst_ex1,bcontr%rst_cnt,
     &       bcontr%merge_op1(1:len_m1),
     &       bcontr%merge_op1op2(1:len_m12),
     &       bcontr%merge_op2op1(1:len_m21)
        if (nj_cnt.gt.0.and.n_operands.eq.2)
     &       write(lu,err=200)
     &       bcontr%occ_ex2,bcontr%rst_ex2,bcontr%merge_op2(1:len_m2)

      end if

      return

      ! error handling
      ! unexpected EOF on read
 100  call quit(1,'rw_bcontr_kernel',
     &     'unexpected EOF while reading formula file')

      ! write errors:
 200  call quit(1,'rw_bcontr_kernel','error writing contraction')

      end
