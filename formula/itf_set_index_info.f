      subroutine itf_set_index_info(itf_index_info,
     &     contr,contr_red,last_contr,
     &     isvtx1,isvtx2,ivtxres,nj_res)
*     condense the information about index labels (on contraction info contr and contr_red)
*     to an index array of the form
*     itf_index_info(<len1>,<len2>,<len3>,<idx11>,<idx12>,...<idx1(len1)>,<idx21>,...,<idx31>,...)
*     the indices contain orbital space (hpvx) and index (idx) as: 1000*hpvx+idx

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(out) ::
     &     itf_index_info(*)

      type(contraction), intent(in) ::
     &     contr, contr_red
      logical, intent(in) ::
     &     last_contr
      integer, intent(in) ::
     &     isvtx1, isvtx2, nj_res, ivtxres(nj_res)

      integer ::
     &     ivtxop(contr%nvtx)
      integer ::
     &     nvtxop, ivtx, ii, icnt, idx, idx1, idx2, idx3
      logical ::
     &     on_list

      nvtxop = 0
      do ivtx = 1, contr%nvtx
        if (contr%svertex(ivtx)==isvtx1) then
          nvtxop = nvtxop+1
          ivtxop(nvtxop) = ivtx
        end if
      end do

      ! keep first 3 entry for dimension info
      ii = 3
      icnt = 0
      do idx = 1, contr%nidx
        on_list = .false.
        do ivtx = 1, nvtxop
          on_list = on_list.or.contr%contr_string(idx)%vtx==ivtxop(ivtx)
        end do
        if (.not.on_list) cycle
        ii = ii+1
        icnt = icnt+1
        itf_index_info(ii)=contr%contr_string(idx)%hpvx*1000
     &       +contr%contr_string(idx)%idx
      end do
      itf_index_info(1) = icnt

      nvtxop = 0
      icnt = 0
      if (isvtx1.ne.isvtx2) then
        do ivtx = 1, contr%nvtx
          if (contr%svertex(ivtx)==isvtx2) then
            nvtxop = nvtxop+1
            ivtxop(nvtxop) = ivtx
          end if
        end do

        icnt = 0
        do idx = 1, contr%nidx
          on_list = .false.
          do ivtx = 1, nvtxop
            on_list = on_list.or.
     &           contr%contr_string(idx)%vtx==ivtxop(ivtx)
          end do
          if (.not.on_list) cycle
          ii = ii+1
          icnt = icnt+1
          itf_index_info(ii)=contr%contr_string(idx)%hpvx*1000
     &         +contr%contr_string(idx)%idx
        end do
      end if
      itf_index_info(2) = icnt

      icnt = 0
      if (.not.last_contr) then
        do idx = 1, contr_red%nidx
          on_list = .false.
          do ivtx = 1, nj_res
            on_list = on_list.or.
     &           contr_red%contr_string(idx)%vtx==ivtxres(ivtx)
          end do
          if (.not.on_list) cycle
          ii = ii+1
          icnt = icnt+1
          itf_index_info(ii)=contr_red%contr_string(idx)%hpvx*1000
     &         +contr_red%contr_string(idx)%idx
        end do
      else
        do idx = 1, contr_red%nxidx
          ii = ii+1
          icnt = icnt+1
          itf_index_info(ii)=contr_red%result_string(idx)%hpvx*1000
     &         +contr_red%result_string(idx)%idx
        end do
      end if
      itf_index_info(3) = icnt

      if (ntest.ge.100) then
        write(lulog,*) 'contents of itf_index_info:'
        write(lulog,'(1x,a,3i4)') 'info(1:3): ',itf_index_info(1:3)
        idx1 = itf_index_info(1)
        idx2 = itf_index_info(2)
        idx3 = itf_index_info(3)
        write(lulog,'(1x,a,10i6)') 'block1: ',
     &       itf_index_info(3+1:3+1+idx1-1)
        write(lulog,'(1x,a,10i6)') 'block2: ',
     &       itf_index_info(3+1+idx1:3+1+idx1+idx2-1)
        write(lulog,'(1x,a,10i6)') 'block3: ',
     &       itf_index_info(3+1+idx1+idx2:3+1+idx1+idx2+idx3-1)
      end if

      end subroutine
