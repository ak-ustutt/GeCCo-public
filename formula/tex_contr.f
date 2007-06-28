*----------------------------------------------------------------------*
      subroutine tex_contr(luout,first,contr,op_info)
*----------------------------------------------------------------------*
*     write info on contraction in TeX style onto unit luout
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'

      integer, intent(in) ::
     &     luout
      logical, intent(in) ::
     &     first
      type(contraction), intent(in) ::
     &     contr
      type(operator_info), intent(in) ::
     &     op_info

      type(operator), pointer ::
     &     op
      
      integer ::
     &     ipos, nocc, ioff, ivtx, iarc, ifac, p, q
      integer, pointer ::
     &     occ(:,:), occ_dst(:,:,:), idxset(:)
      character ::
     &     str*512

      op => op_info%op_arr(contr%idx_res)%op
      occ => op%ihpvca_occ(1:ngastp,1:2,contr%iblk_res)

      if (first) then
        call tex_op(str,op%name,1,occ,1)
        ipos = len_trim(str)+1
        write(str(ipos:),'("+=")')
      else
        write(str,'("+")')
      end if

      if (abs(contr%fac-1d0).gt.1d-12) then
        call real2rat(p,q,contr%fac)
        ipos = len_trim(str)+1
        write(str(ipos:),'("\frac{",i6,"}{",i6,"}")') p,q
      end if

      allocate(occ_dst(ngastp,2,contr%narc+1),idxset(contr%narc+1))
      do ivtx = 1, contr%nvtx
        op => op_info%op_arr(contr%vertex(ivtx)%idx_op)%op
        if (op%dagger) then
          occ_dst(1:ngastp,1:2,1)
     &       = iocc_dagger(op%ihpvca_occ(1:ngastp,1:2,
     &         contr%vertex(ivtx)%iblk_op))
        else
          occ_dst(1:ngastp,1:2,1)
     &       = op%ihpvca_occ(1:ngastp,1:2,
     &         contr%vertex(ivtx)%iblk_op)
        end if
        nocc = 1
        ioff = 0
        idxset(1) = 1
        do iarc = 1, contr%narc
          if (contr%arc(iarc)%link(1).eq.ivtx.or.
     &        contr%arc(iarc)%link(2).eq.ivtx) then
            nocc = nocc+1
            if (contr%arc(iarc)%link(1).eq.ivtx) then
              occ_dst(1:ngastp,1:2,nocc) =
     &             contr%arc(iarc)%occ_cnt
              ifac = 1
            else
              occ_dst(1:ngastp,1:2,nocc) =
     &             iocc_dagger(contr%arc(iarc)%occ_cnt)
              ifac = -1
            end if
            occ_dst(1:ngastp,1:2,1) =
     &           occ_dst(1:ngastp,1:2,1) - occ_dst(1:ngastp,1:2,nocc)
            ! generate 2 3 5 6 8 9 ... etc.
            idxset(nocc) = (iarc + (iarc+1)/2)*ifac
          end if
        end do
        if (.not.iocc_nonzero(occ_dst(1:ngastp,1:2,1))) then
          nocc = nocc-1
          ioff = 1
        end if
        ipos = len_trim(str)+1
        call tex_op(str(ipos:),op%name,nocc,occ_dst(1,1,ioff+1),
     &       idxset(ioff+1))

      end do

      write(luout,'(a)') trim(str)

      deallocate(occ_dst,idxset)

      return
      end
