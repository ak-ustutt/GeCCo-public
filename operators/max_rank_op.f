      integer function max_rank_op(mode_str,op,skip_formal)

      implicit none

      include 'opdim.h'
      include 'def_operator.h'

      logical, intent(in) ::
     &     skip_formal
      type(operator), intent(in) ::
     &     op
      character(*), intent(in) ::
     &     mode_str

      integer ::
     &     nblk, njoined, iblk, ica, icast, icand, idx, occsum
      logical, pointer ::
     &     formal(:)
      integer, pointer ::
     &     ca_occ(:,:), hpvx_occ(:,:,:)

      njoined = op%njoined
      nblk = op%n_occ_cls
      ca_occ => op%ica_occ
      hpvx_occ => op%ihpvca_occ
      max_rank_op = -huge(max_rank_op)
      formal => op%formal_blk
      select case(trim(mode_str))
      case ('C-A')
        if (njoined.ne.1)
     &       call quit(1,'max_rank_op',
     &       'option '//trim(mode_str)//' is not adapted for njoined>1')
        do iblk = 1, nblk
          if (skip_formal.and.formal(iblk)) cycle
          max_rank_op = max(max_rank_op,ca_occ(1,iblk)-ca_occ(2,iblk))
        end do
      case ('A-C')
        if (njoined.ne.1)
     &       call quit(1,'max_rank_op',
     &       'option '//trim(mode_str)//' is not adapted for njoined>1')
        do iblk = 1, nblk
          if (skip_formal.and.formal(iblk)) cycle
          max_rank_op = max(max_rank_op,ca_occ(2,iblk)-ca_occ(1,iblk))
        end do
      case ('C')
        if (njoined.ne.1)
     &       call quit(1,'max_rank_op',
     &       'option '//trim(mode_str)//' is not adapted for njoined>1')
        do iblk = 1, nblk
          if (skip_formal.and.formal(iblk)) cycle
          max_rank_op = max(max_rank_op,ca_occ(1,iblk))
        end do
      case ('A')
        if (njoined.ne.1)
     &       call quit(1,'max_rank_op',
     &       'option '//trim(mode_str)//' is not adapted for njoined>1')
        do iblk = 1, nblk
          if (skip_formal.and.formal(iblk)) cycle
          max_rank_op = max(max_rank_op,ca_occ(2,iblk))
        end do
      case ('H','CH','AH')
        icast = 1
        icand = 2
        if (trim(mode_str).eq.'CH') icand = 1
        if (trim(mode_str).eq.'AH') icast = 2
        do iblk = 1, nblk
          if (skip_formal.and.formal(iblk)) cycle
          do ica = icast, icand
            occsum = 0
            do idx = (iblk-1)*njoined+1, iblk*njoined
              occsum = occsum + hpvx_occ(IHOLE,ica,idx)
            end do
            max_rank_op = max(max_rank_op,occsum)
          end do
        end do
      case ('P','CP','AP')
        icast = 1
        icand = 2
        if (trim(mode_str).eq.'CP') icand = 1
        if (trim(mode_str).eq.'AP') icast = 2
        do iblk = 1, nblk
          if (skip_formal.and.formal(iblk)) cycle
          do ica = icast, icand
            occsum = 0
            do idx = (iblk-1)*njoined+1, iblk*njoined
              occsum = occsum + hpvx_occ(IPART,ica,idx)
            end do
            max_rank_op = max(max_rank_op,occsum)
          end do
        end do
      case ('X','CX','AX')
        icast = 1
        icand = 2
        if (trim(mode_str).eq.'CX') icand = 1
        if (trim(mode_str).eq.'AX') icast = 2
        do iblk = 1, nblk
          if (skip_formal.and.formal(iblk)) cycle
          do ica = icast, icand
            occsum = 0
            do idx = (iblk-1)*njoined+1, iblk*njoined
              occsum = occsum + hpvx_occ(IEXTR,ica,idx)
            end do
            max_rank_op = max(max_rank_op,occsum)
          end do
        end do
      case default
        call quit(1,'rank_occ','unknown mode_str: ',trim(mode_str))
      end select

      return
      end
