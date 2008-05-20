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
     &     nblk, iblk, ica, icast, icand, njoined
      logical, pointer ::
     &     formal(:)
      integer, pointer ::
     &     ca_occ(:,:), hpvx_occ(:,:,:)

      nblk = op%n_occ_cls
      njoined = op%njoined
      ca_occ => op%ica_occ
      hpvx_occ => op%ihpvca_occ
      max_rank_op = -huge(max_rank_op)
      formal => op%formal_blk
      select case(trim(mode_str))
      case ('C-A')
        do iblk = 1, nblk
          if (skip_formal.and.formal(iblk)) cycle
          max_rank_op = max(max_rank_op,ca_occ(1,iblk)-ca_occ(2,iblk))
        end do
      case ('A-C')
        do iblk = 1, nblk
          if (skip_formal.and.formal(iblk)) cycle
          max_rank_op = max(max_rank_op,ca_occ(2,iblk)-ca_occ(1,iblk))
        end do
      case ('C')
        do iblk = 1, nblk
          if (skip_formal.and.formal(iblk)) cycle
          max_rank_op = max(max_rank_op,ca_occ(1,iblk))
        end do
      case ('A')
        do iblk = 1, nblk
          if (skip_formal.and.formal(iblk)) cycle
          max_rank_op = max(max_rank_op,ca_occ(2,iblk))
        end do
      case ('H','CH','AH')
        icast = 1
        icand = 2
        if (trim(mode_str).eq.'CH') icand = 1
        if (trim(mode_str).eq.'AH') icast = 2
        do iblk = 1, nblk*njoined
          if (skip_formal.and.formal(iblk)) cycle
          do ica = icast, icand
            max_rank_op = max(max_rank_op,hpvx_occ(IHOLE,ica,iblk))
          end do
        end do
      case ('P','CP','AP')
        icast = 1
        icand = 2
        if (trim(mode_str).eq.'CP') icand = 1
        if (trim(mode_str).eq.'AP') icast = 2
        do iblk = 1, nblk*njoined
          if (skip_formal.and.formal(iblk)) cycle
          do ica = icast, icand
            max_rank_op = max(max_rank_op,hpvx_occ(IPART,ica,iblk))
          end do
        end do
      case ('X','CX','AX')
        icast = 1
        icand = 2
        if (trim(mode_str).eq.'CX') icand = 1
        if (trim(mode_str).eq.'AX') icast = 2
        do iblk = 1, nblk*njoined
          if (skip_formal.and.formal(iblk)) cycle
          do ica = icast, icand
            max_rank_op = max(max_rank_op,hpvx_occ(IEXTR,ica,iblk))
          end do
        end do
      case default
        call quit(1,'rank_occ','unknown mode_str: ',trim(mode_str))
      end select

      return
      end
