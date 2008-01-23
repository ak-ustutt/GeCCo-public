      integer function max_rank_op(mode_str,op)

      implicit none

      include 'opdim.h'
      include 'def_operator.h'

      type(operator), intent(in) ::
     &     op
      character(*), intent(in) ::
     &     mode_str

      integer ::
     &     nblk, iblk
      integer, pointer ::
     &     ca_occ(:,:)

      nblk = op%n_occ_cls
      ca_occ => op%ica_occ
      max_rank_op = -huge(max_rank_op)
      select case(trim(mode_str))
      case ('C-A')
        do iblk = 1, nblk
          max_rank_op = max(max_rank_op,ca_occ(1,iblk)-ca_occ(2,iblk))
        end do
      case ('A-C')
        do iblk = 1, nblk
          max_rank_op = max(max_rank_op,ca_occ(2,iblk)-ca_occ(1,iblk))
        end do
      case ('C')
        do iblk = 1, nblk
          max_rank_op = max(max_rank_op,ca_occ(1,iblk))
        end do
      case ('A')
        do iblk = 1, nblk
          max_rank_op = max(max_rank_op,ca_occ(2,iblk))
        end do
      case default
        call quit(1,'rank_occ','unknown mode_str: ',trim(mode_str))
      end select

      return
      end
