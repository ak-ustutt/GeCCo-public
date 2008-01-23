*----------------------------------------------------------------------*
      subroutine print_target_list(luout,tgt_info)
*----------------------------------------------------------------------*
*     print target list on standard output
*----------------------------------------------------------------------*
      implicit none

      include 'mdef_target_info.h'

      integer, intent(in) ::
     &     luout
      type(target_info), intent(in) ::
     &     tgt_info

      character*(15), parameter ::
     &     name_ttype(0:4) = (/
     &                       'phony          ',
     &                       'operator       ',
     &                       'formula        ',
     &                       'opt. formula   ',
     &                       'operator list  '/)

      integer ::
     &     itgt, idx, jdx
      type(target), pointer ::
     &     tgt

c      if (long) then
        do itgt = 1, tgt_info%ntargets
          tgt => tgt_info%array(itgt)%tgt
          if (tgt%required) then
            write(luout,'(x,a," (",a,") primary")')
     &           trim(tgt%name),trim(name_ttype(tgt%type))
          else
            write(luout,'(x,a," (",a,") ")')
     &           trim(tgt%name),trim(name_ttype(tgt%type))
          end if
          if (tgt%n_joined_with.gt.0)
     &         write(luout,'(4x,a)') 'joined targets:'
          do idx = 1, tgt%n_joined_with
            write(luout,'(6x,a)') trim(tgt%joined_with(idx))
          end do
          if (associated(tgt%idx_joined_with))
     &         write(luout,'(6x,8i4)')
     &         tgt%idx_joined_with(1:tgt%n_joined_with)
          if (tgt%n_depends_on.gt.0)
     &         write(luout,'(4x,a)') 'depends on:'
          do idx = 1, tgt%n_depends_on
            write(luout,'(6x,a)') trim(tgt%depends_on(idx))
          end do
          if (associated(tgt%idx_depends_on))
     &         write(luout,'(6x,8i4)')
     &         tgt%idx_depends_on(1:tgt%n_depends_on)
          write(luout,'(4x,a)') 'rules:'
          do idx = 1, tgt%n_rules
            write(luout,'(6x,a)') trim(tgt%rules(idx)%command)
            do jdx = 1, tgt%rules(idx)%n_labels
              write(luout,'(10x,a)')trim(tgt%rules(idx)%labels(jdx))
            end do
            do jdx = 1, tgt%rules(idx)%n_parameter_strings
              write(luout,'(10x,a)')trim(tgt%rules(idx)%parameters(jdx))
            end do
          end do
          
        end do
c      else
          
c      end if
        
      return
      end
