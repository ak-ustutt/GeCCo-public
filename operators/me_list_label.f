*----------------------------------------------------------------------*
      subroutine me_list_label(label,root,sym,spin,ms,msc,spin_adapt)
*----------------------------------------------------------------------*
*     generate a systematic spin and symmetry including label
*----------------------------------------------------------------------*

      implicit none

      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'

      character(*), intent(out) ::
     &     label
      logical, intent(in) ::
     &     spin_adapt
      integer, intent(in) ::
     &     sym, spin, ms, msc
      character(*), intent(in) ::
     &     root

      if (len_trim(root)+8.gt.mxlen_melabel)
     &     call quit(1,'me_list_label',
     &     'root-name too long: '//trim(root))

      if (.not.spin_adapt.and.msc.eq.0) then
        write(label,'(a,"G",i1,"SxxM",i2.2)') trim(root),sym,ms
      else if (.not.spin_adapt.and.msc.eq.+1) then
        write(label,'(a,"G",i1,"C+1M",i2.2)') trim(root),sym,ms
      else if (.not.spin_adapt.and.msc.eq.-1) then
        write(label,'(a,"G",i1,"C-1M",i2.2)') trim(root),sym,ms
      else
        write(label,'(a,"G",i1,"S",i2,"M",i2.2)') trim(root),sym,spin,ms
      end if

      return
      end
