*----------------------------------------------------------------------*
      subroutine dummy_restr(irestr_out,
     &     iocc_fit,njoined_fit,orb_info)
*----------------------------------------------------------------------*
*     work-around: provide dummy restrictions
*     should work for both reversed and normal hole-ordering
*----------------------------------------------------------------------*
      implicit none
      
      include 'opdim.h'
      include 'stdunit.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest = 00

      type(orbinf), intent(in), target ::
     &     orb_info
      integer, intent(in) ::
     &     njoined_fit,
     &     iocc_fit(ngastp,2,njoined_fit)
      integer, intent(out) ::
     &     irestr_out(2,orb_info%ngas,2,2,njoined_fit)

      integer ::
     &     imask, ica, igas, ityp, irstmax, idiff, ijoin, ngas
      integer ::
     &     iocc_sum(ngastp)
      integer, pointer ::
     &     ihpvgas(:,:), iad_gas(:)

c      call warn('dummy_restr','I should be obsolete soon ...!')
      
      ngas = orb_info%ngas
      ihpvgas => orb_info%ihpvgas
      iad_gas => orb_info%iad_gas

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,' here speaks dummy_restr')
        write(lulog,*) 'input occupation:'
        call wrt_occ_n(lulog,iocc_fit,njoined_fit)
      end if

      do ijoin = 1, njoined_fit
        do imask = 1, 2
          if (imask.eq.2) then
            irestr_out(1:2,1:ngas,1:2,imask,1:njoined_fit) = 0
            exit
          end if
          do ica = 1, 2
            do igas = 1, ngas
              ityp = ihpvgas(igas,1) ! OPEN SHELL: adapt
              irstmax = iocc_fit(ityp,ica,ijoin)
              if (iad_gas(igas).ne.2.and.igas.eq.1) irstmax = 0
              irestr_out(1,igas,ica,imask,ijoin) = irstmax
              irestr_out(2,igas,ica,imask,ijoin) = irstmax
            end do
          end do
        end do
      end do

      if (ntest.ge.100) then
        write(lulog,*) 'output restriction: '
        do ijoin = 1, njoined_fit
          call wrt_rstr(lulog,irestr_out(1,1,1,1,ijoin),ngas)
        end do
      end if

      return
      end
