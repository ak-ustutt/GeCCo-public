*----------------------------------------------------------------------*
      subroutine dummy_restr(irestr_out,
     &     iocc_fit,njoined_fit,ihpvgas,ngas)
*----------------------------------------------------------------------*
*     work-around: provide dummy restrictions
*----------------------------------------------------------------------*
      implicit none
      
      include 'opdim.h'
      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     njoined_fit,
     &     ngas, ihpvgas(ngas),
     &     iocc_fit(ngastp,2,njoined_fit)
      integer, intent(out) ::
     &     irestr_out(2,ngas,2,2,njoined_fit)

      integer ::
     &     imask, ica, igas, ityp, irstmax, idiff, ijoin
      integer ::
     &     iocc_sum(ngastp)
      
      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,' here speaks dummy_restr')
        write(luout,*) 'input occupation:'
        call wrt_occ_n(luout,iocc_fit,njoined_fit)
      end if

      do ijoin = 1, njoined_fit
        do imask = 1, 2
          if (imask.eq.2) then
            irestr_out(1:2,1:ngas,1:2,imask,1:njoined_fit) = 0
            exit
          end if
          do ica = 1, 2
            do igas = 1, ngas
              ityp = ihpvgas(igas)
              irstmax = iocc_fit(ityp,ica,ijoin)
              irestr_out(1,igas,ica,imask,ijoin) = irstmax
              irestr_out(2,igas,ica,imask,ijoin) = irstmax
            end do
          end do
        end do
      end do

      if (ntest.ge.100) then
        write(luout,*) 'output restriction: '
        do ijoin = 1, njoined_fit
          call wrt_rstr(luout,irestr_out(1,1,1,1,ijoin),ngas)
        end do
      end if

      return
      end
