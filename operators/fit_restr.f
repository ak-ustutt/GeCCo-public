*----------------------------------------------------------------------*
      subroutine fit_restr(irestr_out,iocc_fit,irestr_raw,ihpvgas,ngas)
*----------------------------------------------------------------------*
*     given is a restriction (irestr_raw), e.g. on some final result of 
*     a multiple contraction. iocc_fit is the occupation of some inter-
*     mediate. the restriction on this intermediate is given by the
*     raw restriction. for easier comparison in later parts of the 
*     program, it is advantageous to "fit" the restriction to the
*     actual maximum occupations (on iocc_fit)
*----------------------------------------------------------------------*
      implicit none
      
      include 'opdim.h'
      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     ngas, ihpvgas(ngas),
     &     iocc_fit(ngastp,2), irestr_raw(2,ngas,2,2)
      integer, intent(out) ::
     &     irestr_out(2,ngas,2,2)

      integer ::
     &     imask, ica, igas, ityp, irstmax, idiff
      integer ::
     &     iocc_sum(ngastp)
      
      if (ntest.ge.100) then
        write(luout,*) '======================='
        write(luout,*) ' here speaks fit_restr'
        write(luout,*) '======================='
        write(luout,*) 'input restriction: '
        call wrt_rstr(luout,irestr_raw,ngas)
        write(luout,*) 'input occupation:'
        call wrt_occ(luout,iocc_fit)
      end if

      do imask = 1, 2
        do ica = 1, 2
c          ! sum up maximum occupations after space of given type
c          iocc_sum(1:ngastp) = 0
          do igas = 1, ngas
            ityp = ihpvgas(igas)
c            iocc_sum(ityp) = iocc_sum(ityp)+iocc_fit(ityp,ica)
            irstmax = min(iocc_fit(ityp,ica),
     &                    irestr_raw(2,igas,ica,imask))
            idiff = irestr_raw(2,igas,ica,imask)
     &             -irestr_raw(1,igas,ica,imask)
            irestr_out(1,igas,ica,imask) = max(0,irstmax-idiff)
            irestr_out(2,igas,ica,imask) = irstmax
          end do
        end do
      end do

      if (ntest.ge.100) then
        write(luout,*) 'output restriction: '
        call wrt_rstr(luout,irestr_out,ngas)
      end if

      return
      end
