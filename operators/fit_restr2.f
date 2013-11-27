*----------------------------------------------------------------------*
      subroutine fit_restr2(irestr_out,
     &     iocc_fit,njoined_fit,irestr_raw,njoined_raw,ihpvgas,ngas)
*----------------------------------------------------------------------*
*     given is a restriction (irestr_raw), e.g. on some final result of 
*     a multiple contraction. iocc_fit is the occupation of some inter-
*     mediate. the restriction on this intermediate is given by the
*     raw restriction. for easier comparison in later parts of the 
*     program, it is advantageous to "fit" the restriction to the
*     actual maximum occupations (on iocc_fit)
*     version for joined vertices; njoined_raw is either 1 
*     (applied to each component) or must be equal njoined_fit
*----------------------------------------------------------------------*
      implicit none
      
      include 'opdim.h'
      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     njoined_raw, njoined_fit,
     &     ngas, ihpvgas(ngas),
     &     iocc_fit(ngastp,2,njoined_fit),
     &     irestr_raw(2,ngas,2,2,njoined_raw)
      integer, intent(out) ::
     &     irestr_out(2,ngas,2,2,njoined_fit)

      integer ::
     &     imask, ica, igas, ityp, irstmax, idiff, ijoin, ijoin_raw
      integer ::
     &     iocc_sum(ngastp)
      
      if (ntest.ge.100) then
        write(lulog,*) '======================='
        write(lulog,*) ' here speaks fit_restr'
        write(lulog,*) '======================='
        write(lulog,*) 'input restriction: '
        do ijoin = 1, njoined_raw
          call wrt_rstr(lulog,irestr_raw(1,1,1,1,ijoin),ngas)
        end do
        write(lulog,*) 'input occupation:'
        call wrt_occ_n(lulog,iocc_fit,njoined_fit)
      end if

      if (njoined_raw.ne.1.and.njoined_raw.ne.njoined_fit)
     &     call quit(1,'fit_restr2','dimension mismatch!')

      do ijoin = 1, njoined_fit
        ijoin_raw = min(ijoin,njoined_raw)
        do imask = 1, 2
          do ica = 1, 2
c          ! sum up maximum occupations after space of given type
c          iocc_sum(1:ngastp) = 0
            do igas = 1, ngas
              ityp = ihpvgas(igas)
c            iocc_sum(ityp) = iocc_sum(ityp)+iocc_fit(ityp,ica)
              irstmax = min(iocc_fit(ityp,ica,ijoin),
     &                    irestr_raw(2,igas,ica,imask,ijoin_raw))
              idiff = irestr_raw(2,igas,ica,imask,ijoin_raw)
     &               -irestr_raw(1,igas,ica,imask,ijoin_raw)
              irestr_out(1,igas,ica,imask,ijoin) = max(0,irstmax-idiff)
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
