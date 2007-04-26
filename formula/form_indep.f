*----------------------------------------------------------------------*
      subroutine form_indep(ffoutput,ffinput,name_out,
     &                      ncmpnd,idxop,
     &                      ops,nops)
*----------------------------------------------------------------------*
*     collect those contractions on ffoutput that do not depend on
*     any operator on list idxop(1:ncmpnd)
*
*     andreas, april 2007
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_contraction.h'
      include 'def_contraction_list.h'

      integer, parameter ::
     &     ntest = 00

      character, intent(in) ::
     &     name_out*(*)
      type(filinf), intent(in) ::
     &     ffinput, ffoutput
      integer, intent(in) ::
     &     ncmpnd, nops, idxop(ncmpnd)
      type(operator), intent(in) ::
     &     ops(nops)
      
      type(contraction) ::
     &     contr

      logical ::
     &     ok
      integer ::
     &     luinput, luoutput,
     &     nterms, idum, idxinp, idx, icmpnd

      logical, external ::
     &     rd_contr

      if (ntest.ge.100) then
        write(luout,*) '========================'
        write(luout,*) ' here speaks form_indep'
        write(luout,*) '========================'
      end if

      call file_open(ffinput)
      call file_open(ffoutput)
      luinput = ffinput%unit
      luoutput = ffoutput%unit
      rewind luinput
      rewind luoutput

      read(luinput)
      read(luinput) idum,idxinp

      write(luoutput) name_out
      write(luoutput) idum,idxinp

      ! signal, that still nothing is allocated
      contr%mxvtx = 0
      contr%mxarc = 0
      contr%mxfac = 0

      nterms = 0
      do while(rd_contr(luinput,contr,idxinp))
        
        ok = .true.
        cmp_loop: do idx = 1, contr%nvtx
          do icmpnd = 1, ncmpnd
            if (contr%vertex(idx)%idx_op.eq.idxop(icmpnd)) then
              ok = .false.
              exit cmp_loop
            end if
          end do
        end do cmp_loop

        if (ok) then
          nterms = nterms+1
          call wrt_contr(luoutput,contr)
          if (ntest.ge.100) then
            call prt_contr(luout,contr,ops)
          end if
        end if

      end do

      call file_close_keep(ffoutput)
      call file_close_keep(ffinput)

      call dealloc_contr(contr)
      
      return
      end
