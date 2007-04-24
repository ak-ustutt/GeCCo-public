*----------------------------------------------------------------------*
      subroutine form_indep(ffoutput,ffinput,name_out,
     &                      idxop,
     &                      ops,nops)
*----------------------------------------------------------------------*
*     collect those contractions on ffoutput that do not depend on
*     operator idxop
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
     &     ntest = 100

      character, intent(in) ::
     &     name_out*(*)
      type(filinf), intent(in) ::
     &     ffinput, ffoutput
      integer, intent(in) ::
     &     nops, idxop
      type(operator), intent(in) ::
     &     ops(nops)
      
      type(contraction) ::
     &     contr

      logical ::
     &     ok
      integer ::
     &     luinput, luoutput,
     &     nterms, idum, idxinp, idx

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
        do idx = 1, contr%nvtx
          if (contr%vertex(idx)%idx_op.eq.idxop) then
            ok = .false.
            exit
          end if
        end do

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
