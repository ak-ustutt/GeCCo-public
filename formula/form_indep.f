*----------------------------------------------------------------------*
      subroutine form_indep(form_output,form_input,
     &                      label_out,title_out,
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
      include 'def_formula.h'

      integer, parameter ::
     &     ntest = 00

      character, intent(in) ::
     &     label_out*(*), title_out*(*)
      type(formula), intent(inout) ::
     &     form_input, form_output
      integer, intent(in) ::
     &     ncmpnd, nops, idxop(ncmpnd)
      type(operator), intent(in) ::
     &     ops(nops)
      
      type(contraction) ::
     &     contr

      logical ::
     &     ok
      integer ::
     &     luinput, lulogput, len,
     &     nterms, idum, idxinp, idx, icmpnd
      character ::
     &     name*(form_maxlen_label*2)

      logical, external ::
     &     rd_contr

      if (ntest.ge.100) then
        write(lulog,*) '========================'
        write(lulog,*) ' here speaks form_indep'
        write(lulog,*) '========================'
        write(lulog,*) ' new label: "',trim(label_out),'"'
        write(lulog,*) ' new title: "',trim(title_out),'"'
      end if

      write(name,'(a,".fml")') label_out
      call file_init(form_output%fhand,name,ftyp_sq_unf,0)      
      form_output%label = label_out
      form_output%comment = title_out

      call file_open(form_input%fhand)
      call file_open(form_output%fhand)
      luinput = form_input%fhand%unit
      lulogput = form_output%fhand%unit
      rewind luinput
      rewind lulogput

      read(luinput)
      read(luinput) idum,idxinp

      len = len_trim(title_out)
      write(lulogput) len,title_out
      write(lulogput) idum,idxinp

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
          call wrt_contr(lulogput,contr)
          if (ntest.ge.100) then
            call prt_contr(lulog,contr,ops)
          end if
        end if

      end do

      call file_close_keep(form_output%fhand)
      call file_close_keep(form_input%fhand)

      call dealloc_contr(contr)

      if (ntest.ge.10) then
        write(lulog,*) 'generated terms: ',nterms
      end if
      
      return
      end
