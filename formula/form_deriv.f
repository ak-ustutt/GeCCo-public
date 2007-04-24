*----------------------------------------------------------------------*
      subroutine form_deriv(ffderiv,ffinput,name_deriv,
     &                      idxder,idxmlt,idxres,
     &                      ops,nops)
*----------------------------------------------------------------------*
*     get the derivatives of all contraction on ffinput
*
*     idxder is the index of the operator with respect to which the
*     derivative has to be taken
*     idxmlt is the index of the operator which the derivative is 
*     multiplied with (0 if only the derivative is taken)
*     idxres is the index of the resulting operator (0, if scalar)     
*
*     andreas, end of 2006, put to gecco april 2007
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_contraction.h'
      include 'def_contraction_list.h'

      character, intent(in) ::
     &     name_deriv*(*)
      type(filinf), intent(in) ::
     &     ffinput, ffderiv
      integer, intent(in) ::
     &     nops, idxder, idxmlt, idxres
      type(operator), intent(in) ::
     &     ops(nops)
      
      type(contraction) ::
     &     contr
      type(contraction_list), target ::
     &     conder
      type(contraction_list), pointer ::
     &     cur_conder

      integer ::
     &     luinput, luderiv,
     &     nterms, nder, idx, idum, idxinp, len

      logical, external ::
     &     rd_contr

      call file_open(ffinput)
      call file_open(ffderiv)
      luinput = ffinput%unit
      luderiv = ffderiv%unit
      rewind luinput
      rewind luderiv

      read(luinput)
      read(luinput) idum,idxinp

      len = len_trim(name_deriv)
      write(luderiv) len,name_deriv
      write(luderiv) idum,idxinp

      ! signal, that still nothing is allocated
      contr%mxvtx = 0
      contr%mxarc = 0
      contr%mxfac = 0

      nterms = 0
      do while(rd_contr(luinput,contr,idxinp))
        
        call contr_deriv(conder,nder,contr,ops,idxder,idxmlt,idxres)

        cur_conder => conder
        do idx = 1, nder
          nterms = nterms+1
          call wrt_contr(luderiv,cur_conder%contr)
          if (idx.lt.nder) cur_conder => cur_conder%next
        end do

      end do

      call file_close_keep(ffderiv)
      call file_close_keep(ffinput)

      call dealloc_contr(contr)
      call dealloc_contr_list(conder)
      
      return
      end
