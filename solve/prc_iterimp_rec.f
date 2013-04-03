*----------------------------------------------------------------------*
      recursive subroutine prc_iterimp_rec(ilevel,nlevel,fac,
     &     xbuf1,xbuf2,nbuf,iopt,ffdia,me_v,me_aoff,me_w,fmvp,
     &     opti_info,orb_info,op_info,str_info,strmap_info)
*----------------------------------------------------------------------*
*     iterative improvement of the preconditioning step according to
*     (A^-1)^(n) = (A^-1)^(n-1)*[(1+fac)-fac*(A^(0)-Aoff)*(A^-1)^(n-1)]
*     ilevel : current level
*     nlevel : total number of iterations;
*     xbuf1  : initial vector on input
*     xbuf2, ffdia  : contain preconditioner => (A^-1)^(0)
*     me_v   : the current record contains the
*              initial vector on input (except for ilevel=0,nlevel>0),
*              and the preconditioned vector on output
*     me_aoff: minus fac times the off-diagonal part of the matrix A
*     me_w   : file for the input vector of the matrix vector product
*     fmvp   : formula for the matrix vector product -Aoff*v
*
*     matthias, mar 2013
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_optimize_info.h'
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     ilevel, nlevel, nbuf, iopt
      real(8), intent(in) ::
     &     fac
      type(filinf), intent(in) ::
     &     ffdia
      type(me_list), intent(inout) ::
     &     me_v, me_w
      type(me_list), intent(in) ::
     &     me_aoff
      real(8), intent(inout) ::
     &     xbuf1(*), xbuf2(*)
      type(formula_item), intent(in) ::
     &     fmvp

      type(optimize_info), intent(in) ::
     &     opti_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf),intent(in) ::
     &     str_info
      type(strmapinf),intent(in) ::
     &     strmap_info
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     nopt, nsec, isec, stsec, ndsec, irec_cur, ii
      real(8) ::
     &     xdum, prc_impfac

      integer, pointer ::
     &     nwfpsec(:), idstsec(:), nsec_arr(:)

      real(8), pointer ::
     &     signsec(:)

      if (ntest.ge.100)
     &   write(luout,'(a,i4)') 'prc_iterimp entered on level ',ilevel

      nopt = opti_info%nopt
      irec_cur = me_v%fhand%current_record
      nsec_arr => opti_info%nsec(1:nopt)
      nsec = sum(nsec_arr)
      nwfpsec => opti_info%nwfpsec(1:nsec)
      idstsec => opti_info%idstsec(1:nsec)
      signsec => opti_info%signsec(1:nsec)!2(1:nsec)
      stsec = 1
      ndsec = 0
      if (iopt.gt.1) stsec = stsec + nsec_arr(iopt-1)
      ndsec = ndsec + nsec_arr(iopt)

      if (ilevel.le.0) then

        if (ntest.ge.100)
     &     write(luout,*) 'Application of (A^-1)^(0)'

        ! Just perform ordinary preconditioning step
        ! Divide by precond., account for sign changes if necessary
        do isec = stsec, ndsec
          call diavc(xbuf1(idstsec(isec)),xbuf1(idstsec(isec)),
     &               signsec(isec),xbuf2(idstsec(isec)),
     &               0d0,nwfpsec(isec))
        end do

        ! Store vector and return
        call vec_to_da(me_v%fhand,me_v%fhand%current_record,xbuf1,nbuf)
        return

      else if (ilevel.eq.1) then

        ! do the initial preconditioning directly to me_w:
        call prc_iterimp_rec(0,nlevel,fac,
     &                       xbuf1,xbuf2,nbuf,iopt,ffdia,
     &                       me_w,me_aoff,
     &                       me_w,fmvp,opti_info,
     &                       orb_info,op_info,str_info,strmap_info)

      else

        ! copy initial vector to next record and call kernel
        if (nlevel.eq.ilevel.and.irec_cur.ne.1)
     &     call quit(1,'prc_iterimp_rec','not the right record?')
        if (me_v%fhand%active_records(2).lt.irec_cur+1)
     &     call quit(1,'prc_iterimp_rec','not enough records')
        call switch_mel_record(me_v,irec_cur+1)
        call vec_to_da(me_v%fhand,irec_cur+1,xbuf1,nbuf)
        call prc_iterimp_rec(ilevel-1,nlevel,fac,
     &                       xbuf1,xbuf2,nbuf,iopt,ffdia,
     &                       me_v,me_aoff,
     &                       me_w,fmvp,opti_info,
     &                       orb_info,op_info,str_info,strmap_info)
        ! copy result to me_w:
        call vec_to_da(me_w%fhand,1,xbuf1,nbuf)
        ! switch back to current record
        call switch_mel_record(me_v,irec_cur)

      end if

      if (ntest.ge.100)
     &   write(luout,'(a,i2.2,a,i2.2,a)')
     &   'Computing [1-A^(0)*(A^-1)^(',ilevel,
     &   ')-Aoff*(A^-1)^(',ilevel-1,')]*|v>'
 
      ! for ilevel=1 we use that 2 - A^(0)*(A^-1)^(0) = 1
      ! for higher levels we manually compute 2 - A^(0)*(A^-1)^(n-1)
      if (ilevel.gt.1) then
        ! multiply with A^(0) (=undo last preconditioning step)
        do isec = stsec, ndsec
          do ii = idstsec(isec), idstsec(isec)-1+nwfpsec(isec)
            xbuf1(ii) = signsec(isec)*xbuf1(ii)*xbuf2(ii)
          end do
        end do
        ! get original vector times two and substract the above
        call vec_from_da(me_v%fhand,irec_cur,xbuf2,nbuf)
        xbuf2(1:nbuf) = (1d0+fac)*xbuf2(1:nbuf) - fac*xbuf1(1:nbuf)
        call vec_to_da(me_v%fhand,irec_cur,xbuf2,nbuf)
        ! restore preconditioner
        call vec_from_da(ffdia,1,xbuf2,nbuf)
      end if

      ! compute v = v - Aoff*v
      call evaluate2(fmvp,.false.,.true.,
     &               op_info,str_info,strmap_info,orb_info,
     &               xdum,.false.)

      if (ntest.ge.100)
     &   write(luout,'(a,i2.2,a)')
     &   'Now applying (A^-1)^(',ilevel-1,') again from the left'

      ! call kernel a second time
      call vec_from_da(me_v%fhand,irec_cur,xbuf1,nbuf)
      call prc_iterimp_rec(ilevel-1,nlevel,fac,
     &                     xbuf1,xbuf2,nbuf,iopt,ffdia,
     &                     me_v,me_aoff,
     &                     me_w,fmvp,opti_info,
     &                     orb_info,op_info,str_info,strmap_info)

      if (ntest.ge.100)
     &   write(luout,'(a,i4)') 'prc_iterimp exits level      ',ilevel

      return
      end
