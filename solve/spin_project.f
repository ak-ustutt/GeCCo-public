*----------------------------------------------------------------------*
      subroutine spin_project(me_amp,me_special,fspc,
     &     nwfpar,xbuf1,xbuf2,normalize,xret,opti_info,
     &     orb_info,op_info,str_info,strmap_info)
*----------------------------------------------------------------------*
*
*     projects out unwanted spin components
*     fspc         : formula for action of S^2 operator
*     me_special   : result vector
*
*     matthias, Jun 2012
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
     &     ntest = 100

      integer, intent(in) ::
     &     nwfpar

      type(me_list), intent(inout) ::
     &     me_amp, me_special

      real(8), intent(inout) ::
     &     xbuf1(*), xbuf2(*)
      real(8), intent(out) ::
     &     xret

      logical, intent(in) ::
     &     normalize

      type(formula_item), intent(in) ::
     &     fspc

      type(optimize_info), intent(in) ::
     &     opti_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf),intent(in) ::
     &     str_info
      type(strmapinf),intent(in) ::
     &     strmap_info
      type(operator_info), intent(inout) ::
     &     op_info


      integer ::
     &     mult_min, mult, mult_max, ab_sym, imult, iblk

      real(8) ::
     &     ssp1, kkp1, fac1, fac2

      type(filinf), pointer ::
     &     ffamp, ffspc

      real(8), external ::
     &     da_ddot

      ! pointers to file handle
      ffamp => me_amp%fhand
      ffspc => me_special%fhand

      ab_sym = me_amp%absym
      if (ab_sym.ne.me_special%absym)
     &   call quit(1,'spin_project','incompatible absym')
      if (me_amp%op%njoined.ne.1)
     &   call quit(1,'spin_project','not yet adapted for nj>1')
      mult = orb_info%imult
      ssp1 = dble(mult**2-1)/4d0
      mult_min = 2 - mod(mult,2) !1 for odd and 2 for even multipl.
      mult_max = 1
      do iblk = 1, me_amp%op%n_occ_cls
        mult_max = max(mult_max,me_amp%op%ica_occ(1,iblk)+
     &                          me_amp%op%ica_occ(2,iblk)+1)
      end do

      ! loop over spin multiplicities 2K+1 to be projected out
      do imult = mult_min, mult_max, 2 ! only every second can contrib.
        if (imult.eq.mult) cycle ! don't project out the one we want
        if (ab_sym.ne.0.and.mod(imult-mult,4).ne.0) cycle !every 4th
        if (ntest.ge.100)
     &     write(lulog,*) 'Project out component of multiplicity',imult
        kkp1 = dble(imult**2-1)/4d0
        fac1 = -kkp1/(ssp1-kkp1)
        fac2 = 1d0/(ssp1-kkp1)

        ! let S^2 act on amplitudes
        call evaluate2(fspc,.true.,.false.,
     &         op_info,str_info,strmap_info,orb_info,xret,.false.)
        ! Psi -> 1/(S(S+1)-K(K+1)) * (S^2 - K(K+1)) * Psi
        call da_vecsum(ffamp,ffamp%current_record,
     &                 ffamp,ffamp%current_record,fac1,
     &                 ffspc,1,fac2,
     &                 nwfpar,xbuf1,xbuf2,nwfpar)
        ! Pretend that result vector is not up to date
        call reset_file_rec(ffspc)
      end do

      ! normalize
      if (normalize) then
        xret = da_ddot(ffamp,ffamp%current_record,
     &                 ffamp,ffamp%current_record,
     &                 nwfpar,xbuf1,xbuf2,nwfpar)
        if (xret.ge.1d-12) then
          call da_sccpvec(ffamp,ffamp%current_record,
     &                    ffamp,ffamp%current_record,
     &                    1d0/sqrt(xret),nwfpar,xbuf1,nwfpar)
          xret = 1d0
        end if
      end if
c dbg
c      print *,'vector after spin projection:'
c      call vec_from_da(ffamp,ffamp%current_record,
c     &                 xbuf1,nwfpar)
c      do iblk = 1, nwfpar
c        print *,iblk,xbuf1(iblk)
c      end do
c dbgend

      return
      end
