*----------------------------------------------------------------------*
      subroutine set_gen_intermediate(op,name,
     &     defop,ndefop,orb_info)
*----------------------------------------------------------------------*
*     set general intermediate
*     given are ndefop operator definitions on the 
*     pointer array defop
*     the first operator defines the external lines above and below
*     the intermediate
*     the remaining operators define (top to bottom in a diagram,
*     left to right in an operator string) the operators which will be
*     inserted into the intermediate;
*     one may think of the intermediate originating from an expression
*     ---originally yielding operator one as result---after derivatives
*     with respect to operators 2 to ndefop (e.g. operator 1 may be the
*     energy and the remaining operators the CC or R12 coefficients).
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_operator.h'
      include 'def_operator_array.h'
      include 'def_orbinf.h'
      include 'stdunit.h'
      include 'ifc_baserout.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 00

      type(operator), intent(inout) ::
     &     op
      character, intent(in) ::
     &     name*(*)
      integer, intent(in) ::
     &     ndefop
      type(operator_array), intent(in) ::
     &     defop(ndefop)

      type(orbinf), intent(in), target ::
     &     orb_info

      integer ::
     &     iblk, ioffblk, ijoin, idef, njoined
      integer, pointer ::
     &     hpvxgas(:,:), ngas, nspin
      integer ::
     &     iblk_min(ndefop), iblk_max(ndefop), iblk_dis(ndefop)

      if (ntest.ge.100)
     &     call write_title(lulog,wst_dbg_subr,
     &     'set_gen_intermediate speaking')
      
      nspin => orb_info%nspin
      ngas => orb_info%ngas
      hpvxgas => orb_info%ihpvgas

      if (len_trim(name).gt.len_opname)
     &     call quit(1,'set_gen_intermediate',
     &     'name too long: "'//trim(name)//'"')

      op%name = '        '
      op%name = name
      
      op%type = optyp_intermediate
      op%njoined = ndefop
      njoined = ndefop

      if (ndefop.le.0) then
        write(lulog,*) 'ndefop: ', ndefop
        call quit(1,'set_gen_intermediate',
     &       'ndefop has no sensible value')
      end if

      op%dagger = defop(1)%op%dagger
c      dagtotal = op%dagger
      if (op%dagger)
     &     call quit(1,'set_gen_intermediate',
     &     'not yet decided how to handle adjungate operators')

      op%n_occ_cls = defop(1)%op%n_occ_cls
      iblk_max(1)  = defop(1)%op%n_occ_cls
      do idef = 2, ndefop
        op%n_occ_cls = op%n_occ_cls*defop(idef)%op%n_occ_cls
        iblk_max(idef) = defop(idef)%op%n_occ_cls
      end do

      call init_operator(op,orb_info)

      iblk_min(1:ndefop) = 1
      
      iblk_dis = iblk_min

      iblk = 0
      do
        iblk = iblk+1
        if (iblk.gt.op%n_occ_cls)
     &       call quit(1,'set_gen_intermediate','something is weird')

        op%formal_blk(iblk) = .false.

        ioffblk = (iblk-1)*njoined
        
        op%ihpvca_occ(1:ngastp,1:2,ioffblk+1:ioffblk+njoined) = 0

        ! first operator contributes to first and last
        op%ihpvca_occ(1:ngastp,1:2,ioffblk+1) =
     &      iocc_xdn(1,defop(1)%op%ihpvca_occ(1:ngastp,1:2,iblk_dis(1)))
        op%ihpvca_occ(1:ngastp,1:2,ioffblk+njoined) =
     &       op%ihpvca_occ(1:ngastp,1:2,ioffblk+njoined) +
     &      iocc_xdn(2,defop(1)%op%ihpvca_occ(1:ngastp,1:2,iblk_dis(1)))

        do idef = 2, ndefop
          op%ihpvca_occ(1:ngastp,1:2,ioffblk+idef-1) =
     &         op%ihpvca_occ(1:ngastp,1:2,ioffblk+idef-1) +
     &         iocc_dagger(iocc_xdn(1,defop(idef)%op%
     &                    ihpvca_occ(1:ngastp,1:2,iblk_dis(idef))))
          op%ihpvca_occ(1:ngastp,1:2,ioffblk+idef) =
     &         op%ihpvca_occ(1:ngastp,1:2,ioffblk+idef) +
     &         iocc_dagger(iocc_xdn(2,defop(idef)%op%
     &                    ihpvca_occ(1:ngastp,1:2,iblk_dis(idef))))
        end do

        op%ica_occ(1:2,iblk) = 0

        do ijoin = 1, njoined
          op%ica_occ(1,iblk) = op%ica_occ(1,iblk)+
     &         ielsum(op%ihpvca_occ(1:ngastp,1,ioffblk+ijoin),ngastp)
          op%ica_occ(2,iblk) = op%ica_occ(2,iblk)+
     &         ielsum(op%ihpvca_occ(1:ngastp,2,ioffblk+ijoin),ngastp)
        end do

        ! dto. for restrictions
        op%igasca_restr(1:2,1:ngas,1:2,1:2,1:nspin,
     &                  ioffblk+1:ioffblk+njoined) = 0

        op%igasca_restr(1:2,1:ngas,1:2,1:2,1:nspin,ioffblk+1) =
     &       irest_xdn(1,defop(1)%op%
     &             igasca_restr(1:2,1:ngas,1:2,1:2,1:nspin,iblk_dis(1)),
     &                 hpvxgas,ngas,nspin)
        op%igasca_restr(1:2,1:ngas,1:2,1:2,1:nspin,ioffblk+njoined) =
     &      op%igasca_restr(1:2,1:ngas,1:2,1:2,1:nspin,ioffblk+njoined)+
     &       irest_xdn(2,defop(1)%op%
     &             igasca_restr(1:2,1:ngas,1:2,1:2,1:nspin,iblk_dis(1)),
     &                 hpvxgas,ngas,nspin)
        do idef = 2, ndefop
          op%igasca_restr(1:2,1:ngas,2:1:-1,1:2,1:nspin,ioffblk+idef-1)=
     &    op%igasca_restr(1:2,1:ngas,2:1:-1,1:2,1:nspin,ioffblk+idef-1)+
     &       irest_xdn(1,defop(idef)%op%
     &       igasca_restr(1:2,1:ngas,1:2,1:2,   1:nspin,iblk_dis(idef)),
     &                 hpvxgas,ngas,nspin)

          op%igasca_restr(1:2,1:ngas,2:1:-1,1:2,1:nspin,ioffblk+idef) =
     &    op%igasca_restr(1:2,1:ngas,2:1:-1,1:2,1:nspin,ioffblk+idef) +
     &       irest_xdn(2,defop(idef)%op%
     &       igasca_restr(1:2,1:ngas,1:2,1:2,   1:nspin,iblk_dis(idef)),
     &                 hpvxgas,ngas,nspin)
        end do

        if (.not.next_dist2(iblk_dis,ndefop,iblk_min,iblk_max,1)) exit
      end do

      if (ntest.ge.100) then
        write(lulog,*) 'Intermediate: ',trim(op%name)
        write(lulog,*) 'generated ',op%n_occ_cls,' blocks'
        do iblk = 1, op%n_occ_cls
          write(lulog,'(/x,a,i4)') 'Occupation Nr. ',iblk
          ioffblk = (iblk-1)*njoined
          call wrt_occ_n(lulog,op%ihpvca_occ(1,1,ioffblk+1),njoined)
          do ijoin = 1, njoined
            call wrt_rstr(lulog,
     &           op%igasca_restr(1,1,1,1,1,ioffblk+ijoin),ngas)
          end do
        end do
      end if

      return
      end
