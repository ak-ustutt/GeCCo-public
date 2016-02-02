*----------------------------------------------------------------------*
      subroutine set_user_op2(op,name,type,
     &     occ_def,nblk,njoined,irestr,freeze,min_formal,orb_info)
*----------------------------------------------------------------------*
*     set up occupations for a general operator described by
*     occ_def:   user provided occupations
*     irestr:    restriction on subspaces 
*                min, max. number of operators after completion of
*                subspace within H/P/V/X, for C/A
*     freeze(1): use settings on iadgas to enforce frozen C occupations
*     freeze(2): use settings on iadgas to enforce frozen A occupations
*       new: the array freeze(1:2,1:njoined) contains the max number
*            of active electrons in frozen shells (0 = all frozen)
*            per vertex -> allow for semi-frozen cores
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'def_operator.h'
      include 'def_orbinf.h'
      include 'stdunit.h'
      include 'ifc_baserout.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 00

      type(operator), intent(inout) ::
     &     op
      character(len=*), intent(in) ::
     &     name
      integer, intent(in) ::
     &     type, nblk, njoined, min_formal
      type(orbinf), intent(in), target ::
     &     orb_info
      integer, intent(in) ::
     &     occ_def(ngastp,2,nblk*njoined), irestr(2,orb_info%ngas,2,2)
c      logical, intent(in) ::
c     &     freeze(2)
      integer, intent(in) ::
     &     freeze(2,njoined)

      logical, parameter ::
     &     inv_hole = .false.

      integer ::
     &     idx, irank, na, nc, ica, igas, igasl, idiff, imaxr, nact, 
     &     iocc, igastp, iprint, nx, iblk, ijoin, ioff_blk, gasst,gasnd,
     &     ispin
      integer ::
     &     hpvxprint(orb_info%ngas)
      integer, pointer ::
     &     nspin, ngas, iad_gas(:), hpvxgas(:,:)

      iprint = max(iprlvl,ntest)

      if (iprint.ge.100) then
        call write_title(lulog,wst_dbg_subr,'set_user_op')
        write(lulog,*) 'nblk,njoined: ',nblk,njoined
        call wrt_occ_n(lulog,occ_def,nblk*njoined)
        call wrt_rstr(lulog,irestr,orb_info%ngas)
      end if

      if (len_trim(name).gt.len_opname)
     &    call quit(1,'set_user_op','name too long: "'//trim(name)//'"')

      if (type.ne.optyp_operator.and.type.ne.optyp_density) then
        if (type.eq.optyp_intermediate) then
          write(lulog,*)
     &         'use set_gen_intermediate to define intermediates'
        else
          write(lulog,*) 'type: ',type,' ?'
        end if
        call quit(1,'set_user_op','illegal type specification')
      end if

      if (nblk.le.0) then
        write(lulog,*) 'nblk = ',nblk,' ???'
        call quit(1,'set_user_op2',
     &    'Must at least define one block for "'//trim(name)//'"')
      end if

      nspin => orb_info%nspin
      ngas => orb_info%ngas
      iad_gas => orb_info%iad_gas
      hpvxgas => orb_info%ihpvgas

      ! basic settings:
      op%name = '        '
      op%name = trim(name)

      op%type = type
      op%njoined = njoined

      op%dagger = .false.
      op%formal=.false.

      op%n_occ_cls = nblk

      call init_operator(op,orb_info)

      op%formal_blk(1:nblk) = .false.
      if (min_formal.gt.0.and.min_formal.le.nblk)
     &   op%formal_blk(min_formal:nblk) = .true.

      ! this is basically all:
      op%ihpvca_occ = occ_def
        
      ! well, a few loops for restrictions
      do iblk = 1, nblk
        ioff_blk = (iblk-1)*njoined

        op%ica_occ(1:2,iblk) = 0
        do ijoin = 1, njoined
          op%ica_occ(1,iblk) = op%ica_occ(1,iblk) +
     &         ielsum(op%ihpvca_occ(1:,1,ioff_blk+ijoin),ngastp)
          op%ica_occ(2,iblk) = op%ica_occ(2,iblk) +
     &         ielsum(op%ihpvca_occ(1:,2,ioff_blk+ijoin),ngastp)
        end do

        ! set restrictions
        do ijoin = 1, njoined
         do ica = 1, 2
           do igastp = 1, ngastp
             gasst = orb_info%idx_gas(igastp)
             gasnd = gasst-1+orb_info%ngas_hpv(igastp)
             do igas = gasst, gasnd
               ! set a/c rank as upper bound
               idiff = -irestr(1,igas,ica,1)+irestr(2,igas,ica,1)
               imaxr = min(irestr(2,igas,ica,1),
     &           op%ihpvca_occ(hpvxgas(igas,1),ica,ioff_blk+ijoin))
               ! not sure whether this will work in all cases:
               if (igas.lt.gasnd) then
                 op%igasca_restr(1,igas,ica,1,1:nspin,ioff_blk+ijoin) =
     &                max(0,imaxr - idiff)
               else
                 op%igasca_restr(1,igas,ica,1,1:nspin,ioff_blk+ijoin)=
     &                imaxr
               end if
               op%igasca_restr(2,igas,ica,1,1:nspin,ioff_blk+ijoin) =
     &              imaxr                    
             end do
           end do
         end do
        end do
c        ! set restrictions (OLD)
c        do ijoin = 1, njoined
c         do ica = 1, 2
c          do igas = 1, ngas
c            ! set a/c rank as upper bound
c            idiff = - irestr(1,igas,ica,1)+irestr(2,igas,ica,1)
c            imaxr = min(irestr(2,igas,ica,1),
c     &           op%ihpvca_occ(hpvxgas(igas,1),ica,ioff_blk+ijoin))
c            ! not sure whether this will work in all cases:
c            if (igas.lt.ngas) then
c              op%igasca_restr(1,igas,ica,1,1:nspin,ioff_blk+ijoin) =
c     &             max(0,imaxr - idiff)
c            else
c              op%igasca_restr(1,igas,ica,1,1:nspin,ioff_blk+ijoin)=imaxr
c            end if
cc very quick fix (HM, A COPIED ONE, ACTUALLY):
c            op%igasca_restr(1:2,igas,ica,1,1:nspin,ioff_blk+ijoin) =
c     &           imaxr                    
c          end do
c         end do
c        end do
        ! post-processing for frozen shells
        do ijoin = 1, njoined
         do ica = 1, 2
          nact = freeze(ica,ijoin)
c          if (.not.freeze(ica)) cycle
          do igas = 1, ngas
            if (iad_gas(igas).ne.2) then
              if (igas.eq.1) then
                do ispin = 1, nspin
                  imaxr = min(nact,
     &          op%igasca_restr(2,igas,ica,1,ispin,ioff_blk+ijoin))
                  op%igasca_restr(1:2,igas,ica,1,ispin,ioff_blk+ijoin)=
     &               (/0,imaxr/)
                end do
              else
C                call quit(1,'set_user_op2','ever accessed this part?')
                op%igasca_restr(1,igas,ica,1,1:nspin,ioff_blk+ijoin) =
     &            op%igasca_restr(2,igas,ica,1,1:nspin,ioff_blk+ijoin) 
                op%igasca_restr(1,igas-1,ica,1,1:nspin,ioff_blk+ijoin) =
     &            op%igasca_restr(2,igas-1,ica,1,1:nspin,ioff_blk+ijoin) 
              end if
            end if
          end do
         end do
        end do
        ! mask restriction currently unused
        op%igasca_restr(1:2,1:ngas,1:2,2,1:nspin,
     &                  ioff_blk+1:ioff_blk+njoined) = 0
        ! distinguish similar blocks by increasing version number
        do idx = 1, iblk-1
          if (all(op%ihpvca_occ(:,:,(idx-1)*njoined+1:idx*njoined)-
     &        op%ihpvca_occ(:,:,(iblk-1)*njoined+1:iblk*njoined).eq.0))
     &        op%blk_version(iblk) = op%blk_version(iblk) + 1
        end do

      end do

      if (iprint.ge.2)
     &       write(lulog,'(x,3a,i4)')
     &       'Number of occupation classes for ',
     &       trim(name),': ',op%n_occ_cls

      if (iprint.ge.5) then
        write(lulog,*) 'According to your wishes, I set the following:'
        call print_op_occ(lulog,op)
      end if
      
      return
      end
*----------------------------------------------------------------------*
