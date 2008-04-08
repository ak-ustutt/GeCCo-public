*----------------------------------------------------------------------*
      subroutine set_user_op2(op,name,type,
     &     occ_def,nblk,njoined,irestr,orb_info)
*----------------------------------------------------------------------*
*     set up occupations for a general operator described by
*     occ_def:   user provided occupations
*     irestr:    restriction on subspaces 
*                min, max. number of operators after completion of
*                subspace within H/P/V/X, for C/A
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
      character, intent(in) ::
     &     name*(*)
      integer, intent(in) ::
     &     type, nblk, njoined
      type(orbinf), intent(in), target ::
     &     orb_info
      integer, intent(in) ::
     &     occ_def(ngastp,2,nblk*njoined), irestr(2,orb_info%ngas,2,2)

      logical, parameter ::
     &     inv_hole = .false.

      integer ::
     &     ifree, ipass, irank, na, nc, ica, igas, igasl, idiff, imaxr,
     &     iocc, igastp, iprint, nx, iblk, ijoin, ioff_blk
      integer ::
     &     hpvxprint(ngastp)
      integer, pointer ::
     &     nspin, ngas, iad_gas(:), hpvxgas(:,:)

      iprint = max(iprlvl,ntest)

      if (iprint.ge.100) then
        call write_title(luout,wst_dbg_subr,'set_user_op')
        call wrt_occ_n(luout,occ_def,nblk*njoined)
        call wrt_rstr(luout,irestr,orb_info%ngas)
      end if

      if (len_trim(name).gt.len_opname)
     &    call quit(1,'set_user_op','name too long: "'//trim(name)//'"')

      if (type.ne.optyp_operator.and.type.ne.optyp_density) then
        if (type.eq.optyp_intermediate) then
          write(luout,*)
     &         'use set_gen_intermediate to define intermediates'
        else
          write(luout,*) 'type: ',type,' ?'
        end if
        call quit(1,'set_user_op','illegal type specification')
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

      ! this is basically all:
      op%ihpvca_occ = occ_def
        
      ! well, a few loops for restrictions
      do iblk = 1, nblk
        ioff_blk = (iblk-1)*njoined

        op%ica_occ(1:2,iblk) = 0
        do ijoin = 1, njoined
          op%ica_occ(1,iblk) = op%ica_occ(1,iblk) +
     &         ielsum(op%ihpvca_occ(1,1,ioff_blk+ijoin),ngastp)
          op%ica_occ(2,iblk) = op%ica_occ(2,iblk) +
     &         ielsum(op%ihpvca_occ(1,2,ioff_blk+ijoin),ngastp)
        end do

        ! set restrictions
        do ijoin = 1, njoined
         do ica = 1, 2
          do igas = 1, ngas
            ! set a/c rank as upper bound
            idiff = - irestr(1,igas,ica,1)+irestr(2,igas,ica,1)
            imaxr = min(irestr(2,igas,ica,1),
     &           op%ihpvca_occ(hpvxgas(igas,1),ica,ioff_blk+ijoin))
            ! not sure whether this will work in all cases:
            if (igas.lt.ngas) then
              op%igasca_restr(1,igas,ica,1,1:nspin,ioff_blk+ijoin) =
     &             max(0,imaxr - idiff)
            else
              op%igasca_restr(1,igas,ica,1,1:nspin,ioff_blk+ijoin)=imaxr
            end if
c very quick fix (HM, A COPIED ONE, ACTUALLY):
            op%igasca_restr(1:2,igas,ica,1,1:nspin,ioff_blk+ijoin) =
     &           imaxr                    
          end do
         end do
        end do
        ! post-processing for frozen shells
        do ijoin = 1, njoined
         do ica = 1, 2
          do igas = 1, ngas
            if (iad_gas(igas).ne.2) then
              if (igas.eq.1) then                        
                op%igasca_restr(1:2,igas,ica,1,1:nspin,ioff_blk+ijoin)=0
              else
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

      end do

      if (iprint.ge.2)
     &       write(luout,'(x,3a,i4)')
     &       'Number of occupation classes for ',
     &       trim(name),': ',op%n_occ_cls

      if (iprint.ge.20) then
        write(luout,*) 'According to your wishes, I set the following:'
        call print_op_occ(luout,op)
      end if
      
      return
      end
*----------------------------------------------------------------------*
