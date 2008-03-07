*----------------------------------------------------------------------*
      subroutine set_user_op(op,name,type,
     &     dagger,
     &     occ_def,nblk,irestr,orb_info)
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
      logical, intent(in) ::
     &     dagger
      integer, intent(in) ::
     &     type, nblk
      type(orbinf), intent(in), target ::
     &     orb_info
      integer, intent(in) ::
     &     occ_def(ngastp,2,nblk), irestr(2,orb_info%ngas,2,2)

      logical, parameter ::
     &     inv_hole = .true.

      integer ::
     &     ifree, ipass, irank, na, nc, ica, igas, igasl, idiff, imaxr,
     &     iocc, igastp, iprint, nx, iblk
      integer ::
     &     hpvxprint(ngastp)
      integer, pointer ::
     &     nspin, ngas, iad_gas(:), hpvxgas(:,:)

      iprint = max(iprlvl,ntest)

      if (iprint.ge.100) then
        call write_title(luout,wst_dbg_subr,'set_user_op')
        call wrt_occ_n(luout,occ_def,nblk)
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
      op%njoined = 1  ! always for operators and densities

      op%dagger = dagger
      op%formal=.true.

      op%n_occ_cls = nblk

      call init_operator(op,orb_info)

      op%formal_blk(1:nblk) = .true.

      ! this is basically all:
      op%ihpvca_occ = occ_def
        
      ! well, a few loops for restrictions
      do iblk = 1, nblk

        op%ica_occ(1,iblk) = ielsum(op%ihpvca_occ(1,1,iblk),ngastp)
        op%ica_occ(2,iblk) = ielsum(op%ihpvca_occ(1,2,iblk),ngastp)

        ! set restrictions
        do ica = 1, 2
          do igas = 1, ngas
            ! set a/c rank as upper bound
            idiff = - irestr(1,igas,ica,1)+irestr(2,igas,ica,1)
            imaxr = min(irestr(2,igas,ica,1),
     &           op%ihpvca_occ(hpvxgas(igas,1),ica,iblk))
            ! not sure whether this will work in all cases:
            if (igas.lt.ngas) then
              op%igasca_restr(1,igas,ica,1,1:nspin,iblk) =
     &             max(0,imaxr - idiff)
            else
              op%igasca_restr(1,igas,ica,1,1:nspin,iblk) = imaxr
            end if
c very quick fix:                    
            op%igasca_restr(1:2,igas,ica,1,1:nspin,iblk) =
     &           imaxr                    
          end do
        end do
        ! post-processing for frozen shells
        do ica = 1, 2
          do igas = 1, ngas
            if (iad_gas(igas).ne.2) then
              if (igas.eq.1) then                        
                op%igasca_restr(1:2,igas,ica,1,1:nspin,iblk) = 0
              else
                op%igasca_restr(1,igas,ica,1,1:nspin,iblk) =
     &               op%igasca_restr(2,igas,ica,1,1:nspin,iblk) 
                op%igasca_restr(1,igas-1,ica,1,1:nspin,iblk) =
     &               op%igasca_restr(2,igas-1,ica,1,1:nspin,iblk) 
              end if
            end if
          end do
        end do
        ! mask restriction currently unused
        op%igasca_restr(1:2,1:ngas,1:2,2,1:nspin,iblk) = 0


      end do

      if (iprint.ge.2)
     &       write(luout,'(x,3a,i4)')
     &       'Number of occupation classes for ',
     &       trim(name),': ',op%n_occ_cls

      if (iprint.ge.20) then

        do igas = 1, ngas
          hpvxprint(igas) = igas
        end do
        if (inv_hole) then
          igasl = 0
          do igas = ngas, 1, -1
            if (hpvxgas(igas,1).eq.1) then
              igasl = igasl+1
              hpvxprint(igas) = igasl
            end if
          end do
        end if
        do iocc = 1, op%n_occ_cls
          write(luout,'(/x,a,i4)') 'Occupation Nr. ',iocc
          call wrt_occ(luout,op%ihpvca_occ(1,1,iocc))
          write(luout,'(/4x,6(2x,i2,x))') hpvxprint(1:ngas)
          call wrt_rstr(luout,op%igasca_restr(1,1,1,1,1,iocc),ngas)
        end do
      end if
      
      return
      end
*----------------------------------------------------------------------*
