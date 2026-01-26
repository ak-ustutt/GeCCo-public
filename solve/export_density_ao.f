*----------------------------------------------------------------------*
      subroutine export_density_ao(ndens,rank,label_den,name_export,
     &     trplt,add_ref,
     &     env_type,op_info,str_info,orb_info)
*----------------------------------------------------------------------*
*
*     for a given list of densities (all have rank "rank"): export
*     all of them in ao basis
*
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'ioparam.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'
      include 'mdef_operator_info.h'
      include 'multd2h.h'

      integer, intent(in) ::
     &     ndens, rank
      logical, intent(in) ::
     &     trplt, add_ref
      character(*), intent(in) ::
     &     label_den(ndens), name_export
      character(*), intent(in) ::
     &     env_type
      type(operator_info) ::
     &     op_info
      type(strinf) ::
     &     str_info  
      type(orbinf) ::
     &     orb_info

      real(8), pointer ::
     &     dao(:)
      integer ::
     &     cmo_type, idens, idxden, ifree, nblkd, isym, jsym
      character ::
     &     mat_name*256
      type(filinf) ::
     &     ffcmo, ffdao, ffout

      integer, external ::
     &     idx_mel_list

      if (env_type(1:6).ne.'MOLPRO') then
          call quit(1,"export_density_ao","only possible for molpro")
      end if

      ! get MO-AO trafo from environment
      call file_init(ffcmo,'CMO',ftyp_da_unf,lblk_da)
      cmo_type = -1
      call import_cmo(ffcmo,cmo_type,env_type,orb_info)

      do idens = 1, ndens

        idxden = idx_mel_list(label_den(idens),op_info)
        if (idxden.le.0)
     &       call quit(1,'prop_evaluate',
     &       'label not found: '//trim(label_den(idens)))

        ! back-transform densities
        if (rank.eq.1) then
          call file_init(ffdao,'DAOtmp',ftyp_da_unf,lblk_da)
          call btran_one(ffdao,ffcmo,trplt,add_ref,
     &         op_info%mel_arr(idxden)%mel,orb_info,str_info)
        else
          call quit(1,'export_density_ao','only rank==1 supported')
        end if

        nblkd = 0
        do isym = 1, orb_info%nsym
           jsym = multd2h(isym,op_info%mel_arr(idxden)%mel%gamt)
           nblkd = nblkd+orb_info%nbas(isym)*orb_info%nbas(jsym)
        end do
      
        ifree = mem_alloc_real(dao,nblkd,'dao')

        call file_open(ffdao)

        call get_vec(ffdao,dao,1,nblkd)

        call file_close_delete(ffdao)

        call file_init(ffout,trim(name_export),ftyp_sq_frm,0)

        mat_name = op_info%mel_arr(idxden)%mel%label
        call write_matrix_molpro(dao,ffout,trim(mat_name),
     &                        "DENSITY CHARGE",
     &                        idens==1,op_info%mel_arr(idxden)%mel%gamt,
     &                        orb_info%nbas,orb_info%nsym)

        ifree =  mem_dealloc('dao')

      end do

      call file_delete(ffcmo)


      return
      end
