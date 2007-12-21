*----------------------------------------------------------------------*
      subroutine prop_evaluate(ndens,rank,label_den,
     &     env_type,op_info,str_info,orb_info)
*----------------------------------------------------------------------*
*
*     for a given list of densities (all have rank "rank") evaluate
*     all properties available in the environment
*
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'ioparam.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      
      integer, intent(in) ::
     &     ndens, rank
      character(*), intent(in) ::
     &     label_den(ndens)
      character(*), intent(in) ::
     &     env_type
      type(operator_info) ::
     &     op_info
      type(strinf) ::
     &     str_info  
      type(orbinf) ::
     &     orb_info

      integer ::
     &     cmo_type, idens, idxden
      character ::
     &     label*8
      type(filinf) ::
     &     ffcmo, ffdao, ffprop

      integer, external ::
     &     idx_mel_list

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
          call file_init(ffdao,'DAO',ftyp_da_unf,lblk_da)
          call btran_one(ffdao,ffcmo,.true.,
     &         op_info%mel_arr(idxden)%mel,orb_info)
        else
          call quit(1,'prop_evaluate','only rank==1 supported')
        end if

        ! calculate trace with one-electron integrals
        ! provided by environment
        call oneprop_ao(ffdao,op_info%mel_arr(idxden)%mel,
     &       env_type,orb_info)

      end do

      call file_delete(ffcmo)


      return
      end
