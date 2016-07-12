*----------------------------------------------------------------------*
      subroutine import_propint_dalton(plist,prop_type,psym,trplt,
     &                                 str_info,orb_info)
*----------------------------------------------------------------------*
*     import property integrals from DALTON environment
*     we need:
*      AOPROPER
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'ioparam.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      
      type(me_list), intent(inout) ::
     &     plist
      character(len=*), intent(in) ::
     &     prop_type
      integer, intent(in) ::
     &     psym
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info
      logical, intent(in) ::
     &     trplt

      type(filinf) ::
     &     ffcmo, ffao

      ! obtain cmo coefficients on file ffcmo
      call file_init(ffcmo,'CMO.da',ftyp_da_unf,lblk_da)
      call file_open(ffcmo)
      call import_cmo_dalton(ffcmo,orb_info)
      if (orb_info%caborb.gt.0)
     &     call import_cmox_dalton_special(ffcmo,orb_info)

      ! read property integrals in SAO basis
      call file_init(ffao,'PRAO.da',ftyp_da_unf,lblk_da)
      call file_open(ffao)
      call import_propao_dalton(ffao,prop_type,plist%gamt,psym,orb_info)

      ! transfrom to MO basis
      call tran_one(plist,ffao,ffcmo,trplt,orb_info)


      call file_delete(ffcmo)
      call file_delete(ffao)
      
      return
      end
