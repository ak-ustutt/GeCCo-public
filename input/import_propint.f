*----------------------------------------------------------------------*
      subroutine import_propint(plist,prop_type,env_type,
     &                             str_info,orb_info)
*----------------------------------------------------------------------*
*     import one electron integrals from MOLPRO or DALTON environment
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
     &     env_type, prop_type
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      type(filinf) ::
     &     ffcmo, ffao
      integer ::
     &     psym
      logical ::
     &     trplt
      character(len=30) ::
     &     labelAll

      ! obtain cmo coefficients on file ffcmo
      call file_init(ffcmo,'CMO.da',ftyp_da_unf,lblk_da)
      call file_open(ffcmo)
      
      select case(env_type(1:1))
      case('m','M')
        call import_cmo_molpro(ffcmo,orb_info)
      case ('d','D')
        call import_cmo_dalton(ffcmo,orb_info)
        if (orb_info%caborb.gt.0)
     &       call import_cmox_dalton_special(ffcmo,orb_info)
      end select

      ! operator information
      call get_oper_info(env_type(1:1),prop_type,labelAll,psym,trplt)

      ! read property integrals in SAO basis
      call file_init(ffao,'PRAO.da',ftyp_da_unf,lblk_da)
      call file_open(ffao)

      select case(env_type(1:1))
      case('m','M')
        call import_propao_molpro(ffao,labelAll,plist%gamt,orb_info)
      case('d','D')
        call import_propao_dalton(ffao,labelAll,plist%gamt,
     &       psym,orb_info)
        end select

      ! transfrom to MO basis
      call tran_one(plist,ffao,ffcmo,trplt,orb_info)

      call file_delete(ffcmo)
      call file_delete(ffao)
      
      return
      end subroutine import_propint
