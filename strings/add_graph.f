*----------------------------------------------------------------------*
      subroutine add_graph(ihpv,nocc,ica,op_restr,str_info,orb_info)
*----------------------------------------------------------------------*
*     append a new slot for a graph in str_info
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'

      type(strinf), intent(inout), target ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info
      integer, intent(in) ::
     &     ihpv, nocc, ica, op_restr(2,orb_info%ngas,2,2)

      integer ::
     &     idx, igtyp, nnew, max_igtyp, jgas, igas, ngas

      type(graph), pointer ::
     &     new_g(:)
      integer, pointer ::
     &     new_spc_typ(:), new_spc_occ(:), new_gas_restr(:,:,:,:)
      integer, pointer ::
     &     ihpvgas(:)

      ngas = orb_info%ngas
      ihpvgas => orb_info%ihpvgas

      str_info%ngraph = str_info%ngraph+1
      nnew = str_info%ngraph

      ! update spc_occ, spc_typ, gas_restr
      allocate(new_spc_typ(nnew),new_spc_occ(nnew),
     &     new_gas_restr(2,ngas,2,nnew))

      do idx = 1, nnew-1
        new_spc_typ(idx) = str_info%ispc_typ(idx)
        new_spc_occ(idx) = str_info%ispc_occ(idx)
        new_gas_restr(1:2,1:ngas,1:2,idx) =
     &       str_info%igas_restr(1:2,1:ngas,1:2,idx)
      end do

      if (nnew.gt.1)
     &     deallocate(str_info%ispc_occ,str_info%ispc_typ,
     &     str_info%igas_restr)

      new_spc_typ(nnew) = ihpv
      new_spc_occ(nnew) = nocc

      jgas = 1
      igas_loop: do igas = 1, ngas
        if (ihpvgas(igas).ne.ihpv) cycle igas_loop
        new_gas_restr(1:2,jgas,1:2,nnew)
     &       = op_restr(1:2,igas,ica,1:2)
        jgas = jgas+1
      end do igas_loop

      str_info%ispc_typ => new_spc_typ
      str_info%ispc_occ => new_spc_occ
      str_info%igas_restr => new_gas_restr

      igtyp = ngastp*(nocc-1) + ihpv

      str_info%max_igtyp = max(igtyp,str_info%max_igtyp)
      max_igtyp = str_info%max_igtyp
      if (nnew.gt.1) deallocate(str_info%gtab)
      allocate(str_info%gtab(ld_gtab,max_igtyp))

      ! set hash table for finding matching graphs
      call set_hash4gtyp(str_info,max_igtyp)

      ! extend graph array
      allocate(new_g(nnew))
      do idx = 1, nnew-1
        new_g(idx) = str_info%g(idx)
      end do

      if (nnew.gt.1) deallocate(str_info%g)
      str_info%g => new_g

      call set_new_graph(nnew,str_info,orb_info)

      return
      end
