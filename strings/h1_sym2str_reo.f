*----------------------------------------------------------------------*
      subroutine h1_sym2str_reo(eref,h1sym,h1str,hop,
     &     str_info,orb_info)
*----------------------------------------------------------------------*
*     reorder 1-electron part 
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_operator.h'
      include 'def_orbinf.h'
      include 'ifc_baserout.h'

      real(8), intent(in) ::
     &     eref, h1sym(*)
      real(8), intent(out) ::
     &     h1str(*)
      type(operator), intent(in) ::
     &     hop
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info

      integer ::
     &     ngas, nsym, igas,
     &     idxstr, idxsym, iocc_cls, ihpv_c, ihpv_a, igc, iga,
     &     ms, idxms, isym, isymoff, imo_off, lena, lenc, igasa, igasc,
     &     mosta, monda, mostc, mondc, imo_a, imo_c, imo, jmo

      integer, pointer ::
     &     mostnd(:,:,:), ntoobs(:), ihpvgas(:), iad_gas(:), ireots(:)

      ! for easy typing (+ efficiency)
      mostnd => orb_info%mostnd
      ntoobs => orb_info%ntoobs
      ihpvgas => orb_info%ihpvgas
      iad_gas => orb_info%iad_gas
      ireots => orb_info%ireots
      ngas = orb_info%ngas
      nsym = orb_info%nsym

      if (ntest.ge.100) then
        write(luout,*) '========================='
        write(luout,*) 'here comes h1_sym2str_reo'
        write(luout,*) '========================='
        write(luout,*) 'input matrix'
        call wr_blkmat(h1sym,ntoobs,ntoobs,nsym,1)
      end if      

      ! ----------------------------
      ! loop over occupation classes
      ! ----------------------------
      occ_cls: do iocc_cls = 1, hop%n_occ_cls

        ! reference energy
        if (max(hop%ica_occ(1,iocc_cls),
     &          hop%ica_occ(2,iocc_cls)).eq.0) then
          idxstr = hop%off_op_occ(iocc_cls)+1
          h1str(idxstr) = eref
          cycle
        end if

        ! else: 1-electron operators only ....
        if (max(hop%ica_occ(1,iocc_cls),hop%ica_occ(2,iocc_cls)).ne.1)
     &       cycle
        ! ... but the only normal ones (i.e. no 
        ! external/auxiliary orbitals)
        if (iextr.gt.0.and.max(hop%ihpvca_occ(iextr,1,iocc_cls),
     &                         hop%ihpvca_occ(iextr,2,iocc_cls)).gt.0)
     &       cycle

        ! offset for reordered operator
        idxstr = hop%off_op_occ(iocc_cls)

        ! get type of space for C/A
        ihpv_c = idxlist(1,hop%ihpvca_occ(1,1,iocc_cls),ngastp,1)
        ihpv_a = idxlist(1,hop%ihpvca_occ(1,2,iocc_cls),ngastp,1)

        ! get indices of graphs
        igc = hop%idx_graph(ihpv_c,1,iocc_cls)
        iga = hop%idx_graph(ihpv_a,2,iocc_cls)
        if (min(iga,igc).le.0)
     &       call quit(1,'h1_sym2str_reo','corrupted idx_graph array')
        
        ! -----------------------------
        ! loop over Ms (actually: 2*Ms)
        ! -----------------------------
        ms_loop: do ms = 1, -1, -2
          
          ! the actual index in arrays:
          idxms = 1
          if (ms.eq.-1) idxms = 2

          ! ----------------
          ! loop over IRREPS
          ! ----------------
          isymoff = 0
          imo_off = 0
          gam_loop: do isym = 1, nsym
            
            ! get the lengthes of strings for current Ms, IRREP
            lena=str_info%g(iga)%lenstr_gm(isym,idxms)
            if (lena.le.0) cycle gam_loop
                            
            lenc=str_info%g(igc)%lenstr_gm(isym,idxms)
            if (lenc.le.0) cycle gam_loop
                          
            ! loop over subspaces which belong to current type
            ! of A space (and which are active)
            do igasa = 1, ngas
              if (ihpvgas(igasa).ne.ihpv_a) cycle
              if (iad_gas(igasa).ne.2) cycle
              mosta = mostnd(1,isym,igasa)
              monda = mostnd(2,isym,igasa)

              do imo_a = mosta, monda

                imo = ireots(imo_a) - imo_off
            
                ! loop over subspaces which belong to current type
                ! of C space (and which are active)
                do igasc = 1, ngas
                  if (ihpvgas(igasc).ne.ihpv_c) cycle
                  if (iad_gas(igasc).ne.2) cycle
                  mostc = mostnd(1,isym,igasc)
                  mondc = mostnd(2,isym,igasc)

                  do imo_c = mostc, mondc
                    jmo = ireots(imo_c) - imo_off

                    idxstr = idxstr+1
                    if (imo.le.jmo) idxsym = isymoff+jmo*(jmo-1)/2+imo
                    if (imo.gt.jmo) idxsym = isymoff+imo*(imo-1)/2+jmo
                    
                    h1str(idxstr) = h1sym(idxsym)

                  end do
                end do
              end do
            end do

            isymoff = isymoff + ntoobs(isym)*(ntoobs(isym)+1)/2
            imo_off = imo_off + ntoobs(isym)
          end do gam_loop
        end do ms_loop
      end do occ_cls

      if (ntest.ge.100) then
        write(luout,*) 'reordered operator:'
        ! caution: blocks 1-5 will not always work
        call wrt_op_buf(luout,5,h1str,hop,1,5,
     &       str_info,orb_info)
      end if

      return
      end
