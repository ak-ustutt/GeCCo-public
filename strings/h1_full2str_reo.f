*----------------------------------------------------------------------*
      subroutine h1_full2str_reo(eref,h1tri,h1str,hlist,
     &     str_info,orb_info)
*----------------------------------------------------------------------*
*     reorder 1-electron part, h1 is given in upper triangular form
*     without symmetry blocking
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_orbinf.h'
      include 'ifc_baserout.h'
      include 'hpvxseq.h'

      real(8), intent(in) ::
     &     eref, h1tri(*)
      real(8), intent(out) ::
     &     h1str(*)
      type(me_list), intent(in) ::
     &     hlist
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info

      integer ::
     &     ngas, nsym, igas, nspin, ispin, nblk, ntoob, ncaborb,
     &     idxstr, idxtri, iblk, ihpv_c, ihpv_a, igc, iga,
     &     ihpv_o, ihpv_i,
     &     ms, idxms, isym, lena, lenc, igas_o, igas_i,
     &     mosto, mondo, mosti, mondi, imo_o, imo_i, imo, jmo,
     &     ica_i, ica_o

      integer, pointer ::
     &     mostnd(:,:,:), hpvx_gas(:,:), iad_gas(:), ireots(:),
     &     restr(:,:,:,:,:,:), iad_gasca(:,:)

      type(operator), pointer ::
     &     hop

      ! for easy typing (+ efficiency)
      hop => hlist%op
      mostnd => orb_info%mostnd
      hpvx_gas => orb_info%ihpvgas
      iad_gas => orb_info%iad_gas
      ireots => orb_info%ireots
      ntoob = orb_info%ntoob
      ncaborb = orb_info%caborb
      ngas = orb_info%ngas
      nspin = orb_info%nspin
      nsym = orb_info%nsym

      restr => hop%igasca_restr

      allocate(iad_gasca(ngas,2))

      if (ntest.ge.100) then
        write(luout,*) '========================='
        write(luout,*) 'here comes h1_full2str_reo'
        write(luout,*) '========================='
        write(luout,*) 'input matrix'
        call prtrlt(h1tri,ntoob+ncaborb)
      end if      

      ! ----------------------------
      ! loop over occupation classes
      ! ----------------------------
      nblk  = 0
      idxstr = 0
      occ_cls: do iblk = 1, hop%n_occ_cls

        if (ntest.ge.100.and.
     &     max(hop%ica_occ(1,iblk),hop%ica_occ(2,iblk)).le.1)
     &                                                             then
          write(luout,*) 'formal? ',hop%formal_blk(iblk)
          call wrt_occ(luout,hop%ihpvca_occ(1,1,iblk))
        end if

        ! reference energy
        if (max(hop%ica_occ(1,iblk),
     &          hop%ica_occ(2,iblk)).eq.0) then
c          idxstr = hlist%off_op_occ(iblk)+1
          idxstr = idxstr + 1
          h1str(idxstr) = eref
          nblk = nblk+1
          cycle
        end if

        if (hop%formal_blk(iblk)) cycle
        ! else: 1-electron operators only ....
        if (max(hop%ica_occ(1,iblk),hop%ica_occ(2,iblk)).ne.1)
     &       cycle
        nblk = nblk+1

        ! set standard frozen core settings:
        do igas = 1, ngas
          iad_gasca(igas,1) = iad_gas(igas)
          iad_gasca(igas,2) = iad_gas(igas)
        end do
        
        ! get type of space for C/A ...
        ihpv_c = idxlist(1,hop%ihpvca_occ(1:,1,iblk),ngastp,1)
        ihpv_a = idxlist(1,hop%ihpvca_occ(1:,2,iblk),ngastp,1)

        ! ... and actual frozen core settings:
        if (ihpv_c.eq.IHOLE) then
          do igas = 1, ngas
            if (hpvx_gas(igas,1).ne.IHOLE) cycle
            if ( restr(2,igas,1,1,1,iblk).eq.0.or.
     &           (igas.gt.1.and.
     &           restr(1,igas-1,1,1,1,iblk).eq.1)) then
              iad_gasca(igas,1) = 1
            else
              iad_gasca(igas,1) = 2
            end if
          end do
        end if
        if (ihpv_a.eq.IHOLE) then
          do igas = 1, ngas
            if (hpvx_gas(igas,1).ne.IHOLE) cycle
            if ( restr(2,igas,2,1,1,iblk).eq.0.or.
     &           (igas.gt.1.and.
     &           restr(1,igas-1,2,1,1,iblk).eq.1)) then
              iad_gasca(igas,2) = 1
            else
              iad_gasca(igas,2) = 2
            end if
          end do
        end if

        ihpv_i = ihpv_c
        ihpv_o = ihpv_a
        ica_i  = 1
        ica_o  = 2
        if (hpvxseq(ihpv_c).gt.hpvxseq(ihpv_a)) then
          ihpv_i = ihpv_a
          ihpv_o = ihpv_c
          ica_i  = 2
          ica_o  = 1
        end if

        ! get indices of graphs
        igc = hlist%idx_graph(ihpv_c,1,iblk)
        iga = hlist%idx_graph(ihpv_a,2,iblk)
        if (min(iga,igc).le.0)
     &       call quit(1,'h1_sym2str_reo','corrupted idx_graph array')
        
        ! -----------------------------
        ! loop over Ms (actually: 2*Ms)
        ! -----------------------------
        ms_loop: do ms = 1, -1, -2
          
          ! the actual index in arrays:
          ispin = 1
          idxms = 1
          if (ms.eq.-1) idxms = 2
          if (ms.eq.-1.and.nspin.eq.2) ispin = 2

          ! ----------------
          ! loop over IRREPS
          ! ----------------
          gam_loop: do isym = 1, nsym
            
            ! get the lengths of strings for current Ms, IRREP
            lena=str_info%g(iga)%lenstr_gm(isym,idxms)
            lenc=str_info%g(igc)%lenstr_gm(isym,idxms)
            if (lena.le.0.or.lenc.le.0) cycle gam_loop
                          
            ! loop over subspaces which belong to current type
            ! of A space (and which are active)
            do igas_o = 1, ngas
              if (hpvx_gas(igas_o,ispin).ne.ihpv_o) cycle
              if (iad_gasca(igas_o,ica_o).ne.2) cycle
              mosto = mostnd(1,isym,igas_o)
              mondo = mostnd(2,isym,igas_o)

              do imo_o = mosto, mondo

                imo = ireots(imo_o)
            
                ! loop over subspaces which belong to current type
                ! of C space (and which are active)
                do igas_i = 1, ngas
                  if (hpvx_gas(igas_i,ispin).ne.ihpv_i) cycle
                  if (iad_gasca(igas_i,ica_i).ne.2) cycle
                  mosti = mostnd(1,isym,igas_i)
                  mondi = mostnd(2,isym,igas_i)

                  do imo_i = mosti, mondi
                    jmo = ireots(imo_i)

                    idxstr = idxstr+1
                    if (imo.le.jmo) idxtri = jmo*(jmo-1)/2+imo
                    if (imo.gt.jmo) idxtri = imo*(imo-1)/2+jmo
                    
                    h1str(idxstr) = h1tri(idxtri)

                  end do
                end do
              end do
            end do

          end do gam_loop
        end do ms_loop
      end do occ_cls

      if (ntest.ge.100) then
        write(luout,*) 'reordered operator:'
        call wrt_mel_buf(luout,5,h1str,hlist,1,nblk,
     &       str_info,orb_info)
      end if

      deallocate(iad_gasca)

      return
      end
