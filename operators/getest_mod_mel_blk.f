*----------------------------------------------------------------------*      
      subroutine getest_mod_mel_blk(buffer,mel,
     &     iRdef,norb,icase,
     &     iblk,igam,idxms,idxdis,nel,str_info,orb_info)
*----------------------------------------------------------------------*      
*
*  core routine for getest_mod_mel
*
*  andreas, feb 2012
*----------------------------------------------------------------------*      

      implicit none

      include 'opdim.h'
      include 'hpvxseq.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'multd2h.h'

      integer, intent(in) ::
     &     iblk, igam, idxms, idxdis, nel, norb, iRdef(norb),
     &     icase
      type(me_list), intent(inout) ::
     &     mel
      real(8), intent(inout) ::
     &     buffer(*)
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      logical ::
     &     first
      integer ::
     &     did, idxstr, iel, ioff, msa, ihpv,
     &     idx_occ, njoined, idxst, idxnd, iadd, ijoin,
     &     iorb, nmatch
      integer ::
     &     msdst(ngastp,2,mel%op%njoined),
     &     igamdst(ngastp,2,mel%op%njoined),
     &     lexlscr(nel,3),idorb(nel), idspn(nel), idspc(nel),
     &     nelc(mel%op%njoined), nela(mel%op%njoined)

      type(operator), pointer ::
     &     op

      logical, external ::
     &     next_tupel_ca

      njoined = mel%op%njoined
      idx_occ = (iblk-1)*njoined+1
      op => mel%op

      !ioff = mel%off_op_gmox(iblk)%d_gam_ms(idxdis,igam,idxms)
      ! get Ms and IRREP distribution info from
      ! distribution ID
      did = mel%off_op_gmox(iblk)%did(idxdis,igam,idxms)
      call did2msgm(msdst,igamdst,did,
     &               op%ihpvca_occ(1,1,idx_occ),orb_info%nsym,njoined)
        
      ! loop over all possible index tuples:
      first = .true.
      idxstr = 0
      do while(next_tupel_ca(idorb,idspn,idspc,
     &     nel,njoined,op%ihpvca_occ(1,1,idx_occ),
     &     mel%idx_graph(1,1,idx_occ),
     &     msdst,igamdst,first,
     &     str_info%igas_restr,
     &     orb_info%mostnd,orb_info%igamorb,
     &     orb_info%nsym,orb_info%ngas,
     &     orb_info%ngas_hpv,orb_info%idx_gas,
     &     hpvxseq,lexlscr))
        first = .false.
        idxstr = idxstr+1

        ! compare indices on idorb with indices of list R
        ! nel indices on list: subsystem R
        ! 0   indices on list: subsystem S
        ! 1-nel-1 indices on list: coupling element 
        nmatch = 0
        do iel = 1, nel
          do iorb = 1, norb
            if (idorb(iel).eq.iRdef(iorb)) nmatch = nmatch+1
          end do
        end do

        if (nmatch.eq.0  .and.icase.ne.2 .or.
     &      nmatch.eq.nel.and.icase.ne.3) cycle

        buffer(idxstr) = 0d0

      end do

      return
      end
