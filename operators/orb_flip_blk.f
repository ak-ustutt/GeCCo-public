*----------------------------------------------------------------------*      
      subroutine orb_flip_blk(buffer,
     &     mel,iblk,igam,idxms,idxdis,nel,flipsigns,str_info,orb_info)
*----------------------------------------------------------------------*      
*     flip the signs of the ME-lists' elements according to
*     the flipsign assigned to each orbital.
*
*     Matthias, Nov 2013
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

      integer, intent(in) ::
     &     iblk, igam, idxms, idxdis, nel
      type(me_list), intent(in) ::
     &     mel
      real(8), intent(inout) ::
     &     buffer(*)
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info
      logical, intent(in) ::
     &     flipsigns(orb_info%ntoob)

      logical ::
     &     first
      integer ::
     &     did, idxstr, ioff, idx_occ, njoined, ijoin
      integer ::
     &     msdst(ngastp,2,mel%op%njoined),
     &     igamdst(ngastp,2,mel%op%njoined),
     &     lexlscr(nel,3),idorb(nel), idspn(nel), idspc(nel)

      type(operator), pointer ::
     &     op

      logical, external ::
     &     next_tupel_ca

      njoined = mel%op%njoined
      idx_occ = (iblk-1)*njoined+1
      op => mel%op

      ioff = mel%off_op_gmox(iblk)%d_gam_ms(idxdis,igam,idxms)
      ! get Ms and IRREP distribution info from
      ! distribution ID
      did = mel%off_op_gmox(iblk)%did(idxdis,igam,idxms)
      call did2msgm(msdst,igamdst,did,
     &             op%ihpvca_occ(1,1,idx_occ),orb_info%nsym,njoined)
      
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

        ! flip sign if odd number of flipped orbital indices
c dbg
c        print *,'idxstr+ioff,flsgns,mod:',idxstr+ioff,
c     &    flipsigns(idorb(1:nel)),mod(count(flipsigns(idorb(1:nel))),2)
c        print *,'buffer:',buffer(idxstr)
c dbgend
        if (mod(count(flipsigns(idorb(1:nel))),2).eq.1)
     &     buffer(idxstr) = -buffer(idxstr)

      end do

      return
      end
