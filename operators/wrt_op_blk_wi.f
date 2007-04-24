*----------------------------------------------------------------------*      
      subroutine wrt_op_blk_wi(luout,buffer,
     &     op,iblk,igam,idxms,idxdis,nel,str_info,orb_info)
*----------------------------------------------------------------------*      
*     write opertator block to unit luout with complete index for
*     each element
*----------------------------------------------------------------------*      

      implicit none

      include 'opdim.h'
      include 'hpvxseq.h'
      include 'def_operator.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'multd2h.h'

      integer, intent(in) ::
     &     luout, iblk, igam, idxms, idxdis, nel
      type(operator), intent(in) ::
     &     op
      real(8), intent(in) ::
     &     buffer(*)
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      logical ::
     &     first
      integer ::
     &     did, idxstr, iel, nelc, nela, ioff, msa, ihpv
      integer ::
     &     msdst(ngastp,2), igamdst(ngastp,2), lexlscr(nel,3),
     &     idorb(nel), idspn(nel), idspc(nel)
      character ::
     &     spnstr*(nel+1),fmtstr*56

      logical, external ::
     &     next_tupel_ca

      nelc = op%ica_occ(1,iblk)
      nela = op%ica_occ(2,iblk)
      if (op%off_op_gmox(iblk)%maxd.eq.1) then
        ! we cannot be sure that the extended information on
        ! gmox is set for this operator, so ...
        ioff = op%off_op_gmo(iblk)%gam_ms(igam,idxms)
        ! ... and we set the Ms and IRREP distributions
        ! ourselves .... 
        msdst(1:ngastp,1:2) = 0
        igamdst(1:ngastp,1:2) = 0
        msa = nela - (idxms-1)*2
        do ihpv = 1, ngastp
          if (op%ihpvca_occ(ihpv,1,iblk).gt.0) then
            msdst(ihpv,1) = msa + op%mst
            igamdst(ihpv,1) = multd2h(igam,op%gamt)
          end if
          if (op%ihpvca_occ(ihpv,2,iblk).gt.0) then
            msdst(ihpv,2) = msa
            igamdst(ihpv,2) = igam
          end if
        end do
      else
        ioff = op%off_op_gmox(iblk)%d_gam_ms(idxdis,igam,idxms)
        ! get Ms and IRREP distribution info from
        ! distribution ID
        did = op%off_op_gmox(iblk)%did(idxdis,igam,idxms)
        call did2msgm(msdst,igamdst,did,
     &               op%ihpvca_occ(1,1,iblk),orb_info%nsym)
      end if
        

      write(fmtstr,'("(x,i5,2x,",i2,"i3,x,",i2,"i3,2x,a,2x,",'//
     &                       'i2,"i1,x,",i2,"i1,x,g20.10)")')
     &     nelc,nela,nelc,nela
      
      spnstr(1:nel+1) = ' '  

      ! loop over all possible index tuples:
      first = .true.
      idxstr = 0
      do while(next_tupel_ca(idorb,idspn,idspc,
     &     nel,op%ihpvca_occ(1,1,iblk),
     &     op%idx_graph(1,1,iblk),
     &     msdst,igamdst,first,
     &     str_info%igas_restr,
     &     orb_info%mostnd,orb_info%igamorb,
     &     orb_info%nsym,orb_info%ngas,
     &     orb_info%ngas_hpv,orb_info%idx_gas,
     &     hpvxseq,lexlscr))
        first = .false.
        idxstr = idxstr+1
      
        do iel = 1, nelc
          if (idspn(iel).eq.1) write(spnstr(iel:iel),'("+")')
          if (idspn(iel).eq.-1) write(spnstr(iel:iel),'("-")')
          if (idspn(iel).eq.2) write(spnstr(iel:iel),'("2")')
        end do
        do iel = nelc+1, nelc+nela
          if (idspn(iel).eq.1) write(spnstr(iel+1:iel+1),'("+")')
          if (idspn(iel).eq.-1) write(spnstr(iel+1:iel+1),'("-")')
          if (idspn(iel).eq.2) write(spnstr(iel+1:iel+1),'("2")')
        end do
        write(luout,fmtstr) idxstr+ioff,
     &       idorb(1:nel), spnstr, idspc(1:nel),
     &       buffer(idxstr)

      end do

      return
      end
