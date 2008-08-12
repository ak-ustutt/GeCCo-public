*----------------------------------------------------------------------*
      integer function idx4sg(nel,idspc,idorb,idspn,idgam,
     &     iy4sg,iy_info,iyssg,iwssg,ioff_dgm,ndis,mostnd,
     &     nelmax,ngam,nspc)
*----------------------------------------------------------------------*
*     get index of string of length nel, with orbitals idorb(1:n),
*     spins idspn(1:n), orbital IRREPS idgam(1:n), and subspace distribution
*     idspc(1:n)
*----------------------------------------------------------------------*

      implicit none

      integer, intent(in) ::
     &     nel, idspc(nel), idorb(nel), idspn(nel), idgam(nel),
     &     ngam, nspc, nelmax,ndis,
     &     iy4sg(*), iy_info(3,*),
     &     iyssg(nelmax,nspc), iwssg(0:nelmax,nspc),
     &     mostnd(2,ngam,nspc),
     &     ioff_dgm(ndis,ngam,*)

      integer ::
     &     idx, idx_sg, nn, nm, ig, len, lenprev, ispc,
     &     ioff, iooff, nelmax_sg, msmax_sg

      integer, external ::
     &     idx4sg_sb, idxssg, idxssd, lensubspc
*----------------------------------------------------------------------*

      nm = 0
      ig = 1
      nn = 0
      idx = 0

      lenprev = 0
      do ispc = 1, nspc
        if (ispc.eq.1) then
          idx_sg = 1
        else
          ! index of relevant subspace graph
          idx_sg = idxssg(lenprev,idspc,iyssg,iwssg,nelmax,nspc)
        end if
        ! offset of relevant subspace graph
        ioff = iy_info(1,idx_sg)
        ! max. number of electrons considered in that subspace graph
        nelmax_sg = iy_info(2,idx_sg)
        ! max. ms considered in that subspace graph
        msmax_sg = iy_info(3,idx_sg)

        ! length of current subspace string
        len = lensubspc(ispc,idspc(nn+1),nel-nn)
        iooff = mostnd(1,1,ispc)-1
        ! increment index using chosen subspace graph
        idx = idx4sg_sb(len,idx,nm,ig,
     &       idorb(nn+1),idspn(nn+1),idgam(nn+1),iy4sg(ioff),
     &       iooff,nelmax_sg,msmax_sg,ngam)
        nn = nn+len
        lenprev = len
        if (nn.eq.nel) exit
      end do

      idx_sg = idxssd(idspc,iyssg,nel,nspc)
      idx4sg = idx + ioff_dgm(idx_sg,ig,(nelmax-nm)/2+1)

      return
      end
