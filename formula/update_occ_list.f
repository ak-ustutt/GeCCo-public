      subroutine update_occ_list(olist,vtx1,vtx2,occ_raw,nraw)

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_occ_list.h'

      integer, parameter ::
     &     ntest = 00
      character(len=15) ::
     &     i_am = 'update_occ_list'

      type(occ_list), intent(inout) ::
     &     olist
      integer, intent(in) ::
     &     vtx1, vtx2, nraw, occ_raw(ngastp,2,nraw)

      integer, pointer ::
     &     vtx_inf(:,:), occ(:,:,:), vtx_inf_new(:,:), occ_new(:,:,:)
      logical, pointer ::
     &     occ_is_new(:)
      integer ::
     &     n_vtx_inf, n_occ, idx, jdx, ioff, idxpair, n_occ_new,
     &     fac
      logical ::
     &     newpair, new

      logical, external ::
     &     list_cmp

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,i_am)
        write(luout,*) 'on entry:'
        call print_occ_list(luout,olist)
      end if

      n_vtx_inf  = olist%n_vtx_inf
      n_occ = olist%n_occ

      vtx_inf => olist%vtx_inf
      occ     => olist%occ

      ! check whether (vtx1,vtx2) is already on vtx_inf
      idxpair = 0
      do idx = 1, n_vtx_inf
        if (vtx1.eq.vtx_inf(1,idx).and.vtx2.eq.vtx_inf(2,idx)) then
          idxpair = idx
          exit
        end if
      end do

      newpair = idxpair.eq.0

      ! check whether some occupations already exist
      if (.not.newpair) then
        allocate (occ_is_new(nraw))
        n_occ_new = 0
        do idx = 1, nraw
          new = .true.
          do jdx = 1, n_occ
            new = new.and..not.list_cmp(occ    (:,:,jdx),
     &                                  occ_raw(:,:,idx),ngastp*2)
          end do
          if (new) n_occ_new = n_occ_new + 1
          occ_is_new(idx) = new
        end do
      else
        n_occ_new = nraw
      end if

      ! re-allocate, if necessary
      if (newpair.and.n_vtx_inf+1.gt.olist%max_vtx_inf) then
        fac = olist%max_vtx_inf / olist_vtx_inf_inc
        olist%max_vtx_inf = (fac+1)*olist_vtx_inf_inc
        allocate (vtx_inf_new(4,(fac+1)*olist_vtx_inf_inc))
        vtx_inf_new(1:4,1:n_vtx_inf) = vtx_inf(1:4,1:n_vtx_inf)
        deallocate(olist%vtx_inf)
        olist%vtx_inf => vtx_inf_new
        vtx_inf       => vtx_inf_new
      end if

      if (n_occ+n_occ_new.gt.olist%max_occ) then
        fac = olist%max_occ / olist_occ_inc
        allocate (occ_new(ngastp,2,(fac+1)*olist_occ_inc))
        olist%max_occ = (fac+1)*olist_occ_inc
        occ_new(1:ngastp,1:2,1:n_occ) =
     &         occ(1:ngastp,1:2,1:n_occ)
        deallocate(olist%occ)
        olist%occ => occ_new
        occ       => occ_new
      end if

      if (newpair) then

        vtx_inf(1,n_vtx_inf+1) = vtx1
        vtx_inf(2,n_vtx_inf+1) = vtx2
        vtx_inf(3,n_vtx_inf+1) = n_occ+1
        vtx_inf(4,n_vtx_inf+1) = n_occ+n_occ_new
        occ(1:ngastp,1:2,n_occ+1:n_occ+n_occ_new) =
     &       occ_raw(1:ngastp,1:2,1:n_occ_new)
        olist%n_vtx_inf  = olist%n_vtx_inf+1
        olist%n_occ = olist%n_occ+n_occ_new

      else

        ! shift occupations of pairs with higher index
        if (idxpair.lt.n_vtx_inf) then
          do idx = n_occ, vtx_inf(3,idxpair+1), -1
            occ(1:ngastp,1:2,idx+n_occ_new) =
     &           occ(1:ngastp,1:2,idx)
          end do
        end if

        ! insert new occupations
        jdx = 0
        ioff = vtx_inf(4,idxpair)
        do idx = 1, nraw
          if (.not.occ_is_new(idx)) cycle
          jdx = jdx+1
          occ(1:ngastp,1:2,ioff+jdx) = occ_raw(1:ngastp,1:2,idx)
        end do
        olist%n_occ = olist%n_occ+n_occ_new

        ! update vtx_inf
        vtx_inf(4,idxpair) = vtx_inf(4,idxpair) + n_occ_new
        do idx = idxpair+1, n_vtx_inf
          vtx_inf(3,idx) = vtx_inf(3,idx) + n_occ_new
          vtx_inf(4,idx) = vtx_inf(4,idx) + n_occ_new
        end do

      end if

      if (.not.newpair) deallocate (occ_is_new)

      if (ntest.ge.100) then
        write(luout,*) 'on exit:'
        call print_occ_list(luout,olist)
      end if

      return
      end
