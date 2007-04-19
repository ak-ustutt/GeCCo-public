      logical function next_ssd(idss,nel,nelmax,nspc,mnmxspc)

      implicit none
      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nelmax, nspc,
     &     mnmxspc(2,*)

      integer, intent(inout) ::
     &     idss(1:nelmax), nel

      logical ::
     &     found
      integer ::
     &     idmin(nelmax), idmax(nelmax)

      logical, external ::
     &     next_rvlex, allow_sbsp_dis


      if (ntest.gt.0) then
        write(luout,*) '================'
        write(luout,*) 'Entered next_ssd'
        write(luout,*) '================'
        write(luout,*) ' number of subspaces to be considered: ',nspc
      end if

      ! quick return if no subspace distributions are possible
      if (nspc.le.0) then
        next_ssd = .false.
        return
      end if

      if (ntest.ge.100) then
        write(luout,*) ' initial subspace distribution:'
        write(luout,*) '  nel = ',nel
        write(luout,*) '  dis = (',idss(1:nel),')'
      end if

      ! 1) try to get next distribution with same number of electrons
      found = .false.
      if (nel.gt.0) then
        idmin(1:nel) = 1
        idmax(1:nel) = nspc
        rvlex_loop: do
          ! get next distribution, if any
          if (.not.next_rvlex(nel,idss,idmin,idmax)) exit rvlex_loop
          ! check against constraints on mnmxspc
          if (allow_sbsp_dis(idss,nel,nspc,mnmxspc)) then
            found = .true.
            exit rvlex_loop
          end if
        end do rvlex_loop
      end if

      ! 2) increase number of electrons and return lexically lowest distr.
      if (.not.found.and.nel.lt.nelmax) then
        nel = nel+1
        idss(1:nel) = 1  ! lowest possible string
        idmin(1:nel) = 1
        idmax(1:nel) = nspc
        rvlex_loop2: do
          ! check against constraints on mnmxspc
          if (allow_sbsp_dis(idss,nel,nspc,mnmxspc)) then
            found = .true.
            exit rvlex_loop2
          end if
          ! get next distribution, if any
          if (.not.next_rvlex(nel,idss,idmin,idmax)) exit rvlex_loop2
        end do rvlex_loop2

      end if

      if (ntest.ge.100) then
        if (found) then
          write(luout,*) ' final subspace distribution:'
          write(luout,*) '  nel = ',nel
          write(luout,*) '  dis = (',idss(1:nel),')'
        else
          write(luout,*) ' no further subspace distribution found' 
        end if
      end if

      next_ssd = found

      return
      end
