*----------------------------------------------------------------------*
      logical function check_occ_partition(p_tgt,n_tgt,p_chk,n_chk)
*----------------------------------------------------------------------*
*     check that p_tgt = ([OT1] [OT2] ... [OTn_tgt]) or a sub-partion
*     thereof is contained in p_chk = ([OC1] [OC2] ... [OCn_chk])
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'ifc_operators.h'
      include 'stdunit.h'

      integer, intent(in) ::
     &     n_tgt, n_chk,
     &     p_tgt(ngastp,2,n_tgt), p_chk(ngastp,2,n_chk)

      logical ::
     &     init, init2, found, assigned(n_chk)
      integer ::
     &     maxpart, npart, npart1, npart2, idx,
     &     ichk
      integer ::
     &     zocc(ngastp,2)
      integer, pointer ::
     &     occ(:,:,:)

      logical, external ::
     &     next_part_connection

c dbg
c      print *,'in check_occ_partition: ',n_chk, n_tgt
c dbg
      ! quick answer
      check_occ_partition = .false.
      if (n_chk.lt.n_tgt) return

      ! zero occupation
      zocc = 0
c dbg
c      print *,'p_tgt = '
c      call wrt_occ_n(6,p_tgt,n_tgt)
c      print *,'p_chk = '
c      call wrt_occ_n(6,p_chk,n_chk)
c dbg

      if (n_tgt.eq.1) then
c dbg
c        print *,'N==1 part entered'
c dbg
        
        ! maximum partition grade
        maxpart = min(sum(p_tgt(1:ngastp,1:2,1)),n_chk)
        allocate(occ(ngastp,2,maxpart))
        ! loop over partition grades
        pgrade: do npart = 1, maxpart
          ! generate the partitions of p_tgt
          init = .true.
          do while(next_part_connection(occ,
     &         init,p_tgt,npart,
     &         zocc,p_tgt))
            init = .false.
            ! array for assigning p_chk
            assigned(1:n_chk) = .false.
            ! try to match each element of occ with those
            do idx = 1, npart 
              found = .false.
              do ichk = 1, n_chk
                if (assigned(ichk)) cycle
                if (iocc_equal(occ(1:ngastp,1:2,idx),.false.,
     &                       p_chk(1:ngastp,1:2,ichk),.false.) ) then
                  found = .true.
                  assigned(ichk) = .true.
                  exit
                end if
              end do
              ! if nothing was found, no need to look for the 
              ! other elements
              if (.not.found) exit
            end do
            ! if the partition matches the p_chk, we are done
            if (found) exit pgrade
          end do

        end do pgrade

        check_occ_partition = found

        deallocate(occ)

      else if (n_tgt.eq.2) then
c dbg
c        print *,'N==2 part entered'
c dbg

c        call quit(1,'check_occ_partition','case n_tgt.eq.2 not tested')

        ! maximum partition grade
        maxpart = min(sum(p_tgt(1:ngastp,1:2,1:2)),n_chk)
        allocate(occ(ngastp,2,maxpart))
        ! loop over partition grades
        pgrade2: do npart = 1, maxpart
          ! distribute the partition grade between the two p_tgt's
          do npart1 = 1, maxpart-1
            npart2 = maxpart-npart1
            ! generate the partitions of p_tgt
            init = .true.
            do while(next_part_connection(occ,
     &         init,p_tgt,npart1,
     &         zocc,p_tgt))
              init = .false.
              init2 = .true.
              do while(next_part_connection(occ(1,1,npart1+1),
     &             init2,p_tgt(1,1,2),npart2,
     &             zocc,p_tgt(1,1,2)))
                init2 = .false.

                ! array for assigning p_chk
                assigned(1:n_chk) = .false.
                ! try to match each element of occ with those
                do idx = 1, npart 
                  found = .false.
                  do ichk = 1, n_chk
                    if (assigned(ichk)) cycle
                    if (iocc_equal(occ(1:ngastp,1:2,idx),.false.,
     &                       p_chk(1:ngastp,1:2,ichk),.false.) ) then
                      found = .true.
                      assigned(ichk) = .true.
                      exit
                    end if
                  end do
                  ! if nothing was found, no need to look for the 
                  ! other elements
                  if (.not.found) exit
                end do
                ! if the partition matches the p_chk, we are done
                if (found) exit pgrade2
              end do
            end do

          end do
            
        end do pgrade2

        check_occ_partition = found

        deallocate(occ)

      else if (n_tgt.gt.2) then
        call quit(1,'check_occ_partition',
     &     'only implemented for max. two targets')
      end if
      
c dbg
c      print *,'at the end: ',check_occ_partition
c dbg

      return
      end
