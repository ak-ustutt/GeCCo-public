*----------------------------------------------------------------------*
      subroutine sort_string(sign,string,nidx)
*----------------------------------------------------------------------*
*
*     sort string by insertion sort according to the principles in the
*     function check_sequence(), see below
*
*     sign: updated according to the number of transpositions (-1)^ntransp
*     string: input/output string
*     nidx: length of string
*
*     andreas, dec 2020
*
*----------------------------------------------------------------------*

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'

      integer, parameter ::
     &     ntest = 00
      
      integer, intent(in) ::
     &     nidx
      integer, intent(inout) ::
     &     sign
      type(string_element), intent(inout) ::
     &     string(nidx)

      type(string_element) ::
     &     tmp
      integer ::
     &     ii, jj, ica, ihpvx, idx, jca, jhpvx, jdx,
     &     ntransp

      if (ntest.ge.100) then
        write(lulog,*) 'sort: info on entry'
        write(lulog,'(1x,a,15i4)') 'CA:   ',string(1:nidx)%ca
        write(lulog,'(1x,a,15i4)') 'HPVX: ',string(1:nidx)%hpvx
        write(lulog,'(1x,a,15i4)') 'IDX:  ',string(1:nidx)%idx
      end if
      
!     insertion sort of string elements
      do ii = 2, nidx
        ! take current index
        ica = string(ii)%ca
        ihpvx = string(ii)%hpvx
        idx = string(ii)%idx
        tmp = string(ii)
!     compare to all lower ones:
        ntransp = 0
        transp_loop: do jj = ii-1, 1, -1
          jca = string(jj)%ca
          jhpvx = string(jj)%hpvx
          jdx = string(jj)%idx
          if (check_sequence(ica,ihpvx,idx,jca,jhpvx,jdx)) then
            ntransp = ntransp+1
            string(jj+1) = string(jj)
          else
            exit transp_loop
          end if
        end do transp_loop
        write(lulog,*) ' sort: ',ii,ntransp
        string(ii-ntransp) = tmp
        if (mod(ntransp,2)>0) sign = -sign
      end do

      if (ntest.ge.100) then
        write(lulog,*) 'sort: info on exit'
        write(lulog,'(1x,a,15i4)') 'CA:   ',string(1:nidx)%ca
        write(lulog,'(1x,a,15i4)') 'HPVX: ',string(1:nidx)%hpvx
        write(lulog,'(1x,a,15i4)') 'IDX:  ',string(1:nidx)%idx
      end if


      return
      
      contains

      pure function check_sequence(ica,ihpvx,idx,jca,jhpvx,jdx)
     &                             result(is_ok)
!     compare two string entries and return, whether they are in the correct
!     sequence (if i... is before j...)
      integer, intent(in) :: ica, ihpvx, idx, jca, jhpvx, jdx

      logical :: is_ok
      
!     C before A
      if (ica.lt.jca) then
        is_ok = .true.
        return
      end if
      if (ica.gt.jca) then
        is_ok = .false.
        return
      end if
!     if ica == jca: check index sequence
!     ascending for C and descending for A
!     HPVX for C and HXVP for A (to heuristically support paired C A[PV] or C[PV] AH indices)
      if (ica==1) then
        if (ihpvx.lt.jhpvx) then
          is_ok = .true.
          return
        end if
        if (ihpvx.gt.jhpvx) then
          is_ok = .false.
          return
        end if
        ! if HPVX equal, the index decides
        is_ok = idx.le.jdx   ! "==" should in principle never happen
        return
      else
        ! we shift HPVX (=1,2,3,4) by MOD(<val>+2,4) -> (3,0,1,2)
        if (mod(ihpvx+2,4).gt.mod(jhpvx+2,4)) then
          is_ok = .true.
          return
        end if
        if (mod(ihpvx+2,4).lt.mod(jhpvx+2,4)) then
          is_ok = .false.
          return
        end if
! index decides finally (see above)
        is_ok = idx.ge.jdx
        return
      end if

      end function check_sequence

      
      end
