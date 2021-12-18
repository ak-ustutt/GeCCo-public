      subroutine reorder_string_info(contr,
     &     from_vtx,to_vtx,cnt,occ_shift,ext)

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'

      integer, parameter ::
     &     ntest = 00
      
      type(contraction), intent(inout) ::
     &     contr
      integer, intent(in) ::
     &     from_vtx, to_vtx, cnt, occ_shift(ngastp,2)
      logical, intent(in) ::
     &     ext

      integer ::
     &     idx, idxst, idxnd, inc, ivtx, ica, ihpvx, jdx, nidx, sign
      integer ::
     &     occ_count(ngastp,2)
      type(string_element) ::
     &     tmp

      if (ntest.ge.100) then
        write(lulog,*) 'reorder: string on entry'
        call print_string(contr%contr_string,contr%nidx,contr%nvtx)
        call print_string_cnt(contr%contr_string,contr%nidx,contr%nvtx)
        call print_string_idx(contr%contr_string,contr%nidx,contr%nvtx)
        write(lulog,*) 'shift from, to: ',from_vtx,to_vtx
        call wrt_occ(lulog,occ_shift)
        write(lulog,*) 'external = ',ext
      end if
        
      occ_count = 0             ! counter for shifted vertices
      
      if (from_vtx.lt.to_vtx) then 
        idxst = contr%nidx      ! shift rightmost elements to the right
        idxnd = 1
        inc = -1
      else
        idxst = 1               ! shift leftmost elements to the left
        idxnd = contr%nidx
        inc = 1
      end if
c      write(lulog,'(1x,a,16i4)') 'vtx before: ',
c     &     contr%contr_string(1:contr%nidx)%vtx
      do idx = idxst, idxnd, inc
        if ((contr%contr_string(idx)%vtx==from_vtx).and.
     &       (contr%contr_string(idx)%ext.eqv.ext)) then
          if (.not.ext.and.cnt.ne.contr%contr_string(idx)%cnt) cycle
          ica =  contr%contr_string(idx)%ca
          ihpvx =  contr%contr_string(idx)%hpvx
          if (occ_count(ihpvx,ica).lt.occ_shift(ihpvx,ica)) then
c            write(lulog,*) ' > ',idx,ica,ihpvx
            occ_count(ihpvx,ica) = occ_count(ihpvx,ica)+1
            contr%contr_string(idx)%vtx = to_vtx ! first round: just rename vertex
          end if
        end if
      end do
c      write(lulog,'(1x,a,16i4)') 'vtx after:  ',
c     &     contr%contr_string(1:contr%nidx)%vtx
 
!     need to sort all entries such that vtx are ascending
      do idx = 2, contr%nidx
        tmp = contr%contr_string(idx)
        jdx = idx-1
        do while (contr%contr_string(jdx)%vtx.gt.tmp%vtx)
          contr%contr_string(jdx+1) = contr%contr_string(jdx)
          jdx = jdx-1
          if (jdx==0) exit
        end do
        contr%contr_string(jdx+1) = tmp
      end do

!     we sort all vertices to standard order, no implications for overall sign of diagram
      idx = 1
      sign = 0
      do ivtx = 1, contr%nvtx
        if (idx>nidx) exit
        nidx = 0
        jdx = idx
        do while ( contr%contr_string(idx)%vtx==ivtx )
          nidx = nidx+1
          idx = idx+1
          if (idx>contr%nidx) exit
        end do
c     dbg
c        write(lulog,*) 'ivtx,idx,jdx,nidx: ',ivtx,idx,jdx,nidx
c     dbg
        call sort_string(sign,contr%contr_string(jdx),nidx)
      end do

      if (ntest.ge.100) then
        write(lulog,*) 'reorder: string on exit'
        call print_string(contr%contr_string,contr%nidx,contr%nvtx)
        call print_string_cnt(contr%contr_string,contr%nidx,contr%nvtx)
        call print_string_idx(contr%contr_string,contr%nidx,contr%nvtx)
      end if
      
      end subroutine
