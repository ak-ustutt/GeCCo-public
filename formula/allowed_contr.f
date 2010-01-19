*----------------------------------------------------------------------*
      logical function allowed_contr(contr,arclist,lenlist)
*----------------------------------------------------------------------*
*     given a contraction and a list of arcs that have to be
*     contracted over in a single step, determines if that step
*     is currently possible: no overlap between two involved arcs
*     at the same vertex.
*
*     matthias, jan. 2010
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'ifc_operators.h'
      include 'def_contraction.h'

      integer, intent(in) ::
     &     lenlist, arclist(lenlist)
      type(contraction), intent(in) ::
     &     contr
      integer ::
     &     icnt(ngastp,2), ivtx, ilist, icnt_arc(ngastp,2)
      logical ::
     &     dag

      allowed_contr = .true.
c dbg
c      print *,'output from allowed_contr:'
c dbg

      ! for each vertex ...
      do ivtx = 1, contr%nvtx
c dbg
c      print *,'ivtx: ',ivtx
c dbg
        icnt = 0
        ! ... loop over arcs involved in contraction ...
        do ilist = 1, lenlist
          if (contr%arc(arclist(ilist))%link(1).eq.ivtx) then
            dag = .false.
          else if (contr%arc(arclist(ilist))%link(2).eq.ivtx) then
            dag = .true.
          else
            cycle
          end if
c dbg
c      print *,'iarc, dag:',arclist(ilist), dag
c dbg
          ! ... and find out if there are overlapping arcs
          icnt_arc = contr%arc(arclist(ilist))%occ_cnt
          if (iocc_zero(iocc_overlap(icnt,.false.,
     &        icnt_arc,dag))) then
            icnt = iocc_add(1,icnt,.false.,1,icnt_arc,dag)
c dbg
c            print *,'ok'
c dbg
          else
            allowed_contr = .false.
c dbg
c            print *,'overlap!'
c dbg
            return
          end if
        end do
      end do

      return
      end
