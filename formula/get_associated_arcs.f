*----------------------------------------------------------------------*
      subroutine get_associated_arcs(arc_list,len_list,iarc0,contr)
*----------------------------------------------------------------------*
*     given a spcific arc of a contraction (iarc0, contr) find all
*     arcs which link the same two super-vertices
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'

      integer, parameter ::
     &     ntest = 00

c dbg fix by mh
c original line      integer, intent(out) ::
c original line     &     len_list, arc_list(len_list)
c dbg original

      type(contraction), intent(in), target ::
     &     contr
      integer, intent(in) ::
     &     iarc0

c dbg resume fix
      integer, intent(out) ::
     &     len_list, arc_list(contr%narc)
c dbg end fix      
      integer ::
     &     iarc, ivtx1, ivtx2, isupvtx1, isupvtx2
      integer, pointer ::
     &     narc, joined(:,:), svertex(:)
      type(cntr_arc), pointer ::
     &     arc(:)

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'get_associated_arcs')
      end if

      narc => contr%narc

      arc => contr%arc
      svertex => contr%svertex
      joined => contr%joined
      len_list = 0

      ! super vertex? find further primitive vertices ...
      ! and further contractions between the supervertices
      ivtx1 = arc(iarc0)%link(1)
      ivtx2 = arc(iarc0)%link(2)      
c dbg
c      print *,'joined(0,ivtx1),joined(0,ivtx2)',
c     &         joined(0,ivtx1),joined(0,ivtx2)
c dbg
      if (joined(0,svertex(ivtx1)).gt.1.or.
     &    joined(0,svertex(ivtx2)).gt.1) then
c dbg
c        print *, 'supervertex detected'
c dbg
        isupvtx1 = svertex(ivtx1)
        isupvtx2 = svertex(ivtx2)
        len_list = 0
        do iarc = 1, narc
          ivtx1 = arc(iarc)%link(1)
          ivtx2 = arc(iarc)%link(2)
          if ((svertex(ivtx1).eq.isupvtx1.and.
     &         svertex(ivtx2).eq.isupvtx2)  .or.
     &        (svertex(ivtx1).eq.isupvtx2.and.
     &         svertex(ivtx2).eq.isupvtx1)) then
            len_list = len_list+1
            arc_list(len_list) = iarc
          end if
        end do
      else
        len_list = 1
        arc_list(1) = iarc0
      end if

      if (ntest.ge.100) then
        write(luout,*) 'len_list = ',len_list
        write(luout,*) 'arc_list = ',arc_list(1:len_list)
      end if

      return
      end
