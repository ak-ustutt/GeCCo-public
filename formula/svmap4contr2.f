*----------------------------------------------------------------------*
      subroutine svmap4contr2(svmap,contr,unique)
*----------------------------------------------------------------------*
*     assign each vertex in contr a number which tells us to which
*     super-vertex node this vertex contributes an external line
*     (0 if completely contracted).
*     version for contr with xarc info -> much easier
*     unique is true unless some vertex contributes to more than one
*     super-vertex node
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 00

      type(contraction), intent(in), target ::
     &     contr
      integer, intent(out) ::
     &     svmap(contr%nvtx)
      logical, intent(out) ::
     &     unique

      integer ::
     &     nvtx, nxarc, ivtx, isvtx, ixarc
      type(cntr_arc), pointer ::
     &     xarc(:)
      
      if (ntest.ge.100)
     &     call write_title(luout,wst_dbg_subr,'svmap4contr2')

      nvtx = contr%nvtx
      nxarc = contr%nxarc
      xarc => contr%xarc

      unique = .true.
      svmap(1:nvtx) = 0
      do ivtx = 1, nvtx
        do ixarc = 1, nxarc
          if (xarc(ixarc)%link(1).eq.ivtx) then
            isvtx = xarc(ixarc)%link(2)
            if (svmap(ivtx).gt.0.and.isvtx.ne.svmap(ivtx)) then
              unique = .false.
              exit
cmh              call prt_contr3(luout,contr,-1)
cmh              call quit(1,'svmap4contr2',
cmh     &                       'inconsistent xarc detected!')
            end if
            svmap(ivtx) = isvtx
          end if
        end do
      end do

      if (ntest.ge.100)
     &     write(luout,'(x,a,10i5)') 'svmap: ',svmap(1:nvtx)

      return
      end
