*----------------------------------------------------------------------*
      subroutine set_strmap_kernel(strmap,
     &     iocc1,irestr1,idxms1,igam1,
     &     iocc2,irestr2,idxms2,igam2,
     &     g12y4sg,g12yinf,
     &     g12yssg,g12wssg,
     &     g_off_dgm, g_ndis,
     &     mostnd_cur,nsym,ngas_cur,igamorb)
*----------------------------------------------------------------------*
*     core routine for setting up the (str1,str2)->(str12) mapping
*
*     very first version
*     we could improve by running over batches of str1
*     and generating strings for graph 1 for this batch
*     outside loop over str2
*
*     andreas, may 2007
*      
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 000

      integer, intent(out) ::
     &     strmap(*)
      integer, intent(in) ::
     &     iocc1, irestr1(2,ngas_cur,2), idxms1, igam1,
     &     iocc2, irestr2(2,ngas_cur,2), idxms2, igam2,
     &     nsym, ngas_cur,
     &     mostnd_cur(2,nsym),igamorb(*),
     &     g12y4sg(*),g12yinf(*),
     &     g12yssg(*),g12wssg(*),
     &     g_off_dgm(*), g_ndis

      logical ::
     &     first1, first2
      integer ::
     &     iocc12, ms1, ms2, idxmap, isgn12
      integer ::
     &     idorb1(iocc1),idspn1(iocc1),idgam1(iocc1),idss1(iocc1),
     &     idorb2(iocc2),idspn2(iocc2),idgam2(iocc2),idss2(iocc2),
     &     idorb12(iocc1+iocc2),idspn12(iocc1+iocc2),
     &     idgam12(iocc1+iocc2),idss12(iocc1+iocc2),
     &     istr1, istr2, nstr1, nstr2, idx

      integer, external ::
     &     iordstr2, idx4sg
      logical, external ::
     &     next_string

      ms1 = iocc1-(idxms1-1)*2
      ms2 = iocc2-(idxms2-1)*2
      iocc12 = iocc1+iocc2

      idxmap = 0

      if (ntest.ge.1000) nstr2 = 0

      first2 = .true.
      ! loop over strings from graph 2
      do while(next_string(idorb2,idspn2,idss2,
     &     iocc2,ms2,igam2,first2,
     &     irestr2,
     &     mostnd_cur,igamorb,
     &     nsym,ngas_cur))

        if (ntest.ge.1000) nstr2 = nstr2+1
        if (ntest.ge.1000) nstr1 = 0

        first2 = .false.
        first1 = .true.
        ! loop over strings from graph 1
        do while(next_string(idorb1,idspn1,idss1,
     &       iocc1,ms1,igam1,first1,
     &       irestr1,
     &       mostnd_cur,igamorb,
     &       nsym,ngas_cur))
          first1 = .false.

          if (ntest.ge.1000) nstr1 = nstr1+1

          isgn12 = iordstr2(idorb12,idspn12,idss12,
     &         idorb1,idspn1,idss1,iocc1,
     &         idorb2,idspn2,idss2,iocc2)
c dbg
c          print *,'idorb1:  ',idorb1(1:iocc1)
c          print *,'idorb2:  ',idorb2(1:iocc2)
c          print *,'idorb12: ',isgn12,idorb12(1:iocc1+iocc2)
c dbg

          idxmap = idxmap+1
          strmap(idxmap) = isgn12
          if (isgn12.eq.0) cycle

          do idx = 1, iocc12
            idgam12(idx) = igamorb(idorb12(idx))
          end do
          
          strmap(idxmap) = isgn12*
     &        (idx4sg(iocc12,idss12,idorb12,idspn12,idgam12,
     &                g12y4sg,g12yinf,
     &                g12yssg,g12wssg,
     &                g_off_dgm,g_ndis,
     &                mostnd_cur,
     &                iocc12,nsym,ngas_cur)+1)

        end do
      end do

      if (ntest.ge.1000) then
        write(luout,*) 'set_strmap_kernel: the result of my work is'
        idxmap = 0
        do istr2 = 1, nstr2
          do istr1 = 1, nstr1
            idxmap = idxmap+1
            write(luout,'(x,2i5,"->",i5)') istr1,istr2,strmap(idxmap)
          end do
        end do
      end if

      return
      end
