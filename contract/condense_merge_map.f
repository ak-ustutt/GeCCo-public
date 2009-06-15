*----------------------------------------------------------------------*
      subroutine condense_merge_map(mmap_out,
     &     mmap_in,ld_mmap,nvtx_tgt,reverse)
*----------------------------------------------------------------------*
*     on input: merge map in long form as
*
*          map(1:n(i),1:2,target_vertex) = vertices(i)
*
*     on output: merge map in short form as
*
*          map(:) = 
*            (...,n(i),vertex(1),vertex(2),...,vertex(n(i)),n(i+1),..)
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     ld_mmap, nvtx_tgt,
     &     mmap_in(ld_mmap,2,nvtx_tgt)
      integer, intent(out) ::
     &     mmap_out(*)
      logical, intent(in) ::
     &     reverse

      integer ::
     &     ivtx_tgt, ii, jj, ivtx, idx, iist, iind, iinc

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'condense_merge_map')
        write(luout,*) 'nvtx_tgt, ld_mmap: ',nvtx_tgt, ld_mmap
        if (nvtx_tgt.gt.1000 .or. ld_mmap.gt.1000) stop 'hilfe'
        write(luout,*) 'input map'
        do ivtx_tgt = 1, nvtx_tgt
          write(luout,'(x,i3,"<- op1 vtx ",10i4)')
     &         ivtx_tgt, mmap_in(1:ld_mmap,1,ivtx_tgt)
          write(luout,'(6x," op2 vtx ",10i4)')
     &                   mmap_in(1:ld_mmap,2,ivtx_tgt)
        end do
        if (reverse) write(luout,*) 'reverting Op1 and Op2'
      end if

      if (.not.reverse) then
        iist  = 1
        iind = 2
        iinc  = 1
      else
        iist  = 2
        iind = 1
        iinc  = -1
      end if

      idx = 1
      do ivtx_tgt = 1, nvtx_tgt
        do ii = iist, iind, iinc
          mmap_out(idx) = 0
          do jj = 1, ld_mmap
            ivtx = mmap_in(jj,ii,ivtx_tgt)
            if (ivtx.eq.0) exit
            mmap_out(idx) = mmap_out(idx) + 1
            mmap_out(idx+jj) = ivtx
          end do
          idx = idx+mmap_out(idx)+1
        end do
      end do

      if (ntest.ge.100) then
        write(luout,*) 'output map:'
        idx = 1
        do ivtx_tgt = 1, nvtx_tgt
          ivtx = mmap_out(idx)
          if (ivtx.gt.0) then
            write(luout,'(3x,i3,"<-",i3," from op1: ",10i3)')
     &                        ivtx_tgt, ivtx,
     &                        mmap_out(idx+1:idx+ivtx)
          else
            write(luout,'(3x,i3,"<-  0 from op1")') ivtx_tgt
          end if
          idx = idx + ivtx + 1
          ivtx = mmap_out(idx)
          if (ivtx.gt.0) then
            write(luout,'(3x,i3,"<-",i3," from op2: ",10i3)')
     &                        ivtx_tgt, ivtx,
     &                        mmap_out(idx+1:idx+ivtx)
          else
            write(luout,'(3x,i3,"<-  0 from op2")') ivtx_tgt
          end if
          idx = idx + ivtx + 1
        end do
c dbg
c          print *,'raw map looks like: '
c          print *,mmap_out(1:idx-1)
c dbg
        
      end if

      return
      end
