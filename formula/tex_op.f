*----------------------------------------------------------------------*
      subroutine tex_op(str,name,dagger,nset,typset,occset,idx0set,nj)
*----------------------------------------------------------------------*
*     print operator in TeX style on string str
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'par_opnames_gen.h'

      character, parameter ::
     &     idxchar(4*2*3) = (/'i','k','l',
     &                        'a','c','d',
     &                        'u','w','x',
     &                        'P','R','S',
     &                        'j','l','k',
     &                        'b','d','c',
     &                        'v','x','w',
     &                        'Q','S','R'/)

      character, intent(out) ::
     &     str*(*)
      character, intent(in) ::
     &     name*(*)
      logical, intent(in) ::
     &     dagger
      integer, intent(in) ::
     &     nset, nj, typset(nset), occset(ngastp,2,nj,nset),
     &                             idx0set(ngastp,2,nj,nset)

      character ::
     &     charidx
      integer ::
     &     ij, ica, iset, ipos, idx, iidx, hpvx, nsub, nidx,
     &     idx0

      ipos = 1
      if (dagger) then
        write(str,'("[")')
        ipos = ipos+1
      end if

      write(str(ipos:),'("\op",a)') trim(name)

      ! a preliminary translation table:
c      if (len_trim(name).gt.1) then
c        select case(trim(name))
c        case(op_tbar) 
c          write(str(ipos:),'("\bar{T}")')
c        case(op_dia)
c          write(str(ipos:),'("D")')
c        case(op_hhat)
c          write(str(ipos:),'("\hat{H}")')
c        case(op_omg)
c          write(str(ipos:),'("\Omega")')
c        case(op_ccen)
c          write(str(ipos:),'("E_{\text{CC}}")')
c        case(op_mpr12en) 
c          write(str(ipos:),'("E(\text{MP-R12})")')
c        case(op_mpr12lg)
c          write(str(ipos:),'("L(\text{MP2})")')
c        case(op_c12)
c          write(str(ipos:),'("C2")')
c        case(op_cba)
c          write(str(ipos:),'("\bar{C2}")')
c        case(op_cex)
c          write(str(ipos:),'("C1")')
c        case(op_cexbar)
c          write(str(ipos:),'("\bar{C1}")')
c        case(op_r12,op_rint)
c          write(str(ipos:),'("R")')
c        case(op_omgcex)
c          write(str(ipos:),'("\Omega\,(C1)")')
c        case(op_bh_inter)
c          write(str(ipos:),'("{B^h}")')
c        case(op_c_inter)
c          write(str(ipos:),'("C")')
c        case(op_p_inter)
c          write(str(ipos:),'("P")')
c        case(op_z_inter)
c          write(str(ipos:),'("Z")')
c        case default
c          write(str(ipos:),'(a)') trim(name)
cc          call quit(1,'tex_op','adapt for operator "'//trim(name)//'"')
c        end select
c      else
c        write(str(ipos:),'(a)') trim(name)
c      end if

      ipos = len_trim(str)+1
      if (dagger) then
        write(str(ipos:),'("^\dagger]")')
      end if

      do ica = 1, 2
        ipos = len_trim(str)+1
        if (ica.eq.1) write(str(ipos:),'("_{")')
        if (ica.eq.2) write(str(ipos:),'("^{")')
        do ij = 1, nj
c          do iset = 1, nset
c            do hpvx = 1, ngastp
          do hpvx = 1, ngastp
            do iset = 1, nset
              nidx = occset(hpvx,ica,ij,iset)
              if (nidx.eq.0) cycle
              idx0 = idx0set(hpvx,ica,ij,iset)
              charidx = idxchar((ica-1)*ngastp*3 +
     &                          (hpvx-1)*3 + typset(iset))

              do iidx = idx0+1, idx0+nidx
                ipos = len_trim(str)+1
                write(str(ipos:),'(a,"_{",i2,"}")') charidx,iidx
              end do
            end do
          end do
        end do
        ipos = len_trim(str)+1
        write(str(ipos:),'("}")')
      end do

      return
      end
