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

      ! a preliminary translation table:
      if (len_trim(name).gt.1) then
        select case(trim(name))
        case(op_tbar) 
          write(str(ipos:),'("\bar{T}")')
        case(op_dia)
          write(str(ipos:),'("D")')
        case(op_hhat)
          write(str(ipos:),'("\hat{H}")')
        case(op_omg)
          write(str(ipos:),'("\Omega")')
        case(op_ccen)
          write(str(ipos:),'("E_{\text{CC}}")')
        case(op_mpr12en) 
          write(str(ipos:),'("E(\text{MP-R12})")')
        case(op_mpr12lg)
          write(str(ipos:),'("L(\text{MP2})")')
        case(op_c12)
          write(str(ipos:),'("C2")')
        case(op_cba)
          write(str(ipos:),'("\bar{C2}")')
        case(op_cex)
          write(str(ipos:),'("C1")')
        case(op_cexbar)
          write(str(ipos:),'("\bar{C1}")')
        case(op_r12,op_rint)
          write(str(ipos:),'("R")')
        case(op_omgcex)
          write(str(ipos:),'("\Omega\,(C1)")')
        case(op_bh_inter)
          write(str(ipos:),'("{B^h}")')
        case(op_c_inter)
          write(str(ipos:),'("C")')
        case(op_p_inter)
          write(str(ipos:),'("P")')
        case(op_z_inter)
          write(str(ipos:),'("Z")')
        case default
          write(str(ipos:),'(a)') trim(name)
c          call quit(1,'tex_op','adapt for operator "'//trim(name)//'"')
        end select
      else
        write(str(ipos:),'(a)') trim(name)
      end if

      ipos = len_trim(str)+1
      if (dagger) then
        write(str(ipos:),'("^\dagger]")')
      end if

      do ica = 1, 2
        ipos = len_trim(str)+1
        if (ica.eq.1) write(str(ipos:),'("_{")')
        if (ica.eq.2) write(str(ipos:),'("^{")')
        do ij = 1, nj
          do iset = 1, nset
            do hpvx = 1, ngastp
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
