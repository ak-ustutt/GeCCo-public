*----------------------------------------------------------------------*
      subroutine tex_op(str,name,nocc,occ,idxset)
*----------------------------------------------------------------------*
*     print operator in TeX style on string str
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'par_opnames_gen.h'

      character, parameter ::
     &     idxchar(6*4) = (/'J','I','L','K','N','M',
     &                      'A','B','C','D','E','F',
     &                      'U','V','W','X','Y','Z',
     &                      'i','j','k','l','m','n'/)

      character, intent(out) ::
     &     str*(*)
      character, intent(in) ::
     &     name*(*)
      integer, intent(in) ::
     &     nocc, occ(ngastp,2,nocc), idxset(nocc)

      integer ::
     &     ica, ipos, idx, hpvx, nsub, nidx, nprm, iprm

      ! a preliminary translation table:
      if (len_trim(name).gt.1) then
        select case(trim(name))
c        case('CBAR') 
c          write(str,'("\bar{C}")')
        case(op_tbar) 
          write(str,'("\bar{T}")')
        case(op_dia) 
          write(str,'("D")')
        case(op_hhat) 
          write(str,'("\hat{H}")')
        case(op_omg) 
          write(str,'("\Omega")')
        case(op_ccen) 
          write(str,'("E_{\text{CC}}")')
        case(op_mpr12en) 
          write(str,'("E(\text{MP-R12})")')
        case(op_mpr12lg)
          write(str,'("L(\text{MP2})")')
        case(op_c12)
          write(str,'("C2")')
        case(op_cba)
          write(str,'("\bar{C2}")')
        case(op_cex)
          write(str,'("C1")')
        case(op_cexbar)
          write(str,'("\bar{C1}")')
        case(op_r12)
          write(str,'("R")')
        case(op_omgcex)
          write(str,'("\Omega\,(C1)")')
        case default
          call quit(1,'tex_op','adapt for operator "'//trim(name)//'"')
        end select
      else
        write(str,'(a)') trim(name)
      end if
      do ica = 1, 2
        ipos = len_trim(str)+1
        if (ica.eq.1) write(str(ipos:),'("_{")')
        if (ica.eq.2) write(str(ipos:),'("^{")')
        do idx = 1, nocc
          do hpvx = 1, ngastp
            nsub = occ(hpvx,ica,idx)
            if (nsub.gt.0) then
              if (idxset(idx).gt.0) then
                nidx = mod((idxset(idx)-1)*2+ica-1,6)+1
                nprm = ((idxset(idx)-1)*2+ica-1)/6
              else
                nidx = mod((-idxset(idx)-1)*2+3-ica-1,6)+1
                nprm = ((-idxset(idx)-1)*2+3-ica-1)/6
              end if
              ipos = len_trim(str)+1
              write(str(ipos:),'(a)')
     &             idxchar(nidx+(hpvx-1)*6)
              if (nprm.gt.0)
     &             write(str(ipos+1:),'("^{")')
              ipos = len_trim(str)+1
              do iprm = 1, nprm
                write(str(ipos:),'("\prime")')
                ipos = ipos+6
              end do
              if (nprm.gt.0)
     &             write(str(ipos:),'("}")')
              ipos = len_trim(str)+1
              write(str(ipos:),'("_{",i3,"}")') nsub
            end if
          end do
        end do
        ipos = len_trim(str)+1
        write(str(ipos:),'("}")')        
      end do

      return
      end
