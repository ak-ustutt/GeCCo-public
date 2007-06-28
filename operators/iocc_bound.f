*----------------------------------------------------------------------*
      logical function iocc_bound(cbound,iocc,dagi,jocc,dagj)
*----------------------------------------------------------------------*
*     compare occupations element-wise with relation cbound
*     ('>','>=','<','<='). .TRUE. only if .TRUE. for all elements
*----------------------------------------------------------------------*
      
      implicit none
      include 'opdim.h'
      
      character, intent(in) ::
     &     cbound*(*)
      logical, intent(in) ::
     &     dagi, dagj
      integer, intent(in) ::
     &     iocc(ngastp,2), jocc(ngastp,2)
      
      integer ::
     &     ica, ica_i, ica_j, ihpv

      iocc_bound = .true.
      select case(trim(cbound))
      case('<','LT','lt')
        do ica = 1,2
          ica_i = ica
          ica_j = ica
          if (dagi) ica_i = 3-ica
          if (dagj) ica_j = 3-ica
          do ihpv = 1,ngastp
            iocc_bound = iocc_bound.and.
     &           iocc(ihpv,ica_i).lt.jocc(ihpv,ica_j)
          end do
        end do
      case('<=','LE','le')
        do ica = 1,2
          ica_i = ica
          ica_j = ica
          if (dagi) ica_i = 3-ica
          if (dagj) ica_j = 3-ica
          do ihpv = 1,ngastp
            iocc_bound = iocc_bound.and.
     &           iocc(ihpv,ica_i).le.jocc(ihpv,ica_j)
          end do
        end do
      case('>','GT','gt')
        do ica = 1,2
          ica_i = ica
          ica_j = ica
          if (dagi) ica_i = 3-ica
          if (dagj) ica_j = 3-ica
          do ihpv = 1,ngastp
            iocc_bound = iocc_bound.and.
     &           iocc(ihpv,ica_i).gt.jocc(ihpv,ica_j)
          end do
        end do
      case('>=','GE','ge')
        do ica = 1,2
          ica_i = ica
          ica_j = ica
          if (dagi) ica_i = 3-ica
          if (dagj) ica_j = 3-ica
          do ihpv = 1,ngastp
            iocc_bound = iocc_bound.and.
     &           iocc(ihpv,ica_i).ge.jocc(ihpv,ica_j)
          end do
        end do
      case default
        call quit(1,'iocc_bound','Unknow relation: "'//cbound//'"')
      end select
        
      return
      end
