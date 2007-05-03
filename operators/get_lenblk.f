*----------------------------------------------------------------------*
      subroutine get_lenblk(lenc,lena,
     &     iocc,irst,msdist,igamdist,
     &     str_info,ihpvgas,ngas)
*----------------------------------------------------------------------*
*     get length of current block of operator
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'stdunit.h'
      include 'def_graph.h'
      include 'def_strinf.h'

      type(strinf), intent(in) ::
     &     str_info
      integer, intent(in) ::
     &     ngas, ihpvgas(ngas),
     &     iocc(ngastp,2), irst(2,ngas,2,2),
     &     msdist(ngastp,2), igamdist(ngastp,2)

      integer, intent(out) ::
     &     lenc, lena

      integer ::
     &     ica, igastp, igtyp, ngr4typ, igr4typ,
     &     idxgraph, idxms
      logical, external ::
     &     restr_cmp

      lenc = 1
      lena = 1
      do ica = 1, 2
        do igastp = 1, ngastp
          ! empty spaces have length 1
          if (iocc(igastp,ica).eq.0) cycle
          ! get graph type
          igtyp = ngastp*(iocc(igastp,ica)-1) + igastp
          ! number of graphs with same type
          ngr4typ = str_info%gtab(1,igtyp)
          ! check restrictions
          gr4typ: do igr4typ = 1, ngr4typ
            ! actual index of graph
            idxgraph = str_info%gtab(1+igr4typ,igtyp)
            if (restr_cmp(irst,str_info%igas_restr(1,1,1,idxgraph),
     &                       ica,igastp,ihpvgas,ngas)) exit gr4typ
            idxgraph = -idxgraph ! indicate that this was not what
                                 ! we wanted
          end do gr4typ
          if (idxgraph.le.0) then
            write(luout,*) 'ERROR: string not in list'
            write(luout,*) 'Operator was'
            call wrt_occ(luout,iocc)
            call wrt_rstr(luout,irst,ngas)
            write(luout,*) 'C/A, GAS-TYP: ',ica,igastp
            call quit(1,'get_lenblk','string not in list')
          end if

          idxms = (iocc(igastp,ica)-msdist(igastp,ica))/2+1
          !
          if (ica.eq.1)
     &         lenc = lenc*str_info%g(idxgraph)%
     &         lenstr_gm(igamdist(igastp,ica),idxms)
          if (ica.eq.2)
     &         lena = lena*str_info%g(idxgraph)%
     &         lenstr_gm(igamdist(igastp,ica),idxms)

        end do
      end do
      return
      end
