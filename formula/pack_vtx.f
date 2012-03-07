      integer(8) function pack_vtx(idx_op,iblk_op,dag)

      implicit none

      include 'stdunit.h'
      include 'opdim.h'

      integer, intent(in) ::
     &     idx_op, iblk_op
      logical, intent(in) ::
     &     dag

      integer(8) ::
     &     iadj, pb8, idx8

      pb8 = pack_base
      idx8 = idx_op
      iadj = 0
      if (dag) iadj = 1
      pack_vtx = sign(iadj     *(pb8**6)
     &               +abs(idx8)*(pb8**2)
     &               +iblk_op ,
     &                           idx8)

      if (iblk_op.gt.pb8*pb8) then
        write(luout,*) 'pack_base    = ',pb8*pb8
        write(luout,*) 'iblk_op = ',iblk_op
        call quit(1,'pack_vtx',
     &       'incredibly large iblk_op encountered')
      end if

      if (idx_op.gt.pb8**4) then
        write(luout,*) 'pack_base    = ',pb8**4
        write(luout,*) 'idx_op = ',idx_op
        call quit(1,'pack_vtx',
     &       'large number of operators. increase pack_base')
      end if

      return
      end

