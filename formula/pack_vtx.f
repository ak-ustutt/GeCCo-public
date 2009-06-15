      integer(8) function pack_vtx(idx_op,iblk_op,dag)

      implicit none

      include 'stdunit.h'
      include 'opdim.h'

      integer, intent(in) ::
     &     idx_op, iblk_op
      logical, intent(in) ::
     &     dag

      integer ::
     &     iadj

      iadj = 0
      if (dag) iadj = 1
      pack_vtx = sign(iadj       *(pack_base**6)
     &               +abs(idx_op)*(pack_base**2)
     &               +iblk_op ,
     &                           idx_op)

      if (iblk_op.gt.pack_base*pack_base) then
        write(luout,*) 'pack_base    = ',pack_base*pack_base
        write(luout,*) 'iblk_op = ',iblk_op
        call quit(1,'pack_vtx',
     &       'incredibly large iblk_op encountered')
      end if

      if (idx_op.gt.pack_base**4) then
        write(luout,*) 'pack_base    = ',pack_base**4
        write(luout,*) 'idx_op = ',idx_op
        call quit(1,'pack_vtx',
     &       'large number of operators. increase pack_base')
      end if

      return
      end

