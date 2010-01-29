*----------------------------------------------------------------------*
      logical function opblk_restr_cmp(op1,iblk1,dag1,op2,iblk2,dag2)
*----------------------------------------------------------------------*
*
*     compare restrictions on operator blocks 
*
*     .true. if equal 
*
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'def_operator.h'

      type(operator), intent(in) ::
     &     op1, op2
      integer, intent(in) ::
     &     iblk1, iblk2
      logical, intent(in) ::
     &     dag1, dag2

      logical ::
     &     same, reverse
      integer ::
     &     nspin, ngas, nj, idx1, idx2,
     &     ij, igas, ica, ica2, imsk, ispin

      integer, pointer ::
     &     restr1(:,:,:,:,:,:),
     &     restr2(:,:,:,:,:,:)

      opblk_restr_cmp = .false.

      reverse = dag1.xor.dag2

      nspin = op1%nspin
      ngas  = op1%ngas
      nj    = op1%njoined

      if (nspin.ne.op2%nspin.or.ngas.ne.op2%ngas) return   ! or error?

      if (nj.ne.op2%njoined) return

      restr1 => op1%igasca_restr
      restr2 => op2%igasca_restr

      ! start within restr-array
      idx1 = (iblk1-1)*nj
      idx2 = (iblk2-1)*nj
      if (reverse) idx2 = iblk2*nj+1

      same = .true.
      cmp_loop: do ij = 1, nj
        idx1 = idx1+1
        if (.not.reverse) then
          idx2 = idx2+1
        else
          idx2 = idx2-1
        end if
        do ispin = 1, nspin
          do imsk = 1, 2
            do ica = 1, 2
              ica2 = ica
              if (reverse) ica2 = 3-ica
              do igas = 1, ngas
                same = same.and.
     &             (restr1(1,igas,ica ,imsk,ispin,idx1).eq.
     &              restr2(1,igas,ica2,imsk,ispin,idx2)).and.
     &             (restr1(2,igas,ica ,imsk,ispin,idx1).eq.
     &              restr2(2,igas,ica2,imsk,ispin,idx2))
              end do
            end do
            if (.not.same) exit cmp_loop
          end do
        end do
      end do cmp_loop

      opblk_restr_cmp = same

      return
      end

