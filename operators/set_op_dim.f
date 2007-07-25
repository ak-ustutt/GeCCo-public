*----------------------------------------------------------------------*
      subroutine set_op_dim(ipass,reset_on_rank,op,str_info,ngam)
*----------------------------------------------------------------------*
*
*     set up the dimension arrays for operator op
*
*     the operator has a total Ms and IRREP, as given in "op".
*     it is stored as
*
*          Op(Cp,Ch,Cv,Ap,Ah,Av)
*
*     where C are indices associated with creation strings and
*           A are indices associated with annihilation strings
*     in (p)article, (h)ole, and (v)alence spaces
*     
*     in this routine we define the sequence in which the different blocks
*     of the operator are to be stored; with "different blocks" we mean
*     all the different cases of Ms and symmetry labels distributed among
*     the strings (Cp,Ch,Cv,Ap,Ah,Av) such that they add up to the 
*     required total Ms and IRREP
*
*     ipass = 1   set len_op_gmo, len_op_occ, 
*                     off_op_gmo, off_op_occ, len_op
*                 and get dimension info for off_op_gmox
*
*     ipass = 2   set off_op_gmox
*
*     the len_op_gmox is only needed for operators, where more than
*     one class of spaces (I mean H/P/V) is occupied for C or for A
*     usual H->P excitation operators (e.g. T-operators of SR-CC) do
*     not fall into this group and ipass=2 may be skipped
*
*     reset_on_rank: reset counter when going to new operator rank
*           needed for Hamiltonian, as we treat 1- and 2-particle
*           operators separately (obsolete for gecco-version)
*
*     andreas, end of 2006
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_operator.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'multd2h.h'

      integer, parameter ::
     &     ntest = 100
      
      integer, intent(in) ::
     &     ipass, ngam
      logical, intent(in) ::
     &     reset_on_rank
      type(strinf), intent(in) ::
     &     str_info
      type(operator), intent(inout) ::
     &     op

      logical ::
     &     first
      integer ::
     &     idxstr, idxstr_tot, idxdis, iocc_cls, nexc,
     &     msmax, msmax_last,
     &     msa, msc, msah, msap, msav, msch, mscp, mscv,
     &     idxmsa, idxms, msahmax, msapmax, msavmax,
     &     mschmax, mscpmax, mscvmax,
     &     igch, igcp, igcv, igah, igap, igav,
     &     lench, lencp, lencv, lenah, lenap, lenav,
     &     igamah, igamap, igamav, igamch, igamcp, igamcv, igama, igamc,
     &     did, iexc, igam
      integer ::
     &     msd(ngastp,2), igamd(ngastp,2)

      logical, external ::
     &     next_msgamdist
      integer, external ::
     &     msgmdid
      
      if (ntest.gt.5) then
        write(luout,*) '============'
        write(luout,*) ' set_op_dim'
        write(luout,*) '============'
        write(luout,*) ' ipass = ',ipass
        write(luout,*) ' operator ',trim(op%name)
      end if

      idxstr = 0
      idxstr_tot = 0
      msmax_last = 0
      occ_cls: do iocc_cls = 1, op%n_occ_cls

        if(op%formal_blk(iocc_cls))cycle

        if (ntest.ge.150) then
          write(luout,*) 'class: ',iocc_cls
          call wrt_occ(luout,op%ihpvca_occ(1,1,iocc_cls))
        end if

        msmax = min(op%ica_occ(1,iocc_cls),op%ica_occ(2,iocc_cls))

        ! we take msmax as an indication that the operator rank
        ! has changed.
        ! if requested, count offset arrays from new:
        if (reset_on_rank.and.msmax_last.ne.msmax) then
          idxstr = 0 
        end if
        msmax_last = msmax

        op%off_op_occ(iocc_cls) = idxstr

        if (ipass.eq.1) then
          op%off_op_gmox(iocc_cls)%maxd = 0
        else
          op%off_op_gmox(iocc_cls)%
     &       d_gam_ms(1:op%off_op_gmox(iocc_cls)%maxd,1:ngam,1:msmax)=-1
          op%off_op_gmox(iocc_cls)%
     &       did(1:op%off_op_gmox(iocc_cls)%maxd,1:ngam,1:msmax) = 0
          op%off_op_gmox(iocc_cls)%
     &       ndis(1:ngam,1:msmax) = 0
        end if

        idxmsa = 0
        msa_loop: do msa = msmax, -msmax, -2

          idxmsa = idxmsa+1

          ! C <-> A  means alpha <-> beta !!
          msc = msa + op%mst
          ! to make this point clear:
          ! usually, we have mst=0 operators, which implies
          ! Ms(C) == Ms(A)
          igama_loop: do igama = 1, ngam
            igamc = multd2h(igama,op%gamt)
            
            ! store the current position in offset array
            op%off_op_gmo(iocc_cls)%gam_ms(igama,idxmsa) = idxstr

            ! now, to be general, we have to loop over all
            ! possible MS and IRREP distributions over H/P/V spaces
            ! for both C and A strings

            idxdis = 0

            ! retrieve info on graph (for legibility of code)
            igch = op%idx_graph(ihole,1,iocc_cls)
            igcp = op%idx_graph(ipart,1,iocc_cls)            
            igcv = 0
            if (ivale.gt.0) igcv = op%idx_graph(ivale,1,iocc_cls)
            igah = op%idx_graph(ihole,2,iocc_cls)
            igap = op%idx_graph(ipart,2,iocc_cls)
            igav = 0
            if (ivale.gt.0) igav = op%idx_graph(ivale,2,iocc_cls)
            
            ! MS(A)/IRREP(A)/MS(C)/IRREP(C) distributions:
            msahmax = op%ihpvca_occ(ihole,2,iocc_cls)
            msapmax = op%ihpvca_occ(ipart,2,iocc_cls)
            msavmax = 0
            if (ivale.gt.0)
     &           msavmax = op%ihpvca_occ(ivale,2,iocc_cls)
            mschmax = op%ihpvca_occ(ihole,1,iocc_cls)
            mscpmax = op%ihpvca_occ(ipart,1,iocc_cls)
            mscvmax = 0
            if (ivale.gt.0)
     &           mscvmax = op%ihpvca_occ(ivale,1,iocc_cls)

            first = .true.
            distr_loop: do
                  
              if (.not.next_msgamdist(first,
     &           msa, msc, igama, igamc, op%ihpvca_occ(1,1,iocc_cls),
     &           ngam,msd,igamd)) exit distr_loop
              first = .false.

              if (igah.eq.0) then
                lenah = 1
              else
                idxms = (msahmax-msd(1,2))/2+1
                lenah = str_info%g(igah)%lenstr_gm(igamd(1,2),idxms)
                if (lenah.eq.0) cycle distr_loop
              end if
              if (igap.eq.0) then
                lenap = 1
              else
                idxms = (msapmax-msd(2,2))/2+1
                lenap = str_info%g(igap)%lenstr_gm(igamd(2,2),idxms)
                if (lenap.eq.0) cycle distr_loop
              end if
              if (igav.eq.0) then
                lenav = 1
              else
                idxms = (msavmax-msd(3,2))/2+1
                lenav = str_info%g(igav)%lenstr_gm(igamd(3,2),idxms)
                if (lenav.eq.0) cycle distr_loop
              end if
              
              if (igch.eq.0) then
                lench = 1
              else
                idxms = (mschmax-msd(1,1))/2+1
                lench = str_info%g(igch)%lenstr_gm(igamd(1,1),idxms)
                if (lench.eq.0) cycle distr_loop
              end if
              if (igcp.eq.0) then
                lencp = 1
              else
                idxms = (mscpmax-msd(2,1))/2+1
                lencp = str_info%g(igcp)%lenstr_gm(igamd(2,1),idxms)
                if (lencp.eq.0) cycle distr_loop
              end if
              if (igcv.eq.0) then
                lencv = 1
              else
                idxms = (mscvmax-msd(3,1))/2+1
                lencv = str_info%g(igcv)%lenstr_gm(igamd(3,1),idxms)
                if (lencv.eq.0) cycle distr_loop
              end if
              
              ! increment distribution index
              idxdis = idxdis+1
              
              if (ntest.ge.150) then
                write(luout,*) 'current MS and IRREP distr:'
                call wrt_occ(luout,msd)
                call wrt_occ(luout,igamd)
                write(luout,*) 'accepted => current idxdis = ',idxdis
              end if

              ! save current offset
              if (ipass.eq.2) then
                op%off_op_gmox(iocc_cls)%
     &               d_gam_ms(idxdis,igama,idxmsa)=idxstr
                ! get ID of current distr
                did = msgmdid(op%ihpvca_occ(1,1,iocc_cls),
     &                        msd,igamd,ngam)
                ! save ID of current distr
                op%off_op_gmox(iocc_cls)%
     &               did(idxdis,igama,idxmsa) = did
              end if 

              ! increment string element index
              idxstr = idxstr+lencp*lenap*
     &                        lench*lenah*
     &                        lencv*lenav
              idxstr_tot = idxstr_tot+lencp*lenap*
     &                            lench*lenah*
     &                            lencv*lenav

              if (ipass.eq.2) then
                op%len_op_gmox(iocc_cls)%
     &               d_gam_ms(idxdis,igama,idxmsa) =
     &                 lencp*lenap*
     &                 lench*lenah*
     &                 lencv*lenav
              end if

              if (ntest.ge.150) then
                write(luout,*) 'current block: ',lencp,lenap,
     &                        lench,lenah,
     &                        lencv,lenav,' => ',lencp*lenap*
     &                        lench*lenah*
     &                        lencv*lenav
              end if
              
            end do distr_loop

            op%len_op_gmo(iocc_cls)%gam_ms(igama,idxmsa) = idxstr -
     &           op%off_op_gmo(iocc_cls)%gam_ms(igama,idxmsa)

            op%off_op_gmox(iocc_cls)%maxd =
     &           max(op%off_op_gmox(iocc_cls)%maxd,idxdis)

            if (ipass.eq.2) then
              op%off_op_gmox(iocc_cls)%ndis(igama,idxmsa) = idxdis
            end if

          end do igama_loop

        end do msa_loop

        op%len_op_occ(iocc_cls) = idxstr - op%off_op_occ(iocc_cls)

      end do occ_cls

      op%len_op = idxstr_tot

      if (ntest.ge.100) then
        if (ipass.eq.1) then
          write(luout,*) 'total number of operator elements: ',op%len_op
          write(luout,*) 'length per occupation class:'
          call iwrtma(op%len_op_occ,op%n_occ_cls,1,op%n_occ_cls,1)
          write(luout,*) 'offsets per occupation class:'
          call iwrtma(op%off_op_occ,op%n_occ_cls,1,op%n_occ_cls,1)
          write(luout,*) 'info per occupation class, IRREP, MS:'
          do iocc_cls = 1, op%n_occ_cls
            if(op%formal_blk(iocc_cls))cycle
            nexc = min(op%ica_occ(1,iocc_cls),
     &                 op%ica_occ(2,iocc_cls))
            write(luout,*) 'occ-class: ',iocc_cls
            write(luout,*) 'lengths:'
            call iwrtma(op%len_op_gmo(iocc_cls)%gam_ms,
     &           ngam,nexc+1,ngam,nexc+1)
            write(luout,*) 'offsets:'
            call iwrtma(op%off_op_gmo(iocc_cls)%gam_ms,
     &           ngam,nexc+1,ngam,nexc+1)
          end do
        else
          write(luout,*) 'info per occupation class, DISTR, IRREP, MS:'
          write(luout,*) 'offsets:'
          do iocc_cls = 1, op%n_occ_cls
            if(op%formal_blk(iocc_cls))cycle
            nexc = min(op%ica_occ(1,iocc_cls),
     &                 op%ica_occ(2,iocc_cls))
            write(luout,*) 'occ-class: ',iocc_cls
            call wrt_occ(luout,op%ihpvca_occ(1,1,iocc_cls))
            do iexc = 1, nexc+1
              do igam = 1, ngam
                if (op%off_op_gmox(iocc_cls)%ndis(igam,iexc).eq.0) cycle
                write(luout,*) iexc,igam,' -> ',
     &               op%off_op_gmox(iocc_cls)%
     &               d_gam_ms(1:op%off_op_gmox(iocc_cls)%
     &               ndis(igam,iexc),igam,iexc)
              end do
            end do            
          end do
          write(luout,*) 'distribution IDs:'
          do iocc_cls = 1, op%n_occ_cls
            if(op%formal_blk(iocc_cls))cycle
            nexc = min(op%ica_occ(1,iocc_cls),
     &                 op%ica_occ(2,iocc_cls))
            write(luout,*) 'occ-class: ',iocc_cls
            call wrt_occ(luout,op%ihpvca_occ(1,1,iocc_cls))
            do iexc = 1, nexc+1
              do igam = 1, ngam
                if (op%off_op_gmox(iocc_cls)%ndis(igam,iexc).eq.0) cycle
                write(luout,*) iexc,igam,' -> ',
     &               op%off_op_gmox(iocc_cls)%
     &               did(1:op%off_op_gmox(iocc_cls)%
     &               ndis(igam,iexc),igam,iexc)
              end do
            end do
          end do

        end if
        
      end if

      return
      end
*----------------------------------------------------------------------*
