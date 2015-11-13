      program assemble

      implicit none

      include "include/parameter.dc"
      include "include/coord.h"
      include "include/Mmol.h"
      include "include/cbkv_gap.h"
      include "include/sequ.h"
      include "include/symbols.h"
      include "include/cutoffs.h"

      integer i
      character*(1024) inppdb
      real*8 E_tot

      inppdb="/home/vesso/monte/opt.pdb"

      cutpp=40.0d0

      call check_env(1)
      call pretopq(inppdb)
      call topoin(.true.)
      write(*,*) npro ! Number of protein atoms
      write(*,*) bkvcode(34) ! Atom name

!      write(*,*) nnr

      ! The atom IDs supplied as integers per residue
      ! The residue atom ID starts from nnr(1)+1 and
      ! ends at nnr(2).
      do i=nnr(1)+1,nnr(2)
         print *,i
      end do

      do i=1,npro
         print *,i," ",iac_name(iac(i))," ",bkvcode(i)," ",iac_num(iac(i))
      end do


      write(*,*) x(1),x(2),x(3)
      write(*,*) x(4),x(5),x(6)

      print *,natom
      print *,nphi

      call set_exclude(.false.)
      call find_center
      call nblist(1)

      call enz_default_p(0)
      call system_e(E_tot)
      print *,"Energy before perturbation:",E_tot

      x(4)=x(4)-0.02
      call system_e(E_tot)
      print *,"Energy after perturbation: ",E_tot


      call init_cg

      end

