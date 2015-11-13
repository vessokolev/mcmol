      program assemble

      implicit none

      include "include/parameter.dc"
      include "include/coord.h"
      include "include/Mmol.h"
      include "include/cbkv_gap.h"
      include "include/sequ.h"
      include "include/symbols.h"
      include "include/cutoffs.h"
      include "include/Mdyn.h"
      include "include/cg_ion_hyd_pol_constants.h"

      integer i
      character*(1024) inppdb
      real*8 E_tot

      call phy_cnst

      inppdb="/home/vesso/cg/sso_100ps_relax.pdb"
!      inppdb="/home/vesso/monte/opt.pdb"

      cutpp=100.0d0
      c_scale=10.0

      call check_env(1)

      call pretopq(inppdb)
!      call topoin(.true.)
!      write(*,*) npro ! Number of protein atoms
!      write(*,*) bkvcode(34) ! Atom name



!      call task_analyze -> we don't need this since it is an
!                           interactive shell only. Instead
!                           allres need to be called.
      call allres
      call resang(1)

      call set_exclude(.false.)
      call find_center
      call nblist(1)

      default_cg_energy_c= .true.
      call init_cg(0)

      call write_pdb("/home/vesso/cg/test.pdb")

      call enz_default_p(0)

      print *,"here"

      do i=1,1000

         call simple_neighbor_self_e
         print *,"Self E ",t_sim_ngb
      end do


      call system_e(E_tot)
      print *,"Energy before perturbation:",E_tot

!      x(4)=x(4)-0.02
!      call system_e(E_tot)
!      print *,"Energy after perturbation: ",E_tot


!      call init_cg

      end

