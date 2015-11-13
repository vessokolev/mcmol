c===========================================================================
c     ztchu  10/30/2015  very simple version used for Vesolin's APP
c     called each time the total system energy needs to be calculated
c----------------------------------------------------------------------------
      subroutine system_e(E_tot)

      implicit none
      include 'include/parameter.dc'
      include 'include/Mmol.h'
      include 'include/Mdyn.h'
      include 'include/cparam_batyp.h'
      include 'include/coord.h'

c:::  input variables
      real*8 E_tot

c:::  local variables
      real*8 unit,ebondp,ethetp,ephi,eitor,evdwp,emumup,ehbpp,ene


         E_tot=0.0d0
c============================================================================
         unit=1.0

chu      [classical bond energy for protein]
         call bondp(nbp,ibtp,jbtp,icbp,unit,ebondp,.true.,1)
         print *,"ebondp",ebondp

chu      [classical angle energy for protein]
         call thetap(ntp,ittp,jttp,kttp,ictp,unit,ethetp,1,.true.)
         print *,"ethetp",ethetp

chu      [classical torsion energy for protein]
         l_t_it=.true.
         call  phip(nphi,ipt,jpt,kpt,lpt,icp,cpn,cpd,cp,unit,ephi,.true.)
         call phipsi(ene)
         ephi=ephi+ene
         print *,"ephi",ephi

chu      [classical improper torsion energy for protein]
         l_t_it=.false.
         call phip(nitor,iit,jit,kit,lit,ici,cpni,cpdi,cpi,unit,eitor,.true.)
         print *,"itor",eitor

chu      [nonbonded electrostatic and VDW energy between protein atoms]
         call nonbond(lnon,inon,jnon,iac,crg,unit,evdwp,emumup,ehbpp,1,.true.)
         print *,"nb",evdwp,emumup,ehbpp

         call simple_vdw

         print *,"simple vdw",sim_vdw

         E_tot=ebondp+ethetp+ephi+eitor+evdwp+emumup+ehbpp

         return

         end
c==============================================================================


