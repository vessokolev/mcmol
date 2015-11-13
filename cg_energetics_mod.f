c   This file, with all the subroutines, was written originally by 

c   Dr. Z.T.Chu, under the guidance of Dr. Warshel
c   It was originally written for the research purposes of Dr. Ben Messer
c   Since 2009, it was heavily modified by Spyridon Vicatos, Anna Rychkova
c   Anatoly Dryga, and Arieh Warshel, for research purposes
c   Most comments, written by Spyridon Vicatos
c
c   I apologize in advance for some inelegant comments, and some 
c   innefficient subroutines, but due to the lack of time
c   I had to do my best with coding, while at the same time the code had
c   to be as understandable as possible.
C   Spyridon Vicatos, Jan 2011


 
c--------------------------------------------------------------------------
      subroutine init_cg(hit_lib)

      implicit none

      include 'include/parameter.dc'
      include 'include/Mdyn.h'
      include 'include/seq.h'
      include 'include/cbkv_gap.h'
      include 'include/sequ.h'
      include 'include/coord.h'
      include 'include/dgm_pt.h'
      include 'include/cg_ion_hyd_pol_constants.h'


c     local var:
      integer i,ii,j,n_sim,i_sim,n_sim_type,i_polr(nsimpres),hit,hit_lib
      character*3 semonte_predict_charge,inpq_n,seq_n

c     This common is used ONLY in this subroutine!
      common/init_pol_data/i_polr

      write(6,*)"INIT CG IS RUNNING!!!!!"

c     The DATA for the simplified residues!!!!
c      include 'include/simp_data.h'
       call simp_data

       if(default_cg_energy_c) then
          call cg_energy_c_data
       endif

chu    get the nonbond list for each hydrophobic CG residue
chu    looks you guys copyed and renamed cg_nonbond_list to cb_nonbond_list1, 
chu    then we can removed cg_nonbond_list
chu        call cg_nonbond_list
c.......................................................................
c'

      n_sim_type=nsimpres    ! Set the number of simp types as nsimpres
c                            ! which is a parameter set into the parameter.dc file

c    Provide the library, where you set all the charges for the
c    Ionizable residues (aka Arg, Lys, Glu, Asp, His)

      if(.not.lsimple_crg) then     ! If the simple_charge file does not exist

c       Set all the "types" of residues to zero!
c       Set all "simple charges" to zero!

         do i=1,numres
            is_sim(i)=0
            i_sim_type(i)=0
            sim_q(i)=0.0
         enddo


c        assign the correct dielectrics 
         eff_diel(1) = 40.0
         eff_diel(2) = 10.0

c        Try to initialize the charges and the "polarities" correctly      

         do i=1,numres
            do j=1,5
               if(resseq(i).eq.typ_nam(j)) then
                  sim_q(i) = typ_crg(j)
                  is_sim(i)=1 !
                  goto 944
               endif
            enddo
 944     enddo

      endif

c     read in simple_side_lib, if user provided, otherwise use the above
c     default data 

      if(hit_lib.eq.1) then

         read(89,*)n_sim_type
         if(n_sim_type.gt.50) then
            write(6,'('' ERROR: # of simplified residue types > 50,'',
     $           '' there are no > 50 of folding types'',/,
     $           ''        in the amino.lib'')')
            goto 666
         endif

c       Read the VDW parameters and the "type" of the residue
         do i=1,n_sim_type
            read(89,*)seq_sim1(i),sim_a(i),sim_c(i),i_polr(i)
         enddo

c         if(n_sim_type.lt.38) then
c            write(6,'(/,'' ERROR: missing residues in your'',
c     $                  '' simple_side_lib'',/)')
c            goto 666
c         endif

      endif

      do i=1,numres

         if(seq(i).eq.'MEB')  then
            is_sim(i)=4      ! Set the "type" as "membrane"

         else

            do j=1,n_sim_type

               if(l_exp_simp) then
                  seq_n=ex_n(j)
               else
                  seq_n=seq_sim1(j)
               endif

               if(seq(i).eq.seq_n)then
                  i_sim_type(i)=j

                  if(lsimple_crg.and.is_sim(i).eq.1) goto 10

                  if(is_sim(i).eq.0)then

                     if(i_polr(j).eq.1) then
                        is_sim(i)=2 !polar

                     else
                        is_sim(i)=3 !nonpolar

                     endif
                  endif
                  goto 10
               endif
            enddo

         endif
        
 10   enddo

c     set charge to 0.0 for polar and nonpolar sidechain united atoms
      do i=1,numres
         if(is_sim(i).eq.2.or.is_sim(i).eq.3) then
            do ii=np(i)+1,np(i+1)
               if(bkvcode(ii)(1:2).eq.'CB')crg(ii)=0.0d0
            enddo
         endif
      enddo

      write(6,'(/,'' polarity type of the system:'')')
      write(6,'('' res#    resname   polarity (1: ionized, 2: polar,'',
     $            '' 3: nonpolr, 0: N/A)'')')
      do i=1,numres
         write(6,'(i5,6x,a3,8x,i2)')i,seq(i),is_sim(i)
      enddo
      print *


      close(89)
c.......................................................................
      return 
 666  call killme ('init_polarity')
      end
c-----------------------------------------------------------------------
      subroutine simple_vdw
      implicit none
      include 'include/parameter.dc'
      include 'include/coord.h'

      if (lsimple) then
          call coarseg_vdw ! Coarse Grained
      elseif (lsimple_alfa) then
          call simp_ca_vdw ! Simple CA
      else 
        write(6,*)"ERROR, the type of Simp Model was not specified correctly"
        stop " STOPPED in simple_vdw subroutine, due to error!"
      endif

      return
      end



c-----------------------------------------------------------------------
      subroutine coarseg_vdw
c


      implicit none

      include 'include/parameter.dc'
      include 'include/coord.h'
      include 'include/seq.h'
      include 'include/cbkv_gap.h'
      include 'include/Mdyn.h'
      include 'include/task.h'
      include 'include/derivative.h'
      include 'include/stepsize.h'
      include 'include/symbols.h'
      include 'include/sequ.h'

c     local var:
      integer i,j,ii,jj,k,ihit,jhit,i3,j3,m,iaci
      real*8 dx1,dx2,dx3,r1,r2,r3,r6,r8,a,b,c,df,df_vdw,c_factor,
     $       r_eq,r_eq2,r_eq4,r_eq6,r_eq8,r_eq12
      real*8 c_opt             ! scale constant for the collection of E_opt

      real*8  vdw_cut_param, r6old   ! vdw between protein and membrane (Anna Apr 2011)

c......................................................................
      sim_vdw=0.0
      pseudo_sim_vdw = 0.0d0
      sim_po_vdw=0.0
      e_opt = 0.0d0

c      write(6,*)" Coarse Grained VDW IS RUNNING"

      do i=1,numres-1
         ihit=0
         do j=i+1,numres
            if(seq(i).eq.'MEB'.and.seq(j).eq.'MEB') goto 30
            jhit=0
            if((is_sim(i).eq.1.and.is_sim(j).eq.2).or.     !ionic/polar
     $         (is_sim(j).eq.1.and.is_sim(i).eq.2).or.     !polar/ionic
     $         (is_sim(i).eq.3.and.is_sim(j).eq.3).or.     !nonpolr/nonpolr
     $         (is_sim(i).eq.1.and.is_sim(j).eq.1).or.     !ion/ion
     $         (is_sim(i).eq.2.and.is_sim(j).eq.2).or.     !polar/polar
     $         (is_sim(i).eq.1.and.is_sim(j).eq.3).or.     !ionic/nonpolar
     $         (is_sim(i).eq.3.and.is_sim(j).eq.1).or.     !nonpolar/ionic
     $         (is_sim(i).eq.2.and.is_sim(j).eq.3).or.     !polar/nonpolar
     $         (is_sim(i).eq.3.and.is_sim(j).eq.2)) then   !nonpolar/polar

               do ii=np(i)+1,np(i+1)
                  if(poly_is(ii).ne.0) then
                     iaci = iacpa(poly_is(ii))                  
                     if(iac_name(iaci).eq.'DY') goto 40
                     iaci = iacpb(poly_is(ii))                  
                     if(iac_name(iaci).eq.'DY') goto 40
                  endif
                  if(bkvcode(ii)(1:2).eq.'CB') then
                     ihit=ii
                     goto 10
                  endif
               enddo
               if(ihit.eq.0) goto 40

 10            do jj=np(j)+1,np(j+1)
                  if(poly_is(jj).ne.0) then
                     iaci = iacpa(poly_is(jj))                  
                     if(iac_name(iaci).eq.'DY') goto 30
                     iaci = iacpb(poly_is(jj))                  
                     if(iac_name(iaci).eq.'DY') goto 30
                  endif

                  if(bkvcode(jj)(1:2).eq.'CB') then
                     jhit=jj
                     goto 20
                  endif
               enddo
               if(jhit.eq.0) goto 30

 20            i3   = ihit*3-3
               j3   = jhit*3-3
               dx1  = x(i3+1)-x(j3+1)
               dx2  = x(i3+2)-x(j3+2)
               dx3  = x(i3+3)-x(j3+3)
               r2   = dx1*dx1+dx2*dx2+dx3*dx3
               r2   = 1.d0/r2
               r1   = sqrt(r2)            
               r3   = r1*r2
               r6   = r3*r3
               r8   = r2*r6

c       Declare some important scale constants
c########################################################################
c  
               c_factor=1.0
               c_opt = 1.30d0

               if((is_sim(i).eq.3.and.(is_sim(j).eq.1.or.is_sim(j).eq.2)).or.
     $            (is_sim(j).eq.3.and.(is_sim(i).eq.1.or.is_sim(i).eq.2))) then
                  c_factor=c_scale
               endif

               c=sqrt(sim_c(i_sim_type(i))*sim_c(i_sim_type(j)))/c_factor

               r_eq=sqrt(sim_a(i_sim_type(i))*sim_a(i_sim_type(j)))

c               write(6,*)"sim_c i        sim_c  j    e_factor    c"
c               write(6,'(f8.2,3x,f8.2,3x,f8.2,3x,f8.2)')
c     $          sim_c(i_sim_type(i)), sim_c(i_sim_type(j)), c_factor, c


c               write(6,*)
c               write(6,*)"sim_a i        sim_a  j    r_eq"
c               write(6,'(f8.2,3x,f8.2,3x,f8.2)')
c     $          sim_a(i_sim_type(i)), sim_a(i_sim_type(j)), r_eq


               r_eq2=r_eq*r_eq
               r_eq4=r_eq2*r_eq2
               r_eq6=r_eq4*r_eq2
               r_eq8=r_eq4*r_eq4



               a = 3.0*r_eq8*r8
               b = -4.0*r_eq6*r6

c               write(6,*)
c               write(6,*)"a         b "
c               write(6,'(f8.2,3x,f8.2)')
c     $          a,b





               df_vdw=8.0*a+6.0*b
               df=-r2*c*df_vdw

               do m=1,npoly
                  if(ihit.eq.ipoly(m).or.jhit.eq.ipoly(m))then
                     sim_po_vdw=sim_po_vdw+c*(a+b)
                  endif
               enddo

               do m=1,n_mut_gly
                  if(ihit.eq.i_mut_gly(m).or.jhit.eq.i_mut_gly(m))then
                     sim_po_vdw=sim_po_vdw+c*(a+b)
                     df=df*formwgt
                  endif
               enddo

c
c
c
c
c   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   
c   Calculate van der waals interactions
c   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c          |||||||||||||||||
c               |||||||
c                 |||
c                  V


               sim_vdw=sim_vdw+c*(a+b)

                if(1.0d0/r1.lt.r_eq) then
                    pseudo_sim_vdw = pseudo_sim_vdw +c*(3.0*r_eq8/r_eq8  -4.0*r_eq6/r_eq6  )
                else
                    pseudo_sim_vdw = pseudo_sim_vdw +c*(a+b)
                endif 

                
c              write(6,*)"PSEUDO VDW is ", pseudo_sim_vdw

c              write(6,*)"Contribution ",i,j,c*(a+b), sim_vdw 
c              write(6,*)
c               write(6,*)" sim_vdw   ihit    jhit  r1"
c               write(6,*) sim_vdw, ihit, jhit, 1.0/r1
c               write(6,*)
c
c
c
c


c###############################################################
c###############################################################
C    CHECKING AND COLLECTING values for the E_opt
c


               if(1.0d0/r1.gt.r_eq*c_opt) then
                e_opt = e_opt + c*(a+b)
c            write(6,*)'DEBUG!!!! Residue, Residue, E_OPT '
c            write(6,*)i,j,e_opt,c*(a+b)               
               else
 
                e_opt = e_opt - c*c_factor 
c           write(6,*)'DEBUG!!!! Residue, Residue, E_OPT '               
c           write(6,*)i,j,e_opt,c*c_factor               
               endif


c       if((i.ge.63).and.(i.le.67)) Then

c       write(6,*)'DEBUG!!!! Residue, Residue, r1, r_eq VDW  c*factor '
c'
c       write(6,*)i,j,1.0d0/r1, r_eq, c*(a+b),c*c_factor

c       endif
         
c      write(6,*)'e_hydro',i,e_hydro


c###############################################################




               if(check_sim.eq.1) then
                  write(6,'('' simple polarity for residues:'',2i4,
     $                      ''   side atoms:'',2i6)')i,j,ihit,jhit
                  write(6,'('' r_eq:'',2f10.3,'' C_factor:'',2f10.3)')
     $                      sim_a(i_sim_type(i)),sim_a(i_sim_type(j)),
     $                      sim_c(i_sim_type(i)),sim_c(i_sim_type(j))
                  write(6,'('' c_scaling factor:'',f8.2)')c_factor
                  write(6,'('' rij :'',f10.3,11x,'' sim_vdw:'',f10.3,
     $                      ''  sum:'',f10.3,/)')
     $                      1.0/r1,c*(a+b),sim_vdw
               endif

 35            continue

               d(i3+1)=d(i3+1)+df*dx1
               d(i3+2)=d(i3+2)+df*dx2
               d(i3+3)=d(i3+3)+df*dx3
               d(j3+1)=d(j3+1)-df*dx1
               d(j3+2)=d(j3+2)-df*dx2
               d(j3+3)=d(j3+3)-df*dx3
            endif

c#################################################################
c    vdw between protein and membrane atoms [Anna, Apr 2011]

            if((is_sim(i).eq.4.and.is_sim(j).ne.4).or. !CB of MEB with CB of protein
     $         (is_sim(j).eq.4.and.is_sim(i).ne.4)) then

               do ii=np(i)+1,np(i+1)
                  if(poly_is(ii).ne.0) then
                     iaci = iacpa(poly_is(ii))
                     if(iac_name(iaci).eq.'DY') goto 40
                     iaci = iacpb(poly_is(ii))
                     if(iac_name(iaci).eq.'DY') goto 40
                  endif
                  if(bkvcode(ii)(1:2).eq.'CB') then
                     ihit=ii
                     goto 50
                  endif
               enddo
               if(ihit.eq.0) goto 40

 50            do jj=np(j)+1,np(j+1)
                  if(poly_is(jj).ne.0) then
                     iaci = iacpa(poly_is(jj))
                     if(iac_name(iaci).eq.'DY') goto 30
                     iaci = iacpb(poly_is(jj))
                     if(iac_name(iaci).eq.'DY') goto 30
                  endif
                  if(bkvcode(jj)(1:2).eq.'CB') then
                     jhit=jj
                     goto 60
                  endif
               enddo
               if(jhit.eq.0) goto 30

 60            i3   = ihit*3-3
               j3   = jhit*3-3
               dx1  = x(i3+1)-x(j3+1)
               dx2  = x(i3+2)-x(j3+2)
               dx3  = x(i3+3)-x(j3+3)
               r2   = dx1*dx1+dx2*dx2+dx3*dx3
               r2   = 1.d0/r2
               r1   = sqrt(r2)
               r3   = r1*r2
               r6   = r3*r3

               vdw_cut_param = 7452.75

c   vdw_cut_param was calculated as an average of all the solutions {V[r]=0,r=0} for all
c   possible interactions between CB of 19 residues and CB of MEB.
c   Equilibrium diameters and well depths for CB atoms were taken from the simp_data.f. [Anna]

               r6 = 1.0d0/r6
               r6old = r6        ! Actual r^6
               r6 = r6 + vdw_cut_param
               r6 = 1.0d0/r6     ! 1 / (vdw_cut_param + r^6)

               c=sqrt(sim_c(i_sim_type(i))*sim_c(i_sim_type(j)))
               r_eq=0.5d0*(sim_a(i_sim_type(i))+sim_a(i_sim_type(j)))

               r_eq2=r_eq*r_eq
               r_eq4=r_eq2*r_eq2
               r_eq6=r_eq4*r_eq2
               r_eq12=r_eq6*r_eq6

               a = 4.0d0*c*r_eq12
               b = -4.0d0*c*r_eq6

               sim_vdw=sim_vdw+(a*r6*r6+b*r6)
               pseudo_sim_vdw = pseudo_sim_vdw + (a*r6*r6+b*r6)

c               write(6,*)'DEBUG: MEB-prot'
c               write(6,'(''i j'',2i3)')i,j
c               write(6,'(''i_sim_type(j) i_sim_type(j)'',2i2)')i_sim_type(i),i_sim_type(j)
c               write(6,'(''well_depth eq_radius for i'',2f10.2)')sim_c(i_sim_type(i)),sim_a(i_sim_type(i))
c               write(6,'(''well_depth eq_radius for j'',2f10.2)')sim_c(i_sim_type(j)),sim_a(i_sim_type(j))
c               write(6,'(''c r_eq ''2f10.2)'),c,r_eq
c               write(6,'(''dist'',f10.2)')1/r1
c               write(6,'(''sim_vdw'',f20.2)')a*r6*r6+b*r6
 
               df_vdw = (12*a*r6 + 6*b)*r6old*r6*r6
               df = -r2*df_vdw

               d(i3+1)=d(i3+1)+df*dx1
               d(i3+2)=d(i3+2)+df*dx2
               d(i3+3)=d(i3+3)+df*dx3
               d(j3+1)=d(j3+1)-df*dx1
               d(j3+2)=d(j3+2)-df*dx2
               d(j3+3)=d(j3+3)-df*dx3

            endif

c end of the vdw between membrane and protein
c##################################################################



 30      enddo




 40   enddo
c...........................................................................

c      write(6,*)"COARSE VDW RESULT IS ", sim_vdw
      return
      end
c---------------------------------------------------------------------------



  
c-----------------------------------------------------------------------
      subroutine simp_ca_vdw
c


      implicit none

      include 'include/parameter.dc'
      include 'include/coord.h'
      include 'include/seq.h'
      include 'include/cbkv_gap.h'
      include 'include/Mdyn.h'
      include 'include/task.h'
      include 'include/derivative.h'
      include 'include/stepsize.h'
      include 'include/symbols.h'
      include 'include/sequ.h'

c     local var:
      integer i,j,ii,jj,k,ihit,jhit,i3,j3,m,iaci
      real*8 dx1,dx2,dx3,r1,r2,r3,r6,r8,a,b,c,df,df_vdw,c_factor,
     $       r_eq,r_eq2,r_eq4,r_eq6,r_eq8
      real*8 r_i, r_j, epsilon_i, epsilon_j,alanine_r,alanine_eps

      real*8 c_opt             ! scale constant for the collection of E_opt

c......................................................................
      sim_vdw=0.0
      sim_po_vdw=0.0
      e_opt = 0.0d0

c      write(6,*)" Simple C_Alpha VDW IS RUNNING"

C     parameters for the CA
      alanine_r = 4.24
      alanine_eps = 0.050d0

      do i=1,numres-1
         do j=i+1,numres

            if(seq(i).eq.'MEB'.and.seq(j).eq.'MEB') goto 30  ! get out of the j loop

            if((is_sim(i).eq.1.and.is_sim(j).eq.2).or.     !ionic/polar
     $         (is_sim(j).eq.1.and.is_sim(i).eq.2).or.     !polar/ionic
     $         (is_sim(i).eq.3.and.is_sim(j).eq.3).or.     !nonpolr/nonpolr
     $         (is_sim(i).eq.1.and.is_sim(j).eq.1).or.     !ion/ion
     $         (is_sim(i).eq.2.and.is_sim(j).eq.2).or.     !polar/polar
     $         (is_sim(i).eq.1.and.is_sim(j).eq.3).or.     !ionic/nonpolar
     $         (is_sim(i).eq.3.and.is_sim(j).eq.1).or.     !nonpolar/ionic
     $         (is_sim(i).eq.2.and.is_sim(j).eq.3).or.     !polar/nonpolar
     $         (is_sim(i).eq.3.and.is_sim(j).eq.2)) then   !nonpolar/polar

               do ii=np(i)+1,np(i+1)

C                define the vdw variables for residue i
                  if(bkvcode(ii)(1:2).eq.'CB') then
                   r_i = sim_a(i_sim_type(i))
                   epsilon_i =  sim_c(i_sim_type(i))
                  else if(bkvcode(ii)(1:2).eq.'CA') then
                   r_i = alanine_r
                   epsilon_i = alanine_eps
                  endif


                do jj=np(j)+1,np(j+1)
                 
C                define the vdw variables for residue j
                  if(bkvcode(jj)(1:2).eq.'CB') then
                   r_j = sim_a(i_sim_type(j))
                   epsilon_j =  sim_c(i_sim_type(j))
                  else if(bkvcode(ii)(1:2).eq.'CA') then
                   r_j = alanine_r
                   epsilon_j = alanine_eps
                  endif

c              Calculate vdw values
c 

               i3   = ii*3-3
               j3   = jj*3-3
               dx1  = x(i3+1)-x(j3+1)
               dx2  = x(i3+2)-x(j3+2)
               dx3  = x(i3+3)-x(j3+3)
               r2   = dx1*dx1+dx2*dx2+dx3*dx3
               r2   = 1.d0/r2
               r1   = sqrt(r2)            
               r3   = r1*r2
               r6   = r3*r3
               r8   = r2*r6


c       Declare some important scale constants
c########################################################################
c  
               c_factor=1.0
               c_opt = 1.30d0


               if((is_sim(i).eq.3.and.(is_sim(j).eq.1.or.is_sim(j).eq.2)).or.
     $            (is_sim(j).eq.3.and.(is_sim(i).eq.1.or.is_sim(i).eq.2))) then

                   c_factor=c_scale

               endif

               c=sqrt(epsilon_i*epsilon_j)/c_factor
               r_eq=sqrt(r_i*r_j)


               r_eq2=r_eq*r_eq
               r_eq4=r_eq2*r_eq2
               r_eq6=r_eq4*r_eq2
               r_eq8=r_eq4*r_eq4



               a = 3.0*r_eq8*r8
               b = -4.0*r_eq6*r6


               df_vdw=8.0*a+6.0*b
               df=-r2*c*df_vdw

c        NOTE
c        Some parts of the original vdw subroutine that addreses the
c        Calculations of sim_po_vdw are ommited!!!!

c
c
c
c
c   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   
c   Calculate van der waals interactions
c   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c          |||||||||||||||||
c               |||||||
c                 |||
c                  V


               sim_vdw=sim_vdw+c*(a+b)

c            write(6,*)"SIM CA "
c            write(6,*)"R =",1/r1
c            write(6,*)"a and b =",a,b
c            write(6,*)"c =",c
c            write(6,*)"contr = ",c*(a+b)
c            write(6,*)"VDW TOTAL IS ",sim_vdw
c            write(6,*)
c


c###############################################################
c###############################################################
C    CHECKING AND COLLECTING values for the E_opt
c


               if(1.0d0/r1.gt.r_eq*c_opt) then
                e_opt = e_opt + c*(a+b)              
               else
 
                e_opt = e_opt - c*c_factor         
               endif


c###############################################################
c         UPDATE THE FORCES!

               d(i3+1)=d(i3+1)+df*dx1
               d(i3+2)=d(i3+2)+df*dx2
               d(i3+3)=d(i3+3)+df*dx3
               d(j3+1)=d(j3+1)-df*dx1
               d(j3+2)=d(j3+2)-df*dx2
               d(j3+3)=d(j3+3)-df*dx3


         enddo   ! End of do jj=np(j)+1,np(j+1)
         enddo   ! End of do ii=np(i)+1,np(i+1)



         endif    ! End of the long if statement where (
c                      is_sim(i).eq.1.and.is_sim(j).eq.2).or.     !ionic/polar
c                      etc!!!!




 30      enddo
      enddo
c...........................................................................


      return
      end
c---------------------------------------------------------------------------




      subroutine simple_neighbor_self_e
c     calculate the interaction of ion and nonpolar residues with polar and nonpolar neighbors 
c     and membrane atoms
 
      implicit none

      include 'include/parameter.dc'
      include 'include/coord.h'
      include 'include/seq.h'
      include 'include/cbkv_gap.h'
      include 'include/Mdyn.h'
      include 'include/task.h'
      include 'include/derivative.h'
      include 'include/symbols.h'
      include 'include/stepsize.h'
      include 'include/lib.h'
      include 'include/cg_ion_hyd_pol_constants.h'
      include 'include/dgm_pt.h'
      include 'include/rWatProtein.h'
      include 'include/filename.h'

c     local var:>
      integer i,j,ii,jj,k,m,ihit,jhit,i3,j3,irun,iw,
     $        hit,hit_gly,option,
     $        icount,par_derivative_type, res_type,res_type0, get
      real*8 dx1,dx2,dx3,r1,r2,r3,df,u_temp,temp_u,e_self_np,e_self_p,e_self_mem,DM,
     $       temp_u_np,temp_u_mem, temp_u1, temp_u2

      real*8 c_polar, c_np, c_mem,temp_energy
      real*8 max_polar, max_nonpolar,  max_mem, r_min_cote1, r_min_cote2

      real*8 c_pol(19), c_nonpol(19), c_membrane(19), coefficient_np(19),
     $ coefficient_p(19), coeff_polar
      real*8 bins_polar(19,10), bins_nonpolar(19,20)
      integer bin

      integer current_i_type
      real*8 fact_u1, fact_u2

      real*8 temp_nonpolar, temp_polar, temp_mem, temp_general
 
      real*8 sum_t_np_hydro, sum_t_p_hydro, sum_t_mem_hydro
      real*8 sum_t_np_polar, sum_t_p_polar, sum_t_mem_polar

      integer Nw(numres), Nw_ring(numres), Nw_water(20)
      integer ng_range_old, start_old, end_old

c     The integers where the number of the following residues are stored
      integer n_ala, n_leu, n_ile, n_val, n_prol, n_met, n_phe, n_trp
      integer n_ser, n_thr, n_tyr, n_cys, n_asn, n_gln

c     Langevin grid variable (part are obsolete)
      real*8 r_ld_min0, alpha_min,x1_1,y1,y2, cote

c     all variables describing neighbours are real now (so we can have 1.67 neighbours)
c     (required for correct usage of F(rij) = exp (-6(r-rnp)^2))
      real*8 t_mem_ngb, t_p_ngb, t_np_ngb,  c_mem_factor

c     F(rij) values (see Rychkova, Vicatos, Warshel PNAS)
      real*8 mem_sim_ngb, p_sim_ngb, np_sim_ngb
      
c     cutoff radius for mem neighbours
c     will be assigned to mem_spacing*f (f=2.05)
c      real*8 r_mem_ngb_cutoff

c     number of neighbours for ion in the membrane
c     should not strongly depend on membrane spacing
c     mem in two solvation shells are included
c     within r<mem_spacing*f
      real*8 max_mem_ngb
 

c    closest distance to water from participated residue
      real*8 distanceToWater 

      character seq_tmp*3, type*5
      logical debug, print_option
 
c      real*8 E_main_CG  ! main chain solvatation penalty
c      real*8 E_hydro_CG ! Hbond in CG model
      real*8 MainChainSolvation  ! function for calculation
      real*8 HydrogenBond_CG  ! function for calculation
c......................................................................


      if(run_old_self_e) then
          call old_self_energy

          return
      endif

      get = 1

      max_polar = 6.0d0     ! At 6/25/2012 was 10.0 . At 6/27/2012 as a test, became 6
      max_nonpolar = 15.0d0
      max_mem_ngb = 28.0d0  !number is approximate smth between 25 and 30 will probably work

 
C     The two factors for the hydrophobic energy
      fact_u1 = 1.0d0                   ! for Nw
      fact_u2 = 0.0d0         ! for N neighbor

      r_min_cote1 = 18.0d0
      r_min_cote2 = 12.0d0

      icount=icount+1


      debug=.true.
      debug=.false.

      print_option=.false.

      if(lcheck_self) debug=.true.

c     is_sim()=1   ion
c     is_sim()=2   polar
c     is_sim()=3   nonpolar
c     is_sim()=4   membrane
      

c     Keep this for development 
c      write(6,*)"Success!!!! simple neighbor is running"


      t_sim_ngb=0.0


c      DEFAULT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c      Nw_water(1) = 26    ! ALA
c      Nw_water(2) = 26    ! VAL
c      Nw_water(3) = 38    ! LEU
c      Nw_water(4) = 54    ! ILE
c      Nw_water(5) = 26    ! PRO






c      Nw_water(6) = 100    ! LYS
c      Nw_water(7) = 100   ! HIS
c      Nw_water(8) = 54    ! PHE
c      Nw_water(9) = 100    ! TYR
c      Nw_water(10) = 90    ! TRP    ! could be 90 as well, its under debate
c      Nw_water(11) = 100    ! ASN
c      Nw_water(12) = 100    ! GLN






c      Nw_water(13) = 100   ! SER
c      Nw_water(14) = 100   ! THR
c      Nw_water(15) = 100    ! ARG
c      Nw_water(16) = 100    ! ASP
c      Nw_water(17) = 100    ! GLU
c      Nw_water(18) = 100    ! CYS
c      Nw_water(19) = 26     ! MET




C     Set some fictitious values for Nw_water

c      'A*2','V*2','L*2','I*2','P*2','K*2','H*2','F*2','Y*2','W*2',
c      'N*2','Q*2','S*2','T*2','R*2','D*2','E*2','C*2','M*2','H*2',
c      'A*3','T*3','C*3','G*3' 
      Nw_water(1) = 60    ! ALA
      Nw_water(2) = 110    ! VAL
      Nw_water(3) = 115    ! LEU
      Nw_water(4) = 120    ! ILE
      Nw_water(5) = 50    ! PRO
      Nw_water(6) = 100    ! LYS
      Nw_water(7) = 100   ! HIS
      Nw_water(8) = 130    ! PHE
      Nw_water(9) = 100    ! TYR
      Nw_water(10) = 140    ! TRP    ! could be 90 as well, its under debate
      Nw_water(11) = 100    ! ASN
      Nw_water(12) = 100    ! GLN
      Nw_water(13) = 100   ! SER
      Nw_water(14) = 100   ! THR
      Nw_water(15) = 100    ! ARG
      Nw_water(16) = 100    ! ASP
      Nw_water(17) = 100    ! GLU
      Nw_water(18) = 100    ! CYS
      Nw_water(19) = 110     ! MET



c     variables for the gyration calculation
c      ng_range_old = ng_range
c      start_old = ig_start(1)
c      end_old = ig_end(1)

c      ng_range = 1
c      ig_start(ng_range) = 1
c      ig_end(ng_range) = numres
            
c      call get_gyration_new
 
c     Now revert to the old values, so tha the gyration subroutine be normaly used
c     in other levels of MOLARIS
 
c      ng_range = ng_range_old
c      ig_start(1) = start_old
c      ig_end(1) = end_old
                 



C     Reset the constants for all residues
      do i=1,19
          c_pol(i)=0.0d0
          c_nonpol(i) = 0.0d0
          c_membrane(i) = 0.0d0
          coefficient_np(i) = 0.0d0
          coefficient_p(i) = 0.0d0
          do j=1,10
              bins_polar(i,j) = 0.0d0
          enddo 
          do j=1,20
               bins_nonpolar(i,j) = 0.0d0
          enddo 


      enddo

      coeff_polar = 0.0d0

C     Call subroutine that calculates and assigns ALL the constants!
      call assign_cg_constants(c_pol,c_nonpol,c_membrane)


c Reassign c_membrane constants for np resi to get correct membrane contribution
c Increse by 4 for debugging purposes

c      c_mem_factor=2.15

c      c_membrane(6) =  -0.323*c_mem_factor
c      c_membrane(7) =  -0.965*c_mem_factor
c      c_membrane(8) =  -0.643*c_mem_factor
c      c_membrane(9) =  -0.698*c_mem_factor
c      c_membrane(10) =  -1.028*c_mem_factor
c      c_membrane(11) =  -0.643*c_mem_factor
c      c_membrane(12) =  -1.066*c_mem_factor
c      c_membrane(13) =  -1.397*c_mem_factor




c     Constant for the LD grid calculations, IF the LD is used! (obsolete)
      alpha_min = 0.50d0


c     Initializing 
      e_hydro = 0.0d0
      e_hydro_np = 0.0d0
      e_hydro_mem = 0.0d0
      e_hydro_p = 0.0d0

      e_polar = 0.0d0
      e_polar_np = 0.0d0
      e_polar_mem = 0.0d0
      e_polar_p   = 0.0d0


c     This is not very necessary, but it was important to check
c     these values with the e_hydro etc values
C     Eventualy they will be removed
      sum_t_np_hydro = 0.0d0
      sum_t_p_hydro = 0.0d0
      sum_t_mem_hydro = 0.0d0

      sum_t_np_polar = 0.0d0
      sum_t_p_polar = 0.0d0
      sum_t_mem_polar = 0.0d0

c     Initialize and set the N-neighbors to zero
      do i =1,numres
          Npol_CB(i) = 0.0d0 !  Polar neighbors of residue i
          Nnonpol_CB(i) = 0.0d0 ! NON  Polar (aka hydrophobic) neighbors of residue i
          Nmem_CB(i) = 0.0d0 !  Membrane neighbors of residue i
      enddo

c     Set the partial derivatives of the "function" f to zero
      do i=1,5
         pder_p(i) = 0.0d0
         pder_np(i) = 0.0d0
      enddo

C     Reset the residue type!
      par_derivative_type = 0
      
      if (has_membrane) then
          DM = mem_spacing 
          r_mem_ngb_cutoff = DM*2.05 ! only neighbours in the "two solvation shells" counted as 1
                                    ! at further distance exp(...) function is used
c          write(6,'('' membrane spacing =''f8.3)') DM
          if ( (DM.le.0.5).or.(DM.ge.5.0) ) then
              write (*,*) 
     $             'WARNING: DM - membrane distance is unphysical.'
              write (*,*) 'WARNING: DM = ', DM
              write (*,*) 'Program stopped in cg_energetics.f'
              stop
          endif
      endif

c     Write the "label"  for the Table!
c     ################

      if(print_option) then
          option = 0
          call print_self_energy(option,i,seq_tmp, c_polar, c_np, c_mem,
     $     temp_polar, temp_nonpolar, temp_mem, min_wat_distance(1), temp_energy,
     $     t_p_ngb,t_np_ngb,t_mem_ngb)


              call print_hydro_energy(option,i,seq_tmp, c_polar, c_np, c_mem,
     $             temp_polar, temp_u1, temp_u2, temp_mem, min_wat_distance(1), temp_energy,
     $             t_p_ngb, t_np_ngb,temp_general,t_mem_ngb,fact_u1,fact_u2)



      endif


c     Call subroutine that calculates Nw (Number of waters)
c     if the residue is hydrophobic
c


      call find_Nw(Nw, Nw_ring, print_option)


      total_p_res = 0
      total_np_res = 0
      total_ion_res = 0

      n_ala = 0
      n_leu = 0
      n_ile = 0
      n_val = 0
      n_prol = 0
      n_met = 0
      n_phe = 0
      n_trp = 0

      n_ser = 0
      n_thr = 0
      n_tyr = 0
      n_cys = 0
      n_asn = 0
      n_gln = 0

c#########################################################################################
c#########################################################################################
c     STARTING OP THE LOOP to calculate the energy cotribution for EACH residue
C     IF residue (i) is ionizable, the contribution would be the self energy (i)
C     IF residue (i) is hydrophobic, the contribution would the hydrophobic contribution
C     IF residue (i) is polar, the contribution would the polar contribution
      do 40 i=1,numres ! loop over all ions and nonpolar residues 


          if(is_sim(i).eq.2) total_p_res = total_p_res +1
          if(is_sim(i).eq.3) total_np_res = total_np_res +1
          if(is_sim(i).eq.1) total_ion_res = total_ion_res +1


          if(is_sim(i).eq.4) goto 40   !go to the next residue  if current is MEB
 
              res_type = 0
C         Call the subroutine which simply "recognises" the residue
c         and gives the appropriate value of the "type"
          call find_res_type(get,seq_sim1(i_sim_type(i)), res_type,res_type0)

c         Set the values for residue i (based on the "type" of i)
          c_polar = c_pol(res_type)
          c_np = c_nonpol(res_type)
          c_mem = c_membrane(res_type)

          current_i_type = res_type

c         SET important variables to zero!
c         ###################################################
          t_mem_ngb=0.0  !number of membrane neighbours
          t_p_ngb=0.0    ! .....    polar ....
          t_np_ngb=0.0   ! ....     nonpolar ...
          e_self_mem=0.0  !contib to self_e from membrane
          e_self_np=0.0  !contib to self_e from nonpolar
          e_self_p=0.0  !contib to self_e from polar residues
C         First, set the temp_energy term for residue i to zero
          temp_energy = 0.0d0 !total sel_e fo a given (i-th) residue
          temp_nonpolar = 0.0d0
          temp_polar =0.0d0
          temp_mem = 0.0d0

          ihit=0
          do ii=np(i)+1,np(i+1)
             if(iac_name(iac(ii)).eq.'DY') goto 40
             if(bkvcode(ii)(1:2).eq.'CB') then
                 ihit=ii
             endif
          enddo
          if(ihit.eq.0) goto 40

c         Find the distance of the CB atom and its center!
c         Important for the hydrophobic calculations

          i3   = ihit*3-3
          r3 = (x(i3+1)-sys_xcen(1))**2
     $         + (x(i3+2)-sys_xcen(2))**2
     $         + (x(i3+3)-sys_xcen(3))**2

          r3 = dsqrt(r3)

          if(is_sim(i).eq.3) then      ! if i is hydrophobic
             
              if(res_type.eq.6) n_ala = n_ala+1
              if(res_type.eq.7) n_leu = n_leu+1
              if(res_type.eq.8) n_ile = n_ile+1
              if(res_type.eq.9) n_val = n_val+1
              if(res_type.eq.10) n_prol = n_prol+1
              if(res_type.eq.11) n_met = n_met+1
              if(res_type.eq.12) n_phe = n_phe+1
              if(res_type.eq.13) n_trp = n_trp+1

          elseif (is_sim(i).eq.2) then    ! if i is polar

              if(res_type.eq.14) n_ser = n_ser+1
              if(res_type.eq.15) n_thr = n_thr+1
              if(res_type.eq.16) n_tyr = n_tyr+1
              if(res_type.eq.17) n_cys = n_cys+1
              if(res_type.eq.18) n_asn = n_asn+1
              if(res_type.eq.19) n_gln = n_gln+1
          else
          endif




C        #################      Loop for i = 1 and i = 2
c        during irun = 1 number of neighbors are calculated and then Self energy for a given residue
c        during irun = 2 forces on atoms are calculated
C        #########################################################
         do 50 irun=1,2




c           #########################################################################################
c           #########################################################################################
c           LOOP OVER ALL RESIDUES!
c           Residue (i) will "interact with each residue j, and all contributions will be calculated
c           #########################################################################################

             do 30 j=1,numres                      !loop over all polar or nonpolar, or ion or mem


                 if(i.eq.j) goto 30
                 jhit=0

                 do jj=np(j)+1,np(j+1)
                     if(iac_name(iac(jj)).eq.'DY') goto 30
                     if(bkvcode(jj)(1:2).eq.'CB') then
                         jhit=jj
                     endif
                 enddo

                 if(jhit.eq.0) goto 30

c              inefficient part : distances are calculated two times, for irun = 1,2
C              Eventually this needs to be fixed.
                 i3   = ihit*3-3
                 j3   = jhit*3-3
                 dx1  = x(i3+1)-x(j3+1)
                 dx2  = x(i3+2)-x(j3+2)
                 dx3  = x(i3+3)-x(j3+3)
                 r2   = dx1*dx1+dx2*dx2+dx3*dx3
                 r1   = sqrt(r2)            
               !write (*,*) i, j, x(i3+1),x(j3+1), x(i3+2),x(j3+2),  x(i3+3),x(j3+3)


C              HIT CLYCINE, need more explanation from Dr. Chu!
                 hit_gly=0
                 if(npoly.gt.0) then
                     do m=1,npoly
                         if(ihit.eq.ipoly(m).or.jhit.eq.ipoly(m)) hit_gly=1
                    enddo
                 else
                    do m=1,n_mut_gly
                     if(ihit.eq.i_mut_gly(m).or.jhit.eq.i_mut_gly(m))then
                        hit_gly=1
                     endif
                  enddo
               endif

C 
c             irun = 1  
C             --------------------------------------------------------------------------------------
              if(irun.eq.1)  then ! get total Nnp and Np and Nmem (for future calculation of energy)

                                                             !i is ion or nonpolar or polar
                  if(is_sim(j).eq.2.or.is_sim(j).eq.1) then  !j is polar or ion
                     if(r1.le.r_p) then
                        p_sim_ngb=1.0
                     else
                        p_sim_ngb=exp(-a_p*(r1-r_p)**2)      !need to write the equation in the comment
                     endif
                     t_p_ngb=t_p_ngb+p_sim_ngb
                     if(debug) write(6,'('' ion-pol, residues  r  N_p:'',2i4,
     $                         f12.3,i5)')i,j,r1,p_sim_ngb

c                                                            !i is ion or nonpolar or polar
                  else if(is_sim(j).eq.3) then               !j is nonpolar
                     if(r1.le.r_np) then
                        np_sim_ngb=1.0
                     else
                        np_sim_ngb=exp(-a_np*(r1-r_np)**2)  !need to write the equation in the comment
c                        np_sim_ngb= 0.0d0  !need to write the equation in the comment
                     endif
                     t_np_ngb=t_np_ngb+np_sim_ngb
                     if(debug) write(6,'('' ion-npol,residues  r  N_np'',2i4,
     $                         f12.3,i5)')i,j,r1,np_sim_ngb

c                                                            !i is ion or nonpolar or polar
                  else if(is_sim(j).eq.4) then               !j is membrane
                     ! if below is formula 11 in Rychkova et al PNAS but rcutoff is new
                     ! mem_sim_ngb is F(rij) in article's notation
                     if(r1.le.r_mem_ngb_cutoff) then
                        mem_sim_ngb=1.0
                     else
                        mem_sim_ngb=exp(-a_mem*(r1-r_mem_ngb_cutoff)**2) !need to write the equation in the comment
                     endif
                     ! t_mem_ngb is N(i)mem in article's notation
                     ! t_mem_ngb meaning - number of mem residues around given residue
                     t_mem_ngb=t_mem_ngb+mem_sim_ngb
                     if(debug) write(6,'('' ion-membrane,residues  r  N_np'',
     $                         2i4,f12.3,i5)')i,j,r1,mem_sim_ngb

                  endif     


c             irun = 2  
c             after getting total Nnp (t_np_ngb)  and Np (t_p_ngb), 
c             which are functions of r_ij, calculating forces for atom  j
C             --------------------------------------------------------------------------------------
C             --------------------------------------------------------------------------------------
C             --------------------------------------------------------------------------------------
              else if (irun.eq.2)   then

c
c              Calculate u_temp and eventually df
c              -------------------------------------------------------------------------------------  

c                                                            ! i is ion or nonpolar or polar (ANYTHING)
                  if(is_sim(j).eq.2.or.is_sim(j).eq.1) then    ! if j is ion or polar
                      if(r1.le.r_p) then
                          goto 30
                      else

                          if(t_p_ngb.le.max_polar) then
                              u_temp= c_polar*exp(-a_p_u*(t_p_ngb-max_polar)**2)
                              df=u_temp*(-a_p_u*2.0d0*(t_p_ngb-max_polar))
     $                        *exp(-a_p*(r1-r_p)**2)*2.0d0*(-a_p*(r1-r_p))/r1
                          else
                              goto 30
                          endif
                      endif
                                                             ! i is ion or nonpolar or polar (ANYTHING)
                  else if(is_sim(j).eq.3) then               ! j is nonpolar



c                     Warning, here we need to differenciate for residue i
c                     If residue i is non polar, aka hydrophobic, then the 
c                     formula involves Nw as well!

C                     If residue i is anything else, then we revert back to the
c                     old traditional formula

                      if(is_sim(i).eq.3) then     ! IF RESIDUE (i) is NON POLAR, aka Hydrophobic

c                          u_temp = c_np*(  1.0d0 - exp(  -a_np_u_water*
c     $                     (  Nw_water(res_type0)-Nw_ring(i)  )/
c     $                        Nw_water(res_type0)  )  )
                          goto 30
                           

                      else                        ! IF RESIDUE (i) is anything else 
c                                                   (that is polar/ion or membrane)
                          if(r1.le.r_np) then
                              goto 30
                          else
                              if(t_np_ngb.le.max_nonpolar) then

                                  u_temp=c_np*exp(-a_np_u*(t_np_ngb-max_nonpolar)**2)

                                  df=u_temp*(-a_np_u*2.0d0*(t_np_ngb-max_nonpolar))
     $                            *exp(-a_np*(r1-r_np)**2)*2.0d0*(-a_np*(r1-r_np))/r1
                              else
                                  goto 30
                              endif
                          endif
                      endif
c                                                            !i is ion or nonpolar or polar (ANYTHING)
                  else if(is_sim(j).eq.4) then               !j is membrane
c                 if below and df are not checked 

                     if(r1.le.r_mem_ngb_cutoff) then
                        goto 30
                     else
                        if(t_mem_ngb.le.max_mem_ngb) then

c                          Try to calculate the part "cote", if it exists
C                          This part is an EXACT COPY of the "Mem_energy" calculation,
c                          found later in the code

                           cote = 1.0d0     ! if there is no ld_grid, then keep this as 1
c                          --------------------------------------------------------------
c                          which would be a multiplying factor for the final temp_u
c                          ONLY WHEN AN LD GRID EXISTS!!!!! or ionic grid exist
  
                           if(ld_mem_r.ne.0.0d0.and.ld_mem_space.ne.0.0d0) then
                               cote = exp(-((dis_close_ld(i,1)-r_min_cote1)/r_min_cote2)**2)
                           elseif (l_run_mc_pt) then !ionic grid is switched on
                               !formula below from PNAS Ryshkova, Vicatos, Warshel 2010
                               if (dis_close_ig(i,1).ge.r_min_cote1) then
                                  cote = 1.0d0
                        else
c                                 formula approved by Warshel 
                                  cote = exp(-((dis_close_ig(i,1)-r_min_cote1)/r_min_cote2)**2)
                        endif
                     endif

                           u_temp=c_mem*cote*exp(-a_mem_u*(t_mem_ngb-max_mem_ngb)**2)

                           df=u_temp*(-2.0d0*a_mem_u*(t_mem_ngb-max_mem_ngb))
     $                       *exp(-a_np*(r1-r_mem_ngb_cutoff)**2)*2*(-a_np*(r1-r_mem_ngb_cutoff))/r1
                        else
                           goto 30
                        endif

                  endif    

                  endif    
c
c              End of the Calculation of u_temp and df
c              -------------------------------------------------------------------------------------  

c                 GLY thing. Need more explanation from Dr Chu
                  if(hit_gly.eq.1) df=df*formwgt

C                 Calculate forces!!!!!!
c                 ##############################################################################
                  d(i3+1)=d(i3+1)+df*dx1
                  d(i3+2)=d(i3+2)+df*dx2
                  d(i3+3)=d(i3+3)+df*dx3

                  d(j3+1)=d(j3+1)-df*dx1
                  d(j3+2)=d(j3+2)-df*dx2
                  d(j3+3)=d(j3+3)-df*dx3
c                 ##############################################################################



c                  if(check_num.eq.1) write(6,*) 'self j forces:',
c     $                             jhit,-df*dx1,-df*dx2,-df*dx3




              endif
c             END of irun = 2  
C             --------------------------------------------------------------------------------------
C             --------------------------------------------------------------------------------------
C             --------------------------------------------------------------------------------------
C             --------------------------------------------------------------------------------------


 30         continue  !end of the loop over all polar or nonpolar, or ion or mem residues

c           #########################################################################################
c           #########################################################################################
c           End of j loop OVER ALL RESIDUES!
c           Residue (i) 
c           #########################################################################################




            if(irun.eq.1) then  ! once we know how many neighbours we have let's compute U 
c                                 for a given residue i as well as all sums
c           ###############

c              First we save the neightbors into the global arrays
               Npol_CB(i)    =  t_p_ngb
               Nnonpol_CB(i) =  t_np_ngb
               Nmem_CB(i)    =  t_mem_ngb




c              If residue i is in a STRONGLY NON POLAR (aka Hydrophobic) environment
c              ------------------------------------------------------------------
c              ------------------------------------------------------------------
c              ------------------------------------------------------------------
               if(is_sim(i).eq.3) then   ! If it is hydrophobic, use a different formula
                 

c                 temp_u1 = c_np*(  1.0d0 - exp(  -a_np_u_water*(  Nw_water(res_type0)-Nw_ring(i)  )/Nw_water(res_type0)  ))
                 temp_u1 = c_np* exp(-a_np_u_water*Nw_ring(i)/Nw_water(res_type0)  )

c                temp_u1 = temp_u1 -0.30d0*min_wat_distance(i)

                  if(t_np_ngb.gt.0.01d0.and.t_np_ngb.lt.max_nonpolar) then
                     temp_u2=c_np*exp(-a_np_u*(t_np_ngb-max_nonpolar)**2)
                  else if(t_np_ngb.ge.max_nonpolar)  then
                     temp_u2=c_np
                  else 
                     temp_u2=0.0d0
                  endif

c                  One of many formulas for the free energy
c                  temp_u = fact_u1*temp_u1 + fact_u2*temp_u2

c                  Current formula
                  temp_u = fact_u1*temp_u1 
 
c                  write(6,*)"TEMPU start is ",temp_u
c                  write(6,*)"u1 and u2 are ",temp_u1, temp_u2



C                 additional term, used as a test for Dr Warshel hypothesis

c                  temp_u = 1.0d0*temp_u*(1.0d0-exp(-3.0d0*abs(temp_u)))
c                  temp_u = temp_u*(exp(-1.0d0*abs(temp_u)))

c                  write(6,*)"TEMPU before gyration is ",temp_u
                  temp_u = temp_u*(1.0d0-exp(-3.0d0*(min_wat_distance(i)/7.0d0)))

                  temp_u1 = temp_u


c     $               (  1.0d0 - exp(  -a_np_u_water*(  Nw_water(res_type0)-Nw_ring(i)  )/Nw_water(res_type0)  ))


c                  write(6,'(''TEMPU AFTER gyration is '',f8.2,2x,f8.2,2x,f8.2)')
c     $                temp_u, 1.0d0-exp(-3.0d0*(1.0d0-r3/gyration)-1.6), r3
c                  write(6,'(''coefficient is '',i5,2x,f8.3)')
c     $                i,  (  1.0d0 - exp(  -a_np_u_water*(  Nw_water(res_type0)-Nw_ring(i)  )/Nw_water(res_type0)  ))

c                  write(6,*)"Center "
c                  write(6,*)sys_xcen(1), sys_xcen(2), sys_xcen(3)


               else    !  Else use the same formula as with the rest 

                  if(t_np_ngb.gt.0.01d0.and.t_np_ngb.lt.max_nonpolar) then
                     temp_u=c_np*exp(-a_np_u*(t_np_ngb-max_nonpolar)**2)
                  else if(t_np_ngb.ge.max_nonpolar)  then
                     temp_u=c_np
                  else 
                     temp_u=0.0d0
                  endif
               endif


               temp_u_np=temp_u

               coefficient_np(current_i_type) =   coefficient_np(current_i_type) +
     $         temp_u 


c              ------------------------------------------------------------------
               if(is_sim(i).eq.1) then   ! if residue (i) is ion then....
                   e_self_np=e_self_np+temp_u
                   t_sim_ngb=t_sim_ngb+temp_u
                   if(hit_gly.eq.1) t_sim_ngb_gly=t_sim_ngb_gly+temp_u
                   if(debug) write(6,'('' ion-npol, total_N_np  E  :'',i4,f16.3)')
     $                   t_np_ngb,temp_u
c'
                  temp_energy = temp_energy + temp_u
                  temp_nonpolar = temp_nonpolar +temp_u      ! The non polar contribution
                  pder_np(par_derivative_type) =  pder_np(par_derivative_type) + temp_u/c_np   

c              ------------------------------------------------------------------
               elseif(is_sim(i).eq.3) then                  ! else if residue (i) is nonpolar (a.k.a. hydrophobic)   
                  e_hydro = e_hydro + temp_u
                  e_hydro_np = e_hydro_np + temp_u
                  temp_energy = temp_energy + temp_u
                  temp_nonpolar = temp_nonpolar + temp_u
                  sum_t_np_hydro = sum_t_np_hydro +temp_u
c              ------------------------------------------------------------------

               elseif(is_sim(i).eq.2) then                  ! Last option, residue(i) is polar
                  e_polar = e_polar + temp_u
                  e_polar_np = e_polar_np + temp_u
                  temp_energy = temp_energy + temp_u 
                  temp_nonpolar = temp_nonpolar +temp_u     ! The non polar contribution
                  sum_t_np_polar = sum_t_np_polar +temp_u

              endif








c              If residue i is in a strongly MEMBRANE environment
c              ------------------------------------------------------------------
c              ------------------------------------------------------------------
c              ------------------------------------------------------------------




               if(t_mem_ngb.gt.0.01d0.and.t_mem_ngb.le.max_mem_ngb) then
                  temp_u=c_mem*exp(-a_mem_u*(t_mem_ngb-max_mem_ngb)**2)
               else if (t_mem_ngb.gt.max_mem_ngb) then
                  temp_u=c_mem
c               else if (t_mem_ngb.eq.0.0d0) then
               else
                  temp_u=0.0d0
               endif


               cote = 1.0d0     ! if there is no ld_grid, then keep this as 1
c              --------------------------------------------------------------
c              which would be a multiplying factor for the final temp_u
c              ONLY WHEN AN LD GRID EXISTS!!!!! or ionic grid exist
               if(ld_mem_r.ne.0.0d0.and.ld_mem_space.ne.0.0d0) then

                     ! in this if branch the scling is done
                    ! for ionizable, polar and nonpolar residues
                    ! considered as incorrect by A. Warshel
                    cote = exp(-((dis_close_ld(i,1)-18.0d0)/12)**2)
               elseif (l_run_mc_pt) then !ionic grid is switched on

                   ! introduce scaling only for participated residues 
                   if (res_to_participated(i).ne.0) then   !res is actually participated


                       distanceToWater = rWat(res_to_participated(i)) 
                       !formula below from PNAS Ryshkova, Vicatos, Warshel 2010
                       !but params maybe different


                       if (distanceToWater.ge.18.0) then
                            cote = 1.0d0

                       else

                            cote = exp(-((distanceToWater-18.0d0)/9)**2)
                                   ! formula approved by Warshel 
                                   ! 18 is universal no need for any adjustment
                       endif
c                       write (*,*) ' i=', i,' res = ', seq_sim1(i_sim_type(i)),
c     $                             ' r = ' ,  distanceToWater,'cote = ', cote  
                   endif
               endif
               temp_u=temp_u*cote
               temp_u_mem=temp_u




c              ------------------------------------------------------------------
               if(is_sim(i).eq.1) then                         ! if residue (i) is ion    

                  e_self_mem=e_self_mem+temp_u   
                  t_sim_ngb=t_sim_ngb+temp_u

                  temp_energy =  temp_energy + temp_u
                  temp_mem = temp_mem + temp_u

c                 write statement
                  if(debug) write(6,'('' ion-mem, total_N_mem, total_N_mem_clbr  E  :'',2i4,2f16.3)')
     $                   t_mem_ngb,temp_u,e_self_np


c'              ------------------------------------------------------------------
               elseif(is_sim(i).eq.3) then     !                ! If residue (i) is nonpolar (a.k.a. hydrophobic)   
                  e_hydro = e_hydro + temp_u
                  e_hydro_mem = e_hydro_mem + temp_u
                  temp_energy =  temp_energy + temp_u
                  sum_t_mem_hydro = sum_t_mem_hydro +temp_u
                  temp_mem = temp_mem +temp_u
c              ------------------------------------------------------------------
               elseif(is_sim(i).eq.2) then                      ! Last option, residue(i) is polar
                  e_polar = e_polar + temp_u
                  e_polar_mem = e_polar_mem + temp_u
                  temp_energy = temp_energy + temp_u 
                  sum_t_mem_polar = sum_t_mem_polar +temp_u
                  temp_mem = temp_mem + temp_u
               endif



c              If residue i is in a STRONGLY POLAR (OR ION) environment
c              ------------------------------------------------------------------
c              ------------------------------------------------------------------
c              ------------------------------------------------------------------

               if(t_p_ngb.gt.0.01d0.and.t_p_ngb.le.max_polar) then
                  temp_u=c_polar*exp(-a_p_u*(t_p_ngb-max_polar)**2)
               else if (t_p_ngb.gt.max_polar) then
                  temp_u=c_polar
               else
                  temp_u=0.0d0
               endif


                  coefficient_p(current_i_type) =   coefficient_p(current_i_type) +
     $            temp_u 


c             ------------------------------------------------------------------
              if(is_sim(i).eq.1) then                      ! if residue (i) is ion 

                   e_self_p=e_self_p+temp_u
                   t_sim_ngb=t_sim_ngb+temp_u

                   if(hit_gly.eq.1) t_sim_ngb_gly=t_sim_ngb_gly+temp_u

                   if(debug) write(6,'('' ion-pol,  total_N_p   E  :'',i4,
     $                   f16.3,/)')t_p_ngb,temp_u

                   res_self(i)=e_self_np+e_self_p+e_self_mem    ! Calculate the self energy

                   temp_energy =  temp_energy + temp_u
                   temp_polar = temp_polar + temp_u
                   pder_p(par_derivative_type) =  pder_p(par_derivative_type) + temp_u/c_polar   


                   if(debug.and.istep.eq.nsteps) then
                    write(6,'('' ion_res:'',i5,'', N_np:'',i4,'', E_np:'',f7.3,
     $                      '', N_p:'',i4,'', E_p:'',f8.3,'', E_tot:'',f8.3)')
     $                      i,t_np_ngb,e_self_np,t_p_ngb,e_self_p,
     $                      e_self_np+e_self_p
c'               
                   endif 

c             ------------------------------------------------------------------
              elseif(is_sim(i).eq.3) then                       ! if it is non-polar (a.k.a Hydrophobic)
                   e_hydro = e_hydro + temp_u
                   e_hydro_p = e_hydro_p + temp_u
                   temp_energy =  temp_energy + temp_u
                   sum_t_p_hydro = sum_t_p_hydro +temp_u
                   temp_polar = temp_polar + temp_u
                   coeff_polar = coeff_polar + exp(-a_p_u*(t_p_ngb-max_polar)**2)

c             ------------------------------------------------------------------
              elseif(is_sim(i).eq.2) then                                           ! Last option, residue(i) is polar
                  e_polar = e_polar + temp_u
                  e_polar_p = e_polar_p + temp_u
                  temp_energy = temp_energy + temp_u
                  sum_t_p_polar = sum_t_p_polar + temp_u
                  temp_polar = temp_polar + temp_u 


              endif

            endif     ! END of if i_run = 1
c                                 
c           ###############




 50      continue !over irun=1,2




c            ###########################################################
c            PRINT OUTS!!!!!!!!!!
c            ############################################################


c                 Trying to find out what contributions are important            
c                 for certain residues. Residue i is the important one

         if(print_option) then

                  do iw=1,20
                     if(seq_sim1(i_sim_type(i)).eq.to_amino(iw)) seq_tmp=from_amino(iw)
                  enddo    


C           .... and print the results
            if(is_sim(i).ne.0) then

               if(is_sim(i).ne.3) then
                   call print_self_energy(is_sim(i),i,seq_tmp, c_polar, c_np, c_mem,
     $             temp_polar, temp_nonpolar, temp_mem, min_wat_distance(i),temp_energy,
     $             t_p_ngb,t_np_ngb,t_mem_ngb)
               else

               temp_general = float(Nw_ring(i))

               call print_hydro_energy(is_sim(i),i,seq_tmp, c_polar, c_np, c_mem,
     $          temp_polar, temp_u2, temp_u1, temp_mem,min_wat_distance(i), temp_energy,
     $          t_p_ngb, t_np_ngb,temp_general,t_mem_ngb, fact_u1, fact_u2)


c     subroutine print_hydro_energy    (option,   i,seq_name, c_polar, c_nonpolar, c_mem,
c     $ u_polar, u_nonpolar, u_ring, u_mem, min_wat_d, u_total,
c       n_p, n_np, n_ring, n_mem,factor_u1,factor_u2)



C              Save the values of Np and N_np into bins!

c              Polar bins

               bin = int(t_p_ngb) + 1
               if(bin.gt.10) bin = 10
               bins_polar(current_i_type,bin) = bins_polar(current_i_type,bin) + 1


c              NON POLAR (Hydrophobic) bins

               bin = int(t_np_ngb) + 1
               if(bin.gt.20) bin = 20
               bins_nonpolar(current_i_type,bin) = bins_nonpolar(current_i_type,bin) + 1

               endif            
            endif

         endif    
c            ###########################################################
c            END OF PRINT OUTS
c            ############################################################


 40   continue ! end of the loop over residue (i)



C         END OP THE LOOP to calculate the energy cotribution for EACH residue
Cc
C
c#########################################################################################
c#########################################################################################
c#########################################################################################
c#########################################################################################

c                   MAIN CHAIN AND HYDROGEN BOND SOLVATION!!!!

c#########################################################################################
c#########################################################################################
c#########################################################################################       

c         ! lets' calcualte main chain solvatation and H_bond term

         call  number_of_Calpha_neighbours ! for each C-alpha calulate 
                                             ! how many nonpol and mem neighbour it has
         call  calculate_Ualpha_array ! U_alpha array is created
         E_main_CG = MainChainSolvation()
         E_hb_CG = HydrogenBond_CG()
         if(print_option) then
            write(6,'('' E_main_solv = '',f10.3)'),E_main_CG
            write(6,'('' E_hb = '',f10.3)'), E_hb_CG
         endif

C     AFTERMATH, all calculations are over!!!!!
C#########################################################################################


c     PRINT ALL the necessary parameters into the official log file of molaris
      if(istep.eq.0) then

          get = 2
C     
          do i=1,5
c         call the subroutine to find res_type0
              call find_res_type(get,seq_sim1(i_sim_type(i)), i, res_type0)

c              write(6,'(3x,a3,5x,f8.3,5x,f8.3,8x,f8.3)'),
c     $        from_amino(res_type0), c_pol(i), c_nonpol(i), c_membrane(i)

          enddo

C
          do i=6,13
c             call the subroutine to find res_type0
              call find_res_type(get,seq_sim1(i_sim_type(i)), i, res_type0)
          enddo

          do i=14,19
c             call the subroutine to find res_type0
              call find_res_type(get,seq_sim1(i_sim_type(i)), i, res_type0)
          enddo

c         SOME EXTRA WRITES, before the closing of files 987, 988 and 989
c

          if(print_option) then
              option = 4 
              call print_self_energy(option,i,seq_tmp, c_polar, c_np, c_mem,
     $        e_polar, e_hydro, temp_mem, min_wat_distance(i), temp_energy,t_p_ngb,t_np_ngb,t_mem_ngb)

              call print_hydro_energy(option,i,seq_tmp, c_polar, c_np, c_mem,
     $             temp_polar, temp_u1, temp_u2, temp_mem, min_wat_distance(i), e_hydro,
     $             t_p_ngb, t_np_ngb,temp_general,t_mem_ngb,fact_u1,fact_u2)


C             Print some more stuff in the hydro_energy table
              write(989,*)
              write(989,*)
              write(989,*)" Statistics for the hydrophobic residues"
              write(989,*)"-----------------------------------------------------------------"
              write(989,*)
              write(989,*)" Total number of hydrophobic residues is ",total_np_res
              write(989,*)
              write(989,*)" Total polar contribution is ",e_hydro_p, coeff_polar
              write(989,*)
              write(989,'(''ALA '',i5,2x,''%hydro '',f5.2,2x,''%total '',f5.2,2x,''polar '',f8.3,'' hydro '',f8.3)') 
     $          n_ala, float(n_ala)/float(total_np_res)*100.0d0,
     $          float(n_ala)/float(numres)*100.0d0, coefficient_p(6), coefficient_np(6) 
              write(989,*)  

              if(n_ala.ne.0) then

              do iw =1,10
                 bins_polar(6,iw) = bins_polar(6,iw)/dfloat(n_ala)
              enddo

              do iw =1,20
                 bins_nonpolar(6,iw) = bins_nonpolar(6,iw)/dfloat(n_ala)
              enddo
              write(989,'(''         ALANINE  ALA  '',1x,f7.2,4x,f7.2,4x,f7.2)'),c_pol(6),c_nonpol(6),c_membrane(6)
              write(989,*)
              write(989,*)"  POLAR Bins       HYDROPHOBIC BINS          "
              write(989,*)"-------------------------------------------------------- "

              do iw = 1,10
                  write(989,'(i2,1x,f7.2,7x,f7.2)') iw, bins_polar(6,iw), bins_nonpolar(6,iw)
              enddo
              do iw = 11,20
                  write(989,'(i2,15x,f7.2)') iw, bins_nonpolar(6,iw)
              enddo
              write(989,*)
              write(989,*)
              endif
   


              write(989,'(''LEU '',i5,2x,''%hydro '',f5.2,2x,''%total '',f5.2,2x,''polar '',f8.3,'' hydro '',f8.3)') 
     $          n_leu, float(n_leu)/float(total_np_res)*100.0d0,
     $          float(n_leu)/float(numres)*100.0d0, coefficient_p(7), coefficient_np(7)  
              write(989,*)   

              if(n_leu.ne.0) then

              do iw =1,10
                 bins_polar(7,iw) = bins_polar(7,iw)/dfloat(n_leu)
              enddo

              do iw =1,20
                 bins_nonpolar(7,iw) = bins_nonpolar(7,iw)/dfloat(n_leu)
              enddo

              write(989,'(''         LEUCINE  LEU  '',1x,f7.2,4x,f7.2,4x,f7.2)'),c_pol(7),c_nonpol(7),c_membrane(7)
              write(989,*)
              write(989,*)"  POLAR Bins       HYDROPHOBIC BINS          "
              write(989,*)"-------------------------------------------------------- "

              do iw = 1,10
                  write(989,'(i2,1x,f7.2,7x,f7.2)') iw, bins_polar(7,iw), bins_nonpolar(7,iw)
              enddo
              do iw = 11,20
                  write(989,'(i2,15x,f7.2)') iw, bins_nonpolar(7,iw)
              enddo
              write(989,*)
              write(989,*)
              endif
  
              write(989,'(''ILE '',i5,2x,''%hydro '',f5.2,2x,''%total '',f5.2,2x,''polar '',f8.3,'' hydro '',f8.3)') 
     $          n_ile, float(n_ile)/float(total_np_res)*100.0d0,
     $          float(n_ile)/float(numres)*100.0d0, coefficient_p(8), coefficient_np(8) 
              write(989,*)   

              if(n_ile.ne.0) then

              do iw =1,10
                 bins_polar(8,iw) = bins_polar(8,iw)/dfloat(n_ile)
              enddo

              do iw =1,20
                 bins_nonpolar(8,iw) = bins_nonpolar(8,iw)/dfloat(n_ile)
              enddo

              write(989,'(''         ISOLEUCINE  ILE  '',1x,f7.2,4x,f7.2,4x,f7.2)'),c_pol(8),c_nonpol(8),c_membrane(8)
              write(989,*)
              write(989,*)"  POLAR Bins       HYDROPHOBIC BINS          "
              write(989,*)"-------------------------------------------------------- "

              do iw = 1,10
                  write(989,'(i2,1x,f7.2,7x,f7.2)') iw, bins_polar(8,iw), bins_nonpolar(8,iw)
              enddo
              do iw = 11,20
                  write(989,'(i2,15x,f7.2)') iw, bins_nonpolar(8,iw)
              enddo
              write(989,*)
              write(989,*)
              endif


              write(989,'(''VAL '',i5,2x,''%hydro '',f5.2,2x,''%total '',f5.2,2x,''polar '',f8.3,'' hydro '',f8.3') 
     $          n_val, float(n_val)/float(total_np_res)*100.0d0,
     $          float(n_val)/float(numres)*100.0d0, coefficient_p(9), coefficient_np(9)   
              write(989,*) 
  
             if(n_val.ne.0) then

              do iw =1,10
                 bins_polar(9,iw) = bins_polar(9,iw)/dfloat(n_val)
              enddo

              do iw =1,20
                 bins_nonpolar(9,iw) = bins_nonpolar(9,iw)/dfloat(n_val)
              enddo

              write(989,'(''         VALINE  VAL  '',1x,f7.2,4x,f7.2,4x,f7.2)'),c_pol(9),c_nonpol(9),c_membrane(9)
              write(989,*)
              write(989,*)"  POLAR Bins       HYDROPHOBIC BINS          "
              write(989,*)"-------------------------------------------------------- "

              do iw = 1,10
                  write(989,'(i2,1x,f7.2,7x,f7.2)') iw, bins_polar(9,iw), bins_nonpolar(9,iw)
              enddo
              do iw = 11,20
                  write(989,'(i2,15x,f7.2)') iw, bins_nonpolar(9,iw)
              enddo
              write(989,*)
              write(989,*)
              endif



              write(989,'(''PRO '',i5,2x,''%hydro '',f5.2,2x,''%total '',f5.2,2x,''polar '',f8.3,'' hydro '',f8.3') 
     $          n_prol, float(n_prol)/float(total_np_res)*100.0d0,
     $          float(n_prol)/float(numres)*100.0d0, coefficient_p(10), coefficient_np(10) 
              write(989,*)


            if(n_prol.ne.0) then

              do iw =1,10
                 bins_polar(10,iw) = bins_polar(10,iw)/dfloat(n_prol)
              enddo

              do iw =1,20
                 bins_nonpolar(10,iw) = bins_nonpolar(10,iw)/dfloat(n_prol)
              enddo

              write(989,'(''         PROLINE  PRO  '',1x,f7.2,4x,f7.2,4x,f7.2)'),c_pol(10),c_nonpol(10),c_membrane(10)
              write(989,*)
              write(989,*)"  POLAR Bins       HYDROPHOBIC BINS          "
              write(989,*)"-------------------------------------------------------- "

              do iw = 1,10
                  write(989,'(i2,1x,f7.2,7x,f7.2)') iw, bins_polar(10,iw), bins_nonpolar(10,iw)
              enddo
              do iw = 11,20
                  write(989,'(i2,15x,f7.2)') iw, bins_nonpolar(10,iw)
              enddo

              write(989,*)
              write(989,*)
              endif

   
              write(989,'(''MET '',i5,2x,''%hydro '',f5.2,2x,''%total '',f5.2,2x,''polar '',f8.3,'' hydro '',f8.3)') 
     $          n_met, float(n_met)/float(total_np_res)*100.0d0,
     $          float(n_met)/float(numres)*100.0d0, coefficient_p(11), coefficient_np(11) 
              write(989,*)   


            if(n_met.ne.0) then

              do iw =1,10
                 bins_polar(11,iw) = bins_polar(11,iw)/dfloat(n_met)
              enddo

              do iw =1,20
                 bins_nonpolar(11,iw) = bins_nonpolar(11,iw)/dfloat(n_met)
              enddo

              write(989,'(''         METHIONINE  MET  '',1x,f7.2,4x,f7.2,4x,f7.2)'),c_pol(11),c_nonpol(11),c_membrane(11)
              write(989,*)
              write(989,*)"  POLAR Bins       HYDROPHOBIC BINS          "
              write(989,*)"-------------------------------------------------------- "

              do iw = 1,10
                  write(989,'(i2,1x,f7.2,7x,f7.2)') iw, bins_polar(11,iw), bins_nonpolar(11,iw)
              enddo
              do iw = 11,20
                  write(989,'(i2,15x,f7.2)') iw, bins_nonpolar(11,iw)
              enddo

              write(989,*)
              write(989,*)
              endif

              write(989,'(''PHE '',i5,2x,''%hydro '',f5.2,2x,''%total '',f5.2,2x,''polar '',f8.3,'' hydro '',f8.3)') 
     $          n_phe, float(n_phe)/float(total_np_res)*100.0d0,
     $          float(n_phe)/float(numres)*100.0d0, coefficient_p(12), coefficient_np(12) 
              write(989,*)
   

             if(n_phe.ne.0) then

              do iw =1,10
                 bins_polar(12,iw) = bins_polar(12,iw)/dfloat(n_phe)
              enddo

              do iw =1,20
                 bins_nonpolar(12,iw) = bins_nonpolar(12,iw)/dfloat(n_phe)
              enddo

              write(989,'(''      PHENYLALANINE  PHE  '',1x,f7.2,4x,f7.2,4x,f7.2)'),c_pol(12),c_nonpol(12),c_membrane(12)
              write(989,*)
              write(989,*)"  POLAR Bins       HYDROPHOBIC BINS          "
              write(989,*)"-------------------------------------------------------- "

              do iw = 1,10
                  write(989,'(i2,1x,f7.2,7x,f7.2)') iw, bins_polar(12,iw), bins_nonpolar(12,iw)
              enddo
              do iw = 11,20
                  write(989,'(i2,15x,f7.2)') iw, bins_nonpolar(12,iw)
              enddo

              write(989,*)
              write(989,*)
              endif


              write(989,'(''TRP '',i5,2x,''%hydro '',f5.2,2x,''%total '',f5.2,2x,''polar '',f8.3,'' hydro '',f8.3)') 
     $          n_trp, float(n_trp)/float(total_np_res)*100.0d0,
     $          float(n_trp)/float(numres)*100.0d0, coefficient_p(13), coefficient_np(13) 
              write(989,*)   


           if(n_trp.ne.0) then

              do iw =1,10
                 bins_polar(13,iw) = bins_polar(13,iw)/dfloat(n_trp)
              enddo

              do iw =1,20
                 bins_nonpolar(13,iw) = bins_nonpolar(13,iw)/dfloat(n_trp)
              enddo

              write(989,'(''         TRYPTOPHANE  TRP  '',1x,f7.2,4x,f7.2,4x,f7.2)'),c_pol(13),c_nonpol(13),c_membrane(13)
              write(989,*)
              write(989,*)"  POLAR Bins       HYDROPHOBIC BINS          "
              write(989,*)"-------------------------------------------------------- "

              do iw = 1,10
                  write(989,'(i2,1x,f7.2,7x,f7.2)') iw, bins_polar(13,iw), bins_nonpolar(13,iw)
              enddo
              do iw = 11,20
                  write(989,'(i2,15x,f7.2)') iw, bins_nonpolar(13,iw)
              enddo

              write(989,*)
              write(989,*)
              endif

              write(989,'(''ALA '',i5,2x,''%hydro '',f5.2,2x,''%total '',f5.2,2x,''polar '',f8.3,'' hydro '',f8.3)') 
     $          n_ala, float(n_ala)/float(total_np_res)*100.0d0,
     $          float(n_ala)/float(numres)*100.0d0, coefficient_p(6), coefficient_np(6) 
              write(989,'(''LEU '',i5,2x,''%hydro '',f5.2,2x,''%total '',f5.2,2x,''polar '',f8.3,'' hydro '',f8.3)') 
     $          n_leu, float(n_leu)/float(total_np_res)*100.0d0,
     $          float(n_leu)/float(numres)*100.0d0, coefficient_p(7), coefficient_np(7)  
             write(989,'(''ILE '',i5,2x,''%hydro '',f5.2,2x,''%total '',f5.2,2x,''polar '',f8.3,'' hydro '',f8.3)') 
     $          n_ile, float(n_ile)/float(total_np_res)*100.0d0,
     $          float(n_ile)/float(numres)*100.0d0, coefficient_p(8), coefficient_np(8)  
             write(989,'(''VAL '',i5,2x,''%hydro '',f5.2,2x,''%total '',f5.2,2x,''polar '',f8.3,'' hydro '',f8.3') 
     $          n_val, float(n_val)/float(total_np_res)*100.0d0,
     $          float(n_val)/float(numres)*100.0d0, coefficient_p(9), coefficient_np(9)   
              write(989,'(''PRO '',i5,2x,''%hydro '',f5.2,2x,''%total '',f5.2,2x,''polar '',f8.3,'' hydro '',f8.3') 
     $          n_prol, float(n_prol)/float(total_np_res)*100.0d0,
     $          float(n_prol)/float(numres)*100.0d0, coefficient_p(10), coefficient_np(10) 
              write(989,'(''MET '',i5,2x,''%hydro '',f5.2,2x,''%total '',f5.2,2x,''polar '',f8.3,'' hydro '',f8.3)') 
     $          n_met, float(n_met)/float(total_np_res)*100.0d0,
     $          float(n_met)/float(numres)*100.0d0, coefficient_p(11), coefficient_np(11) 
              write(989,'(''PHE '',i5,2x,''%hydro '',f5.2,2x,''%total '',f5.2,2x,''polar '',f8.3,'' hydro '',f8.3)') 
     $          n_phe, float(n_phe)/float(total_np_res)*100.0d0,
     $          float(n_phe)/float(numres)*100.0d0, coefficient_p(12), coefficient_np(12) 
              write(989,'(''TRP '',i5,2x,''%hydro '',f5.2,2x,''%total '',f5.2,2x,''polar '',f8.3,'' hydro '',f8.3)') 
     $          n_trp, float(n_trp)/float(total_np_res)*100.0d0,
     $          float(n_trp)/float(numres)*100.0d0, coefficient_p(13), coefficient_np(13) 
              write(989,*)   

              write(989,*)"HYD RESIDUES"   
              write(989,'(i5,2x,i5,2x,i5,2x,i5,2x,i5,2x,i5,2x,i5,2x,i5)')
     $        n_ala, n_leu, n_ile, n_val, n_prol, n_met, n_phe, n_trp 

             write(989,*)   
             write(989,*)"----------------------------------------------"
c             write(989,'(''   prints ALA+LEU+VAL '',2x,f8.3)')
c     $         coefficient_np(6)+coefficient_np(7)+coefficient_np(9) 
c
c             write(989,'(''   prints ILE+PRO+MET '',2x,f8.3)')
c     $         coefficient_np(8)+coefficient_np(10)+coefficient_np(11) 
c
c            write(989,'(''   prints PHE+TRP '',2x,f8.3)')
c     $         coefficient_np(12)+coefficient_np(13) 
c
              close(989)

C             Print some more stuff in the polar table
              write(988,*)
              write(988,*)
              write(988,*)" Statistics for the polar residues"
              write(988,*)"-----------------------------------------------------------------"
              write(988,*)
              write(988,*)" Total number of polar residues is ",n_ser+n_thr+n_tyr+n_cys+n_asn+n_gln
              write(988,*)
              write(988,'(''SER '',i5,2x,''%polar '',f5.2,2x,''%total '',f5.2,2x,''polar '',f8.3,'' hydro '',f8.3)') 
     $          n_ser, float(n_ser)/float(n_ser+n_thr+n_tyr+n_cys+n_asn+n_gln)*100.0d0,
     $          float(n_ser)/float(numres)*100.0d0, coefficient_p(14), coefficient_np(14) 
              write(988,*)  


              write(988,'(''         SERINE   SER  '',1x,f7.2,4x,f7.2,4x,f7.2)'),c_pol(14),c_nonpol(14),c_membrane(14)
              write(988,*)

              write(988,*)
              write(988,*)


              write(988,'(''THR '',i5,2x,''%polar '',f5.2,2x,''%total '',f5.2,2x,''polar '',f8.3,'' hydro '',f8.3)') 
     $          n_thr, float(n_thr)/float(n_ser+n_thr+n_tyr+n_cys+n_asn+n_gln)*100.0d0,
     $          float(n_thr)/float(numres)*100.0d0, coefficient_p(15), coefficient_np(15) 
              write(988,*)  


              write(988,'(''        THREONINE THR  '',1x,f7.2,4x,f7.2,4x,f7.2)'),c_pol(15),c_nonpol(15),c_membrane(15)
              write(988,*)

              write(988,*)
              write(988,*)              


              write(988,'(''TYR '',i5,2x,''%polar '',f5.2,2x,''%total '',f5.2,2x,''polar '',f8.3,'' hydro '',f8.3)') 
     $          n_tyr, float(n_tyr)/float(n_ser+n_thr+n_tyr+n_cys+n_asn+n_gln)*100.0d0,
     $          float(n_tyr)/float(numres)*100.0d0, coefficient_p(16), coefficient_np(16) 
              write(988,*)  


              write(988,'(''        TYROSINE  TYR  '',1x,f7.2,4x,f7.2,4x,f7.2)'),c_pol(16),c_nonpol(16),c_membrane(16)
              write(988,*)

              write(988,*)
              write(988,*)              


              write(988,'(''CYS '',i5,2x,''%polar '',f5.2,2x,''%total '',f5.2,2x,''polar '',f8.3,'' hydro '',f8.3)') 
     $          n_cys, float(n_cys)/float(n_ser+n_thr+n_tyr+n_cys+n_asn+n_gln)*100.0d0,
     $          float(n_cys)/float(numres)*100.0d0, coefficient_p(17), coefficient_np(17) 
              write(988,*)  

              write(988,'(''        CYSTEINE  CYS  '',1x,f7.2,4x,f7.2,4x,f7.2)'),c_pol(17),c_nonpol(17),c_membrane(17)
              write(988,*)

              write(988,*)
              write(988,*)              


              write(988,'(''ASN '',i5,2x,''%polar '',f5.2,2x,''%total '',f5.2,2x,''polar '',f8.3,'' hydro '',f8.3)') 
     $          n_asn, float(n_asn)/float(n_ser+n_thr+n_tyr+n_cys+n_asn+n_gln)*100.0d0,
     $          float(n_asn)/float(numres)*100.0d0, coefficient_p(18), coefficient_np(18) 
              write(988,*)  

              write(988,'(''       ASPARAGINE  ASN  '',1x,f7.2,4x,f7.2,4x,f7.2)'),c_pol(18),c_nonpol(18),c_membrane(18)
              write(988,*)

              write(988,*)
              write(988,*)        

              write(988,'(''GLN '',i5,2x,''%polar '',f5.2,2x,''%total '',f5.2,2x,''polar '',f8.3,'' hydro '',f8.3)') 
     $          n_gln, float(n_gln)/float(n_ser+n_thr+n_tyr+n_cys+n_asn+n_gln)*100.0d0,
     $          float(n_gln)/float(numres)*100.0d0, coefficient_p(19), coefficient_np(19) 
              write(988,*)  

              write(988,'(''       GLUTAMINE  GLN  '',1x,f7.2,4x,f7.2,4x,f7.2)'),c_pol(19),c_nonpol(19),c_membrane(19)
              write(988,*)

              write(988,*)
              write(988,*)   

              close(988)





          endif

      endif ! End of the if print is true





c     Variables necessary for the studies of Molecular Dynamics with CG
      sum_e_hydro = sum_e_hydro + e_hydro
      n_occur_hydro = n_occur_hydro +1        ! counter for the hydro energy
c                                               which starts from the beginning of 
c                                               the MD simulation!



c...........................................................................
      return
      end




c                        EEEE   N   N   DD
c                        E      NN  N   D D
c                        EEEE   N N N   D  D
c                        E      N  NN   D D
c                        EEEE   N   N   DD
c                
c###########################################################################################
c###########################################################################################
c###########################################################################################
c###########################################################################################
c###########################################################################################
c'"



      subroutine old_self_energy

c
      implicit none

      include 'include/parameter.dc'
      include 'include/coord.h'
      include 'include/seq.h'
      include 'include/cbkv_gap.h'
      include 'include/Mdyn.h'
      include 'include/task.h'
      include 'include/derivative.h'
      include 'include/symbols.h'
      include 'include/stepsize.h'
      include 'include/lib.h'
      include 'include/cg_ion_hyd_pol_constants.h'
      include 'include/dgm_pt.h'
      include 'include/rWatProtein.h'
      include 'include/filename.h'

c     local var:>
      integer i,j,ii,jj,k,m,ihit,jhit,i3,j3,irun,iw,
     $        hit,hit_gly,option,
     $        icount,par_derivative_type, res_type,res_type0, get
      real*8 dx1,dx2,dx3,r1,r2,df,u_temp,temp_u,e_self_np,e_self_p,e_self_mem,DM,
     $       temp_u_np,temp_u_mem

      real*8 c_polar, c_np, c_mem,temp_energy
      real*8 max_polar, max_nonpolar,  max_mem, r_min_cote1, r_min_cote2

      real*8 c_pol(19), c_nonpol(19), c_membrane(19)

      real*8 temp_nonpolar, temp_polar, temp_mem, temp_general
 
      real*8 sum_t_np_hydro, sum_t_p_hydro, sum_t_mem_hydro
      real*8 sum_t_np_polar, sum_t_p_polar, sum_t_mem_polar

      integer Nw(numres), Nw_ring(numres), Nw_water(20)

c     Langevin grid variable (part are obsolete)
      real*8 r_ld_min0, alpha_min,x1_1,y1,y2, cote

c     all variables describing neighbours are real now (so we can have 1.67 neighbours)
c     (required for correct usage of F(rij) = exp (-6(r-rnp)^2))
      real*8 t_mem_ngb, t_p_ngb, t_np_ngb

c     F(rij) values (see Rychkova, Vicatos, Warshel PNAS)
      real*8 mem_sim_ngb, p_sim_ngb, np_sim_ngb
      
c     cutoff radius for mem neighbours
c     will be assigned to mem_spacing*f (f=2.05)
c      real*8 r_mem_ngb_cutoff

c     number of neighbours for ion in the membrane
c     should not strongly depend on membrane spacing
c     mem in two solvation shells are included
c     within r<mem_spacing*f
      real*8 max_mem_ngb


c    closest distance to water from participated residue
      real*8 distanceToWater 

      character seq_tmp*3, type*5
      logical debug, print_option
 
c      real*8 E_main_CG  ! main chain solvatation penalty
c      real*8 E_hydro_CG ! Hbond in CG model
      real*8 MainChainSolvation  ! function for calculation
      real*8 HydrogenBond_CG  ! function for calculation
c......................................................................


          write(6,*)
          write(6,*)          

          write(6,*)"#########################################################################"
          write(6,*)          
          write(6,*)" THIS IS A RUN WITH THE OLD MODEL PARAMETERS AND EQUATION!"
          write(6,*)          
          write(6,*)"#########################################################################"
 
          write(6,*)

      get = 1

      icount=icount+1


      debug=.true.
      debug=.false.

c     Printout fq is determined by log_write_fq
      print_option=.false.
      if (mod(istep,log_write_fq).eq.0.or.istep.eq.nsteps) then
         print_option=.true.
      endif

      if(lcheck_self) debug=.true.

c     is_sim()=1   ion
c     is_sim()=2   polar
c     is_sim()=3   nonpolar
c     is_sim()=4   membrane
      

c     Keep this for development 
c      write(6,*)"Success!!!! simple neighbor is running"


      t_sim_ngb=0.0


C     Reset the constants for all residues
      do i=1,19
          c_pol(i)=0.0d0
          c_nonpol(i) = 0.0d0
          c_membrane(i) = 0.0d0
      enddo


C     Assign the cg parameters!






C     Parameters which are used to calculate N_p N_np and N_mem

      r_p =  5.0d0
      r_np = 7.0d0

      a_p= 6.0d0
      a_np = 6.0d0
      a_mem = 6.0d0

C     Parameters which are used to calculate the energy contributions

      max_polar = 4.0d0
      max_nonpolar = 6.0d0
      max_mem_ngb = 27.0d0  

      a_p_u = 0.2d0
      a_np_u = 0.2d0
      a_mem_u = 0.2d0

C     These ones are to calculate cote, which is used for the 
c     energy contributions of the membrane
      r_min_cote1 = 18.0d0
      r_min_cote2 = 12.0d0


C     DO the ionizable residues first 
c     Those are the most elaborate ones:

C     Aspartic Acid, ASP:
      c_pol(1)      = -1.3d0
      c_nonpol(1)   = 3.1d0
      c_membrane(1) = 15.0d0

C     Glutamic Acid GLU:
      c_pol(2)      = -1.30d0
      c_nonpol(2)   = 5.1d0
      c_membrane(2) = 15.0d0

C     Lysine LYS:
      c_pol(3)      = -1.0d0
      c_nonpol(3)   = 5.4d0
      c_membrane(3) = 15.0d0

C     Arginine ARG:
      c_pol(4)      = -1.0d0
      c_nonpol(4)   = 4.3d0
      c_membrane(4) = 15.0d0

C     Histidine HIS:
      c_pol(5)      = -1.0d0
      c_nonpol(5)   = 2.3d0
      c_membrane(5) = 15.0d0


C     Start assigning "intermediate constants for the hydrophobic residues)


C     All values taken from run which were in November 2011     
c
c     Alanine  A*2 8
      c_pol(6)      = 0.5d0*1.0d0
      c_nonpol(6)   = -1.22d0*3.52d0/24.0d0
      c_membrane(6) = 15.0d0

c     Leucine  L*2 

      c_pol(7)      = 0.5d0*1.0d0
      c_nonpol(7)   = -1.22d0*5.61d0/24.0d0
      c_membrane(7) = 15.0d0

c     Isoleucine I*2

      c_pol(8)      = 0.5d0*1.0d0
      c_nonpol(8)   = -1.22d0*5.62d0/24.0d0
      c_membrane(8) = 15.0d0

c     Valine V*2

      c_pol(9)      = 0.5d0*1.0d0
      c_nonpol(9)   = -1.22d0*4.95d0/24.0d0
      c_membrane(9) = 15.0d0

c     Proline P*2

      c_pol(10)      = 0.5d0*1.0d0
      c_nonpol(10)   = -1.22d0*4.85d0/24.0d0
      c_membrane(10) = 15.0d0

c     Methionine M*2

      c_pol(11)      = 0.5d0*1.0d0
      c_nonpol(11)   = -1.22d0*3.62d0/24.0d0
      c_membrane(11) = 15.0d0

c     Phenylalanine F*2

      c_pol(12)      = 0.5d0*1.0d0
      c_nonpol(12)   = -1.22d0*5.67d0/24.0d0
      c_membrane(12) = 15.0d0

c     Tryptophane  W*2: f_p = 1  , f_np = 5.90/4 = 1.48

      c_pol(13)      = 0.5d0*1.0d0
      c_nonpol(13)   = -1.22d0*5.9d0/24.0d0
      c_membrane(13) = 15.0d0
c ##########################
      

C    JMB 2001 312, 927-934

c     Serine  S*2

      c_pol(14)      = -0.2d0
      c_nonpol(14)   = 0.2d0
      c_membrane(14) = 15.0d0

c     Threonine T*2

      c_pol(15)      = -0.2d0
      c_nonpol(15)   = 0.2d0
      c_membrane(15) = 15.0d0

c     Tyrosine Y*1

      c_pol(16)      = -0.2d0
      c_nonpol(16)   = 0.2d0
      c_membrane(16) = 15.0d0

c     Cysteine C*2

      c_pol(17)      = -0.2d0
      c_nonpol(17)   = 0.2d0
      c_membrane(17) = 15.0d0
c     Asparagine N*2

      c_pol(18)      = -0.2d0
      c_nonpol(18)   = 0.2d0
      c_membrane(18) = 15.0d0

c     Glutamine  Q*2

      c_pol(19)      = -0.2d0
      c_nonpol(19)   = 0.2d0
      c_membrane(19) = 15.0d0


c
c###################################################################################
c###################################################################################
c###################################################################################


c     Constant for the LD grid calculations, IF the LD is used! (obsolete)
      alpha_min = 0.50d0


c     Initializing 
      e_hydro = 0.0d0
      e_hydro_np = 0.0d0
      e_hydro_mem = 0.0d0
      e_hydro_p = 0.0d0

      e_polar = 0.0d0
      e_polar_np = 0.0d0
      e_polar_mem = 0.0d0
      e_polar_p   = 0.0d0


c     This is not very necessary, but it was important to check
c     these values with the e_hydro etc values
C     Eventualy they will be removed
      sum_t_np_hydro = 0.0d0
      sum_t_p_hydro = 0.0d0
      sum_t_mem_hydro = 0.0d0

      sum_t_np_polar = 0.0d0
      sum_t_p_polar = 0.0d0
      sum_t_mem_polar = 0.0d0

c     Initialize and set the N-neighbors to zero
      do i =1,numres
          Npol_CB(i) = 0.0d0 !  Polar neighbors of residue i
          Nnonpol_CB(i) = 0.0d0 ! NON  Polar (aka hydrophobic) neighbors of residue i
          Nmem_CB(i) = 0.0d0 !  Membrane neighbors of residue i
      enddo

c     Set the partial derivatives of the "function" f to zero
      do i=1,5
         pder_p(i) = 0.0d0
         pder_np(i) = 0.0d0
      enddo

C     Reset the residue type!
      par_derivative_type = 0
      
      if (has_membrane) then
          DM = mem_spacing 
          r_mem_ngb_cutoff = DM*2.05 ! only neighbours in the "two solvation shells" counted as 1
                                    ! at further distance exp(...) function is used
c          write(6,'('' membrane spacing =''f8.3)') DM
          if ( (DM.le.0.5).or.(DM.ge.5.0) ) then
              write (*,*) 
     $             'WARNING: DM - membrane distance is unphysical.'
              write (*,*) 'WARNING: DM = ', DM
              write (*,*) 'Program stopped in cg_energetics.f'
              stop
          endif
      endif


c     Write the "label"  for the Table!
c     ################

      if(print_option) then
          option = 0
          call print_self_energy(option,i,seq_tmp, c_polar, c_np, c_mem,
     $     temp_polar, temp_nonpolar, temp_mem, min_wat_distance(i), temp_energy,
     $     t_p_ngb,t_np_ngb,t_mem_ngb)

      endif



      total_p_res = 0
      total_np_res = 0
      total_ion_res = 0

c#########################################################################################
c#########################################################################################
c     STARTING OP THE LOOP to calculate the energy cotribution for EACH residue
C     IF residue (i) is ionizable, the contribution would be the self energy (i)
C     IF residue (i) is hydrophobic, the contribution would the hydrophobic contribution
C     IF residue (i) is polar, the contribution would the polar contribution
      do 40 i=1,numres ! loop over all ions and nonpolar residues 

      if(is_sim(i).eq.2) total_p_res = total_p_res +1
      if(is_sim(i).eq.3) total_np_res = total_np_res +1
      if(is_sim(i).eq.1) total_ion_res = total_ion_res +1


          if(is_sim(i).eq.4) goto 40   !go to the next residue  if current is MEB
 
              res_type = 0
C         Call the subroutine which simply "recognises" the residue
c         and gives the appropriate value of the "type"
          call find_res_type(get,seq_sim1(i_sim_type(i)), res_type,res_type0)

c         Set the values for residue i (based on the "type" of i)
          c_polar = c_pol(res_type)
          c_np = c_nonpol(res_type)
          c_mem = c_membrane(res_type)


c         SET important variables to zero!
c         ###################################################
          t_mem_ngb=0.0  !number of membrane neighbours
          t_p_ngb=0.0    ! .....    polar ....
          t_np_ngb=0.0   ! ....     nonpolar ...
          e_self_mem=0.0  !contib to self_e from membrane
          e_self_np=0.0  !contib to self_e from nonpolar
          e_self_p=0.0  !contib to self_e from polar residues
C         First, set the temp_energy term for residue i to zero
          temp_energy = 0.0d0 !total sel_e fo a given (i-th) residue
          temp_nonpolar = 0.0d0
          temp_polar =0.0d0
          temp_mem = 0.0d0

          ihit=0
          do ii=np(i)+1,np(i+1)
             if(iac_name(iac(ii)).eq.'DY') goto 40
             if(bkvcode(ii)(1:2).eq.'CB') then
                 ihit=ii
             endif
          enddo
          if(ihit.eq.0) goto 40


C
C        #################      Loop for i = 1 and i = 2
c        during irun = 1 number of neighbors are calculated and then Self energy for a given residue
c        during irun = 2 forces on atoms are calculated
C        #########################################################
         do 50 irun=1,2




c           #########################################################################################
c           #########################################################################################
c           LOOP OVER ALL RESIDUES!
c           Residue (i) will "interact with each residue j, and all contributions will be calculated
c           #########################################################################################

             do 30 j=1,numres                      !loop over all polar or nonpolar, or ion or mem

                 if(i.eq.j) goto 30
                 jhit=0

                 do jj=np(j)+1,np(j+1)
                     if(iac_name(iac(jj)).eq.'DY') goto 30
                     if(bkvcode(jj)(1:2).eq.'CB') then
                         jhit=jj
                     endif
                 enddo

                 if(jhit.eq.0) goto 30

c                inefficient part : distances are calculated two times, for irun = 1,2
C                Eventually this needs to be fixed.
                 i3   = ihit*3-3
                 j3   = jhit*3-3
                 dx1  = x(i3+1)-x(j3+1)
                 dx2  = x(i3+2)-x(j3+2)
                 dx3  = x(i3+3)-x(j3+3)
                 r2   = dx1*dx1+dx2*dx2+dx3*dx3
                 r1   = sqrt(r2)            


C                HIT CLYCINE, need more explanation from Dr. Chu!
                 hit_gly=0
                 if(npoly.gt.0) then
                     do m=1,npoly
                         if(ihit.eq.ipoly(m).or.jhit.eq.ipoly(m)) hit_gly=1
                     enddo
                 else
                     do m=1,n_mut_gly
                         if(ihit.eq.i_mut_gly(m).or.jhit.eq.i_mut_gly(m))then
                             hit_gly=1
                         endif
                     enddo
                 endif

C 
c                irun = 1  
C                 --------------------------------------------------------------------------------------
                 if(irun.eq.1)  then ! get total Nnp and Np and Nmem (for future calculation of energy)

                                                             !i is ion or nonpolar or polar
                     if(is_sim(j).eq.2.or.is_sim(j).eq.1) then  !j is polar or ion
                         if(r1.le.r_p) then
                             p_sim_ngb=1.0
                         else
                             p_sim_ngb=exp(-a_p*(r1-r_p)**2)      !need to write the equation in the comment
                         endif
 
                         t_p_ngb=t_p_ngb+p_sim_ngb
 
                         if(debug) write(6,'('' ion-pol, residues  r  N_p:'',2i4,
     $                         f12.3,i5)')i,j,r1,p_sim_ngb

c                                                            !i is ion or nonpolar or polar
                     else if(is_sim(j).eq.3) then               !j is nonpolar
                         if(r1.le.r_np) then
                             np_sim_ngb=1.0
                         else
                             np_sim_ngb=exp(-a_np*(r1-r_np)**2)  !need to write the equation in the comment
c                             np_sim_ngb= 0.0d0  !need to write the equation in the comment
                         endif

                         t_np_ngb=t_np_ngb+np_sim_ngb
                         if(debug) write(6,'('' ion-npol,residues  r  N_np'',2i4,
     $                         f12.3,i5)')i,j,r1,np_sim_ngb

c                                                            !i is ion or nonpolar or polar
                     else if(is_sim(j).eq.4) then               !j is membrane
                         ! if below is formula 11 in Rychkova et al PNAS but rcutoff is new
                         ! mem_sim_ngb is F(rij) in article's notation
                         if(r1.le.r_mem_ngb_cutoff) then
                              mem_sim_ngb=1.0
                         else
                              mem_sim_ngb=exp(-a_mem*(r1-r_mem_ngb_cutoff)**2) !need to write the equation in the comment
                         endif
  
                         ! t_mem_ngb is N(i)mem in article's notation
                         ! t_mem_ngb meaning - number of mem residues around given residue

                         t_mem_ngb=t_mem_ngb+mem_sim_ngb
 
                         if(debug) write(6,'('' ion-membrane,residues  r  N_np'',
     $                         2i4,f12.3,i5)')i,j,r1,mem_sim_ngb

                     endif     


c                irun = 2  
c                after getting total Nnp (t_np_ngb)  and Np (t_p_ngb), 
c                which are functions of r_ij, calculating forces for atom  j
C                --------------------------------------------------------------------------------------
C                --------------------------------------------------------------------------------------
C                --------------------------------------------------------------------------------------

                 else if (irun.eq.2)   then

c
c                    Calculate u_temp and eventually df
c                    -------------------------------------------------------------------------------------  


c                    ###############    2   ###################    ! i is ion or nonpolar or polar (ANYTHING)
                     if(is_sim(j).eq.2.or.is_sim(j).eq.1) then    ! if j is ion or polar
                         if(r1.le.r_p) then
                             goto 30
                         else

                             if(t_p_ngb.le.max_polar) then
                                  u_temp= c_polar*exp(-a_p_u*(t_p_ngb-max_polar)**2)
                                  df=u_temp*(-a_p_u*2.0d0*(t_p_ngb-max_polar))
     $                            *exp(-a_p*(r1-r_p)**2)*2.0d0*(-a_p*(r1-r_p))/r1
                             else
                                 goto 30
                             endif
                         endif
                              



c                    ###############    3   ###################    ! i is ion or nonpolar or polar (ANYTHING)
                     else if(is_sim(j).eq.3) then               ! j is nonpolar


c                        Warning, here we need to differenciate for residue i
c                        If residue i is non polar, aka hydrophobic, then the 
c                        formula involves Nw as well!

C                        If residue i is anything else, then we revert back to the
c                        old traditional formula

                         if(is_sim(i).eq.3) then     ! IF RESIDUE (i) is NON POLAR, aka Hydrophobic


                             if(r1.le.r_np) then
                                 goto 30
                             else
                                 if(t_np_ngb.le.max_nonpolar) then

                                     u_temp=c_np*exp(-a_np_u*(t_np_ngb-max_nonpolar)**2)

                                     df=u_temp*(-a_np_u*2.0d0*(t_np_ngb-max_nonpolar))
     $                                *exp(-a_np*(r1-r_np)**2)*2.0d0*(-a_np*(r1-r_np))/r1
                                 else
                                     goto 30
                                 endif
                             endif

                         endif




c                        
c                    ###############    4   ###################   !i is ion or nonpolar or polar (ANYTHING)
                     else if(is_sim(j).eq.4) then               !j is membrane


                         if(r1.le.r_mem_ngb_cutoff) then
                             goto 30
                         else
                             if(t_mem_ngb.le.max_mem_ngb) then

c                                Try to calculate the part "cote", if it exists
C                                This part is an EXACT COPY of the "Mem_energy" calculation,
c                                found later in the code

                                 cote = 1.0d0     ! if there is no ld_grid, then keep this as 1
c                                 --------------------------------------------------------------
c                                 which would be a multiplying factor for the final temp_u
c                                 ONLY WHEN AN LD GRID EXISTS!!!!! or ionic grid exist
  
                                  if(ld_mem_r.ne.0.0d0.and.ld_mem_space.ne.0.0d0) then
                                      cote = exp(-((dis_close_ld(i,1)-r_min_cote1)/r_min_cote2)**2)
                                  elseif (l_run_mc_pt) then !ionic grid is switched on
                                        !formula below from PNAS Ryshkova, Vicatos, Warshel 2010
                                      if (dis_close_ig(i,1).ge.r_min_cote1) then
                                           cote = 1.0d0
                                      else
c                                          formula approved by Warshel 
                                           cote = exp(-((dis_close_ig(i,1)-r_min_cote1)/r_min_cote2)**2)
                                      endif
                                  endif

                                  u_temp=c_mem*cote*exp(-a_mem_u*(t_mem_ngb-max_mem_ngb)**2)

                                  df=u_temp*(-2.0d0*a_mem_u*(t_mem_ngb-max_mem_ngb))
     $                             *exp(-a_np*(r1-r_mem_ngb_cutoff)**2)*2*(-a_np*(r1-r_mem_ngb_cutoff))/r1
                             else
                                 goto 30
                             endif

                         endif    

                     endif    
c
c                    End of the Calculation of u_temp and df
c                    ---------------------------------------------------------------------------  

c                    GLY thing. Need more explanation from Dr Chu
                     if(hit_gly.eq.1) df=df*formwgt

C                        Calculate forces!!!!!!
c                        ################################################################
                         d(i3+1)=d(i3+1)+df*dx1
                         d(i3+2)=d(i3+2)+df*dx2
                         d(i3+3)=d(i3+3)+df*dx3

                         d(j3+1)=d(j3+1)-df*dx1
                         d(j3+2)=d(j3+2)-df*dx2
                         d(j3+3)=d(j3+3)-df*dx3
c                        ###############################################################




                 endif      ! End of if irun = 2

C                ---------------------------------------------------------------------
C                ---------------------------------------------------------------------
C                ---------------------------------------------------------------------
C                ---------------------------------------------------------------------


 30          continue  !end of the loop do j =1 over all polar or nonpolar, or ion or mem residues





c            #########################################################################################
c            #########################################################################################
c            #########################################################################################
c                      CALCULATE U FOR A GIVEN RESIDUE
c            ########################################################################################
c            ########################################################################################
c            ########################################################################################



             if(irun.eq.1) then  ! once we know how many neighbours we have let's compute U 
c                                 for a given residue i as well as all sums
c                ###############

c                First we save the neightbors into the global arrays

                 Npol_CB(i)    =  t_p_ngb
                 Nnonpol_CB(i) =  t_np_ngb
                 Nmem_CB(i)    =  t_mem_ngb


c                If residue i is in a STRONGLY NON POLAR (aka Hydrophobic) environment
c                ------------------------------------------------------------------
c                ------------------------------------------------------------------
c                ------------------------------------------------------------------
                 if(is_sim(i).eq.3) then   ! If it is hydrophobic, use a different formula
 
                     if(t_np_ngb.gt.0.01d0.and.t_np_ngb.lt.max_nonpolar) then
                         temp_u=c_np*exp(-a_np_u*(t_np_ngb-max_nonpolar)**2)
                     else if(t_np_ngb.ge.max_nonpolar)  then
                         temp_u=c_np
                     else 
                         temp_u=0.0d0
                     endif
                 endif

                 temp_u_np=temp_u

c                ------------------------------------------------------------------
                 if(is_sim(i).eq.1) then   ! if residue (i) is ion then....
                     e_self_np=e_self_np+temp_u
                     t_sim_ngb=t_sim_ngb+temp_u
                     if(hit_gly.eq.1) t_sim_ngb_gly=t_sim_ngb_gly+temp_u
                     if(debug) write(6,'('' ion-npol, total_N_np  E  :'',i4,f16.3)')
     $                   t_np_ngb,temp_u
c'
                     temp_energy = temp_energy + temp_u
                     temp_nonpolar = temp_nonpolar +temp_u      ! The non polar contribution
                     pder_np(par_derivative_type) =  pder_np(par_derivative_type) + temp_u/c_np   

c                ------------------------------------------------------------------
                 elseif(is_sim(i).eq.3) then                  ! else if residue (i) is nonpolar (a.k.a. hydrophobic)   
                     e_hydro = e_hydro + temp_u
                     e_hydro_np = e_hydro_np + temp_u
                     temp_energy = temp_energy + temp_u
                     temp_nonpolar = temp_nonpolar + temp_u
                     sum_t_np_hydro = sum_t_np_hydro +temp_u
c              ------------------------------------------------------------------

                 elseif(is_sim(i).eq.2) then                  ! Last option, residue(i) is polar
                      e_polar = e_polar + temp_u
                      e_polar_np = e_polar_np + temp_u
                      temp_energy = temp_energy + temp_u 
                      temp_nonpolar = temp_nonpolar +temp_u     ! The non polar contribution
                      sum_t_np_polar = sum_t_np_polar +temp_u

                 endif


c              If residue i is in a strongly MEMBRANE environment
c              ------------------------------------------------------------------
c              ------------------------------------------------------------------
c              ------------------------------------------------------------------
               if(t_mem_ngb.gt.0.01d0.and.t_mem_ngb.le.max_mem_ngb) then
                  temp_u=c_mem*exp(-a_mem_u*(t_mem_ngb-max_mem_ngb)**2)
               else if (t_mem_ngb.gt.max_mem_ngb) then
                  temp_u=c_mem
c               else if (t_mem_ngb.eq.0.0d0) then
               else
                  temp_u=0.0d0
               endif

               cote = 1.0d0     ! if there is no ld_grid, then keep this as 1
c              --------------------------------------------------------------
c              which would be a multiplying factor for the final temp_u
c              ONLY WHEN AN LD GRID EXISTS!!!!! or ionic grid exist
               if(ld_mem_r.ne.0.0d0.and.ld_mem_space.ne.0.0d0) then
                    ! in this if branch the scling is done
                    ! for ionizable, polar and nonpolar residues
                    ! considered as incorrect by A. Warshel
                    cote = exp(-((dis_close_ld(i,1)-18.0d0)/12)**2)
               elseif (l_run_mc_pt) then !ionic grid is switched on
                   ! introduce scaling only for participated residues 
                   if (res_to_participated(i).ne.0) then   !res is actually participated
                       distanceToWater = rWat(res_to_participated(i)) 
                       !formula below from PNAS Ryshkova, Vicatos, Warshel 2010
                       !but params maybe different
                       if (distanceToWater.ge.18.0) then
                            cote = 1.0d0
                       else
                            cote = exp(-((distanceToWater-18.0d0)/9)**2)
                                   ! formula approved by Warshel 
                                   ! 18 is universal no need for any adjustment
                       endif
                       write (*,*) ' i=', i,' res = ', seq_sim1(i_sim_type(i)),
     $                             ' r = ' ,  distanceToWater,'cote = ', cote  
                   endif
               endif
               temp_u=temp_u*cote
               temp_u_mem=temp_u

c              ------------------------------------------------------------------
               if(is_sim(i).eq.1) then                         ! if residue (i) is ion    

                  e_self_mem=e_self_mem+temp_u   
                  t_sim_ngb=t_sim_ngb+temp_u

                  temp_energy =  temp_energy + temp_u
                  temp_mem = temp_mem + temp_u

c                 write statement
                  if(debug) write(6,'('' ion-mem, total_N_mem, total_N_mem_clbr  E  :'',2i4,2f16.3)')
     $                   t_mem_ngb,temp_u,e_self_np


c'              ------------------------------------------------------------------
               elseif(is_sim(i).eq.3) then     !                ! If residue (i) is nonpolar (a.k.a. hydrophobic)   
                  e_hydro = e_hydro + temp_u
                  e_hydro_mem = e_hydro_mem + temp_u
                  temp_energy =  temp_energy + temp_u
                  sum_t_mem_hydro = sum_t_mem_hydro +temp_u
                  temp_mem = temp_mem +temp_u
c              ------------------------------------------------------------------
               elseif(is_sim(i).eq.2) then                      ! Last option, residue(i) is polar
                  e_polar = e_polar + temp_u
                  e_polar_mem = e_polar_mem + temp_u
                  temp_energy = temp_energy + temp_u 
                  sum_t_mem_polar = sum_t_mem_polar +temp_u
                  temp_mem = temp_mem + temp_u
               endif


c              If residue i is in a STRONGLY POLAR (OR ION) environment
c              ------------------------------------------------------------------
c              ------------------------------------------------------------------
c              ------------------------------------------------------------------

               if(t_p_ngb.gt.0.01d0.and.t_p_ngb.le.max_polar) then
                  temp_u=c_polar*exp(-a_p_u*(t_p_ngb-max_polar)**2)
               else if (t_p_ngb.gt.max_polar) then
                  temp_u=c_polar
               else
                  temp_u=0.0d0
               endif


c             ------------------------------------------------------------------
              if(is_sim(i).eq.1) then                      ! if residue (i) is ion 

                   e_self_p=e_self_p+temp_u
                   t_sim_ngb=t_sim_ngb+temp_u

                   if(hit_gly.eq.1) t_sim_ngb_gly=t_sim_ngb_gly+temp_u

                   if(debug) write(6,'('' ion-pol,  total_N_p   E  :'',i4,
     $                   f16.3,/)')t_p_ngb,temp_u

                   res_self(i)=e_self_np+e_self_p+e_self_mem    ! Calculate the self energy

                   temp_energy =  temp_energy + temp_u
                   temp_polar = temp_polar + temp_u
                   pder_p(par_derivative_type) =  pder_p(par_derivative_type) + temp_u/c_polar   


                   if(debug.and.istep.eq.nsteps) then
                    write(6,'('' ion_res:'',i5,'', N_np:'',i4,'', E_np:'',f7.3,
     $                      '', N_p:'',i4,'', E_p:'',f8.3,'', E_tot:'',f8.3)')
     $                      i,t_np_ngb,e_self_np,t_p_ngb,e_self_p,
     $                      e_self_np+e_self_p
c'               
                   endif 

c             ------------------------------------------------------------------
              elseif(is_sim(i).eq.3) then                       ! if it is non-polar (a.k.a Hydrophobic)
                   e_hydro = e_hydro + temp_u
                   e_hydro_p = e_hydro_p + temp_u
                   temp_energy =  temp_energy + temp_u
                   sum_t_p_hydro = sum_t_p_hydro +temp_u
                   temp_polar = temp_polar + temp_u

c             ------------------------------------------------------------------
              elseif(is_sim(i).eq.2) then                                           ! Last option, residue(i) is polar
                  e_polar = e_polar + temp_u
                  e_polar_p = e_polar_p + temp_u
                  temp_energy = temp_energy + temp_u
                  sum_t_p_polar = sum_t_p_polar + temp_u
                  temp_polar = temp_polar + temp_u 


              endif



            endif     ! END of if i_run = 1
c                                 
c           ###############




 50      continue !over irun=1,2




c            ###########################################################
c            PRINT OUTS!!!!!!!!!!
c            ############################################################


c                 Trying to find out what contributions are important            
c                 for certain residues. Residue i is the important one

         if(print_option) then

                  do iw=1,20
                     if(seq_sim1(i_sim_type(i)).eq.to_amino(iw)) seq_tmp=from_amino(iw)
                  enddo    


C           .... and print the results
            if(is_sim(i).ne.0) then

               if(is_sim(i).ne.3) then
                   call print_self_energy(is_sim(i),i,seq_tmp, c_polar, c_np, c_mem,
     $             temp_polar, temp_nonpolar, temp_mem, min_wat_distance(i), temp_energy,
     $             t_p_ngb,t_np_ngb,t_mem_ngb)
               else
               temp_general = float(Nw_ring(i))
                   call print_self_energy(is_sim(i),i,seq_tmp, c_polar, c_np, c_mem,
     $             temp_polar, temp_nonpolar, temp_mem, min_wat_distance(i), temp_energy,
     $             t_p_ngb,temp_general,t_mem_ngb)
               endif            
            endif

         endif    
c            ###########################################################
c            END OF PRINT OUTS
c            ############################################################


 40   continue ! end of the loop over residue (i)




c     PRINT ALL the necessary parameters into the official log file of molaris
      if(istep.eq.0) then

          write(6,*)          

          write(6,*)"#########################################################################"
          write(6,*)          
          write(6,*)" PARAMETERS WHICH WERE USED IN THE SIMPLIFIED MODEL"
          write(6,*)" Pamameters involve the self energy, hydrophobic energy, polar energy"
          write(6,*)" PARAMETERS WHICH WERE USED IN THE SIMPLIFIED MODEL"
          write(6,*)          

          write(6,*)"#########################################################################"
 
          write(6,*)"             Ionizable residues "
          write(6,*)"           -----------------------"
          write(6,*)"Residue     C_polar     C_non_polar     C_membrane"
          write(6,*)"---------------------------------------------------------"


          get = 2
C     
          do i=1,5
c         call the subroutine to find res_type0
              call find_res_type(get,seq_sim1(i_sim_type(i)), i, res_type0)

              write(6,'(3x,a3,5x,f8.3,5x,f8.3,8x,f8.3)'),
     $        from_amino(res_type0), c_pol(i), c_nonpol(i), c_membrane(i)

          enddo

          write(6,*)

          write(6,*)"             Non Polar (Hydrophobic) residues "
          write(6,*)"           --------------------------------------"
          write(6,*)"Residue     C_polar     C_non_polar     C_membrane"
          write(6,*)"---------------------------------------------------------"

C
          do i=6,13
c             call the subroutine to find res_type0
              call find_res_type(get,seq_sim1(i_sim_type(i)), i, res_type0)

              write(6,'(3x,a3,5x,f8.3,5x,f8.3,8x,f8.3)'),
     $        from_amino(res_type0), c_pol(i), c_nonpol(i), c_membrane(i)

          enddo

          write(6,*)
          write(6,*)"                  Polar residues "
          write(6,*)"                -------------------"
          write(6,*)"Residue     C_polar     C_non_polar     C_membrane"
          write(6,*)"---------------------------------------------------------"

C     
          do i=14,19
c             call the subroutine to find res_type0
              call find_res_type(get,seq_sim1(i_sim_type(i)), i, res_type0)

              write(6,'(3x,a3,5x,f8.3,5x,f8.3,8x,f8.3)'),
     $        from_amino(res_type0), c_pol(i), c_nonpol(i), c_membrane(i)

          enddo

          write(6,*)
          write(6,*)"Parameters for Neighbors Np, Nnp, Nmem equations "
          write(6,*)"-------------------------------------------------"
          write(6,*)
          write(6,*)"a_p = ",a_p
          write(6,*)"a_np = ",a_np
          write(6,*)"a_mem = ",a_mem
          write(6,*)
          write(6,*)"r_p = ",r_p
          write(6,*)"r_np = ",r_np
          write(6,*)"r_mem_ngb_cutoff = ",r_mem_ngb_cutoff
          write(6,*)"r_mem is dependent of Membrane Spacing, which is  ", DM
          write(6,*)

          write(6,*)"Parameters for Energy Contributions Up, Unp, Umem equations "
          write(6,*)"-------------------------------------------------"
          write(6,*)
          write(6,*)"a_p_u = ", a_p_u
          write(6,*)"a_np_u = ",a_np_u
          write(6,*)"a_np_u_water = ",a_np_u_water
          write(6,*)"a_mem_u = ",a_mem_u
          write(6,*)
          write(6,*)"max polar = ",max_polar
          write(6,*)"max non polar = ",max_nonpolar
          write(6,*)"max # of membrane atoms in two solvation shells = ",max_mem_ngb
          write(6,*)
          write(6,*)"Additional terms for Umem "
          write(6,*)"-------------------------------------------------"
          write(6,*)
          write(6,*)"r_min_cote1 = ",r_min_cote1
          write(6,*)"r_min_cote2 = ",r_min_cote2
          write(6,*)"div = ",div
          write(6,*)


          write(6,*)" END OF PARAMETERS PRINTS"
          write(6,*)          
          write(6,*)"#########################################################################"
          write(6,*)          



          if(print_option) then
              option = 4 
              call print_self_energy(option,i,seq_tmp, c_polar, c_np, c_mem,
     $        e_polar, e_hydro, temp_mem, min_wat_distance(1), temp_energy,t_p_ngb,t_np_ngb,t_mem_ngb)

          endif

      endif ! End of the if print is true

c     Variables necessary for the studies of Molecular Dynamics with CG
      sum_e_hydro = sum_e_hydro + e_hydro
      n_occur_hydro = n_occur_hydro +1        ! counter for the hydro energy
c                                               which starts from the beginning of 
c                                               the MD simulation!


      print *,"here -->>" 

      if(istep.eq.nsteps) print *
c...........................................................................
      return
      end




c        
c        
c         EEEE   N N N   D  D   of old_self_energy
c        
c        




c###########################################################################################
c###########################################################################################
c###########################################################################################


      subroutine simple_neighbor
c
c     calculate the interaction of ion with polar and nonpolar neighbors

      implicit none

      include 'include/parameter.dc'
      include 'include/coord.h'
      include 'include/seq.h'
      include 'include/cbkv_gap.h'
      include 'include/Mdyn.h'
      include 'include/derivative.h'
      include 'include/cg_ion_hyd_pol_constants.h'
      include 'include/symbols.h'

c     local var:>
      integer i,j,ii,jj,k,ihit,jhit,i3,j3
      real*8 dx1,dx2,dx3,r1,r2,df,e_sim_ngb,e_fp,y
c......................................................................
c     is_sim()=1   ion
c     is_sim()=2   polar
c     is_sim()=3   nonpolar
c     is_sim()=4   membrane

      t_sim_ngb=0.0
      t_sim_ngb_gly=0.0

      do i=1,numres                                     !i is ion
         if(is_sim(i).eq.2.or.is_sim(i).eq.3) goto 40
         ihit=0
         do ii=np(i)+1,np(i+1)
            if(iac_name(iac(ii)).eq.'DY') goto 40
            if(bkvcode(ii)(1:2).eq.'CB') then
               ihit=ii
            endif
         enddo
         if(ihit.eq.0) goto 40

         do j=1,numres                                  !j is polar or nonpolar
            if(is_sim(j).eq.1) goto 30
            jhit=0
            do jj=np(j)+1,np(j+1)
               if(iac_name(iac(jj)).eq.'DY') goto 30
               if(bkvcode(jj)(1:2).eq.'CB') then
                  jhit=jj
               endif
            enddo
            if(jhit.eq.0) goto 30

            i3   = ihit*3-3
            j3   = jhit*3-3
            dx1  = x(i3+1)-x(j3+1)
            dx2  = x(i3+2)-x(j3+2)
            dx3  = x(i3+3)-x(j3+3)
            r2   = dx1*dx1+dx2*dx2+dx3*dx3
            r1   = sqrt(r2)            

            if(is_sim(j).eq.2) then                       !polar
               if(r1.le.r_p) then
                  e_sim_ngb=1.0
               else
                  e_sim_ngb=exp(-a_p*(r1-r_p)**2)
               endif

               e_fp=exp(0.89588*(e_sim_ngb-0.20351))-1.20
               df=-2.0*a_p*e_fp*e_sim_ngb*(r1-r_p)/r1
c               write(6,'('' ion-po  res e_fp r1:'',2i4,f12.3)')i,j,e_fp,r1

            else if(is_sim(j).eq.3) then                  !nonpolar
               if(r1.le.r_np) then
                  e_sim_ngb=1.0
               else
                  e_sim_ngb=exp(-a_np*(r1-r_np)**2)
               endif

               y=atan(4.0*e_sim_ngb-16.0)
               e_fp=0.5+y/3.1416
               e_fp=a_npn*e_fp
               df=-8.0*a_np*e_sim_ngb*(r1-r_np)/(3.1416*r1*(tan(y)**2+1.0))
               df=a_npn*df
c               write(6,'('' ion-np  res e_fp r1:'',2i4,f12.3)')i,j,e_fp,r1
            endif

            t_sim_ngb=t_sim_ngb+e_fp

c           forces only on polar and nonpolar atoms, not on ions
            d(j3+1)=d(j3+1)-df*dx1
            d(j3+2)=d(j3+2)-df*dx2
            d(j3+3)=d(j3+3)-df*dx3

 30      enddo
 40   enddo
c...........................................................................
      return
      end
c----------------------------------------------------------------------







c----------------------------------------------------------------------
      subroutine expl_neighbor_self_e
c
c     calculate the interaction of ion with polar and nonpolar neighbors
c     in explicit forms

      implicit none

      include 'include/parameter.dc'
      include 'include/coord.h'
      include 'include/seq.h'
      include 'include/cbkv_gap.h'
      include 'include/Mdyn.h'
      include 'include/derivative.h'
      include 'include/symbols.h'
      include 'include/cg_ion_hyd_pol_constants.h'
      include 'include/stepsize.h'

c     local var:>
      integer i,j,k,irun,t_p_ngb,t_np_ngb,p_expl_ngb,np_expl_ngb,m,m3,i_cb,
     $        t_hit,n,n3,ni,nj,write_ion,write_np,mem_expl_ngb,t_mem_ngb
      real*8  dx1,dx2,dx3,r1,r2,df,u_temp,temp1_x(3),temp2_x(3),temp_u,
     $        e_self_np,e_self_p,e_self_mem
      logical debug
      data write_np,write_ion/0,0/
      save write_np,write_ion
c......................................................................
      debug=.true.
      debug=.false.

c     iexpl_polr()=1   ion
c     iexpl_polr()=2   polar
c     iexpl_polr()=3   nonpolar
c     iexpl_polr()=4   membrane

      t_expl_ngb=0.0
      a_p=6.0
      a_np=6.0

c     check if ion is defined
      t_hit=0
      do i=1,numres
         if(iexpl_polr(i).eq.1) t_hit=t_hit+1
      enddo
      if(t_hit.eq.0) then
         if(write_ion.eq.0) then
            write(6,'('' WARNING: there is no ion defined for expl_ion_self_e'',/,
     $             '' calculation'')')
            write_ion=1
         endif
         return
      endif

c      if(istep.eq.nsteps) write(6,'(/,'' explicit_ion_self_e for each ion:'')')

      do i=1,numres                                               !i is ion
         if(iexpl_polr(i).eq.0.or.iexpl_polr(i).eq.2.or.
     $         iexpl_polr(i).eq.3.or.iexpl_polr(i).eq.4) goto 40
         do irun=1,2

            if(irun.eq.1) then
               t_p_ngb=0
               t_np_ngb=0
               t_mem_ngb=0
               e_self_np=0.0
               e_self_mem=0.0
               e_self_p=0.0
               
            endif

            do k=1,3
               temp1_x(k)=0.0
            enddo
            ni=0
            do j=np(i)+5,np(i+1)-2
               if(iac_name(iac(j)).eq.'DY') goto 40
               ni=ni+1
               do k=1,3
                  temp1_x(k)=temp1_x(k)+x(3*(j-1)+k)
               enddo
            enddo

            do k=1,3
               temp1_x(k)=temp1_x(k)/ni              !get side center of i_ion
            enddo

            t_hit=0
            do j=1,numres !j is polar or nonpolar, or ion
               if(iexpl_polr(j).eq.0.or.i.eq.j) goto 30
               i_cb=5
               if(resseq(i).eq.'PRO') i_cb=4
               t_hit=t_hit+1

               do k=1,3
                  temp2_x(k)=0.0
               enddo
               nj=0
               do m=np(j)+i_cb,np(j+1)-2
                  if(iac_name(iac(m)).eq.'DY') goto 30
                  nj=nj+1
                  do k=1,3
                     temp2_x(k)=temp2_x(k)+x(3*(m-1)+k)  
                  enddo
               enddo

               do k=1,3
                  temp2_x(k)=temp2_x(k)/nj         !get side center of j_nonion
               enddo

               dx1  = temp1_x(1)-temp2_x(1)
               dx2  = temp1_x(2)-temp2_x(2)
               dx3  = temp1_x(3)-temp2_x(3)
               r2   = dx1*dx1+dx2*dx2+dx3*dx3
               r1   = sqrt(r2)            

               if(irun.eq.1) then     !get total Nnp and Np to calculate energy
                  if(iexpl_polr(j).eq.2.or.iexpl_polr(j).eq.1) then  !polar,ion
                     if(r1.le.r_p_expl) then
                        p_expl_ngb=1
                     else
                        p_expl_ngb=int(exp(-a_p*(r1-r_p)**2))
                     endif
                     t_p_ngb=t_p_ngb+p_expl_ngb

                     if(debug) then
                        write(6,'('' for ion '',i4,''_'',a3,'' residue '',i4,
     $                          ''_'',a3,''    polar r ='',f6.2,'' Np  = '',i3,
     $                          '' total_Np  = '',i3)')i,resseq(i),j,resseq(j),
     $                          r1,p_expl_ngb,t_p_ngb
                     endif

                  else if(iexpl_polr(j).eq.3) then               !j is nonpolar
                     if(r1.le.r_np_expl) then
                        np_expl_ngb=1
                     else
                        np_expl_ngb=int(exp(-a_np*(r1-r_np)**2))
                     endif
                     t_np_ngb=t_np_ngb+np_expl_ngb

                     if(debug) then
                        write(6,'('' for ion '',i4,''_'',a3,'' residue '',i4,
     $                          ''_'',a3,'' nonpolar r ='',f6.2,'' Nnp = '',i3,
     $                          '' total_Nnp = '',i3)')i,resseq(i),j,resseq(j),
     $                          r1,np_expl_ngb,t_np_ngb
                     endif

                  else if(iexpl_polr(j).eq.4) then               !j is membrane
                     if(r1.le.r_np_expl) then
                        mem_expl_ngb=1
                     else
                        mem_expl_ngb=int(exp(-a_np*(r1-r_np)**2))
                     endif
                     t_mem_ngb=t_mem_ngb+mem_expl_ngb

                     if(debug) then
                        write(6,'('' for ion '',i4,''_'',a3,'' residue '',i4,
     $                          ''_'',a3,'' nonpolar r ='',f6.2,'' Nnp = '',i3,
     $                          '' total_Nnp = '',i3)')i,resseq(i),j,resseq(j),
     $                          r1,np_expl_ngb,t_np_ngb
                     endif
                  endif

               else        
c                 after getting total Nnp (t_np_ngb)  and Np (t_p_ngb), 
c                 which are functions of r_ij, calculating forces for atom j

                  if(iexpl_polr(j).eq.2.or.iexpl_polr(j).eq.1) then  !polar,ion
                     if(r1.le.r_p_expl) then
                        goto 30
                     else
                        if(t_p_ngb.le.4) then
                           u_temp=-2*exp(-0.2*(t_p_ngb-4)**2)
                           df=u_temp*(-0.2*2*(t_p_ngb-4))
     $                       *exp(-a_np*(r1-r_p)**2)*2*(-a_np*(r1-r_p))/r1
                        else
                           goto 30
                        endif

                     endif

                  else if(iexpl_polr(j).eq.3) then               !j is nonpolar
                     if(r1.le.r_np_expl) then
                        goto 30
                     else
                        if(t_np_ngb.le.6) then
                           u_temp=4.0*exp(-0.2*(t_np_ngb-6)**2)
                           df=u_temp*(-0.2*2*(t_np_ngb-6))
     $                       *exp(-a_np*(r1-r_np)**2)*2*(-a_np*(r1-r_np))/r1
                        else
                           goto 30
                        endif
                     endif

                  else if(iexpl_polr(j).eq.4) then               !j is membrane
                     if(r1.le.r_np_expl) then
                        goto 30
                     else
                        if(t_np_ngb.le.6) then
                           u_temp=6.0*exp(-0.2*(t_mem_ngb-6)**2)
                           df=u_temp*(-0.2*2*(t_mem_ngb-6))
     $                       *exp(-a_np*(r1-r_np)**2)*2*(-a_np*(r1-r_np))/r1
                        else
                           goto 30
                        endif
                     endif


                  endif

c                 forces only on polar and nonpolar atoms, not on ions
                  do n=np(i)+5,np(i+1)-2
                     n3=3*(n-1)
                     d(n3+1)=d(n3+1)+df*dx1/ni
                     d(n3+2)=d(n3+2)+df*dx2/ni
                     d(n3+3)=d(n3+3)+df*dx3/ni
                  enddo

                  do m=np(j)+i_cb,np(j+1)-2
                     m3=3*(m-1)
                     d(m3+1)=d(m3+1)-df*dx1/nj
                     d(m3+2)=d(m3+2)-df*dx2/nj
                     d(m3+3)=d(m3+3)-df*dx3/nj
                     if(check_num.eq.1) write(6,*) 'self j forces:',
     $                          m,-df*dx1/nj,-df*dx2/nj,-df*dx3/nj
                  enddo

               endif
 30         enddo

            if(t_hit.eq.0) then
               if(write_np.eq.0) then
                  write(6,'('' WARNING: there is no polar or nonpolar residues'',
     $                      '' defined for expl_ion_self_e'',/,
     $                      ''          calculation'')')
                  write_np=1
               endif
               return
            endif

            if(irun.eq.1) then
c              irun=1, calculating the total U_self(t_expl_ngb) of U_np+U_polar
c              for nonpolar
               if(t_np_ngb.le.6) then                     !get U_np
                  temp_u=4.0*exp(-0.2*(t_np_ngb-6)**2)
               else
                  temp_u=4.0
               endif

               e_self_np=e_self_np+temp_u
               t_expl_ngb=t_expl_ngb+temp_u
               if(debug) write(6,'('' nonpolar self_e for residue :'',i4,
     $                                f8.3)')i,temp_u

c              for membrane
               if(t_mem_ngb.le.6) then                     !get U_mem
                  temp_u=6.0*exp(-0.2*(t_mem_ngb-6)**2)
               else
                  temp_u=6.0
               endif

               e_self_mem=e_self_mem+temp_u
               e_self_np=e_self_np+temp_u
               t_expl_ngb=t_expl_ngb+temp_u
               if(debug) write(6,'('' membrane self_e for residue :'',i4,
     $                                f8.3)')i,temp_u

c              for polar
               if(t_p_ngb.le.4) then                      !get U_polar
                  temp_u=-2.0*exp(-0.2*(t_p_ngb-4)**2)
               else
                  temp_u=-2.0
               endif

               if(debug) write(6,'('' polar self_e for residue    :'',i4,
     $                                f8.3)')i,temp_u

               t_expl_ngb=t_expl_ngb+temp_u
               e_self_p=e_self_p+temp_u
               res_self(i)=e_self_np+e_self_p
               if(debug) write(6,'('' total self_e for residue    :'',i4,
     $                                f8.3)')i, res_self(i)

               if(debug.and.istep.eq.nsteps) then
                  write(6,'('' for ion_res:'',i5,8x,''  self_E_tot:'',f8.3)')
     $                      i, t_expl_ngb
               endif
            endif
         enddo
 40   enddo
      if(istep.eq.nsteps) print *
c...........................................................................
      return
      end
c----------------------------------------------------------------------




      subroutine assign_cg_constants(c_pol,c_nonpol,c_membrane)

      implicit none
      include 'include/parameter.dc'
      include 'include/cg_ion_hyd_pol_constants.h'

      real*8 c_pol(19), c_nonpol(19), c_membrane(19)
      real*8 hyd_fctr_p(8), hyd_fctr_np(8), hyd_fctr_mem(8), pol_fctr_p(6), pol_fctr_np(6)
      real*8 r_p_def, r_np_def

      integer i

      r_p_def =  5.0d0
      r_np_def = 7.0d0

      a_p=0.10d0
      a_np = a_p
      a_mem = 6.0d0

      a_p_u = 0.1d0
      a_np_u = 0.02d0
      a_np_u_water = 1.4d0
      a_mem_u = 0.005d0

      if(r_p.eq.0.0d0.or.r_np.eq.0.0d0) Then  ! r_p, r_np are not specified, 
c                                             !then we use the default values
          r_p = r_p_def          ! default is 5 (from above)
          r_np = r_np_def        ! default is 7 (from above)
      else                     ! save the new default values of rp and rnp
          r_p_def =  r_p
          r_np_def = r_np        
      endif

C     DO the ionizable residues first 
c     Those are the most elaborate ones:

C     Aspartic Acid, ASP:
      c_pol(1)      = cg_b_polar_def(1)
      c_nonpol(1)   = cg_b_nonpolar_def(1)
      c_membrane(1) = cg_b_membrane_def(1)

C     Glutamic Acid GLU:
      c_pol(2)      = cg_b_polar_def(2)
      c_nonpol(2)   = cg_b_nonpolar_def(2)
      c_membrane(2) = cg_b_membrane_def(2)

C     Lysine LYS:
      c_pol(3)      = cg_b_polar_def(3)
      c_nonpol(3)   = cg_b_nonpolar_def(3)
      c_membrane(3) = cg_b_membrane_def(3)

C     Arginine ARG:
      c_pol(4)      = cg_b_polar_def(4)
      c_nonpol(4)   = cg_b_nonpolar_def(4)
      c_membrane(4) = cg_b_membrane_def(4)

C     Histidine HIS:
      c_pol(5)      = cg_b_polar_def(5)
      c_nonpol(5)   = cg_b_nonpolar_def(5)
      c_membrane(5) = cg_b_membrane_def(5)


C     Start assigning "intermediate constants for the hydrophobic residues)

c     Alanine  A*2 8
c      hyd_fctr_p(1) = 1.0d0
c      hyd_fctr_np(1) = 3.52/div     ! 3.52
      hyd_fctr_mem(1) = -0.69

c     Leucine  L*2 
c      hyd_fctr_p(2) = 1.0d0
c      hyd_fctr_np(2) = 9.5/div    ! 10.5
      hyd_fctr_mem(2) = -2.07

c     Isoleucine I*2
c      hyd_fctr_p(3) = 1.0d0
c      hyd_fctr_np(3) = 8.0/div     ! 7.0
      hyd_fctr_mem(3) = -1.38

c     Valine V*2
c      hyd_fctr_p(4) = 1.0d0
c      hyd_fctr_np(4) = 8.5/div     ! 7.6
      hyd_fctr_mem(4) = -1.50

c     Proline P*2
c      hyd_fctr_p(5) = 1.0d0
c      hyd_fctr_np(5) = 11.18/div   ! 11.18
      hyd_fctr_mem(5) = -2.21

c     Methionine M*2
c      hyd_fctr_p(6) = 1.0d0
c      hyd_fctr_np(6) = 7.0/div     ! 7.0
      hyd_fctr_mem(6) = -1.38

c     Phenylalanine F*2
c      hyd_fctr_p(7) = 1.0d0
c      hyd_fctr_np(7) = 11.6/div    ! 11.6  
      hyd_fctr_mem(7) = -2.29

c     Tryptophane  W*2: f_p = 1  , f_np = 5.90/4 = 1.48
c      hyd_fctr_p(8) = 1.0d0
c      hyd_fctr_np(8) = 11.20/div   ! 15.2
      hyd_fctr_mem(8) = -3.00

C     TEST PURPUSES
c     -----------------------------------
C     -----------------------------------

c     Alanine  A*2 
      hyd_fctr_p(1) = 1.0d0
      hyd_fctr_np(1) = 1.0/div     ! 3.52

c     Leucine  L*2 
      hyd_fctr_p(2) = 1.0d0
      hyd_fctr_np(2) = 1.0/div    ! 10.5

c     Isoleucine I*2
      hyd_fctr_p(3) = 1.0d0
      hyd_fctr_np(3) = 1.0/div     ! 7.0

c     Valine V*2
      hyd_fctr_p(4) = 1.0d0
      hyd_fctr_np(4) = 1.0/div     ! 7.6

c     Proline P*2
      hyd_fctr_p(5) = 1.0d0
      hyd_fctr_np(5) = 1.0/div   ! 11.18

c     Methionine M*2
      hyd_fctr_p(6) = 1.0d0
      hyd_fctr_np(6) = 1.0/div     ! 7.0

c     Phenylalanine F*2
      hyd_fctr_p(7) = 1.0d0
      hyd_fctr_np(7) = 1.0/div    ! 11.6  

c     Tryptophane  W*2: f_p = 1  , f_np = 5.90/4 = 1.48
      hyd_fctr_p(8) = 1.0d0
      hyd_fctr_np(8) = 1.0/div   ! 15.2




c ##########################
      
      do i=6,13
         c_pol(i)      =  cg_b_polar_def(i)*hyd_fctr_p(i-5)
         c_nonpol(i)   =  cg_b_nonpolar_def(i)*hyd_fctr_np(i-5)
         c_membrane(i) =  cg_b_membrane_def(i)*hyd_fctr_mem(i-5)
      enddo


C     Intermediate values for the polar



c     Serine  S*2
      pol_fctr_p(1)  = 1.0d0
      pol_fctr_np(1) = 1.0d0
c     Threonine T*2
      pol_fctr_p(2)  = 1.0d0 
      pol_fctr_np(2) = 1.0d0
c     Tyrosine Y*1
      pol_fctr_p(3)  = 1.0d0 
      pol_fctr_np(3) = 1.0d0
c     Cysteine C*2
      pol_fctr_p(4)  = 1.0d0 
      pol_fctr_np(4) = 1.0d0
c     Asparagine N*2
      pol_fctr_p(5)  = 1.0d0 
      pol_fctr_np(5) = 1.0d0
c     Glutamine  Q*2
      pol_fctr_p(6)  = 1.0d0 
      pol_fctr_np(6) = 1.0d0



      do i=14,19
         c_pol(i)      =  cg_b_polar_def(i)*pol_fctr_p(i-13)
         c_nonpol(i)   =  cg_b_nonpolar_def(i)*pol_fctr_np(i-13)
         c_membrane(i) =  cg_b_membrane_def(i)
      enddo


      return
      end




      subroutine find_res_type(get,seq_name, res_type, type0)

      implicit none
      character*3 seq_name
      integer res_type, get, type0

      if(get.eq.1) then
      
c     Aspartic acid
      if((seq_name.eq.'D*2').or.(seq_name.eq.'D**').or.(seq_name.eq.'ASP')) then  
          res_type = 1 
          type0 =  16
      endif

C     Glutamic acid
      if((seq_name.eq.'E*2').or.(seq_name.eq.'E**').or.(seq_name.eq.'GLU')) then  
          res_type = 2
          type0 =  17 
      endif

C     Lysine
      if((seq_name.eq.'K*2').or.(seq_name.eq.'K**').or.(seq_name.eq.'LYS')) then  
          res_type = 3
          type0 =  6
      endif

C     Arginine
      if((seq_name.eq.'R*2').or.(seq_name.eq.'R**').or.(seq_name.eq.'ARG')) then  
          res_type = 4 
          type0 =  15 
      endif

C     Histidine
      if((seq_name.eq.'H*2').or.(seq_name.eq.'H**').or.(seq_name.eq.'HIS')) then  
          res_type = 5 
          type0 =  7 
      endif

C     Histidine version 2
      if(seq_name.eq.'HIE') then  
          res_type = 5 
          type0 =  7 
      endif

C     Alanine
      if((seq_name.eq.'A*2').or.(seq_name.eq.'A**').or.(seq_name.eq.'ALA')) then  
          res_type = 6 
          type0 =  1 
      endif

C     Leucine
      if((seq_name.eq.'L*2').or.(seq_name.eq.'L**').or.(seq_name.eq.'LEU')) then  
          res_type = 7 
          type0 =  3 
      endif

C     Isoleucine
      if((seq_name.eq.'I*2').or.(seq_name.eq.'I**').or.(seq_name.eq.'ILE')) then  
          res_type = 8 
          type0 =  4 
      endif

C     Valine
      if((seq_name.eq.'V*2').or.(seq_name.eq.'V**').or.(seq_name.eq.'VAL')) then  
          res_type = 9 
          type0 =  2 
      endif

C     Proline
      if((seq_name.eq.'P*2').or.(seq_name.eq.'P**').or.(seq_name.eq.'PRO')) then  
          res_type = 10
          type0 =  5 
      endif

C     Methiodine
      if((seq_name.eq.'M*2').or.(seq_name.eq.'M**').or.(seq_name.eq.'MET')) then  
          res_type = 11
          type0 =  19 
      endif

C     Phenylalanine
      if((seq_name.eq.'F*2').or.(seq_name.eq.'F**').or.(seq_name.eq.'PHE')) then  
          res_type = 12
          type0 =  8 
      endif

C     Tryptophane
      if((seq_name.eq.'W*2').or.(seq_name.eq.'W**').or.(seq_name.eq.'TRP')) then  
          res_type = 13
          type0 =  10 
      endif


C     Serine
      if((seq_name.eq.'S*2').or.(seq_name.eq.'S**').or.(seq_name.eq.'SER')) then  
          res_type = 14
          type0 =  13 
      endif

C     Threonine
      if((seq_name.eq.'T*2').or.(seq_name.eq.'T**').or.(seq_name.eq.'THR')) then  
          res_type = 15
          type0 =  14 
      endif

C     Tyrosine
      if((seq_name.eq.'Y*2').or.(seq_name.eq.'Y**').or.(seq_name.eq.'TYR')) then  
          res_type = 16
          type0 =  9 
      endif

C     Cysteine
      if((seq_name.eq.'C*2').or.(seq_name.eq.'C**').or.(seq_name.eq.'CYS')) then  
          res_type = 17
          type0 =  18 
      endif

C     Asparagine
      if((seq_name.eq.'N*2').or.(seq_name.eq.'N**').or.(seq_name.eq.'ASN')) then  
          res_type = 18
          type0 =  19 
      endif

C     Glutamine
      if((seq_name.eq.'Q*2').or.(seq_name.eq.'Q**').or.(seq_name.eq.'GLN')) then  
          res_type = 19
          type0 =  12 
      endif

      endif

      if(get.eq.2) then

      if(res_type.eq.1) type0 = 16
      if(res_type.eq.2) type0 = 17
      if(res_type.eq.3) type0 = 6
      if(res_type.eq.4) type0 = 15
      if(res_type.eq.5) type0 = 7
      if(res_type.eq.6) type0 = 1
      if(res_type.eq.7) type0 = 3
      if(res_type.eq.8) type0 = 4
      if(res_type.eq.9) type0 = 2
      if(res_type.eq.10) type0 = 5
      if(res_type.eq.11) type0 = 19
      if(res_type.eq.12) type0 = 8
      if(res_type.eq.13) type0 = 10
      if(res_type.eq.14) type0 = 13
      if(res_type.eq.15) type0 = 14
      if(res_type.eq.16) type0 = 9
      if(res_type.eq.17) type0 = 18
      if(res_type.eq.18) type0 = 11
      if(res_type.eq.19) type0 = 12

      endif


c      write(6,*)"type is", res_type, type0,seq_name
      return
      end




      subroutine print_self_energy(option,i,seq_name, c_polar, c_nonpolar, c_mem,
     $ u_polar, u_nonpolar, u_mem, min_wat_d, u_total, n_p, n_np, n_mem)

      implicit none
      include 'include/parameter.dc'
      include 'include/filename.h'
      character*3 seq_name
      character*13  outfile1, outfile2, outfile3
      integer option, i, p_file
      real*8 c_polar, c_nonpolar, c_mem, u_polar, u_nonpolar, u_mem,
     $ u_total, n_p, n_np, n_mem,min_wat_d
c      real*8   sum_sp(3)


c      do i=1,3
c          sum_sp(i) = 0.0d0
c      enddo

      if(option.eq.0) then
          
          outfile1 = 'table_ion.txt'
          outfile2 = 'table_pol.txt'

          open(987,file=path_of_out(1:ndir)//outfile1,form='formatted',
     $           status='unknown')

          open(988,file=path_of_out(1:ndir)//outfile2,form='formatted',
     $           status='unknown')


          write(987,*)"                                    IONIZABLE RESIDUES    "
          write(987,*)"-----------------------------------------------------------------------------------------"
          write(987,*)'Residue   c      c       c       N       N      N    Min     U        U       U       U  '
          write(987,*)'          pol   non_p   mem             NON          Wat             NON                 '
          write(987,*)'                                Polar   Polar   mem  dist  POLAR    POLAR    MEM    TOTAL'
          write(987,*)"----------------------------------------------------------------------------------------"


          write(988,*)"                                    POLAR RESIDUES    "
          write(988,*)"-----------------------------------------------------------------------------------------"
          write(988,*)'Residue   c      c       c       N       N      N    Min     U        U       U       U  '
          write(988,*)'          pol   non_p   mem             NON          Wat             NON                 '
          write(988,*)'                                Polar   Polar   mem  dist  POLAR    POLAR    MEM    TOTAL'
          write(988,*)"-----------------------------------------------------------------------------------------"
  


  

       endif

       if((option.ne.0).and.(option.ne.4)) then

          if(option.eq.1) p_file = 987     ! Ionizable
          if(option.eq.2) p_file = 988     ! Polar



          write(p_file,'(i4,1x,a3,2x,f5.2,2x,f5.2,2x,f5.2,1x,
     $    f8.3, f8.3, f8.3,1x,
     $    f4.1,
     $    f8.3, f8.3, f8.3, f8.3 )')
     $
     $    i,seq_name,     ! i, seq_name, vdw_residue_i
     $    c_polar, c_nonpolar, c_mem, !  The constants: polar  nonpolar mem 
     $    n_p, n_np, n_mem, ! the neighbors N_nonpolar, N_polar, N_mem
     $    min_wat_d,        ! The minimum water distance
     $    u_polar, u_nonpolar, u_mem, u_total  ! The energy terms

c          sum_sp(option) = sum_sp(option)+ u_total

       endif 


      if(option.eq.4) then

          write(987,*)"-----------------------------------------------------------------------------------------"
          write(987,*)'Residue   c      c       c       N       N      N     Min     U        U       U       U  '
          write(987,*)'          pol   non_p   mem             NON           Wat             NON                 '
          write(987,*)'                                Polar   Polar   mem   dist  POLAR    POLAR    MEM    TOTAL'
          write(987,*)"-----------------------------------------------------------------------------------------" 


          write(988,*)"-----------------------------------------------------------------------------------------"
          write(988,*)'Residue   c      c       c       N       N      N    Min     U        U       U       U  '
          write(988,*)'          pol   non_p   mem             NON          Wat             NON                 '
          write(988,*)'                                Polar   Polar   mem  dist  POLAR    POLAR    MEM    TOTAL'
          write(988,*)"-----------------------------------------------------------------------------------------"
          write(988,*)
          write(988,*)"                                                        Sum of u POLAR = ",u_polar  



          close(987)
c          close(988)

    
      endif


      return
      end

c-------------------------------------------------------------------------- 



c               call print_hydro_energy(is_sim(i),i,seq_tmp, c_polar, c_np, c_mem,
c     $          temp_polar, temp_u2, temp_u1, temp_mem, temp_energy,
c     $          t_p_ngb, t_np_ngb,temp_general,t_mem_ngb, fact_u1, fact_u2)


      subroutine print_hydro_energy(option,i,seq_name, c_polar, c_nonpolar, c_mem,
     $ u_polar, u_nonpolar, u_ring, u_mem, min_wat_d, u_total, n_p, n_np, n_ring, n_mem,factor_u1,factor_u2)

      implicit none
      include 'include/parameter.dc'
      include 'include/filename.h'
      character*3 seq_name
      character*13  outfile1, outfile2, outfile3
      integer option, i, p_file
      real*8 c_polar, c_nonpolar, c_mem, u_polar, u_nonpolar, u_ring, u_mem,
     $ u_total, n_p, n_np, n_ring, n_mem, factor_u1, factor_u2, min_wat_d
c      real*8   sum_sp(3)


c      do i=1,3
c          sum_sp(i) = 0.0d0
c      enddo

      if(option.eq.0) then
          
          outfile3 = 'table_hyd.txt'

 
          open(989,file=path_of_out(1:ndir)//outfile3,form='formatted',
     $           status='unknown')



      write(989,*)"                                    NON POLAR (HYDROPHOBIC) RESIDUES    "
      write(989,*)"---------------------------------------------------------------------------------------------------------------"
      write(989,*)'Residue   c      c       c       N       N       N      N      Min      U     U       U        U       U     U'
      write(989,*)'          pol   non_p   mem             NON    Ring            Wat           NON     Ring     np+               '
      write(989,*)'                                Polar   Polar          mem     dis   POLAR   POLAR           ring     MEM   TOT'
      write(989,*)"---------------------------------------------------------------------------------------------------------------"


       endif

       if((option.ne.0).and.(option.ne.4)) then

          if(option.eq.3) p_file = 989     ! Non Polar (Hydrophobic)


          write(p_file,'(i4,1x,a3,2x,f5.2,2x,f5.2,2x,f5.2,1x,
     $    f8.3, f8.3, f8.3, f8.3,1x,
     $    f4.1,
     $    f8.3, f8.3, f8.3, f8.3, f8.3, f8.3 )')
     $
     $    i,seq_name,     ! i, seq_name, vdw_residue_i
     $    c_polar, c_nonpolar, c_mem, !  The constants: polar  nonpolar mem 
     $    n_p, n_np,  n_ring, n_mem, ! the neighbors N_nonpolar, N_polar, N_mem 
c                                    and the N_ring 
     $    min_wat_d,
     $    u_polar, u_nonpolar, u_ring, factor_u2*u_nonpolar+factor_u1*u_ring,u_mem, u_total  ! The energy terms



       endif 


      if(option.eq.4) then

      write(989,*)"---------------------------------------------------------------------------------------------------------------"
      write(989,*)'Residue   c      c       c       N       N       N      N      Min     U      U       U        U       U    U'
      write(989,*)'          pol   non_p   mem             NON    Ring            Wat           NON     Ring     np+               '
      write(989,*)'                                Polar   Polar          mem     d     POLAR   POLAR            ring    MEM  TOT'
      write(989,*)"---------------------------------------------------------------------------------------------------------------"


          write(989,*)
          write(989,*)"     Sum of u HYDRO = ",u_total
  



    
      endif


      return
      end

c-------------------------------------------------------------------------- 













      subroutine cg_nonbond_list
      implicit none
      include 'include/parameter.dc'
      include 'include/Mdyn.h'
      include 'include/Mmol.h'
      include 'include/cmolec.h'
      include 'include/coord.h'
      include 'include/bdparm.h'
      include 'include/exclude_pair.h'
      include 'include/seq.h'

c:::  local var
      integer i,j,k,kk,m,jb,kb,hit_res(2000),ihit_res(2000),i_hit(2000,200),
     $        i_cb,npair,j_atom(200),ii,n_grid_point
      integer*4 imask,ibit
      real*8 dis

      character*3 nonpolar(9)

      data nonpolar/'A*2','F*2','I*2','L*2','M*2','P*2','V*2','W*2','MEB'/
c.......................................................................
chu   this subroutine is only called once

      n_cb_res=0
      do i=1,9
         do j=1,numres
            if(resseq(j).eq.nonpolar(i)) then
               jb=np(j)+5
               if(resseq(j).eq.'P*2') jb=np(j)+4
               if(resseq(j).eq.'MEB') jb=np(j)+1
               n_cb_res=n_cb_res+1
               cb_non(n_cb_res)=0
               cb_inon(n_cb_res)=jb
               hit_res(n_cb_res)=0
               ihit_res(n_cb_res)=j
               write(6,'(/,'' nonbond pair list for hydrophobic'',
     $                     '' residue '',i4,''_'',a3,'' with radius:'',f8.1,
     $                     ''A'')')j,resseq(j),cb_nlist_r
               write(6,'(''    res_i      pair#_i    CB_i    atom_j'',
     $                   ''        res_j     dist_res_ij'')')
               do k=1,numres
                  if(k.ne.j) then
                     kb=np(k)+5
                     if(resseq(j).eq.'P*2') kb=np(k)+4
                     if(resseq(j).eq.'MEB') kb=np(k)+1
                     dis=0.0
                     call dist_atoms(jb,kb,dis)
                     if(dis.le.cb_nlist_r) then
                        hit_res(n_cb_res)=hit_res(n_cb_res)+1
                        i_hit(n_cb_res,hit_res(n_cb_res))=k
                        do m=np(k)+1,np(k+1)
c                          [ is jb sr excluded to m ? ]
                           if (abs(jb-m).le.nrange) then
                              if (m.lt.jb) then
                                 imask=wnb(m)
                                 ibit=32-(jb-m)
                                 if (btest(imask,ibit)) then
                                    goto 44
                                 endif
                              else
                                 imask=wnb(jb)
                                 ibit=32-(m-jb)
                                 if (btest(imask,ibit)) then
                                    goto 44
                                 endif
                              endif
                           endif

                           cb_non(n_cb_res)=cb_non(n_cb_res)+1

                           if(n_cb_res.gt.2000) then
                              write(6,'('' ERROR: # of CG hydra residues>2000'')')
                              goto 666
                           endif
                           if(cb_non(n_cb_res).gt.100) then
                              write(6,'('' ERROR: # of nonbonded atoms>200'')')
                              goto 666
                           endif

                           cb_jnon(n_cb_res,cb_non(n_cb_res))=m
                           write(6,'(i5,''_'',a3,4i10,''_'',a3,f11.1)')j,
     $                          resseq(j),cb_non(n_cb_res),cb_inon(n_cb_res),
     $                          cb_jnon(n_cb_res,cb_non(n_cb_res)),k,resseq(k),
     $                          dis
 44                     enddo
                     endif
                  endif
               enddo
            endif
        enddo
      enddo

chu   Spyros, you may call cg_grid_points in any other place you want
      call cg_grid_points

c.......................................................................
      return
 666  call killme ('cg_nonbond_list')
      end
c=============================================================================

      subroutine cg_grid_points

      implicit none
      include 'include/parameter.dc'
      include 'include/Mdyn.h'
      include 'include/seq.h'
      include 'include/coord.h'

      integer i,ii,npair,i_cb,j_atom(200),n_grid_point,ires

c...........................................................................
       write(6,'(//,'' ..................................................'',/,
     $              '' generating grids for all the hydrophobic_resiudes:'')')
      do i=1,n_cb_res
         i_cb=cb_inon(i)
         npair=cb_non(i)
         do ii=1,npair
            j_atom(ii)=cb_jnon(i,ii)
         enddo
         call gen_grid_cg_nonbond(i_cb,npair,j_atom,n_grid_point)
         n_cg_grid_points(i)=n_grid_point
      enddo

      write(6,'('' ........................................................'',/
     $          '' hydrophobic_residue    #_of_grid_points'')')

      do i=1,n_cb_res
         ires=atom_ires(cb_inon(i))
         write(6,'(i8,''_'',a3,16x,i5)')ires,resseq(ires),
     $                                  n_cg_grid_points(i)
      enddo

      return

      end

c============================================================================
      subroutine  find_Nw(Nw, Nw_ring, print_option)
      implicit none
      include 'include/parameter.dc'
      include 'include/Mdyn.h'
      include 'include/coord.h'
      include 'include/dgm_pt.h'
      include 'include/seq.h'
      include 'include/exclude.h'
      include 'include/filename.h'

      
      integer i, j, k, j3 
      real*8 xGridTemp(3)
      real*8 gridSpacing     ! for building grid(3, mxlgvn)
      real*8 grid(3,mxlgvn)  ! dense grid wich is has no clash with pro and mem but also 
                             ! not too far from it
      integer Nw(numres), Nw_ring(numres)
      real*8 xC, yC,zC ! current coord of CB of ionized residue
      integer Ngrid, temp_value_n
      real*8 closestDistance   ! function to calc closest dist between current point and grid
      real*8 ProteinStart_ZComponent ! leftmost z coord of the pro
      logical is_Good_Point ! if point close but not too close to protein and the same about mem
      real*8 scaleBox   ! we increase box dimentsions by this val
      real*8 R_wat_cutoff, temp_value_distance
      logical print_option

      scaleBox = 1.3    !box for grid is 30% bigger than real box for the system
      gridSpacing = 3.0
      Ngrid = 0
C      R_wat_cutoff =  7.0d0    Debatable! it should be larger
      R_wat_cutoff =  11.0d0   ! Debatable! it should be larger

      if(print_option) then 
          write (6,'('' center = '',3f10.3)') xs(1), xs(2), xs(3)
          write (6,'('' box    = '',3f10.3)') box(1), box(2), box(3)
      endif

      !let's build rec grid but save only grid points close to the protein
      xGridTemp(1) = xs(1) - box(1)*0.5*scaleBox

c     Initialize the Nw
      do i=1,numres
         Nw(i)=0.0d0
         Nw_ring(i) =0.0d0
      enddo

      do while ( (xGridTemp(1) - xs(1)).le.box(1)*0.5*scaleBox) 
          xGridTemp(2) = xs(2) - box(2)*0.5*scaleBox
          do while ( (xGridTemp(2) - xs(2)).le.box(2)*0.5*scaleBox) 
              xGridTemp(3) = xs(3) - box(3)*0.5*scaleBox
              do while ( (xGridTemp(3) - xs(3)).le.box(3)*0.5*scaleBox) 
                  !let's check if point is worth saving
                  if (is_Good_Point(xGridTemp(1), xGridTemp(2),xGridTemp(3)).eq..true.) then
                      Ngrid = Ngrid + 1
                      if (Ngrid.gt.mxlgvn) then
                          stop 'too many grid points'
                      endif
                      grid(1,Ngrid) = xGridTemp(1)
                      grid(2,Ngrid) = xGridTemp(2)
                      grid(3,Ngrid) = xGridTemp(3)
                  endif 
                  xGridTemp(3) = xGridTemp(3) + gridSpacing
              enddo
              xGridTemp(2) = xGridTemp(2) + gridSpacing
          enddo
          xGridTemp(1) = xGridTemp(1) + gridSpacing
      enddo

      ! at this point relevant grid should be calculated

      open (787, file=path_of_out(1:ndir)//'grid_new_spyros.pdb')
      do i = 1, Ngrid
          write(787,8) i, 'CB  ', 'NNN', i, (grid(k,i),k=1,3)
      enddo
      close (787)
      ! now let's determine closest distance from grid to pro or x_{center}
      do i = 1, numres


          j = np(i) +5
          if(resseq(i).eq.'MEB') j=np(i)+1
          if(resseq(i).eq.'P*2') j=np(i)+4

          j3 = 3*(j-1)
          xC = x(j3+1)
          yC = x(j3+2)
          zC = x(j3+3)

          call n_neighbor_hydrophobic(grid, Ngrid, xC,yC,zC,
     $        R_wat_cutoff, temp_value_n,temp_value_distance)

          Nw(i) = temp_value_n
          min_wat_distance(i) = temp_value_distance

          if(is_sim(i).eq.3) then
              call cb_nonbond_list1(i,Nw_ring(i), print_option)
          endif

      enddo
     
      return
 8    format('ATOM',2x,i5,2x,a4,a3,i6,4x,3f8.3)
      end






      logical function is_Good_Point(x,y,z)
      ! the grid is good if 
      ! 2<r_{pro}^{closest}<10
      ! mem_{spacing}*sqrt(3)/2< r_{mem}^{closest}<10
      implicit none
      include 'include/parameter.dc'
      include 'include/seq.h'
      real*8 rPro ! min distance to protein
      real*8 rMem ! min distance to membrane
      real*8 x,y,z ! current grid position
      real*8 rMinPro, rMaxPro
      real*8 rMinMem, rMaxMem
      logical hasProClash, hasMemClash
      rMinPro = 3.5d0     ! maybe 4.0 but need to test (3.0 gives several points inside the protein)
      rMaxPro = 10.0
      rMinMem = mem_spacing*1.75*0.5
      rMaxMem = 10.0
      call closest_distance_mem(x,y,z,rPro, rMem)
      hasProClash = rPro.le.rMinPro
      hasMemClash = rMem.le.rMinMem
      if (hasProClash.or.hasMemClash) then      !too close to pro or mem
          is_Good_Point = .false.
      else                       ! possibly good point - no clash
           if (rMem.lt.rMaxMem .or. rPro.lt.rMaxPro) then  !close enough and not extremely far but no clash
               is_Good_Point = .true.
           else    ! no clash but too far from pro or mem to be interesting point
               is_Good_Point = .false.
           endif 
      endif
      return
      end



      subroutine n_neighbor_hydrophobic(grid, Ngrid,x,y,z, R_wat_cutoff, n_neighbor, min_distance )
      implicit none
      include 'include/parameter.dc'
      include 'include/seq.h'

      
      integer i, j, k, n_neighbor
      integer Ngrid
      real*8 grid(3,mxlgvn), R_wat_cutoff, dist, min_distance
      real*8 x,y,z


C     Initialize the variables
      n_neighbor = 0
      min_distance = 1000.0d0    ! A starting very big value

      do i = 1, Ngrid
           dist = 0.0 + (x - grid(1,i))**2
           dist = dist   + (y - grid(2,i))**2
           dist = dist   + (z - grid(3,i))**2

           dist = dsqrt(dist)

           if (dist.lt.R_wat_cutoff) then
              n_neighbor = n_neighbor + 1
           endif

           if (dist.lt.min_distance) then
              min_distance = dist
           endif

      enddo


      return
      end


      subroutine closest_distance_mem(xC,yC,zC,rPro, rMem)
!     calculate min distance between point (x,y,z) to pro and mem
      implicit none
      include 'include/parameter.dc'
      include 'include/Mdyn.h'
      include 'include/coord.h'
      include 'include/seq.h'
      integer j,jj
      real*8 xC, yC, zC,rPro, rMem
      real*8 dis_minPRO, dis_minMEM, d_mem, d_pro

      dis_minPRO=1000000000.0
      dis_minMEM=1000000000.0
      do j=1,numres
          if (resseq(j).eq.'MEB') then   ! do things for membrane
              jj = 3*( np(j)+1 - 1)
              d_mem = 0.0   + (xC - x(jj+1))**2
              d_mem = d_mem + (yC - x(jj+2))**2
              d_mem = d_mem  + (zC - x(jj+3))**2
              dis_minMEM=dmin1(d_mem,dis_minMEM)
          else
              do jj=np(j)+1,np(j+1)
                  d_pro = 0.0
                  d_pro = d_pro + (xC-x(3*(jj-1)+1))**2
                  d_pro = d_pro + (yC-x(3*(jj-1)+2))**2
                  d_pro = d_pro + (zC-x(3*(jj-1)+3))**2
                  dis_minPRO=dmin1(d_pro,dis_minPRO)
              enddo
          endif
      enddo !end loop over all residues
      rPro = sqrt(dis_minPRO)
      rMem = sqrt(dis_minMEM)
      return
      end


      subroutine cb_nonbond_list1(j, Nw, print_option)
      implicit none
      include 'include/parameter.dc'
      include 'include/Mdyn.h'
      include 'include/Mmol.h'
      include 'include/cmolec.h'
      include 'include/coord.h'
      include 'include/bdparm.h'
      include 'include/exclude_pair.h'
      include 'include/seq.h'

c:::  local var
      integer j, Nw
      integer i,k,kk,m,jb,kb,hit_res(1000),ihit_res(1000),i_hit(1000,200),
     $        i_cb,npair,j_atom(200),ii
      integer*4 imask,ibit
      real*8 dis
      logical print_option
c.......................................................................
      Nw = 0
      cb_nlist_r = 11.0d0
      n_cb_res=1


      jb=np(j)+5
      if(resseq(j).eq.'P*2') jb=np(j)+4
      if(resseq(j).eq.'MEB') jb=np(j)+1

      cb_non(n_cb_res)=0
      cb_inon(n_cb_res)=jb
      hit_res(n_cb_res)=0
      ihit_res(n_cb_res)=j
c      write(6,'(/,'' nonbond pair list for user-specified type'',
c     $       '' for res_i with radius:'',f8.1,''A'')')cb_nlist_r
c      write(6,'(''    res_i      pair#_i    CB_i    atom_j'',
c     $     ''        res_j     dist_resij'')')

      do k=1,numres
            if(k.ne.j) then
                kb=np(k)+5
                if(resseq(k).eq.'P*2') kb=np(k)+4
                if(resseq(k).eq.'MEB') kb = np(k)+1
                dis=0.0
                call dist_atoms(jb,kb,dis)
                if(dis.le.cb_nlist_r) then
                    hit_res(n_cb_res)=hit_res(n_cb_res)+1
                    i_hit(n_cb_res,hit_res(n_cb_res))=k
                    do m=np(k)+1,np(k+1)
c                      [ is jb sr excluded to m ? ]
                       if (abs(jb-m).le.nrange) then
                          if (m.lt.jb) then
                             imask=wnb(m)
                             ibit=32-(jb-m)
                             if (btest(imask,ibit)) then
                                goto 44
                             endif
                          else
                             imask=wnb(jb)
                             ibit=32-(m-jb)
                             if (btest(imask,ibit)) then
                                 goto 44
                             endif
                          endif
                       endif

                       cb_non(n_cb_res)=cb_non(n_cb_res)+1
                       cb_jnon(n_cb_res,cb_non(n_cb_res))=m
c                       write(6,'(i5,''_'',a3,4i10,''_'',a3,f11.1)')j,
c     $                   resseq(j),cb_non(n_cb_res),cb_inon(n_cb_res),
c     $                   cb_jnon(n_cb_res,cb_non(n_cb_res)),k,resseq(k),
c     $                   dis

 44                 enddo
                endif
            endif
      enddo

      do ii=1,cb_non(n_cb_res)
          j_atom(ii)=cb_jnon(n_cb_res,ii)  ! the neighboring atoms of residue j
      enddo
      i_cb=cb_inon(n_cb_res)    ! The CB atom if residue j
      npair=cb_non(n_cb_res)    ! The total number of OTHER atoms, "close" to residue j
               
      call gen_grid_cg_nonbond1(i_cb,npair,j_atom, Nw, print_option)



c.......................................................................
      return
 666  call killme ('cb_nonbond_list')
      end
c=============================================================================
c-------------------------------------------------------------------------
      subroutine gen_grid_cg_nonbond1(icb,npair,j_atom, n_grid_point, print_option)
      implicit none


      include 'include/parameter.dc'
      include 'include/Mdyn.h'
      include 'include/coord.h'
      include 'include/seq.h'
      include 'include/cbkv_gap.h' 
      include 'include/lib.h'

c...  input vars 
      integer icb,npair,j_atom(npair)

c...  local vars
      integer i,j,k,m,ii,jj,kk,limit,mid,n_grid_point,ityp,ires, ia
      logical lcube(0:MXCUBE,0:MXCUBE,0:MXCUBE)
      pointer (p_lcube,lcube)

      real*8 spacing,rg2,drg2,ri,rj,rk,d2,rvdw,rvdw2,r_hydra(24),r_hydra2(24),r2_1,r2_2

      logical ligoodrange
      logical print_option
      external ligoodrange


c                      ORIGINAL
c      data r_hydra/3.0,3.5,4.0,4.5,3.5,3.0,3.0,5.0,3.0,6.0,
c     $             3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,
c     $             3.0,3.0,3.0,3.5/



C                      NEW TRIALS!!!!!
      data r_hydra/3.0,4.5,5.5,6.0,3.5,3.0,3.0,5.0,3.0,6.0,
     $             3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,6.0,3.0,
     $             3.0,3.0,3.0,3.5/

C     Those are the corresponding to the h_hydra residues!!!!
c      data from_amino/'ALA','VAL','LEU','ILE','PRO','LYS','HIS','PHE','TYR',
c     $    'TRP','ASN','GLN','SER','THR','ARG','ASP','GLU','CYS','MET','HIE',
c     $    'A  ','T  ','C  ','G  '/
c      data to_amino/  'A*2','V*2','L*2','I*2','P*2','K*2','H*2','F*2','Y*2',
c     $    'W*2','N*2','Q*2','S*2','T*2','R*2','D*2','E*2','C*2','M*2','H*2',
C     $    'A*3','T*3','C*3','G*3'/

C                     NEW TRIALS
      data r_hydra2/4.0,4.0,4.0,4.0,4.0,3.0,3.0,4.0,3.0,4.0,
     $             3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,4.0,3.0,
     $             3.0,3.0,3.0,3.0/

c                         Always 2
c      data r_hydra2/2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,
c     $             2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,
c     $             2.0,2.0,2.0,2.0/
c.......................................................................
      spacing=2.5
      rvdw=3.0
      rvdw2=rvdw*rvdw

      if(print_option) then
          write(6,'(/,'' RVDW_cutoff     = '',f7.1)') rvdw
          write(6,'('' CG_Grid spacing = '',f7.2)') spacing
          write(6,'('' CG_Grid radius  = '',f7.2,/)') cb_nlist_r
      endif

      p_lcube=LOC(lcube)
      allocate (lcube)

      limit=cb_nlist_r*2/spacing
      if(mod(limit,2).eq.0) limit=limit+1  !make limit into an odd number
      if(.not.ligoodrange('LIMIT:MXCUBE',limit,0,MXCUBE))goto 666
      mid=int(limit/2.+0.5d0)

      drg2=spacing*spacing

      rg2=cb_nlist_r**2

chu   create the grid sphere with a radius=cb_nlist_r

      do ii=1,limit
         do jj=1,limit
            do kk=1,limit
               ri=ii-mid
               rj=jj-mid
               rk=kk-mid

               d2=(ri*ri+rj*rj+rk*rk)*drg2
               if(d2.le.rg2) then
                  lcube(ii,jj,kk)=.true.
               else 
                  lcube(ii,jj,kk)=.false.
               endif
            enddo
         enddo
      enddo

chu   eleminate the grid points outside the sphere layer between r_hydra and
chu   r_hydra+3

      ires=atom_ires(icb)
      do i=1,24
         if(to_amino(i).eq.resseq(ires)) ityp=i
      enddo

      r2_1=r_hydra(ityp)**2
c      r2_2=(r_hydra(ityp)+2.0)**2
      r2_2=(r_hydra(ityp)+r_hydra2(ityp))**2

      do ii=1,limit
         do jj=1,limit
            do kk=1,limit
               if(lcube(ii,jj,kk)) then
                  ri=ii-mid
                  rj=jj-mid
                  rk=kk-mid
                  d2=(ri*ri+rj*rj+rk*rk)*drg2
                  if(d2.lt.r2_1.or.d2.gt.r2_2) lcube(ii,jj,kk)=.false.
               endif
            enddo
         enddo
      enddo

chu   find if any CG Grid point is too close to any cb_nonbond list atoms to
chu   be eleminated

      do ii=1,limit
         do jj=1,limit
            do kk=1,limit
               if(lcube(ii,jj,kk)) then
                  do i=1,npair
                     j=j_atom(i)


                     ri=x(3*(j-1)+1)-(x(3*(icb-1)+1)+(ii-mid)*spacing)
                     rj=x(3*(j-1)+2)-(x(3*(icb-1)+2)+(jj-mid)*spacing)
                     rk=x(3*(j-1)+3)-(x(3*(icb-1)+3)+(kk-mid)*spacing)
                     d2=ri*ri+rj*rj+rk*rk

                    

                     if(d2.lt.rvdw2) then
                        lcube(ii,jj,kk)=.false.
                        goto 10
                     endif

c                    For CB atoms
                     if(bkvcode(j)(1:2).eq.'CB') then
                        if(d2.lt.16.0d0) then
                            lcube(ii,jj,kk)=.false.
                            goto 10
                        endif
                     endif                          


 10               enddo
               endif
            enddo
         enddo
      enddo


c      open (787, file = 'grid_debug_PHE.pdb')

chu   collect the total existing grid points
      n_grid_point=0
      ia = 0
      do ii=1,limit
         do jj=1,limit
            do kk=1,limit
               if(lcube(ii,jj,kk)) then 
                 n_grid_point=n_grid_point+1

c                 ia = ia + 1
c                 write(787,9) ia, 'CB  ', 'NNN', ia, 
c     $           x(3*(icb-1)+1)+(ii-mid)*spacing, 
c     $           x(3*(icb-1)+2)+(jj-mid)*spacing,
c     $           x(3*(icb-1)+3)+(kk-mid)*spacing  
c 9               format('ATOM',2x,i5,2x,a4,a3,i6,4x,3f8.3)

               endif

            enddo
         enddo
      enddo


c      close (787)

      if(print_option) then 
          write(6,'('' there are '',i5,'' grid points around residue '',i4,
     $              '' within the sphere'',/,
     $              '' layer between the radii'',f4.1,'' and '',f4.1/)')
     $         n_grid_point,atom_ires(icb),r_hydra(ityp), r_hydra(ityp)+2.0
      endif

      deallocate(lcube)


      return
c-------------------------------------------------------------------------

666   call killme ('gen_grid_cg_nonbond')

      end
c=========================================================================





      subroutine calculate_Ualpha_array
c     calulate U_alpha array which is required in CG model  for Main chain solvation
c     and H-bond terms
      implicit none
      include  'include/parameter.dc'
      include  'include/Mdyn.h'
      include  'include/coord.h'
      include  'include/seq.h'
      include 'include/cg_ion_hyd_pol_constants.h'
      real*8 temp_u
      integer i
c=============================================================================================

      do 40 i = 1, numres
          if(is_sim(i).eq.4) goto 40  !go to the next residue  if current is meb
c         theta - percentage of polar residues around CA atom of the residue i
          theta(i)=(Nnonpol_CA(i)+Nmem_CA(i)*(N_p_max/N_mem_max)-N_p_max)/N_p_max    
          if(abs(theta(i)).lt.max_theta) then
             temp_u=exp(-a_theta_u*(abs(theta(i))-max_theta)**2)
          else
             temp_u=1.0
          endif
          u_alpha(i)=temp_u

c          if (i.eq.1) then
c          write(6,'('' For residue: '',i4,''_'',a3)'),i, resseq(i)
c          write(6,'('' Nnp Nmem theta U_alpha '',4f10.6)'),Nnonpol_CA(i),Nmem_CA(i),
c     $                                               theta(i),u_alpha(i)
c          endif

 40   continue
      return
      end

      real*8 function MainChainSolvation()
c     calulate main chain solvatation in CG model
c     based on knowldege of U_alpha
      implicit none
      include  'include/parameter.dc'
      include 'include/coord.h'
      include  'include/seq.h'
      include 'include/cbkv_gap.h'
      include 'include/Mdyn.h'
      include 'include/derivative.h'
      include 'include/symbols.h'
      include  'include/cg_ion_hyd_pol_constants.h'
      real*8 e_main_alpha
      real*8 b_solv
      real*8 r1, df, dx1, dx2, dx3

c     Variables for partial derivatives
      real*8 dE_main_du_alpha, du_alpha_dtheta, dtheta_dr, dtheta_dNnp, dtheta_dNmem,
     $       dNnp_dr, dNmem_dr

      integer i, j, ihit, jhit, ii, jj, i3, j3


c     ################################################################################################
c     Energy calculation
c     ################################################################################################
      e_main_alpha = 0.0

c     solvation for the residue in water
      b_solv = -2.0

      do i = 1, numres
          e_main_alpha = e_main_alpha+b_solv*u_alpha(i)
      enddo
      MainChainSolvation = e_main_alpha

c     ################################################################################################
c     Force calculation
c     ################################################################################################
      df=0.0

      do 40 i=1,numres ! loop over all but membrane residues
          if(is_sim(i).eq.4) goto 40  !go to the next residue  if current is meb
          ! get C-aplha atom of the residue
          ihit=0
          do ii=np(i)+1,np(i+1)
              if(iac_name(iac(ii)).eq.'DY') goto 40
              if(bkvcode(ii)(1:2).eq.'CA') then
                  ihit=ii
              endif
          enddo
          if(ihit.eq.0) goto 40
          do 30 j=1,numres  ! loop over all residue and calulate # heighbours
              if(i.eq.j) goto 30
              ! gre C-alpha atom of the j residue or CB of MEB
              jhit=0
              do jj=np(j)+1,np(j+1)
                 if(iac_name(iac(jj)).eq.'DY') goto 30
                 if(bkvcode(jj)(1:2).eq.'CA'.or.resseq(j).eq.'MEB') then
                    jhit=jj
                 endif
              enddo
              if(jhit.eq.0) goto 30
              r1 = 0.0
              call dist_atoms(ihit,jhit,r1)
              i3   = ihit*3-3
              j3   = jhit*3-3
              dx1  = x(i3+1)-x(j3+1)
              dx2  = x(i3+2)-x(j3+2)
              dx3  = x(i3+3)-x(j3+3)


c             E_main is a function of
c             E_main(u_alpha(theta(Nnp(r),Nmem(r))))
c             So, we need the following partial derivatives:

              dE_main_du_alpha=b_solv

              if (abs(theta(i)).ge.max_theta) then
                 du_alpha_dtheta=0
              else
                 du_alpha_dtheta=-2*a_theta_u*(abs(theta(i))-max_theta)
     $                           *exp(-a_theta_u*(abs(theta(i))-max_theta)**2)
              endif

              dtheta_dNnp=1/N_p_max
              dtheta_dNmem=1/N_mem_max

              dNnp_dr=0.0
              if (is_sim(j).eq.3) then
                 if (r1.le.r_np) then
                    dNnp_dr=0
                 else
                    dNnp_dr=-2*a_np*(r1-r_np)*exp(-a_np*(r1-r_np)**2)
                 endif
              endif

              dNmem_dr=0
              if (is_sim(j).eq.4) then
                 if (r1.le.r_mem_ngb_cutoff) then
                    dNmem_dr=0
                 else
                    dNmem_dr=-2*a_mem*(r1-r_mem_ngb_cutoff)*exp(-a_mem*(r1-r_mem_ngb_cutoff)**2)
                 endif
              endif

              dtheta_dr=dtheta_dNnp*dNnp_dr+dtheta_dNmem*dNmem_dr

              df=dE_main_du_alpha*du_alpha_dtheta*dtheta_dr

              df=-df/r1

              d(i3+1)=d(i3+1)+df*dx1
              d(i3+2)=d(i3+2)+df*dx2
              d(i3+3)=d(i3+3)+df*dx3

              d(j3+1)=d(j3+1)-df*dx1
              d(j3+2)=d(j3+2)-df*dx2
              d(j3+3)=d(j3+3)-df*dx3

 30       continue
 40   continue

      return
      end

      real*8 function HydrogenBond_CG()
c     calulate Hydrogen bond in CG model
c     based on U_alpha array
      implicit none
      include  'include/parameter.dc'
      include  'include/Mdyn.h'
      include  'include/coord.h'
      include  'include/bdparm.h'
      include  'include/derivative.h'
      include  'include/stepsize.h'
      include  'include/seq.h'
      include  'include/cbkv_gap.h' 
      include  'include/cg_ion_hyd_pol_constants.h'
      real*8 factor_wat, factor_mem, hb_m, hb_r
      real*8 e_hb_wat, e_hb_mem, e_hb_pair, ehbpp
      real*8 df, dx1, dx2, dx3, r1

      real*8  x1dis, y1dis, z1dis, dis2, angle
      integer atomz,atom_a1,atom_a2

c     Variables for partial derivatives
      real*8 dE_hb_dE_hb_wat,dE_hb_dE_hb_mem,dE_hb_du_alpha(2),
     $       dE_hb_wat_dE_hbij,dE_hb_mem_dE_hbij,dE_hb_wat_dr,
     $       dE_hbij_dr,du_alpha_dtheta(2),
     $       dtheta_dNnp,dtheta_dNmem,dNnp_dr,dNmem_dr,dtheta_dr

      integer i, i_alpha, j_alpha, ica, jca, i3, j3, k

      ehbpp = 0.0   ! total H_bond energy

c     Scaling factors for HB in water and in membrane are taken form the following equation:
c     factor = (1-0.8*u_alpha(ica)*u_alpha(jca))/4.5
c     For water U_alfa=1 which gives
      factor_wat=0.2/4.5
c     For membrane U_alfa=0 which gives
      factor_mem=1/4.5

c     Parameters for the added Gaussian function (steepness and center)
      hb_m=22.2
      hb_r=2.9

c     For force debugging:
c         do i=1,nat3
c            d(i)=0.0
c         enddo

c     DO cycle for all HB pairs


      do i=1,hb_ij_count
c        Find the CA numbers for corresponding HB pair atoms
         if (resseq(atom_ires(hb_i(i))).eq.'GLY') then
c            write(6,*)"DEBUG: hit gly",hb_i(i)
            i_alpha = hb_i(i)-4
         else
            i_alpha = hb_i(i)-5
         endif
         j_alpha = hb_j(i)+1
c        Now find the residue numbers for corresponding CA's
         ica = atom_ires(i_alpha)
         jca = atom_ires(j_alpha)

c        ################################################################################################
c        Energy calculation
c        ################################################################################################

         r1 = 0.0
         call dist_atoms(hb_i(i),hb_j(i),r1)

c        try to find the Nitrogen atom! Only in CG
         if(bkvcode(hb_i(i))(1:1).eq.'H') then 

             atomz = hb_i(i)-1
             atom_a1 = hb_j(i)
             atom_a2 = hb_i(i)

             x1dis=x(atomz*3-3+1)-x(hb_i(i)*3-3+1)
             y1dis=x(atomz*3-3+2)-x(hb_i(i)*3-3+2)
             z1dis=x(atomz*3-3+3)-x(hb_i(i)*3-3+3)
         endif


         if(bkvcode(hb_j(i))(1:1).eq.'H') then

             atomz = hb_j(i)-1
             atom_a1 = hb_i(i)
             atom_a2 = hb_j(i)

             x1dis=x(atomz*3-3+1)-x(hb_j(i)*3-3+1)
             y1dis=x(atomz*3-3+2)-x(hb_j(i)*3-3+2)
             z1dis=x(atomz*3-3+3)-x(hb_j(i)*3-3+3)
         endif

         dis2 = x1dis**2+y1dis**2+z1dis**2
         dis2 = sqrt(dis2)

         call check_ang(atom_a1,atom_a2,atomz,angle)


c        Function for HB in water, e_hb_wat, has a shallow minimum and a Gaussian centered at 2.9 A.
c        This reflects the fact that we have to spend some energy to break the HB, but once it is broken,
c        the HB with water molecules can be formed.

         e_hb_wat=vdw_hb_ij(i)*factor_wat+exp(-hb_m*(r1-hb_r)**2)
c         e_hb_wat=vdw_hb_ij(i)*factor_wat

c        Finction for HB in membrane, e_hb_mem, has a deeper minimum, since it is more diffcult to break the HB
c        in a hydrophobic environment.
         e_hb_mem=vdw_hb_ij(i)*factor_mem

c        e_hb_pair is defined in such a way, that
c        in the membrane (U_alfa=0) e_hb_pair=e_hb_mem
c        and in water    (U_alfa=1) e_hb_pair=e_hb_wat
         e_hb_pair=e_hb_wat*u_alpha(ica)*u_alpha(jca)+e_hb_mem*(1-u_alpha(ica)*u_alpha(jca))
c         e_hb_pair=exp(-hb_m*(hb_dis(i)-hb_r)**2)

C        Some write statements
c         Write(6,*)
c         write(6,*)"PAIR  ",i,atom_ires(i_alpha),atom_ires(j_alpha)
c         write(6,*)"Atoms ",atom_a1, atom_a2, atomz
c         write(6,'(''Distances and angles '',f8.2,2x,f8.2,2x,f8.2)'),r1,dis2, angle

c         write(6,*)"DEBUG: original hb_ij   ",vdw_hb_ij(i)
c         write(6,*)"DEBUG: fact wat, exp  ",factor_wat, exp(-hb_m*(r1-hb_r)**2)
c         write(6,*)"DEBUG: fact mem   ",factor_mem
c         write(6,*)"DEBUG: e_hb_wat e_hp_mem",e_hb_wat, e_hb_mem
c         write(6,*)"DEBUG: u_alpha(ica)*u_alpha(jca)  ",u_alpha(ica)*u_alpha(jca)
c         write(6,*)"DEBUG: e_hb_pair ",e_hb_pair
c         Write(6,*)
 
c         e_hb_pair=e_hb_pair*1000
c         write(6,*)"DEBUG: e_hb_pair at x",e_hb_pair
c         write(6,*)"DEBUG: e_hb_pair at x+0.001",exp(-hb_m*(hb_dis(i)+0.001-hb_r)**2)
c        The following 2 lines are for debug purposes:
c         e_hb_pair=e_hb_mem
c         e_hb_pair=e_hb_wat*0.5*0.5+e_hb_mem*(1-0.5*0.5)
         ehbpp = ehbpp + e_hb_pair

c         Write(6,*)
c         if(mod(istep,log_write_fq).eq.0.or.istep.eq.nsteps) then
c            write(6,'('' HB O...H'',i14,i24)')hb_i(i),hb_j(i)
c            write(6,'('' CA_i from CA-C=O:'',i5,''  CA_j from H-N-CA:'',i5)'),i_alpha,j_alpha
c         endif


c        ################################################################################################
c        Force calculation
c        ################################################################################################
         df=0.0

         i3   = hb_i(i)*3-3
         j3   = hb_j(i)*3-3
         dx1  = x(i3+1)-x(j3+1)
         dx2  = x(i3+2)-x(j3+2)
         dx3  = x(i3+3)-x(j3+3)

c        E_hb is a function of
c        E_hb(E_hb_wat(E_hbij(r_H-O),r_H-O),E_hb_mem(E_hbij(r_H-O)),
c             u_alpha_i(theta_i(Nnp(r_CA-CA),Nmem(r_CA-CA))),
c             u_alpha_j(theta_j(Nnp(r_CA-CA),Nmem(r_CA-CA))))
c        We need to get all the partial derivatives with respect to r_H-O,
c        because the force is applied on these two atoms.
c        Since r_CA-CA is independent on r_H-O, dr_CA-CA/dr_H-O=0!!!
c        (there is no one r_CA-CA for a particular r_H-O)

c         write(6,*)"DEBUG: hb_pair",i,"dist",r1,hb_dis(i)

         dE_hb_dE_hb_wat=u_alpha(ica)*u_alpha(jca)
         dE_hb_dE_hb_mem=1-u_alpha(ica)*u_alpha(jca)
c         dE_hb_du_alpha(ica)=(e_hb_wat-e_hb_mem)*u_alpha(jca)
c         dE_hb_du_alpha(jca)=(e_hb_wat-e_hb_mem)*u_alpha(ica)

         dE_hb_wat_dE_hbij=factor_wat
         dE_hb_mem_dE_hbij=factor_mem

         dE_hb_wat_dr=-2*hb_m*(r1-hb_r)*exp(-hb_m*(r1-hb_r)**2)
c         write(6,*)"DEBUG: hb_m hb_dis hb_r dE_hb_wat_dr",hb_m,hb_dis(i),hb_r,dE_hb_wat_dr
c         dE_hb_wat_dr=0

c        Use the saved value for the dE_hbij_dr (from nonbond.f), taking 
c        into account that df_vdw_hb_ij was devided by r already.
         dE_hbij_dr=df_vdw_hb_ij(i)*r1

c         if (hb_dis(i).le.hb_r) then
c            dE_hbij_dr=0
c         else
c            dE_hbij_dr=2*hb_gas_a*hb_gas_m*(hb_dis(i)-hb_gas_r)*exp(-hb_gas_m*(hb_dis(i)-hb_gas_r)**2)
c         endif

c         do k=ica,jca
c            if (abs(theta(k)).ge.max_theta) then
c               du_alpha_dtheta(k)=0
c            else
c               du_alpha_dtheta(k)=-2*a_theta_u*(abs(theta(k))-max_theta)
c     $                               *exp(-a_theta_u*(abs(theta(k))-max_theta)**2)
c            endif
c         enddo

c         dtheta_dNnp=1/N_p_max
c         dtheta_dNmem=1/N_mem_max

c         dNnp_dr=0.0
c         if (is_sim(jca).eq.3) then
c            if (r1.le.r_np) then
c               dNnp_dr=0
c            else
c               dNnp_dr=-2*a_np*(r1-r_np)*exp(-a_np*(r1-r_np)**2)
c            endif
c         endif
c
c        dNmem_dr=0
c        if (is_sim(jca).eq.4) then
c           if (r1.le.r_mem_ngb_cutoff) then
c              dNmem_dr=0
c           else
c              dNmem_dr=-2*a_mem*(r1-r_mem_ngb_cutoff)*exp(-a_mem*(r1-r_mem_ngb_cutoff)**2)
c           endif
c        endif

c        dtheta_dr=dtheta_dNnp*dNnp_dr+dtheta_dNmem*dNmem_dr

c       Now let's calculate dE_hb_dr

c         df=dE_hb_dE_hb_wat*(dE_hb_wat_dE_hbij*dE_hbij_dr+dE_hb_wat_dr)
c     $      +dE_hb_dE_hb_mem*dE_hb_mem_dE_hbij*dE_hbij_dr
c     $      +dE_hb_du_alpha(ica)*du_alpha_dtheta(ica)*dtheta_dr
c     $      +dE_hb_du_alpha(jca)*du_alpha_dtheta(jca)*dtheta_dr

         df=dE_hb_dE_hb_wat*(dE_hb_wat_dE_hbij*dE_hbij_dr+dE_hb_wat_dr)
     $      +dE_hb_dE_hb_mem*dE_hb_mem_dE_hbij*dE_hbij_dr
c         df=dE_hb_wat_dr


         df=df/r1
c         df=df*1000

c         write(6,*)"DEBUG: df for HB",df

c         write(6,*)"DEBUG: Applying Ehb force on atoms",hb_i(i),hb_j(i),df

         d(i3+1)=d(i3+1)+df*dx1
         d(i3+2)=d(i3+2)+df*dx2
         d(i3+3)=d(i3+3)+df*dx3

         d(j3+1)=d(j3+1)-df*dx1
         d(j3+2)=d(j3+2)-df*dx2
         d(j3+3)=d(j3+3)-df*dx3

      enddo  ! END of DO cycle for all HB pairs
      HydrogenBond_CG = ehbpp

      return
      end


      subroutine number_of_Calpha_neighbours
      implicit none

      include 'include/parameter.dc'
      include 'include/coord.h'
      include 'include/seq.h'
      include 'include/cbkv_gap.h'
      include 'include/Mdyn.h'
      include 'include/task.h'
      include 'include/derivative.h'
      include 'include/symbols.h'
      include 'include/stepsize.h'
      include 'include/lib.h'
      include 'include/cg_ion_hyd_pol_constants.h'
      include 'include/dgm_pt.h'
      include 'include/filename.h'
      include 'include/Mmol.h'
 
      real*8 r1, t_np_ngb, t_mem_ngb, np_sim_ngb, mem_sim_ngb
      integer i, j, ihit, jhit, ii, jj
   
      do 40 i=1,numres ! loop over all but membrane residues 
          t_np_ngb = 0.0
          t_mem_ngb = 0.0
          if(is_sim(i).eq.4) goto 40  !go to the next residue  if current is meb
          ! get C-aplha atom of the residue
          ihit=0
          do ii=np(i)+1,np(i+1)
              if(iac_name(iac(ii)).eq.'DY') goto 40
              if(bkvcode(ii)(1:2).eq.'CA') then
                  ihit=ii
              endif
          enddo
          if(ihit.eq.0) goto 40
          do 30 j=1,numres  ! loop over all residue and calulate # heighbours
              if(i.eq.j) goto 30
              ! gre C-alpha atom of the j residue or CB of MEB
              jhit=0
              do jj=np(j)+1,np(j+1)
                 if(iac_name(iac(jj)).eq.'DY') goto 30
                 if(bkvcode(jj)(1:2).eq.'CA'.or.resseq(j).eq.'MEB') then
                    jhit=jj
                 endif
              enddo
              if(jhit.eq.0) goto 30
              r1 = 0.0
              call dist_atoms(ihit,jhit,r1)
              if(is_sim(j).eq.3) then                   !j is nonpolar
                  if(r1.le.r_np) then
                     np_sim_ngb=1.0
                  else
                     np_sim_ngb=exp(-a_np*(r1-r_np)**2)
c                      np_sim_ngb=1.0
                  endif
                  t_np_ngb=t_np_ngb+np_sim_ngb
c                  write(6,'(''for pair '',2i5,'' Nnp'',f10.3)'),ihit,jhit,np_sim_ngb
              endif
              if(is_sim(j).eq.4) then               !j is membrane
                  if(r1.le.r_mem_ngb_cutoff) then
                     mem_sim_ngb=1.0
                  else
                     mem_sim_ngb=exp(-a_mem*(r1-r_mem_ngb_cutoff)**2) 
                  endif
                  t_mem_ngb=t_mem_ngb+mem_sim_ngb
              endif
 30       continue
          Nnonpol_CA(i) = t_np_ngb
          Nmem_CA(i)    = t_mem_ngb
 40   continue
      end
