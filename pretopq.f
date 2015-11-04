
c     the Id keyword will be substituted by CVS with:
c     path of RCS file, rev. number, date (UTC), author, state, locker
c
c     $Id: pretop.f,v 1.7 2005/04/04 21:04:28 ztchu Exp $
c

c     this subroutine reads in user's amino_lib, parm_lib and pdb files
c     to make a topology file for the strcuture of the system in pdb file
c     it also present all the options for dealing with the pdb file about
c     the H atoms, trimming for folding, using AMBER parameter set, using
c     non-pdb coordinates files, using alternative HIS (HIE) conformation,
c     etc.

      subroutine pretopq(inppdb)
      implicit none
      include 'include/parameter.dc'
      include 'include/task_general.h'
      include 'include/task.h'
      include 'include/bdparm.h'
      include 'include/seq.h'
      include 'include/coord.h'
      include 'include/lib.h'
      include 'include/struct.h'
      include 'include/filename.h'
      include 'include/trim.h'
      include 'include/pretop.h'
      include 'include/water.h'

      integer i,i1,i2,minus1,hit_h
      integer nword,word_pos(MAXWORD)
      real*8  e_bad
      character*(MXCHARACT) tmpline,tmplineup,xyzfile,topo_filename0,inppdb
      logical find_in_string,lcomment
      data rtrim/30.0/
      data alfa_beta_a, alfa_beta_f/109.5,60.0/
c...........................................................
      debug_lib=.false.
      lcheck_neu=.false.
      lcheck_parm=.false.

 3    continue

      tmpline=inppdb

      call process_line(tmpline,tmplineup,nword,word_pos,lcomment,'PDB')

      if (lcomment) goto 3
      if(tmplineup(word_pos(1):word_pos(2)).eq.'NEXT') then
         write(6,'(/,'' WARNING: You have not entered a coordinate file,'',
     $             '' many keywords may not work'',/)')
         l_pdb=.false.
         open(10,file=path_of_lib,status='old',form='formatted',err=9009)
         write(6,'(/,'' reading library data from: '',/,1x,A)')path_of_lib
         call readlib(0)      
         call parmset
         return
      endif
      l_pdb=.true.

c     for checking amino library problems if user added new entry in it:
      if(tmplineup(word_pos(3):word_pos(4)).eq.'CHECK_LIB') then
         debug_lib=.true.
         i_lib=1
         if(nword.gt.2) read (tmpline(word_pos(5):),*,err=663) i_lib
      endif

c     for checking neutral group problems in user's newly created entry
c     in amino_lib:
      if(tmplineup(word_pos(3):word_pos(4)).eq.'CHECK_NEUTRAL')
     $     lcheck_neu=.true.

c     for checking parameter problems if user edited parm_lib:
      if(tmplineup(word_pos(3):word_pos(4)).eq.'CHECK_PARM')
     $     lcheck_parm=.true.

      if(tmplineup(word_pos(1):word_pos(2)).eq.'QUIT'.or.
     $   tmplineup(word_pos(1):word_pos(2)).eq.'X'.or.
     $   tmplineup(word_pos(1):word_pos(2)).eq.'EXIT') stop ' SEE YOU LATER'

      if(tmplineup(word_pos(1):word_pos(2)).eq.'HELP')then
         if(tmplineup(word_pos(3):word_pos(4)).eq.'LEVEL')then
            call help_level(tmplineup(word_pos(5):word_pos(6)))
         else if(tmplineup(word_pos(3):word_pos(4)).eq.'ALPHABETIC')then
            call alpha_keywords(0)
         else if(tmplineup(word_pos(3):word_pos(4)).eq.'CREATE_KEYWORDS_LIST')then
            call alpha_keywords(1)
            stop
         else if(tmplineup(word_pos(3):word_pos(4)).eq.'HIDDEN')then
            call help_level('HIDDEN')
         else if(tmplineup(word_pos(3):word_pos(4)).eq.'TASK')then
            call help_task
         else
            call help(tmplineup(word_pos(3):word_pos(4)))
         endif
         goto 3
      endif

      if(tmplineup(word_pos(1):word_pos(2)).eq.'ALL_KEYWORDS')then
         call all_keywords
      endif

      pdbfile=tmpline(word_pos(1):word_pos(2))

      lup_big_ion=.true.
      lamber=.false.
      lxyz=.false.
      ikeep_h=0
      iw_from=0
      minus1=0
      ltrim=.false.
      lfold=.false.
      lfold2=.false.
      ltop=.false.
      nkeephis=0
      n_ion_atom=0

      if (nword.gt.1) then

         if(index(tmplineup,'KEEPALLH').ne.0) then
               write (6,*) 'Keeping ALL Hydrogens from the PDB'
               write(6,'('' Note: you need also to use "keepallh" for all the'',
     $                   '' subsequent runs'')')
               ikeep_h=2  
         elseif(index(tmplineup,'KEEPH ').ne.0) then
               write (6,*) 'Keeping Good Hydrogens from the PDB'
               write(6,'('' Note: you need also to use "keeph" for all the'',
     $                   '' subsequent runs'')')
               ikeep_h=1
         elseif(index(tmplineup,'KEEPNOH').ne.0) then
            write (6,*) 'Rejecting ALL Hydrogens from the PDB and create',
     $                  ' them by molaris'
            write(6,'('' Note: you need also to use "keepnoh" for all the'',
     $                   '' subsequent runs'')')
               ikeep_h=-1
         endif

         if(index(tmplineup,'XYZ_FILE').ne.0) then
            xyzfile=tmpline(word_pos(3):word_pos(4))
            write(6,'(/,'' Reading coordinate file (not in pdb format) from:'',
     $                1x,a)')xyzfile 
            open(11,file=xyzfile,status='old',form='formatted',err=9002)
            call read_xyz(xyzfile)
            lxyz=.true.
            ikeep_h=2  
         endif

         if(index(tmplineup,'AMBER').ne.0) then
            write (6,'(/,'' Using AMBER force field parameters'',/)')
            lamber=.true.
            lemp=index(path_of_lib,'.lib')-1
            do i=lemp,1,-1                              
               if(path_of_lib(i:i).eq.'/')then
                  lemp0=i
                  goto 4
               endif
            enddo
 4          path_of_lib=path_of_lib(1:lemp0)//'amber_amino.lib'
         endif

         if(tmplineup(word_pos(3):word_pos(4)).eq.'KEEPHIS') then
            nkeephis=-1               !keep all HIS or HIE from pdb

            if(nword.gt.2)then        !keep specified HIS or HIE in pdb
               nkeephis=0
               do i=5, nword*2, 2
                  nkeephis=nkeephis+1
                  read(tmpline(word_pos(i):),*,err=663)ikeephis(nkeephis)
               enddo
            endif

            if(nkeephis.gt.MXST)then
               write(6,'('' # of residue his to be kept >'',i2)')MXST
               goto 663
            endif
         elseif(tmplineup(word_pos(5):word_pos(6)).eq.'KEEPHIS') then
            nkeephis=-1               !keep all HIS or HIE from pdb

            if(nword.gt.3)then        !keep specified HIS or HIE in pdb
               nkeephis=0
               do i=7, nword*2, 2
                  nkeephis=nkeephis+1
                  read(tmpline(word_pos(i):),*,err=663)ikeephis(nkeephis)
               enddo
            endif

            if(nkeephis.gt.MXST)then
               write(6,'('' # of residue his to be kept >'',i2)')MXST
               goto 663
            endif
         elseif(tmplineup(word_pos(7):word_pos(8)).eq.'KEEPHIS') then
            nkeephis=-1               !keep all HIS or HIE from pdb

            if(nword.gt.4)then        !keep specified HIS or HIE in pdb
               nkeephis=0
               do i=9, nword*2, 2
                  nkeephis=nkeephis+1
                  read(tmpline(word_pos(i):),*,err=663)ikeephis(nkeephis)
               enddo
            endif

            if(nkeephis.gt.MXST)then
               write(6,'('' # of residue his to be kept >'',i2)')MXST
               goto 663
            endif
         endif

         if(nkeephis.ne.0) write(6,'('' Note: you need also to use "keephis"'',
     $                               '' for all the subsequent runs'')')

         if(tmplineup(word_pos(5):word_pos(6)).eq.'TRIM') then 
               ltrim=.true.
               read (tmpline(word_pos(7):),*,err=663)(trimcent(i),i=1,3),rtrim 
               write (6,*) 'The PDB will be trimmed'
         elseif(tmplineup(word_pos(3):word_pos(4)).eq.'TRIM') then 
               ltrim=.true.
               read (tmpline(word_pos(5):),*,err=663)(trimcent(i),i=1,3),rtrim 
               write (6,*) 'The PDB will be trimmed'

         elseif(tmplineup(word_pos(3):word_pos(4)).eq.'KEEPENZW') then
            read (tmpline(word_pos(5):word_pos(6)),*,err=663) iw_from
            write (6,*) 'keep enzymix water from the PDB as solvent'

         elseif(tmplineup(word_pos(3):word_pos(4)).eq.'FOLD_ALFA_BETA') then 
            lfold=.true.
            ltrim=.true.
            read (tmpline(word_pos(5):),*,err=663)(trimcent(i),i=1,3),rtrim
            write(6,'(/'' residues in the above pdb file outside the'',
     $              '' folding radius: '',f5.1,''A will be '',/,
     $              '' trimed to C-alpha for the main chain and united'',
     $              '' C-beta for the side chain'',/)')rtrim

            if(nword.gt.6) then
               read (tmpline(word_pos(13):),*,err=663) alfa_beta_a,alfa_beta_f
            endif

         elseif(tmplineup(word_pos(3):word_pos(4)).eq.'FOLD2') then 
            lfold2=.true.
            ltrim=.true.
            read (tmpline(word_pos(5):),*,err=663)(trimcent(i),i=1,3),rtrim
            write(6,'(/,'' all the sidechains of residues in the above pdb'',
     $              '' file outside the folding'',/,'' radius: '',f5.1,
     $              ''A will be trimed to united C-beta atoms'',/)')rtrim
            write(6,'('' If there are DNA residues, they will be trimed'',
     $                '' to mainchain of -C3*-C4*-P-'',/)')

         elseif(tmplineup(word_pos(1):word_pos(2)).eq.'TOP_IN') then
               topo_filename=tmpline(word_pos(3):word_pos(4))
               call topoin(.true.)
               ltop=.true.
               goto 300
         endif
      endif
c----------------------------------------------------------------
      open(10,file=path_of_lib,status='old',form='formatted',err=9009)
      write(6,'(/,'' reading library data from: '',/,1x,A)')path_of_lib
      call readlib(0)
                                          !check pdb file name extension
      lemp0=1
      lemp=0
      
      lemp=index(pdbfile,'.pdb')-1
c      lemp=index(pdbfile,'.ent')-1

      if (lemp.le.0) then
         write(6,'('' The coordinate file name should have '',
     $             ''".pdb" as extension, e.g., test.pdb '',/)')
         goto 666
      end if

      do i=lemp,1,-1                              
         if(pdbfile(i:i).eq.'/')then
            lemp0=i+1
            goto 5
         endif
      enddo

 5    continue

      pdbfile_head=pdbfile(lemp0:lemp)             !"peeled" name of the pdb 
                                                   !or ent file

c----------------------------------------------------READING PDB FILE
      open(11,file=pdbfile,status='old',form='formatted',err=9000)      

c     check if there is any H atom in the pdb file:
      call check_h_exist(hit_h)
      if(hit_h.eq.0) then
         if(ikeep_h.eq.2) then
            write(6,'('' WARNING: there is no H atoms in the pdb file,'',
     $                '' the keyword "keepallh" is not'',/,
     $                ''          meaningful'')')
            ikeep_h=0
         else if(ikeep_h.eq.1) then
            write(6,'('' WARNING: there is no H atoms in the pdb file,'',
     $                '' the keyword "keeph" is not'',/,
     $                ''          meaningful'')')
            ikeep_h=0
         else if(.not.lxyz) then
            write(6,'('' NOTICE: there is no H atoms in the pdb file,'',
     $             '' molaris will create all the H'',/,
     $             ''         atoms'')')
         endif
      endif

c     if there is any new residue in the pdb file, but no entry in amino lib
c     try automatically to make a new new entry in amino lib, for complex
c     residue (molecule it may fail)
      call make_alib                                      

c     reopen the pdb file, in case after creating the new entry in amino lib
c     and having recreated pdb file:
      open(11,file=pdbfile,status='old',form='formatted') 
                                                          
c     read in the parm lib
      call parmset

c     read in pdb file
      call read_brk   

c     checking if parameter misssing
      call chk_nbparm

c----------------------------------------------------WRITING TOPOLOGY FILE
c     generating the topology structure file:

      topo_filename=path_of_out(1:ndir)//pdbfile_head(1:lemp-lemp0+1)//'_0.top'
      call topology_driver(.false.)
      original_topf=topo_filename

      topo_filename=path_of_out(1:ndir)//pdbfile_head(1:lemp-lemp0+1)//'.top'

      call write_topology(0) 

      write(6,'(/,'' Created a topology file: '',/,1x,a,/)') topo_filename   
      call topoin(.true.)                               !read the topology

c     since the system cp depends on the machine buffer size, so we call 
c     topology_driver the 2nd time to get *_.top and *_0.top be written out
c     right away
c      call system('cp '//topo_filename//topo_filename0)

c      call flush

c-----------------------------------------------------------------------
 300  e_bad=50.0
c      call resbond(0,e_bad)                    !check for bad bonds
c      if(numres.le.300) call check_close_atom  !check for bad nonbond atom pair
      call write_sequence(6)
      call center(0,0)

c     the xyz coordinate file is used for running evb/qmmm only, 
c     so automatically create an evb.dat file
      if(lxyz) call create_evb_dat
      print *
      return
c----------------------------------------------------------------
 9009 write(6,9109)path_of_lib
 9109 format(/,1x,'amino.lib is not found in the specified directory',/,a)
      goto 666

 9000 write(6,9100)
 9100 format(/,1x,
     $     'The above pdb file is not found (or a keyword is missing)',/,
     $     ' Reinput the right pdb file name: ',/$,' ')
      goto 665

 9002 write(6,9200)
 9200 format(/,1x,
     $     'The above xyz file is not found',/)
      goto 665

 9001 write(6,9101)
 9101 format(/,1x,'The .top file is not found ',/)
      goto 665

 663  write (6,*) 'There is an error in the line',tmpline
 665  if (interf.eq.0) goto 3
 666  call killme('pretopq')
      end
c----------------------------------------------------------------
