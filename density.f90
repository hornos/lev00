subroutine density(spin)
!..................................................................
!    It prints out the density from a file along a line or in a
! plane to the file 'out.dat_[i]', where 'i' is a character
! which distinguishes between different plots (if any)
!..................................................................
! DIRC - direct lattice vectors
! RECC - reciprocal lattice vectors (with 2*pi)
! VOLC - unit cell volume
!..................................................................
! BCELL - reciprocal lattice vectors (without 2*pi)
!..................................................................
use param
use atoms
use code
implicit none
real*8, dimension(:,:,:),allocatable :: GRID
real*8, parameter :: pi=3.1415927d0
real*8, parameter :: tiny = 0.00001
integer ijk,iErr,lenght
real*8 Dip(3),Quadr(3,3),totdens
character cha,cha2*2
character filen*12,filen1*12,name*6,item*2
logical Dip_Done, Quadr_Done,spin
!
!........... do BCELL = RECC/(2*pi) - it is used in transforming
!            grid point coordinates
!
   BCELL=RECC/(2*pi)

!........... intitalise the density
!
   allocate(GRID(NGX,NGY,NGZ))
   GRID=0.0
!___________ set core charges to 0
   allocate(Z_atom(NSPEC)) 
   Z_atom=0.0
!
!........... read the density to GRID(iX,iY,iZ).
!       The total density 'totdens' is calculated just to take care of
!       the calculations.
!
   totdens=0.0

!_______________ VASP input
   if(Which_Code.eq.'  VASP') then

      if(.not. spin) then
1        write(*,*) '....... Choose the file to be read in: ......'
         write(*,*)
         write(*,*) '  C. Total electron density CHGCAR'
         write(*,*) '  P. Partial electron density PARCHG'
         write(*,*) '  L. Electrostatic potential LOCPOT'
         write(*,*) '  Q. Quit: do not read any file'
         write(*,*)
         write(*,*)'------------>'
         read(*,'(a)',err=1) item
         if(item.eq.'C') then
            name='CHGCAR'
         else if(item.eq.'P') then
            name='PARCHG'
         else if(item.eq.'L') then
            name='LOCPOT'
         else if(item.eq.'Q') then
            return
         else
            go to 1
         end if
      else
         name='CHGCAR'
      end if
      write(*,*)'Reading in the density from '//name//' ...'
      filen1='  '//name//'.new'
      call vasp_dens(grid,totdens,name,spin,iErr)
      if(iErr.eq.1) go to 150

!_______________ SIESTA input
   else if(Which_Code.eq.'SIESTA') then
      write(*,*) 'Reading in the density from '//trim(seed)//'.RHO ...'
      filen1=' job.RHO.new'
      call siesta_dens(grid,totdens,spin,iErr)
      if(iErr.eq.1) go to 150
   end if
   close (1)
   write(*,*) '.....> Total charge = ',totdens/NPLWV,' <.....'
   write(*,*)'Done!'
!
!.......... At this stage we have grid(iX,iY,iZ) which contains the
! total charge density of the system
!..................................................................
!
!....................................................................
!............ General part: let us plot just once ...................
!....................................................................
!.... ijk - counts different cycles of calculations (not more than 9),
!           i.e. different plots
!
   ijk=1
   Dip_Done=.false.
   Quadr_Done=.false.
!
!............ name of the file for the output
!
2  if(ijk.le.9) then
      write(cha,'(i1)') ijk
      filen='out.dat_'//cha
      lenght=9
   else if(ijk.le.99) then
      write(cha2,'(i2)') ijk
      filen='out.dat_'//cha2
      lenght=10
   else
      write(*,*)'DENSITY: You cannot trial my patience so much!'
      go to 100
   end if
!
!____________ choose between a line, plane or charge
!
   write(*,*)'...... Choose between line or plane ......'
   write(*,*)
   write(*,'(a33,i2)')'     NUMBER OF THE CURRENT PLOT: ',ijk
   write(*,*)' pL. Plot density along a line'
   write(*,*)' pP. Plot density in a plane'
   write(*,*)' CS. Amount of charge inside a sphere'
   write(*,*)' Ex. Exploration of the density'
   if(Dip_Done) then
      write(*,*)' DM. Dipole moment <== DONE!'
   else
      write(*,*)' DM. Dipole moment'
   end if
   if(Quadr_Done) then
      write(*,*)' QM. Quadrupole moment <== DONE!'
   else
      write(*,*)' QM. Quadrupole moment'
   end if
   if(Dip_Done.and.Quadr_Done) write(*,*) &
       ' vP. Get density via point charges: match moments & potential'
   write(*,*)' cA. Cutting atoms out of the density'
   write(*,*)' wD. Write non-zero density as '//filen1
   write(*,*)' gO. Write density in gOpenMol cube format'
   write(*,*)' mD. Get density via point charges: match density'
   write(*,*)' Sf. Transform the charge density for a shifted system'
   write(*,*)' TH. STM image (Tersoff-Hamann)'
   write(*,*)'  Q. Return to the previous menu'
   write(*,*)
   write(*,*)'------------>'
   read(*,'(a)',err=3) item
   if(item.eq.'pL') then
      call line(grid,DIRC,BCELL,VOLC,filen,lenght)
      ijk=ijk+1
   else if(item.eq.'pP') then
      call plane(grid,DIRC,BCELL,VOLC,filen,lenght)
      ijk=ijk+1
   else if(item.eq.'CS') then
      call charge_sph(grid,DIRC,BCELL,VOLC,totdens,filen,lenght)
      ijk=ijk+1
   else if(item.eq.'Ex') then
      call max_dens(grid,DIRC,BCELL,VOLC,totdens)
   else if(item.eq.'DM') then
      call dipole(grid,totdens,filen,lenght,Dip)
      Dip_Done=.true.
   else if(item.eq.'QM') then
      call quadrpl(grid,Quadr)
      Quadr_Done=.true.
   else if(item.eq.'vP'  .and. Dip_Done .and. Quadr_Done) then
      call simulate(grid,Dip,Quadr,filen,lenght)
   else if(item.eq.'cA') then
      call cut_atoms(grid,totdens)
   else if(item.eq.'wD') then
      write(*,*)'Writing non-zero elements to a new density file ...'
      call write_dens(grid)
      ijk=ijk+1
   else if(item.eq.'gO') then
      call for_gOpenMol(grid)
   else if(item.eq.'mD') then
      call simul_box(grid)
   else if(item.eq.'Sf') then
#ifdef HORNOS
#ifdef _OPENMP
      call omp_shift_charge(grid,DIRC,BCELL)
#else
      call shift_charge(grid,DIRC,BCELL)
#endif
#endif
   else if(item.eq.'TH') then
      call stm_TH(grid)
   else if(item.eq.'Q') then
      go to 100
   else
      go to 3
   end if
   go to 2
3  write(*,*) "Incorrect item number! Try again!"
   go to 2
!
!............ finish
!
100 deallocate(GRID)
    deallocate(Z_atom)
    return
!
!........... errors
 150  write(*,*)'FATAL! Error while opening '//filen1(1:8)//' file!'
      deallocate(GRID)
      deallocate(Z_atom)
end subroutine density

subroutine simul_box(grid)
!....................................................................
!   A box is chosen inside the cell and the density in the box is
! distributed with point charges using a small grid in the box
!....................................................................
use param
use atoms
use menu
implicit none
integer, parameter :: NN0=10
real*8, parameter :: tiny=0.0001
real*8 GRID(NGX,NGY,NGZ),R(3),x(3),R1(3),denval,ch
real*8 Center(3),Sides(3,3),Direct(3,3),Face(3),RecipS(3,3)
real*8 vect(3,3),vecB(3,3),RecipC(3,3),corner(3),PosCh(3),R2(3)
real*8 RecipD(3,3),dip(3),VolBox,da,VolCell,TotCh,chdens,Qlarge,Qsmall
real*8 Pot_Diff,dipm,Vol,rCharge,dV,factor,den,rCh,dQ,Qleft,Qright,ar,a
integer ngrid(3),NRs(3),iQuit,i,j,nChrg,iPnt,Nchrg1,NN2,NN,k
integer ijk,j0,ix,iy,iz,n1,n2,n3,k1,k2,k3,jj,Ncell,NunitCell,nat,i1,i2,i3
real*8, dimension(:,:,:),allocatable :: Poten,Poten1
character iask,cha1,cha2*2
logical Yes_Do,Yes_Pot,Yes_Comp, Yes_Dip
data Center/3*0.0/,ngrid/3*1/,NunitCell/1/
data Face/3*1.0/,TotCh/0.0/,NRs/3*10/,NN2/3/
!
!________ default for the directions of the box sides: along lattice vectors
      do i=1,3
         Direct(i,1:3)=DIRC(i,1:3)
         call normalize(Direct(i,1),Direct(i,2),Direct(i,3))
      end do
      allocate(Poten(-NN0:NN0,-NN0:NN0,-NN0:NN0))
      allocate(Poten1(-NN0:NN0,-NN0:NN0,-NN0:NN0))
!
!................. start the main menu
!
      Yes_Do=.false.
      Yes_Pot=.false.
      Yes_Comp=.false.
      Yes_Dip=.false.
1     iQuit=0
      write(*,*)'............MENU for SIMULATE in the BOX .............'
      write(*,*)'......... Change these parameters if necessary:.......'
      write(*,*)
      write(*,'(a)') '>>>>> Representation of results: through number of electrons'
      write(*,'(a)') '>>>>> Algorithm for the charge integration: <nonconserving>'
      
      write(*,'(a,3(f10.5,a))') '   1. The box center is at: (', &
           Center(1),',',Center(2),',',Center(3),')'

      write(*,'(a)')'   2. Directions of the box sides are along:'
      write(*,'(10x,a,3(f10.5,a))') &
           '1   ',Direct(1,1),',',Direct(1,2),',',Direct(1,3),')'
      write(*,'(10x,a,3(f10.5,a))') &
           '2   ',Direct(2,1),',',Direct(2,2),',',Direct(2,3),')'
      write(*,'(10x,a,3(f10.5,a))') &
           '3   ',Direct(3,1),',',Direct(3,2),',',Direct(3,3),')'

      do i=1,3
         do j=1,3
            Sides(i,j)=Direct(i,j)*Face(i)
         end do
      end do
      call BASTR(Sides,RecipS,VolBox,0)
      write(*,'(a,3(1x,f10.5))') &
          '   3. Lengths of the box sides (in A) are: ',(Face(i),i=1,3)
      corner(1)=Center(1)-0.5*(Sides(1,1)+Sides(2,1)+Sides(3,1))
      corner(2)=Center(2)-0.5*(Sides(1,2)+Sides(2,2)+Sides(3,2))
      corner(3)=Center(3)-0.5*(Sides(1,3)+Sides(2,3)+Sides(3,3))
      write(*,'(5x,a,3(1x,f10.5))') &
           '>>>> corner of the box is at ',(corner(i),i=1,3)

      nChrg=ngrid(1)*ngrid(2)*ngrid(3)
      write(*,'(a,3(1x,i3))') &
           '   4. The number of charges in each direction: ',(ngrid(i),i=1,3)
      write(*,'(5x,a,i5)')'>>>> total number of charges = ',nChrg
      do i=1,3
         da=Face(i)/ngrid(i)
         do j=1,3
            vect(i,j)=Direct(i,j)*da
         end do
      end do
      call BASTR(vect,RecipC,VolCell,0)
      write(*,'(5x,a,f10.5)')'>>>> the cell volume = ',VolCell

      write(*,'(a,3(1x,i3))') &
       '   5. The integration grid in each cell of the box: ', &
                                                   (NRs(i),i=1,3)
      do i=1,3
         da=Face(i)/(NRs(i)*ngrid(i))
         do j=1,3
            vecB(i,j)=Direct(i,j)*da
         end do
      end do

      write(*,'(a)') &
        '   6. Scan the box and integrate the charge (reference only):'
      if(TotCh.ne.0.0) then
         write(*,'(5x,a,f10.5)') &
              '>>>> total charge in the box = ',TotCh
         write(*,'(5x,a,i10)') &
              '>>>> total grid points in the box = ',iPnt
      end if
      if(Yes_Do) then
         write(*,'(a)') &
              '   7. Get point charges using the grid specified <= DONE!'
         write(*,'(5x,a,i5)')'>>>> total number of charges = ',Nchrg1
         write(*,'(5x,a,f10.5)')'>>>> total charge = ',chdens
         write(*,'(5x,2(a,f10.5))') &
              '>>>> charges between ',Qsmall,' and ',Qlarge
         write(*,'(a)') '   8. Show point charges'
         write(*,'(a)') '  88. Visualise point charges'
         write(*,'(a,i3)') &
              '   9. Number of unit cells to test the potential: ',NunitCell
         NN2=nint( (real(NunitCell)**(0.3333333)-1.)/2. )
         NN=2*NN2+1
         NunitCell=NN**3
         write(*,'(5x,3(a,i3))') &
              '>>>> the test box defined as : ',NN,'x',NN,'x',NN
         if(Yes_Pot) then
            write(*,'(a)') &
                 '  10. Compare potential with the previous one'
            if(Yes_Comp)write(*,'(5x,a,e12.6)')'>>>> error = ',Pot_Diff
         else
            write(*,'(a)') &
                 '  10. Calculate the potential at a set of points'
         end if
         if(Yes_Dip) then
            write(*,'(a,f10.5,a)') &
                 '  11. Calculated dipole moment = ',dipm,' e.A'
         else
            write(*,'(a,f10.5)') '  11. Calculate the dipole moment'
         end if
      else
         write(*,'(a)') &
              '   7. Get point charges using the grid specified'
      end if
      
      write(*,'(a)') '-------  G e n e r a l  s e t t i n g s ---------'
      write(*,'(a)')'  An. Coordinates are specified in: '//angstr
      write(*,'(a)') &
           '  Co. Show current atomic positions in fractional/Cartesian'
      write(*,'(a)')'------ L e a v e   t h e   m e n u -------------'
      write(*,'(a)')'   Q. Return to the previous menu'
      write(*,*)
      write(*,*)'------> Choose the item and press ENTER:'
      read (*,'(a)',err=1) cha2
!
![An]__________ choose the way how the coordinates are given
!
      IF(trim(cha2).eq.'An') THEN
         if(angstr.eq.'<Fractional>') then
            angstr='<Angstroms> '
         else if(angstr.eq.'<Angstroms> ') then
            angstr='<AtomNumber>'
         else if(angstr.eq.'<AtomNumber>') then
            angstr='<Fractional>'
         end if
!
![1]__________ give central point on the plane
!
      ELSE IF(trim(cha2).eq.'1') THEN
         WRITE(*,*)'Give the center of your box'
         call givepoint(Center(1),Center(2),Center(3),angstr)
         Yes_Do=.false.
         TotCh=0.0
!
![2]__________ change directions of the box sides
!
      ELSE IF(trim(cha2).eq.'2') THEN
 19      do i=1,3
 11        WRITE(*,'(a,i1,a)') &
                      'Give direction ',i,' for the box as X,Y,Z'
           read(*,*,err=11) (Direct(i,j),j=1,3)
           call normalize(Direct(i,1),Direct(i,2),Direct(i,3))
         end do
         call BASTR(Direct,RecipD,Vol,0)
         if(Vol.lt.tiny) then
           write(*,*)'ERROR! Collinear directions!'
           go to 19
         end if
         Yes_Do=.false.
         TotCh=0.0
!
![3]__________ choose lengths of the box sides
!
      ELSE IF(trim(cha2).eq.'3') THEN
 12      WRITE(*,'(a)')'Give lengths of the box sides (in A)'
         read(*,*,err=12) (Face(j),j=1,3)
         do j=1,3
           Face(j)=abs(Face(j))
           if(Face(j).lt.tiny) go to 12
         end do
         Yes_Do=.false.
         TotCh=0.0
!
![4]__________ choose the grid in the box: this grid will
!           determine the distribution of point charges by putting
!           by one charge in every of its cells
!
      ELSE IF(trim(cha2).eq.'4') THEN
 13      WRITE(*,'(a)') &
           'Specify the number of charges in each direction:'
         read(*,*,err=13) (ngrid(i),i=1,3)
         do i=1,3
           if(ngrid(i).lt.1) go to 13
         end do
         Yes_Do=.false.
!
![5]__________ number of grid points for each small cell in the box
!           (is used for integration)
!
      ELSE IF(trim(cha2).eq.'5') THEN
 996     write(*,*)'Give the integration grid for every CELL:'
         read(*,*,err=996) (NRs(i),i=1,3)
         do i=1,3
           if(NRs(i).lt.2) go to 996
         end do
         Yes_Do=.false.
!
![6]__________ integrate the charge in the box (for reference only)
!
      ELSE IF(trim(cha2).eq.'6') THEN
         write(*,*)'Using conserving algorithm ...'
!
!____________ for statistics (10%, 20%, ...)
         ijk=0
         j0=NPLWV/10
!
         iPnt=0
         rCharge=0.0
         do iZ=0,NGZ-1
            do iY=0,NGY-1
               do iX=0,NGX-1
!_____________ ask whether the point (iX,iY,iZ) is inside the box
                  call ask_box(iX,iY,iZ,Center,Sides,RecipS,DIRC,iask)
                  if(iask.eq.'y') then
                     iPnt=iPnt+1
                     rCharge = rCharge + grid(iX+1,iY+1,iZ+1)
                  end if
!______________ statistics
                  ijk=ijk+1
                  if(ijk/j0*j0.eq.ijk) &
                       write(*,'(a,i3,a)') '... done ',ijk/j0*10,' %'
               end do
            end do
         end do
         TotCh=rCharge/NPLWV
!
![7]_________ in the loop over all cells in the box: use the fine grid NRs()
!           to integrate the charge and find its position in the cell
!
      ELSE IF(trim(cha2).eq.'7') THEN
         write(*,'(a)') &
              'Opening the file <charges.sim> for simulation charges...'
         open(31,file='charges.sim',form='formatted',status='unknown')
!
!____________ position of the box in the original unit cell (its corner)
!
         dV=VolCell/(NRs(1)*NRs(2)*NRs(3))
         factor=dV/VOLC
         write(*,*)'Using nonconserving algorithm ...'
         chdens=0.0
         j=0
         Qsmall=10000.0
         Qlarge=0.0
         DO N1=0,ngrid(1)-1
            DO N2=0,ngrid(2)-1
               DO N3=0,ngrid(3)-1
!____________________ position of the cell in the box
                  x(1)= N1*vect(1,1)+N2*vect(2,1)+N3*vect(3,1)
                  x(2)= N1*vect(1,2)+N2*vect(2,2)+N3*vect(3,2)
                  x(3)= N1*vect(1,3)+N2*vect(2,3)+N3*vect(3,3)
!____________________ position of the cell in the original unit cell
                  R(1)=x(1)+corner(1)
                  R(2)=x(2)+corner(2)
                  R(3)=x(3)+corner(3)
!____________________ integrate the charge in the current cell (rCharge)
!                     and determine the position of the charge (PosCh)
                  rCharge=0.0
                  PosCh(1)=0.0
                  PosCh(2)=0.0
                  PosCh(3)=0.0
                  do k1=0,NRs(1)-1
                     do k2=0,NRs(2)-1
                        do k3=0,NRs(3)-1
                           R2(1)= k1*vecB(1,1)+k2*vecB(2,1)+k3*vecB(3,1)
                           R2(2)= k1*vecB(1,2)+k2*vecB(2,2)+k3*vecB(3,2)
                           R2(3)= k1*vecB(1,3)+k2*vecB(2,3)+k3*vecB(3,3)
                           R1(1)= R2(1)+R(1)
                           R1(2)= R2(2)+R(2)
                           R1(3)= R2(3)+R(3)
                           call reducn(R1,DIRC,BCELL)
                           call interpolate(R1,BCELL,denval,grid)
                           den=denval*factor
                           rCharge = rCharge + den
                           chdens=chdens + den
                           PosCh(1)=PosCh(1)+R2(1)*den
                           PosCh(2)=PosCh(2)+R2(2)*den
                           PosCh(3)=PosCh(3)+R2(3)*den
                        end do
                     end do
                  end do
!______________________ write the simulating charges to the file;
!                       find the smallest and the largest charge

                  if(rCharge.gt.0.0001) then
                     j=j+1
                     PosCh(1)=PosCh(1)/rCharge+R(1)
                     PosCh(2)=PosCh(2)/rCharge+R(2)
                     PosCh(3)=PosCh(3)/rCharge+R(3)
                     write(31,'(i5,5x,3(f10.5,x),5x,e16.10)') &
                          j,(PosCh(i),i=1,3), -rCharge
                     if(rCharge.ge.Qlarge) Qlarge=rCharge
                     if(rCharge.le.Qsmall) Qsmall=rCharge
                  end if
               END DO
            END DO
         END DO
         close (31)
         WRITE(*,*)'Done! The file <charges.sim> created!'
         Nchrg1=j
         Yes_Do=.true.
!
![8]__________ show point charges by 15 at a time
!
      ELSE IF(trim(cha2).eq.'8' .and. Yes_Do) THEN
         write(*,'(a)') &
              'Opening the file <charges.sim> for simulation charges...'
         open(31,file='charges.sim',form='formatted',status='old')
         j=0
 46      read(31,*,err=49,end=49) jj, (x(i),i=1,3),rCh
         j=j+1
         if(j.eq.16) then
            j=0
            write(*,*)'Press ENTER when ready ...'
            read(*,*)
         end if
         write(*,'(a,i5,a,f10.5,a,3(1x,f10.5))') &
              'ch(',jj,')= ',rCh, ' at ',(x(i),i=1,3)
         go to 46
 49      close (31)
         write(*,*)'Press ENTER when ready ...'
         read(*,*)
!
![88]__________ visualise point charges
!
      ELSE IF(trim(cha2).eq.'88' .and. Yes_Do) THEN

         dQ=(Qlarge-Qsmall)/11.0
         write(*,*)'Creating a clever input for Xmol ...'
         write(*,'(a)') &
              'Opening the file <charges.sim> for simulation charges...'
         open(31,file='charges.sim',form='formatted',status='old')
         write(*,'(a)') &
             'Creating the file <charges.xyz> for simulation charges...'
         open(33,file='charges.xyz',form='formatted',status='unknown')
         write(33,'(i10/)') Nchrg1
 86      read(31,*,err=89,end=89) j, (x(i),i=1,3),rCh
         do i=1,11
            Qleft=Qsmall+dQ*(i-1)
            Qright=Qleft+dQ
            if(i.eq.11) Qright=Qlarge+0.00001
            if(i.eq.1) Qleft=Qsmall-0.00001
            if(-rCh.ge.Qleft.and.-rCh.lt.Qright) then
               if(i.le.9) then
                  write(cha1,'(i1)') i
                  cha2=cha1//' '
               else
                  write(cha2,'(i2)') i
               endif
               write(33,'(a,5x,3(x,f10.5))')'LV'//cha2,(x(j),j=1,3)
               go to 86
            end if
         end do
         write(*,'(a,i5,x,f10.5)') 'ERROR! charge j=',j,rCh
         go to 86
 89      close (31)
         close (33)
!
![9]__________ the number of unit cells where the potential will be calculated
!           to assess the convergence
!
      ELSE IF(trim(cha2).eq.'9' .and. Yes_Do) THEN
         Ncell=(2*NN0+1)**3
 43      WRITE(*,'(a,i5,a)') 'Number of unit cell to check the '// &
                  'potential in (between 27 and ',Ncell,'):'
         read(*,*,err=43) NunitCell
         if(NunitCell.lt.27) go to 43
         NN2=nint( (real(NunitCell)**(0.3333333)-1.)/2. )
         if(NN2.gt.NN0) go to 43
         Yes_Pot=.false.
!
![10]__________ calculate/compare potential at the center of the neighbouring cells
!
      ELSE IF(trim(cha2).eq.'10' .and. Yes_Do) THEN
         if(Yes_Pot) then
            Poten1=0.0
         else
            Poten=0.0
         end if

         write(*,'(a)') 'Opening the file <charges.sim> for simulation charges...'
         open(31,file='charges.sim',form='formatted',status='old')
 56      read(31,*,err=59,end=59) j, (x(i),i=1,3),rCh

         do i1=-NN2,NN2
            do i2=-NN2,NN2
               do 37 i3=-NN2,NN2
                  R(1)=i1*DIRC(1,1)+i2*DIRC(2,1)+i3*DIRC(3,1)-x(1)
                  R(2)=i1*DIRC(1,2)+i2*DIRC(2,2)+i3*DIRC(3,2)-x(2)
                  R(3)=i1*DIRC(1,3)+i2*DIRC(2,3)+i3*DIRC(3,3)-x(3)
                  if(i1.eq.0.and.i2.eq.0.and.i3.eq.0) go to 37
                  ar=sqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
                  if(Yes_Pot) then
                     Poten1(i1,i2,i3)=Poten1(i1,i2,i3)+rCh/ar
                  else
                     Poten(i1,i2,i3)=Poten(i1,i2,i3)+rCh/ar
                  end if
37             end do
            end do
         end do

         go to 56
59       close (31)
         if(Yes_Pot) then
            Pot_Diff=0.0
            do i1=-NN2,NN2
               do i2=-NN2,NN2
                  do i3=-NN2,NN2
                     a=Poten(i1,i2,i3)-Poten1(i1,i2,i3)
                     Pot_Diff=Pot_Diff+abs(a)
                     Poten(i1,i2,i3)=Poten1(i1,i2,i3)
                  end do
               end do
            end do
            Yes_Comp=.true.
         end if
         Yes_Pot=.true.
!
![11]_________ calculate the dipole moment of the simulated charge distribution
!
      ELSE IF(trim(cha2).eq.'11' .and. Yes_Do) THEN
         Dip=0.0
!______________ useful input for the Madelung code (separate)
         nat=0
         do k=1,NSPEC
            do j=1,NspN(k)
               nat=nat+1
            end do
         end do
         open(41,file='mad.inp',form='formatted',status='unknown')
         write(41,*) Nchrg1+nat
!_________________ nuclear part first
 77      write(*,*)'Specify nucleii charges in the order of species:'
         read(*,*,err=77) (Z_atom(i),i=1,NSPEC)
         nat=0
         do k=1,NSPEC
            do j=1,NspN(k)
               nat=nat+1
!___________________(a) find the right image of this nuclei which is
!                       inside the box
               do n1=-1,1
                  do n2=-1,1
                     do n3=-1,1
                        x(1)=n1*DIRC(1,1)+n2*DIRC(2,1)+n3*DIRC(3,1)+TI(1,nat)
                        x(2)=n1*DIRC(1,2)+n2*DIRC(2,2)+n3*DIRC(3,2)+TI(2,nat)
                        x(3)=n1*DIRC(1,3)+n2*DIRC(2,3)+n3*DIRC(3,3)+TI(3,nat)
                        call ask_box2(x,Center,Sides,RecipS,iask)
                        if(iask.eq.'y') go to 33
                     end do
                  end do
               end do
               write(*,'(a,i5,a)') &
                    'ERROR: atom nat=',nat,' is not inside the box!'
               write(*,'(a)') 'Make sure the box is large enough!'
               go to 1
!___________________(b) calculate the contribution to the dipole moment
33             Dip(1)=Dip(1)+Z_atom(k)*(x(1)-Center(1))
               Dip(2)=Dip(2)+Z_atom(k)*(x(2)-Center(2))
               Dip(3)=Dip(3)+Z_atom(k)*(x(3)-Center(3))
               write(*,'(a,2i4,a,i4,a,3(x,f10.5))')'Atom ',k,j, &
                ' charge= ',Z_atom(k),' position ',(x(i),i=1,3)
               write(41,'(3(f10.5,x),i5)') (x(i),i=1,3),Z_atom(k)
            end do
         end do
!_________________ electronic part second
         write(*,'(a)')'Opening the file <charges.sim> for simulation charges...'
         open(31,file='charges.sim',form='formatted',status='old')
         ch=0.0
 76      read(31,*,err=79,end=79) j, (x(i),i=1,3),rCh
         write(41,'(5x,4(f13.8,x))') (x(i),i=1,3),rCh
         do i=1,3
            Dip(i)=Dip(i)+rCh*(x(i)-Center(i))
         end do
         ch=ch+rCh
         go to 76
 79      close (31)
         close (41)
         jj=j
         write(*,'(a,i10)') 'Number of electronic charges found =',jj
         write(*,'(a,f10.5)')'Electronic charge found = ',ch
         dipm=sqrt(Dip(1)*Dip(1)+Dip(2)*Dip(2)+Dip(3)*Dip(3))
         Yes_Dip=.true.
!
![Co].... display atomic positions
!
      ELSE IF(trim(cha2).eq.'Co') THEN
         call show_atoms()
!
!__________ return to the previous menu
!
      ELSE IF(trim(cha2).eq.'Q') THEN
         deallocate(Poten)
         deallocate(Poten1)
         return
      ELSE
         write(*,*)'ERROR! Try again!'
      END IF
      go to 1
end subroutine simul_box

subroutine cut_atoms(grid,totdens)
!....................................................................
!    The density corresponding to specified atoms will be cut out of
! the density by assigning zeros to the corresponding grid points.
!....................................................................
use param
use atoms
use menu
implicit none
real*8,parameter :: tiny=0.01
real*8 GRID(NGX,NGY,NGZ),R(3),totdens,drad,rCharge1,rCharge2,rCharge3,charg
real*8 rad,factor,dx,dv,c1,rad1,ar,rad2,c2,rCharge4,charge_tot,rCharge,denval
real*8,dimension(:),allocatable :: RadCut,NumE_asked,NumE_nonconserv,NumE_conserv
integer,dimension(:),allocatable :: NumAt
character iask,line*40,answer,cha2*2
logical Yes_Spec,Yes_Rad,Yes_Cut
integer iQuit,NumAtCut,jj,i,im,j,iPnt,iPnt3,i0,i1,n1,n2,jj0,nat,ii,iz,iy,ix
integer k1,k2,k3,isp,ijk,j0
real*8 :: TinyCh=0.0
!......................................................................
!_____ choose the starting point, the smallest and the largest radii,
!      the number of points in between and the grid inside the sphere.
! iQuit = 0 - not quit, proceed with plotting in the parent program;
!         1 - quit, do not proceed with plotting.
! method='nonconserv' - for a "non-conserving" algorithm when we scan the
!                       sphere rather than the UC so that each point
!                       may enter several times.
!......................................................................
      allocate(RadCut(NIONS))
      allocate(NumAt(NIONS))
      allocate(NumE_asked(NIONS))
      allocate(NumE_nonconserv(NIONS))
      allocate(NumE_conserv(NIONS))
!
      Yes_Spec=.false.
      Yes_Rad=.false.
      Yes_Cut=.false.
1     iQuit=0
      write(*,*)'..............MENU for CUT ...........................'
      write(*,*)'......... Change these parameters if necessary:.......'
      write(*,*)
      write(*,'(a)') '>>>>> Representation of results: through number of electrons'
      write(*,'(a)') '>>>>> Algorithm for the charge integration: <nonconserving>'

      if(.not.Yes_Spec) then
         iQuit=1
         NumAtCut=0
         jj=0
         write(*,'(a)') &
              '   1. Specify atoms/electrons to be cut out of the density:'
      else
         jj=NumAtCut
         write(*,'(a)')'   1. Atoms to be cut out of the density:'
         write(*,'(a,i5,a)') &
            '      To be cut ',NumAtCut,' atoms with numbers(electrons):'
         do i=1,NumAtCut,6
            im=i+5
            if(NumAtCut-i.lt.6) im=NumAtCut
            write(*,'(5x,15(i3,a,f6.2,a,x))') (NumAt(j),'(',NumE_asked(j),')',j=i,im)
         end do
      end if
      
      write(*,'(a39,f10.5)') '   2. The smallest radius (Angstroms): ', RadiusS
      write(*,'(a38,f10.5)') '   3. The largest radius (Angstroms): ',RadiusL

      if(Nrad.lt.3) then
         iQuit=1
         write(*,'(a)') '   4. The number of points between '// &
                    'these radii: ... undefined ...'
      else
         write(*,'(a48,i5)') &
              '   4. The number of points between these radii: ',Nrad
         dRad=(RadiusL-RadiusS)/(Nrad-1)
      end if

      write(*,'(a48,i5)') '   5. X,Y,Z integration grid inside the sphere: ',NRESOLs
      if(NRESOLs.le.1) iQuit=1

      if(Yes_Rad) then
         write(*,'(a)') '   6. Scan atoms to obtain the radii using charge to cut <= DONE!'
         write(*,'(a)') '   7. Show the list of atoms and their radii + '// &
              ' exact charge to be cut out'
         write(*,'(a,e12.6)') '   8. The threshhold: the smallest density allowed: ',TinyCh
         write(*,'(a)') '   9. Cut atoms out (current density (in memory) is destroyed!)'
         if(Yes_Cut) then
            write(*,'(5x,a,f10.5)') &
                 'Out of the total original density             = ',rCharge1
            write(*,'(5x,a,i10)') &
                 '... grid points removed due to radii          = ',iPnt
            write(*,'(5x,a,f10.5)') &
                 '... with the density cut out                  = ',rCharge2
            write(*,'(5x,a,i10)') &
                 '... add. grid points removed due to threshold = ',iPnt3
            write(*,'(5x,a,e12.6)') &
                 '... with the density cut out                  = ',rCharge4
            write(*,'(5x,a,f10.5)') &
                 '... so that the TOTAL density left is         = ',rCharge3
            write(*,'(a)') '  10. Write charges cut out + core charges into a file'
         end if
      else
         write(*,'(a)') '   6. Scan atoms to obtain the radii using charge to cut'
      end if

      write(*,'(a)') '------- A t o m i c  p o s i t i o n s  ---------'
      write(*,'(a)') '  Co. Show current atomic positions in fractional/Cartesian'
      write(*,'(a)') '------ L e a v e   t h e   m e n u -------------'
      write(*,'(a)')'   Q. Return to the previous menu'
      write(*,*)
      write(*,*)'------> Choose the item and press ENTER:'
      read (*,'(a)',err=1) cha2
!
![1]__________ specify atoms/electrons to be cut of the density
!
      IF(trim(cha2).eq.'1') THEN
 40      write(*,*)'Specify/edit the group of atoms as #1 - #2 '
         write(*,*)'associated with the same number of electrons:'
         read(*,'(a)') line
         do i=1,40
            if(line(i:i).eq.'-') go to 41
         end do
         write(*,*)'ERROR! Dash need to be specified explicitly!'
         go to 40
 41      i0=i-1
         i1=i+1
         read(line(:i0),*,err=40) n1
         if(n1.lt.1 .or. n1.gt.NIONS) then
            write(*,*)'ERROR in the first number!'
            go to 40
         end if
         read(line(i1:),*,err=40) n2
         if(n2.lt.n1 .or. n2.gt.NIONS) then
            write(*,*)'ERROR in the second number!'
            go to 40
         end if
 42      write(*,*)'Specify electronic charge to be associated' &
             //' with these atoms (f10.5)'
         read(*,*,err=42) charg
         jj0=jj
         do 45 i=n1,n2
            if(NumAtCut.ne.0) then
               do j=1,jj0
                  if(NumAt(j).eq.i) then
                     NumE_asked(j)=charg
                     write(*,'(a,i5,a)') &
                          '... the target charge on atom',i,' is changed'
                     go to 45
                  end if
               end do
            end if
            jj=jj+1
            NumAt(jj)=i
            NumE_asked(jj)=charg
45       end do
         NumAtCut=jj
         Yes_Rad=.false.
         Yes_Spec=.true.
!
![2,3]__________ give radii
!
      ELSE IF(trim(cha2).eq.'2') THEN
 10      write(*,*) 'Enter the smallest radius > 0.0 (in Angstroms):'
         read(*,*,err=10) RadiusS
         if(RadiusS.lt.tiny) go to 10
         Yes_Rad=.false.

      ELSE IF(trim(cha2).eq.'3') THEN
 11      write(*,*) 'Enter the largest radius (in Angstroms):'
         read(*,*,err=11) RadiusL
         if(RadiusL.lt.0.0) go to 11
         Yes_Rad=.false.
!
![4]__________ number of points between RadiusS and RadiusL
!
      ELSE IF(trim(cha2).eq.'4') THEN
997      write(*,'(a31,f10.5,a4,f10.5)') 'Give the number of '// &
                  'different radii: '
         read(*,*,err=997) Nrad
         if(Nrad.lt.1) go to 997
         Yes_Rad=.false.
!
![5]__________ number of grid points inside the sphere
!
      ELSE IF(trim(cha2).eq.'5') THEN
996      write(*,*)'Give this number:'
         read(*,*,err=996) NRESOLs
         if(NRESOLs.lt.2) go to 996
         Yes_Rad=.false.
!
![6]__________ calculation: scan all atoms specified in NatAt();
!   for every atom loop over radii from RadiusS till RadiusL to
!   integrate the charge (possible overalp of spheres is ignored),
!   and then interpolate the radius to match the charge in NumE_asked().
!   The result is a vector of radii RadCut(i), i=1,..,NumAtCut
!
      ELSE IF(trim(cha2).eq.'6') THEN
         if(iQuit.ne.0) then
            write(*,*)'ERROR! You still have undefined parameters!'
            go to 1
         end if
         DO nat=1,NumAtCut
            ii=NumAt(nat)
            write(*,'(a,i3)')'Working on the charge for atom ',nat
!__________ get the charge versus radius for this atom.
!__________ make the fitting: get the radius RadCut(nat) to
!           match the charge NumE_asked(nat) (linear interpolation)
!
            rad1=0.0
            c1=0.0
            do i=1,Nrad
               Rad = dRad * (i-1) + RadiusS
               dX=Rad/(NRESOLs/2)
               dV=dX*dX*dX
               factor=dV/VOLC
!
!_________  charge "non-conserving" algorithm:
!     Scan a net of points inside the sphere of Radius using NRESOLs
!     and calculate the amount of charge inside Radius
!
               rCharge=0.0
               do k1=-NRESOLs/2,NRESOLs/2
                  do k2=-NRESOLs/2,NRESOLs/2
                     do k3=-NRESOLs/2,NRESOLs/2
                        R(1)= dX*k1
                        R(2)= dX*k2
                        R(3)= dX*k3
                        aR=sqrt(R(1)*R(1)+R(2)*R(2)+R(3)*R(3))

                        if(aR.le.Rad) then
                           R(1)=R(1)+TI(1,ii)
                           R(2)=R(2)+TI(2,ii)
                           R(3)=R(3)+TI(3,ii)
                           call reducn(R,DIRC,BCELL)
                           call interpolate(R,BCELL,denval,grid)
                           rCharge = rCharge + denval
                        end if

                     end do
                  end do
               end do
               charg = rCharge*factor
               write(*,'(2(a,f10.5))') '>>> Rad= ',Rad,' charge= ',charg
!______________ check for the interpolation interval
               if(NumE_asked(nat).ge.c1 .and. &
                                       NumE_asked(nat).le.charg) then
                  rad2=Rad
                  c2=charg
                  go to 80
               else
                  rad1=Rad
                  c1=charg
               end if
            end do
!______________ error: RadiusL is probably too small
            write(*,'(a,i3)') &
                 'ERROR! Cannot find the fit for the atom nat=',nat
            go to 1

!_____________ fit the radius
 80         RadCut(nat)=rad1+(NumE_asked(nat)-c1)/(c2-c1)*(rad2-rad1)
!
!_____________  calculate the charge using "non-conserving" algorithm again
!               for the interpolated radius
!
            rCharge=0.0
            do k1=-NRESOLs/2,NRESOLs/2
               do k2=-NRESOLs/2,NRESOLs/2
                  do k3=-NRESOLs/2,NRESOLs/2
                     R(1)= dX*k1
                     R(2)= dX*k2
                     R(3)= dX*k3
                     aR=sqrt(R(1)*R(1)+R(2)*R(2)+R(3)*R(3))

                     if(aR.le.RadCut(nat)) then
                        R(1)=R(1)+TI(1,ii)
                        R(2)=R(2)+TI(2,ii)
                        R(3)=R(3)+TI(3,ii)
                        call reducn(R,DIRC,BCELL)
                        call interpolate(R,BCELL,denval,grid)
                        rCharge = rCharge + denval
                     end if

                  end do
               end do
            end do
            NumE_nonconserv(nat)=rCharge*factor
            write(*,'(a,i3,a,f10.5,a,f10.5)') '>>> Fit found: RadCut(', &
                 nat,')= ',RadCut(nat),' for charge = ',NumE_nonconserv(nat)
         END DO
         Yes_Rad=.true.
!
![7]__________ preview all radii for every atom chosen to be cut out
!           and calculate the actual electronic charges to be cut out
!           using the conserving algorithm
!
      ELSE IF(Yes_Rad .and. trim(cha2).eq.'7') THEN
!___________________ first scan the original grid and calculate the
!                    actual charges NumE_conserv(nat) to be cut out
         write(*,*)'Please, wait, scanning the charge density again ...'
         do nat=1,NumAtCut
            NumE_conserv(nat)=0.0d0
         end do
         do iZ=0,NGZ-1
            do iY=0,NGY-1
               do 29 iX=0,NGX-1

                  do nat=1,NumAtCut
                     ii=NumAt(nat)
                     call ask(iX,iY,iZ,TI(1,ii),RadCut(nat),DIRC,iask)
                     if(iask.eq.'y') then
                        NumE_conserv(nat)=NumE_conserv(nat) + grid(iX+1,iY+1,iZ+1)
                        go to 29
                     end if
                  end do

29             end do
            end do
         end do
         do nat=1,NumAtCut
            NumE_conserv(nat)=NumE_conserv(nat)/NPLWV
         end do

!________________ print out the whole bulk of information prior to destroying
!                 the density
         write(*,'(/a)') '>>>>>>>> Radii + charges for chosen atoms <<<<<<<<'
         write(*,'(a)')  '    #  atom        radius   '// &
              '   charge(asked)     charge(estimate)     charge(actual)'
         do nat=1,NumAtCut
            write(*,'(i5,1x,i5,5x,f10.5,x,3(4x,f15.8))') nat,NumAt(nat), &
                 RadCut(nat),NumE_asked(nat),NumE_nonconserv(nat), &
                 NumE_conserv(nat)
         end do
         write(*,'(a)')'Press ENTER when ready ...'
         read(*,*)
!
![8]__________ give the threshold for the density (lower allowed boundary)
!
      ELSE IF(trim(cha2).eq.'8') THEN
 71      write(*,*) 'The smallest density (in electrons) allowed:'
         read(*,*,err=71) TinyCh
         if(TinyCh.lt.0.0) go to 71
         Yes_Cut=.false.
!
![9]__________ cut atoms out into the active GRID(ix,iy,iz) density
!              so that the original one will be destroyed
!
      ELSE IF(Yes_Rad .and. trim(cha2).eq.'9') THEN
         do nat=1,NumAtCut
            NumE_conserv(nat)=0.0d0
         end do
!____________ for statistics (5%, 10%, 15%, 20%, ...)
         ijk=0
         j0=NPLWV/10

         iPnt=0
         iPnt3=0
         rCharge1=0.0
         rCharge2=0.0
         rCharge3=0.0
         rCharge4=0.0
         do iZ=0,NGZ-1
            do iY=0,NGY-1
               do iX=0,NGX-1
                  rCharge1 = rCharge1 + grid(iX+1,iY+1,iZ+1)

!_____________ ask whether the point (iX,iY,iZ) is inside either
!              of the spheres associated with the atoms chosen;
!       item=9: if inside, set it to 0.0; if outside, set it
!               to 0.0 only if grid() < TinyCh, ignore otherwise
!
                  do nat=1,NumAtCut
                     ii=NumAt(nat)
                     call ask(iX,iY,iZ,TI(1,ii),RadCut(nat),DIRC,iask)
                     if(iask.eq.'y') then
                        iPnt=iPnt+1
                        rCharge2 = rCharge2 + grid(iX+1,iY+1,iZ+1)
                        NumE_conserv(nat)=NumE_conserv(nat) + grid(iX+1,iY+1,iZ+1)
                        grid(iX+1,iY+1,iZ+1)=0.0d0
                        go to 21
                     end if
                  end do
                  if(abs(grid(iX+1,iY+1,iZ+1)).gt.TinyCh) then
                     rCharge3 = rCharge3 + grid(iX+1,iY+1,iZ+1)
                  else
                     iPnt3=iPnt3+1
                     rCharge4 = rCharge4 + grid(iX+1,iY+1,iZ+1)
                     grid(iX+1,iY+1,iZ+1)=0.0d0
                  end if
!______________ statistics
21                ijk=ijk+1
                  if(ijk/j0*j0.eq.ijk) write(*,'(a,i3,a)') '... done ',ijk/j0*10,' %'

               end do
            end do
         end do
         rCharge1=rCharge1/NPLWV
         rCharge2=rCharge2/NPLWV
         rCharge3=rCharge3/NPLWV
         rCharge4=rCharge4/NPLWV
         do nat=1,NumAtCut
            NumE_conserv(nat)=NumE_conserv(nat)/NPLWV
         end do
         Yes_Cut=.true.
!
![10]__________ write charges cut out around atoms together with the corresponding
!           core charges into a file 'charges.cut'
!
      ELSE IF(trim(cha2).eq.'10' .and. Yes_Cut) THEN
 77      write(*,*)'Specify nuclei charges for all species:'
         read(*,*,err=77) (Z_atom(i),i=1,NSPEC)

         write(*,*)'Writing total charges into <charges.cut> file ...'
         open(13,file='charges.cut',form='formatted',status='unknown')

         do nat=1,NumAtCut
            ii=NumAt(nat)
            isp=atom_species(ii)
            charge_tot=real(Z_atom(isp))-NumE_conserv(nat)
            write(13,'(3(f10.5,x),e16.10)') TI(1,ii),TI(2,ii),TI(3,ii),charge_tot
         end do
         write(*,*)'Add other core charges (y/other)?'
         read(*,'(a)') answer
         if(answer.eq.'y') then
            do 78 ii=1,NIONS
               do nat=1,NumAtCut
                  if(NumAt(nat).eq.ii) go to 78
               end do
               isp=atom_species(ii)
               write(13,'(3(f10.5,x),e16.10)') &
                    TI(1,ii),TI(2,ii),TI(3,ii),real(Z_atom(isp))
78          end do
         end if
         close (13)
         write(*,*)'Done!'
!
![Co].... display atomic positions
!
      ELSE IF(trim(cha2).eq.'Co') THEN
         call show_atoms()
!
!__________ return to the previous menu
!
      ELSE IF(trim(cha2).eq.'Q') THEN
         deallocate(RadCut)
         deallocate(NumAt)
         deallocate(NumE_asked)
         deallocate(NumE_nonconserv)
         deallocate(NumE_conserv)
         return
      ELSE
         write(*,*)'ERROR! Try again!'
      END IF
      go to 1
end subroutine cut_atoms

subroutine for_gOpenMol(grid)
!....................................................................
!   A box is chosen inside the cell and the density in the box is
! written into a file in the cube gOpenMol format.
!....................................................................
use param
use atoms
use menu
use box
implicit none
real*8,parameter :: tiny=0.0001,A_au=1.889725989
real*8 GRID(NGX,NGY,NGZ),R(3),x(3),R1(3),R2(3)
real*8,dimension(3) :: Center=(/0.0,0.0,0.0/),Box_ext=(/0.0,0.0,0.0/)
real*8,dimension(3) :: Face=(/1.0,1.0,1.0/)
real*8 Sides(3,3),Direct(3,3),RecipS(3,3),face1
real*8 Sides1(3,3),Recip1(3,3),VolBox,den,denval,charge,VolBox1
real*8,dimension(:,:),allocatable  :: Pos
real*8,dimension(:),allocatable  :: Zd_tmp
integer,dimension(:),allocatable :: iType
integer,dimension(3) :: ngrid=(/1,1,1/)
real*8 :: vect(3,3),corner(3),da(3),TotCh=0.0,dV,rCharge
integer i,j,iPnt,ijk,j0,iX,iY,iZ,nt,nat,k,nChrg,i1,i2,i3
character iask,line*40,cha2*2
logical Yes_Do,Yes_xyz,Err
!
!________ default for the directions of the box sides:
!         along the Cartesian axes
!         (Oleg: otheriwse gOpenMol does not understand the file)
      allocate(Pos(3,N_at))
      allocate(Zd_tmp(0:N_z))
      allocate(iType(N_at))
      Direct=0.0
      do i=1,3
        Direct(i,i)=1.0
      end do
!
!................. start the main menu
!
      Yes_Do=.false.
      Yes_xyz=.true.
1     write(*,*)'............MENU for gOpenMol writer .................'
      write(*,*)'......... Change these parameters if necessary:.......'
      write(*,*)

      write(*,'(a,3(f10.5,a))') '   1. The box center is at: (', &
                Center(1),',',Center(2),',',Center(3),')'

      write(*,'(a)')'      Directions of the box sides are fixed along:'
      write(*,'(10x,a,3(f10.5,a))') &
          '1   <',Direct(1,1),',',Direct(1,2),',',Direct(1,3),'  >'
      write(*,'(10x,a,3(f10.5,a))') &
          '2   <',Direct(2,1),',',Direct(2,2),',',Direct(2,3),'  >'
      write(*,'(10x,a,3(f10.5,a))') &
          '3   <',Direct(3,1),',',Direct(3,2),',',Direct(3,3),'  >'

      write(*,'(a,3(1x,f10.5))') &
          '   2. Lengths of the box sides (in A) are: ',(Face(i),i=1,3)

      do i=1,3
         do j=1,3
            Sides(i,j)=Direct(i,j)*Face(i)
         end do
      end do
      call BASTR(Sides,RecipS,VolBox,0)
      corner(1)=Center(1)-0.5*(Sides(1,1)+Sides(2,1)+Sides(3,1))
      corner(2)=Center(2)-0.5*(Sides(1,2)+Sides(2,2)+Sides(3,2))
      corner(3)=Center(3)-0.5*(Sides(1,3)+Sides(2,3)+Sides(3,3))
      write(*,'(5x,a,3(1x,f10.5))') &
          '>>>> the starting point (corner) of the box is at ', &
           (corner(i),i=1,3)

      if(Yes_xyz) then
         write(*,'(a,3(x,f10.5),a)') &
      '   3. Box extension for xyz-file: [',(Box_ext(i),i=1,3),']'
         do i=1,3
            face1=Face(i)+Box_ext(i)
            do j=1,3
               Sides1(i,j)=Direct(i,j)*face1
            end do
         end do
         call BASTR(Sides1,Recip1,VolBox1,0)
      end if

      write(*,'(a,3(1x,i3))') &
      '   4. The grid  in the box in 3 directions: ',(ngrid(i),i=1,3)
      nChrg=ngrid(1)*ngrid(2)*ngrid(3)
      write(*,'(5x,a,i10)')'>>>> total number of grid points = ',nChrg
      do i=1,3
         da(i)=Face(i)/ngrid(i)
         do j=1,3
            vect(i,j)=Direct(i,j)*da(i)
         end do
      end do
      write(*,'(5x,a,3(x,f10.5))')'>>>> steps along 3 directions are ', &
           (da(i),i=1,3)

      write(*,'(a)') &
        '   5. Scan the box and integrate the charge (reference only):'
      if(TotCh.ne.0.0) then
         write(*,'(5x,a,f10.5)') '>>>> total charge in the box = ',TotCh
         write(*,'(5x,a,i10)')   '>>>> total grid points in the box = ',iPnt
      end if

      if(Yes_Do) then
         write(*,'(a)') '   W. Write the cube-format file for gOpenMol  <-- DONE!'
         write(*,'(5x,a,f10.5)') &
                '>>>> check-sum for the density in the box: ',charge*dV
      else
         write(*,'(a)') '   W. Write the cube-format file for gOpenMol'
      end if

      write(*,'(a)') '-------  G e n e r a l  s e t t i n g s ---------'
      write(*,'(a)')'  An. Coordinates are specified in: '//angstr
      write(*,'(a)') '  Co. Show current atomic positions in fractional/Cartesian'
      if(Yes_xyz) then
         write(*,'(a)') &
          '  XY. Atoms in the XYZ box: constrained by the density box'
      else
         write(*,'(a)')   '  XY. Atoms in the XYZ box: all cell atoms'
      end if
      write(*,'(a)') '------------ C a l c u l a t o r ----------------'
      write(*,'(a)') '  cT. Centre of a triangle of points/atoms'
      write(*,'(a)') '  mD. Middle distance between two points/atoms'
      write(*,'(a)') '------ L e a v e   t h e   m e n u --------------'
      write(*,'(a)') '   Q. Return to the previous menu'
      write(*,*)
      write(*,*)'------> Choose the item and press ENTER:'
      read (*,'(a)',err=1) cha2
!
![An]__________ choose the way how the coordinates are given
!
      IF(trim(cha2).eq.'An') THEN
         if(angstr.eq.'<Fractional>') then
            angstr='<Angstroms> '
         else if(angstr.eq.'<Angstroms> ') then
            angstr='<AtomNumber>'
         else if(angstr.eq.'<AtomNumber>') then
            angstr='<Fractional>'
         end if
!
![1]__________ give the central point
!
      ELSE IF(trim(cha2).eq.'1') THEN
         WRITE(*,*)'Give the center of your box'
         call givepoint(Center(1),Center(2),Center(3),angstr)
         Yes_Do=.false.
         TotCh=0.0
!
![2]__________ choose lengths of the box sides
!
      ELSE IF(trim(cha2).eq.'2') THEN
 12      WRITE(*,'(a)')'Give lengths of the box sides (in A)'
         read(*,*,err=12) (Face(j),j=1,3)
         do j=1,3
           Face(j)=abs(Face(j))
           if(Face(j).lt.tiny) go to 12
         end do
         Yes_Do=.false.
         TotCh=0.0
!
![3]__________ choose the extention for the xyz box with respect to the
!              density box (in each direction)
!
      ELSE IF(trim(cha2).eq.'3' .and. Yes_xyz) THEN
 11      write(*,'(a/a)') &
        '====> The box for the atomic coordinates is larger than <=====' &
        ,'   the density box in each direction by (please, specify!):'
         read(*,*,err=11) (Box_ext(i),i=1,3)
!
![4]__________ choose the grid in the box
!
      ELSE IF(trim(cha2).eq.'4') THEN
 13      write(*,'(a)') 'Specify the number of grid points in all 3 directions:'
         read(*,*,err=13) (ngrid(i),i=1,3)
         do i=1,3
           if(ngrid(i).lt.1) go to 13
         end do
         if(ngrid(3).gt.N_z) then
            write(*,*)'ERROR! Cannot be larger than N_z = ',N_z
            go to 13
         end if
         Yes_Do=.false.
!
![5]__________ integrate the charge in the box (for reference only)
!
      ELSE IF(trim(cha2).eq.'5') THEN
         write(*,*)'Using conserving algorithm ...'
!
!____________ for statistics (10%, 20%, ...)
         ijk=0
         j0=NPLWV/10
!
         iPnt=0
         rCharge=0.0
         do iZ=0,NGZ-1
            do iY=0,NGY-1
               do iX=0,NGX-1
!_____________ ask whether the point (iX,iY,iZ) is inside the box
                  call ask_box(iX,iY,iZ,Center,Sides,RecipS,DIRC,iask)
                  if(iask.eq.'y') then
                     iPnt=iPnt+1
                     rCharge = rCharge + grid(iX+1,iY+1,iZ+1)
                  end if
!______________ statistics
                  ijk=ijk+1
                  if(ijk/j0*j0.eq.ijk) &
                       write(*,'(a,i3,a)') '... done ',ijk/j0*10,' %'
               end do
            end do
         end do
         TotCh=rCharge/NPLWV
!
![W]_________ write the gOpenMol cube-format file with the density
!             In addition, an xyz file is also written.
!
      ELSE IF(trim(cha2).eq.'W') THEN

!___________________ write an xyz file with atomic positions
         write(*,'(a)') &
              'Opening the file <gOpenMol.xyz> with the coordinates ...'
         open(37,file='gOpenMol.xyz',form='formatted',status='unknown')
         if(Yes_xyz) then
            call at_in_box(Center,Sides1,Recip1,Pos,iType,N_at,nt,Err)
            if(Err) go to 1
            write(37,'(i5/)') nt
            do i=1,nt
               write(37,'(a2,3(x,f12.6))') Species(iType(i)),(Pos(k,i),k=1,3)
            end do
         else
            write(37,'(i5/)') NIONS
            nat=0
            do i=1,NSPEC
               do j=1,NspN(i)
                  nat=nat+1
                  write(37,'(a,3(x,f12.6))') Species(i),(TI(k,nat),k=1,3)
               end do
            end do
         end if
         close(37)
         WRITE(*,*)'Done! The file <gOpenMol.xyz> created!'

!___________________ open the file and start writing the heading
         call system('date > date.tmp')
         open(31,file='date.tmp',form='formatted',status='old')
         read(31,'(a)') line
         close(31,status='delete')
         write(*,'(a)') 'Opening the file <gOpenMol.cube> with the density ...'
         open(31,file='gOpenMol.cube',form='formatted',status='unknown')
         write(31,*)'cube format gOpenMol input file for the density'
         write(31,'(2a)')'created by lev00 on: ',line
         write(31,'(i5,3f12.6)') NIONS,(corner(i)*A_au,i=1,3)
         do i=1,3
            write(31,'(i5,3f12.6)') ngrid(i),(vect(i,j)*A_au,j=1,3)
         end do
!____________________ write atomic positions in a.u.
         nat=0
         do i=1,NSPEC
            do j=1,NspN(i)
               nat=nat+1
               charge=0.0
               write(31,'(i5,4f12.6)') Mendel(i),charge,(TI(k,nat)*A_au,k=1,3)
            end do
         end do
!____________________ write the density in the box
         charge=0.0
         dV=VolBox/nChrg
         do i1=0,ngrid(1)-1
            do i2=0,ngrid(2)-1
               do i3=0,ngrid(3)-1
!_________________________ Cartesian position of the grid point in the box
                  x(1)= i1*vect(1,1)+i2*vect(2,1)+i3*vect(3,1)+corner(1)
                  x(2)= i1*vect(1,2)+i2*vect(2,2)+i3*vect(3,2)+corner(2)
                  x(3)= i1*vect(1,3)+i2*vect(2,3)+i3*vect(3,3)+corner(3)
!_________________________ interpolate the density at x()
                  call reducn(x,DIRC,BCELL)
                  call interpolate(x,BCELL,denval,grid)
                  den=denval/VOLC
                  Zd_tmp(i3)=den
                  charge = charge + Zd_tmp(i3)
               end do
               write(31,'(6e13.5)') (Zd_tmp(i3),i3=0,ngrid(3)-1)
            end do
         end do
         close (31)
         WRITE(*,*)'Done! The file <gOpenMol.cube> created!'

         Yes_Do=.true.
!
![Co].... display atomic positions
!
      ELSE IF(trim(cha2).eq.'Co') THEN
         call show_atoms()
!
![cT]__________ centre of a triangle
!
      ELSE IF(trim(cha2).eq.'cT') THEN
         WRITE(*,*)'Give the 1st point/atom:'
         call givepoint(R(1),R(2),R(3),angstr)
         WRITE(*,*)'Give the 2nd point/atom:'
         call givepoint(R1(1),R1(2),R1(3),angstr)
         WRITE(*,*)'Give the 3rd point/atom:'
         call givepoint(R2(1),R2(2),R2(3),angstr)
         x=(R+R1+R2)/3.0
         write(*,*)'The centre of the triangle is at:',(x(i),i=1,3)

!
![mD]__________ middle distance between two points
!
      ELSE IF(trim(cha2).eq.'mD') THEN
         WRITE(*,*)'Give the 1st point/atom:'
         call givepoint(R(1),R(2),R(3),angstr)
         WRITE(*,*)'Give the 2nd point/atom:'
         call givepoint(R1(1),R1(2),R1(3),angstr)
         x=(R+R1)/2.0
         write(*,*)'The middle distance between the points is at:', &
              (x(i),i=1,3)
!
![XY]__________ specify if all atoms or only atoms inside the density box
!               should be printed into the xyz file
!
      ELSE IF(trim(cha2).eq.'XY') THEN
         if(Yes_xyz) then
            Yes_xyz=.false.
         else
            Yes_xyz=.true.
         end if
!
![Q]__________ return to the previous menu
!
      ELSE IF(trim(cha2).eq.'Q') THEN
         deallocate(Pos)
         deallocate(Zd_tmp)
         deallocate(iType)
         return
      ELSE
         write(*,*)'ERROR! Try again!'
      END IF
      go to 1
end subroutine for_gOpenMol

subroutine shift_charge(grid,DIRC,BCELL)
!..........................................................................
! the density is recalculated for a shifted system
!..........................................................................
use param
implicit none
real*8 GRID(NGX,NGY,NGZ),DIRC(3,3),BCELL(3,3)
real*8 shift(3),x(3),totdens,denval
real*8,dimension(:,:,:),allocatable :: GRIDn
integer ix,iy,iz

!........ save the current density in gridn()
      allocate(GRIDn(NGX,NGY,NGZ))
      gridn=grid

!........ ask the shift vector

 1    write(*,*)'Give a new origin of the coordinate system as X,Y,Z:'
      read(*,*,err=1) Shift(1),Shift(2),Shift(3)
!
!..................................................................
!............... TRANSFORM the CHARGE DENSITY .....................
!..................................................................
!
      write(*,*)'Transforming the density ...'
      totdens=0.0
      do iZ=0,NGZ-1
         do iY=0,NGY-1
            do iX=0,NGX-1
               x(1)=iX*DIRC(1,1)/NGX+iY*DIRC(2,1)/NGY+ &
                    iZ*DIRC(3,1)/NGZ+Shift(1)
               x(2)=iX*DIRC(1,2)/NGX+iY*DIRC(2,2)/NGY+ &
                    iZ*DIRC(3,2)/NGZ+Shift(2)
               x(3)=iX*DIRC(1,3)/NGX+iY*DIRC(2,3)/NGY+ &
                    iZ*DIRC(3,3)/NGZ+Shift(3)
               call reducn(x,DIRC,BCELL)
               call interpolate(x,BCELL,denval,gridn)
               grid(iX+1,iY+1,iZ+1)=denval
               totdens=totdens+denval
            end do
         end do
      end do
      write(*,*) '.....> New total charge = ',totdens/(NPLWV),' <.....'
      write(*,*)'Done!'
      deallocate(GRIDn)
end subroutine shift_charge

#ifdef HORNOS
#ifdef _OPENMP
subroutine omp_shift_charge(grid,DIRC,BCELL)
!..........................................................................
! the density is recalculated for a shifted system
!..........................................................................
use param
implicit none
real*8 GRID(NGX,NGY,NGZ),DIRC(3,3),BCELL(3,3)
real*8 shift(3),x(3),totdens,denval
real*8,dimension(:,:,:),allocatable :: GRIDn
integer ix,iy,iz

! real*8 xytotdens(NGZ),ztotdens

!........ save the current density in gridn()
      allocate(GRIDn(NGX,NGY,NGZ))
      gridn=grid

!........ ask the shift vector

 1    write(*,*)'Give a new origin of the coordinate system as X,Y,Z:'
      read(*,*,err=1) Shift(1),Shift(2),Shift(3)
!
!..................................................................
!............... TRANSFORM the CHARGE DENSITY .....................
!..................................................................
!
      write(*,*)'OMP: Transforming the density ...'

      totdens=0.0
!      ztotdens=0.0
!      xytotdens=0.0

! parallel do loop for iZ
! slow index: iZ
! fast index: iX, iY
!$OMP PARALLEL DO PRIVATE(iX,iY,x,denval)
      do iZ=0,NGZ-1
         do iY=0,NGY-1
            do iX=0,NGX-1
               x(1)=iX*DIRC(1,1)/NGX+iY*DIRC(2,1)/NGY+ &
                    iZ*DIRC(3,1)/NGZ+Shift(1)
               x(2)=iX*DIRC(1,2)/NGX+iY*DIRC(2,2)/NGY+ &
                    iZ*DIRC(3,2)/NGZ+Shift(2)
               x(3)=iX*DIRC(1,3)/NGX+iY*DIRC(2,3)/NGY+ &
                    iZ*DIRC(3,3)/NGZ+Shift(3)
               ! returns x
               call reducn(x,DIRC,BCELL)
               !
               call interpolate(x,BCELL,denval,gridn)
               grid(iX+1,iY+1,iZ+1)=denval
!$OMP CRITICAL
               totdens=totdens+denval
!$OMP END CRITICAL
!               xytotdens(iZ)=xytotdens(iZ)+denval
            end do
         end do
! !$OMP CRITICAL
!          ztotdens=ztotdens+xytotdens(iZ)
! !$OMP END CRITICAL
      end do
!$OMP END PARALLEL DO

      write(*,*) '.....> New total charge   = ',totdens/(NPLWV),' <.....'
!      write(*,*) '.....> New total charge Z = ',ztotdens/(NPLWV),' <.....'
      write(*,*)'Done!'
      deallocate(GRIDn)
end subroutine omp_shift_charge
#endif
#endif
