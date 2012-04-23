subroutine dipole(grid,totdens,filen,lenght,Dip)
!....................................................................
!   Dipole moment of the charge density grid(i,j,k) is calculated
! inside a sequence of spheres of different radii.
!  nfile - unit number for the file filen(1:lenght) with output data.
!....................................................................
use param
use menu
use atoms
implicit none
real*8 GRID(NGX,NGY,NGZ),R(3),x(3),Dip(3),totdens
character iask,cha
character filen*12,Title*50,title_pl*7
data Title/'                                                  '/
character(len=10) :: method1='nonconserv'
real*8,parameter :: dzero=0.0
integer Nradd,iQuit,item,k1,k2,k3,j,i,k,nat,n1,n2,n3,lenght
data title_pl/'       '/
logical Yes_Do,Ncharge
real*8 fCENTX,fCENTY,fCENTZ,drad,units,dX,dV,factor,ar,ax,dip1,rad,denval

      Point(1:3,1)=dzero
      Nradd=1
      way_res='D'
!......................................................................
!_____ choose the starting point, the smallest and the largest radii,
!      the number of points in between and the grid inside the sphere.
! iQuit = 0 - not quit, proceed with plotting in the parent program;
!         1 - quit, do not proceed with plotting.
! method='nonconserv' - for a "non-conserving" algorithm when we scan the
!                       sphere rather than the UC so that each point
!                       may enter several times.
!......................................................................
      Ncharge=.false.
      Yes_Do=.false.
1     iQuit=0
      write(*,*)'..............MENU for DIPOLE ........................'
      write(*,*)'......... Change these parameters if necessary:.......'
      write(*,*)

      if(.not.Ncharge) then
        write(*,'(a)') '   2. Charges of nucleii: undefined'
        iQuit=1
      else
        write(*,'(a)') '   2. Charges of nucleii by species: '
        write(*,'(6x,20(i2,1x))') (Z_atom(i),i=1,NSPEC)
      end if

      write(*,'(a)')'   3. The center of your sphere:'
      fCENTX=BCELL(1,1)*Point(1,1)+BCELL(1,2)*Point(2,1)+ &
             BCELL(1,3)*Point(3,1)
      fCENTY=BCELL(2,1)*Point(1,1)+BCELL(2,2)*Point(2,1)+ &
             BCELL(2,3)*Point(3,1)
      fCENTZ=BCELL(3,1)*Point(1,1)+BCELL(3,2)*Point(2,1)+ &
             BCELL(3,3)*Point(3,1)
      write(*,'(a,f8.3,2(a1,f8.3),a9,f8.3,2(a1,f8.3),a1)') &
        ' A=> (',Point(1,1),',',Point(2,1),',',Point(3,1), &
        '), fr=> (', fCENTX,',',fCENTY,',',fCENTZ,')'

      if(RadiusS.le.dzero) then
        write(*,'(a)')'   4. The smallest radius (Angstroms): undefined'
        iQuit=1
      else
        write(*,'(a,f10.5)') &
           '   4. The smallest radius (Angstroms): ', RadiusS
      end if

      if(RadiusL.le.dzero .or. RadiusL.lt.RadiusS) then
        write(*,'(a)')'   5. The largest radius (Angstroms): undefined'
        if(RadiusL.lt.RadiusS) write(*,'(a)') &
          '      ERROR: it is less than the smallest one!'
        iQuit=1
      else
        write(*,'(a39,f10.5)') &
           '   5. The largest radius (Angstroms):  ',RadiusL
      end if

      if(Nradd.eq.0) then
        iQuit=1
        write(*,'(a)') '   6. The number of points between '// &
                    'these radii: ... undefined ...'
      else
        write(*,'(a48,i5)') &
         '   6. The number of points between these radii: ',Nradd
        if(Nradd.eq.1) then
          dRad=RadiusS
        else
          dRad=(RadiusL-RadiusS)/(Nradd-1)
        end if
      end if

      if(Yes_Do) then
        write(*,'(a)') &
         '   9. Calculate dipole moment; the file for plotting is ' &
                                            //filen//' <= DONE!'
        if(Nradd.eq.1) write(*,'(6x,a,3(f10.5,a))') &
         'Dipole moment is: (', &
           Dip(1)*units,',',Dip(2)*units,',',Dip(3)*units,') '//cha
      else
        write(*,'(a)') &
         '   9. Calculate dipole moment; the file for plotting is ' &
                                            //filen
      end if

      write(*,'(a)') &
          '  10. Preview the dependence of dipole moment versus Radius'
      write(*,'(a)') &
          '  11. Create a PostScript file '//filen(:lenght)//'.ps'

      write(*,'(a)') &
                    '-------  G e n e r a l  s e t t i n g s ---------'

      write(*,'(a)')'   0. Coordinates are specified in: '//angstr
      if(way_res.eq.'D') then
          write(*,'(a)') '   1. Units: Debye'
          cha='D'
          units=4.80320817
      else
          write(*,'(a)') '   1. Units: electrons*Angstrem'
          cha='e'
          units=1.0
      end if
      write(*,'(a)') '   7. Algorithm for charge integration is FIXED to: <'// &
                         method1//'>'

      write(*,'(a48,i5)') &
       '   8. X,Y,Z integration grid inside the sphere: ',NRESOLd
      if(NRESOLd.le.1) iQuit=1

      write(*,'(a)')'------ L e a v e   t h e   m e n u -------------'
      write(*,'(a)')'  12. Return to the previous menu'
      write(*,*)
      write(*,*)'------> Choose the item and press ENTER:'
      read(*,*,err=100) item
!
!__________ choose the way how the coordinates are given
!
      IF(item.eq.0) THEN
         if(angstr.eq.'<Fractional>') then
            angstr='<Angstroms> '
         else if(angstr.eq.'<Angstroms> ') then
            angstr='<AtomNumber>'
         else if(angstr.eq.'<AtomNumber>') then
            angstr='<Fractional>'
         end if
!
!__________ choose a way to present the results
!
      ELSE IF(item.eq.1) THEN
         if(way_res.eq.'D') then
            way_res='e'
         else
            way_res='D'
         end if
!
!__________ specify charges on nucleii by species
!
      ELSE IF(item.eq.2) THEN
 12      write(*,*)'Specify nucleii charges in the order of species:'
         write(*,'(10(a2,x))') (Species(i),i=1,NSPEC)
         read(*,*,err=12) (Z_atom(i),i=1,NSPEC)
         Ncharge=.true.
         Yes_Do=.false.
!
!__________ give sphere position
!
      ELSE IF(item.eq.3) THEN
         WRITE(*,*)'Give position of your sphere'
         call givepoint(Point(1,1),Point(2,1),Point(3,1),angstr)
         Yes_Do=.false.
!
!__________ give radii
!
      ELSE IF(item.eq.4) THEN
 10      WRITE(*,*) 'Enter the smallest radius (in Angstroms):'
         READ(*,*,err=10) RadiusS
         if(RadiusS.le.dzero) go to 10
         Yes_Do=.false.

      ELSE IF(item.eq.5) THEN
 11      WRITE(*,*) 'Enter the largest radius (in Angstroms):'
         READ(*,*,err=11) RadiusL
         if(RadiusL.le.dzero) go to 11
         Yes_Do=.false.
!
!__________ number of points between RadiusS and RadiusL
!
      ELSE IF(item.eq.6) THEN
997      write(*,'(a)') 'Give the number of different radii: '
         read(*,*,err=997) Nradd
         if(Nradd.lt.1) go to 997
         Yes_Do=.false.
!
!__________ algorithm: charge conserving or not
!
      ELSE IF(item.eq.7) THEN
         write(*,*)'WARNING! You cannot change the method!'
!
!__________ number of grid points inside the sphere
!
      ELSE IF(item.eq.8) THEN
 996      write(*,*)'Specify the grid:'
          read(*,*,err=996) NRESOLd
          if(NRESOLd.lt.2) go to 996
          Yes_Do=.false.
!
!__________ calculation
!   (loop over values of Radius from RadiusS till RadiusL)
!
      ELSE IF(item.eq.9) THEN
        if(iQuit.ne.0) then
           write(*,*)'ERROR! You still have undefined parameters!'
           go to 1
        end if
        open(31,file=filen(:lenght),status='unknown',form='formatted')
        write(*,*)'The file '//filen//' has been opened for the PLOT.'
        write(*,*)'Working on the dipole: writing to '//filen(1:lenght)//' ...'
        DO i=1,Nradd
           Rad = dRad * (i-1) + RadiusS
           write(*,*)'..... Radius = ',Rad,' ......'
           do j=1,3
              Dip(j)=dzero
           end do
!
!____________((((((((((((((  electronic part first  ))))))))))))))
!
!          tot=0.0
!
!_________   "non-conserving" algorithm: since each point of the
!            UC may be met more than once, the charge is not properly
!            normalised. Scan the net of points inside the sphere of
!            Radius using NRESOLd and calculate the amount of charge
!            inside Rad and the dipole wrt THE SPHERE CENTER
!
           dX=Rad/(NRESOLd/2)
           dV=dX*dX*dX
           factor=dV/VOLC

           do k1=-NRESOLd/2,NRESOLd/2
              do k2=-NRESOLd/2,NRESOLd/2
                 do k3=-NRESOLd/2,NRESOLd/2
                    R(1)= dX*k1
                    R(2)= dX*k2
                    R(3)= dX*k3
                    aR=sqrt(R(1)*R(1)+R(2)*R(2)+R(3)*R(3))
                    if(aR.le.Rad) then
                       x(1)=R(1)+Point(1,1)
                       x(2)=R(2)+Point(2,1)
                       x(3)=R(3)+Point(3,1)
                       call reducn(x,DIRC,BCELL)
                       call interpolate(x,BCELL,denval,grid)
                       do j=1,3
                          Dip(j)=Dip(j)-denval*R(j)*factor
                       end do
                    end if
                 end do
              end do
           end do
           
!
!____________((((((((((((((  nuclear part second  ))))))))))))))
!
           nat=0
           do k=1,NSPEC
              do j=1,NspN(k)
                 nat=nat+1
!___________________(a) find the right image of this nuclei which is
!                       inside the sphere
                 do n1=-1,1
                    do n2=-1,1
                       do n3=-1,1
                          x(1)=n1*DIRC(1,1)+n2*DIRC(2,1)+n3*DIRC(3,1)+TI(1,nat)- &
                               Point(1,1)
                          x(2)=n1*DIRC(1,2)+n2*DIRC(2,2)+n3*DIRC(3,2)+TI(2,nat)- &
                               Point(2,1)
                          x(3)=n1*DIRC(1,3)+n2*DIRC(2,3)+n3*DIRC(3,3)+TI(3,nat)- &
                               Point(3,1)
                          ax=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
                          if(ax.lt.Rad) go to 33
                       end do
                    end do
                 end do
                 write(*,'(a,i5,a)') &
                      'ERROR: atom nat=',nat,' is not inside the sphere!'
                 write(*,'(a)') 'Make sure the position and the radius are correct!'
                 go to 1
!___________________(b) calculate the contribution to the dipole moment
33               Dip(1)=Dip(1)+Z_atom(k)*x(1)
                 Dip(2)=Dip(2)+Z_atom(k)*x(2)
                 Dip(3)=Dip(3)+Z_atom(k)*x(3)
              end do
           end do
           dip1=sqrt(Dip(1)*Dip(1)+Dip(2)*Dip(2)+Dip(3)*Dip(3))
           write(*,'(5(f10.5,5x))') Rad, dip1*units,(Dip(j)*units,j=1,3)
           write(31,'(5(f10.5,5x))') Rad,dip1*units,(Dip(j)*units,j=1,3)
        END DO
        close (31)
        write(*,*)'.... File '//filen(1:lenght)//' has been created! ....'
        Yes_Do=.true.
!
!__________ preview the file just created
!
      ELSE IF(item.eq.10) THEN
         if(Yes_Do) then
            call Plot1(filen,lenght,Title,title_pl, &
                 'Sphere radius (A)   ', &
                 'Dipole module   ('//cha//') ',  'Screen', 33,0, &
                 'N',.false.,dzero,dzero)
         else
            write(*,*) 'IGNORED! You have to accomplish the option 9 first!'
         end if
!
!__________ create a PostScript file of the plot
!
      ELSE IF(item.eq.11) THEN
         if(Yes_Do) then
            write(*,*)'Give the title:'
            read(*,'(a)') Title
            call Plot1(filen,lenght,Title,title_pl, &
              'Sphere radius (A)   ', &
              'Dipole module   ('//cha//') ', 'Postsc', 33,0, &
                 'N',.false.,dzero,dzero)
         else
            write(*,*) 'IGNORED! You have to accomplish the option 9 first!'
         end if
!
!__________ return to the previous menu
!
      ELSE IF(item.eq.12) THEN
         return
      ELSE
         go to 100
      END IF
      go to 1
!
!....... error
100   write(*,*)'Incorrect item number! Try again!'
      go to 1
end subroutine dipole

subroutine quadrpl(grid,Quadr)
!....................................................................
!   Quadrupole tensor of the charge density grid(i,j,k) is calculated
! inside a sequence of spheres of different radii.
!  nfile - unit number for the file filen(1:lenght) with output data.
!....................................................................
use param
use menu
use atoms
implicit none
character(len=10) :: method1='nonconserv'
real*8,parameter :: dzero=0.0
real*8 GRID(NGX,NGY,NGZ),R(3),x(3),Dip(3)
real*8 Quadr(3,3),A(3),E(3),quad(3,3)
character iask,cha
logical Yes_Do,Ncharge,Yes_Axes
real*8 fCENTX,fCENTY,fCENTZ,Rad,dX,dV,factor,aR,aR2,ab,denval,ax
integer item,i,iQuit,k1,k2,k3,j,nat,n1,n2,n3,k,i1

      Point(1:3,1)=dzero
!......................................................................
!_____ choose the starting point, the smallest and the largest radii,
!      the number of points in between and the grid inside the sphere.
! iQuit = 0 - not quit, proceed with plotting in the parent program;
!         1 - quit, do not proceed with plotting.
! method1='nonconserv' - for a "non-conserving" algorithm when we scan the
!                       sphere rather than the UC so that each point
!                       may enter several times.
!......................................................................
      Yes_Do=.false.
      Yes_Axes=.false.
      do i=1,NSPEC
        if(Z_atom(i).gt.0) then
          Ncharge=.true.
          go to 1
        end if
      end do
      Ncharge=.false.
1     iQuit=0
      write(*,*)'..............MENU for QUADRUPOLE ....................'
      write(*,*)'......... Change these parameters if necessary:.......'
      write(*,*)

      if(.not.Ncharge) then
        write(*,'(a)') '   1. Charges of nucleii: undefined'
        iQuit=1
      else
        write(*,'(a)') '   1. Charges of nucleii by species: '
        write(*,'(6x,20(i2,1x))') (Z_atom(i),i=1,NSPEC)
      end if

      write(*,'(a)')'   2. The center of your sphere:'
      fCENTX=BCELL(1,1)*Point(1,1)+BCELL(1,2)*Point(2,1)+ &
             BCELL(1,3)*Point(3,1)
      fCENTY=BCELL(2,1)*Point(1,1)+BCELL(2,2)*Point(2,1)+ &
             BCELL(2,3)*Point(3,1)
      fCENTZ=BCELL(3,1)*Point(1,1)+BCELL(3,2)*Point(2,1)+ &
             BCELL(3,3)*Point(3,1)
      write(*,'(a,f8.3,2(a1,f8.3),a9,f8.3,2(a1,f8.3),a1)') &
        ' A=> (',Point(1,1),',',Point(2,1),',',Point(3,1), &
        '), fr=> (', fCENTX,',',fCENTY,',',fCENTZ,')'

      if(RadiusS.le.dzero) then
         write(*,'(a)')'   3. The sphere radius (Angstroms): undefined'
         iQuit=1
      else
         write(*,'(a,f10.5)') '   3. The sphere radius (Angstroms): ', RadiusS
      end if

      if(Yes_Do) then
         write(*,'(a)') '   6. Calculate quadrupole tensor.  <= DONE!'
         write(*,'(6x,a,3(f10.5,a))') &
              '      Quadrupole tensor is [in Electron*Angstrem^2]:'
         do i=1,3
            write(*,'(10x,3(1x,f10.5))') (Quadr(i,j),j=1,3)
         end do
         write(*,'(a)') '   7. Calculate principal axes of the quadrupole tensor'
         if(Yes_Axes) then
            write(*,'(7x,a,3(f10.5,1x))') 'Eigenvalues:  ',(E(j),j=1,3)
            write(*,'(7x,a,3(f10.5,1x))') 'Eigenvectors: ',(quad(1,j),j=1,3)
            write(*,'(21x,3(f10.5,1x))') (quad(2,j),j=1,3)
            write(*,'(21x,3(f10.5,1x))') (quad(3,j),j=1,3)
         end if
      else
         write(*,'(a)') '   6. Calculate quadrupole tensor.'
      end if

      write(*,'(a)') '-------  G e n e r a l  s e t t i n g s ---------'
      write(*,'(a)') '   0. Coordinates are specified in: '//angstr
      write(*,'(a)') &
       '   4. Algorithm for charge integration is FIXED to: <'// &
           method1//'>'

      write(*,'(a48,i5)') &
       '   5. X,Y,Z integration grid inside the sphere: ',NRESOLd
      if(NRESOLd.le.1) iQuit=1


      write(*,'(a)')'------ L e a v e   t h e   m e n u -------------'
      write(*,'(a)')'   8. Return to the previous menu'
      write(*,*)
      write(*,*)'------> Choose the item and press ENTER:'
      read(*,*,err=100) item
!
!__________ choose the way how the coordinates are given
!
      IF(item.eq.0) THEN
         if(angstr.eq.'<Fractional>') then
            angstr='<Angstroms> '
         else if(angstr.eq.'<Angstroms> ') then
            angstr='<AtomNumber>'
         else if(angstr.eq.'<AtomNumber>') then
            angstr='<Fractional>'
         end if
!
!__________ specify charges on nucleii by species
!
      ELSE IF(item.eq.1) THEN
 12      write(*,*)'Specify nucleii charges in the order of species:'
         read(*,*,err=12) (Z_atom(i),i=1,NSPEC)
         Ncharge=.true.
         Yes_Do=.false.
!
!__________ give sphere position
!
      ELSE IF(item.eq.2) THEN
         WRITE(*,*)'Give position of your sphere'
         call givepoint(Point(1,1),Point(2,1),Point(3,1),angstr)
         Yes_Do=.false.
!
!__________ give radius
!
      ELSE IF(item.eq.3) THEN
 10      WRITE(*,*) 'Enter the sphere radius (in Angstroms):'
         READ(*,*,err=10) RadiusS
         if(RadiusS.le.dzero) go to 10
         Yes_Do=.false.
!
!__________ algorithm: charge conserving or not
!
      ELSE IF(item.eq.4) THEN
         write(*,*)'WARNING! You cannot change the method!'
!
!__________ number of grid points inside the sphere
!
      ELSE IF(item.eq.5) THEN
 996       write(*,*)'Give this number:'
           read(*,*,err=996) NRESOLd
           if(NRESOLd.lt.2) go to 996
           Yes_Do=.false.
!
!__________ calculation
!
      ELSE IF(item.eq.6) THEN
         if(iQuit.ne.0) then
            write(*,*)'ERROR! You still have undefined parameters!'
            go to 1
         end if
         write(*,*)'Working on the quadrupole tensor ...'

         Rad = RadiusS
         Quadr=dzero
!
!____________((((((((((((((  electronic part first  ))))))))))))))
!
!          tot=0.0
!
!_________   "non-conserving" algorithm: since each point of the
!            UC may be met more than once, the charge is not properly
!            normalised. Scan the net of points inside the sphere of
!            Radius using NRESOLd and calculate the amount of charge
!            inside Rad and the quadrupole wrt THE SPHERE CENTER
!
         dX=Rad/(NRESOLd/2)
         dV=dX*dX*dX
         factor=dV/VOLC

         do k1=-NRESOLd/2,NRESOLd/2
            do k2=-NRESOLd/2,NRESOLd/2
               do k3=-NRESOLd/2,NRESOLd/2
                  R(1)= dX*k1
                  R(2)= dX*k2
                  R(3)= dX*k3
                  aR=sqrt(R(1)*R(1)+R(2)*R(2)+R(3)*R(3))
                  if(aR.le.Rad) then
                     x(1)=R(1)+Point(1,1)
                     x(2)=R(2)+Point(2,1)
                     x(3)=R(3)+Point(3,1)
                     call reducn(x,DIRC,BCELL)
                     call interpolate(x,BCELL,denval,grid)
                     aR2=aR*aR
                     do i=1,3
                        do j=1,3
                           ab=3*R(i)*R(j)
                           if(i.eq.j) ab=ab-aR2
                           Quadr(i,j)=Quadr(i,j)-denval*ab*factor
                        end do
                     end do
!                      tot=tot+denval*factor
                  end if
               end do
            end do
         end do
         
!
!____________((((((((((((((  nuclear part second  ))))))))))))))
!
         nat=0
         do k=1,NSPEC
            do j=1,NspN(k)
               nat=nat+1
!___________________(a) find the right image of this nuclei which is
!                       inside the sphere
               do n1=-1,1
                  do n2=-1,1
                     do n3=-1,1
                        x(1)=n1*DIRC(1,1)+n2*DIRC(2,1)+n3*DIRC(3,1)+TI(1,nat)- &
                             Point(1,1)
                        x(2)=n1*DIRC(1,2)+n2*DIRC(2,2)+n3*DIRC(3,2)+TI(2,nat)- &
                             Point(2,1)
                        x(3)=n1*DIRC(1,3)+n2*DIRC(2,3)+n3*DIRC(3,3)+TI(3,nat)- &
                             Point(3,1)
                        ax=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
                        if(ax.lt.Rad) go to 33
                     end do
                  end do
               end do
               write(*,'(a,i5,a)') &
                    'ERROR: atom nat=',nat,' is not inside the sphere!'
               write(*,'(a)') 'Make sure the position and the radius are correct!'
               go to 1
!___________________(b) calculate the contribution to the quadrupole moment
33             aR2=ax*ax
               do i=1,3
                  do i1=1,3
                     ab=3*x(i)*x(i1)
                     if(i.eq.i1) ab=ab-aR2
                     Quadr(i,i1)=Quadr(i,i1)+Z_atom(k)*ab
                  end do
               end do
            end do
         end do
         
         Yes_Do=.true.
!
!_____________ principal axes of the quadrupole tensor
!
      else if(item.eq.7 .and. Yes_Do) then
         quad=Quadr
         call diag(quad,A,E,3,.false.)
         Yes_Axes=.true.
!
!__________ return to the previous menu
!
      ELSE IF(item.eq.8) THEN
         return
      ELSE
         go to 100
      END IF
      go to 1
!
!....... error
 100  write(*,*)'Incorrect item number! Try again!'
      go to 1
end subroutine quadrpl






