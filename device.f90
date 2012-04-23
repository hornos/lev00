subroutine check_atoms(nat,good)
use atoms
implicit none
!.............................................................
! checks if atom nat is equivalent to any of k=1,...,nat-1
!.............................................................
  real*8 r(3)
  integer nat,k
  logical good,ask_equiv
  if(nat.gt.1) then
     do k=1,nat-1
        r(1)=TI(1,k)-TI(1,nat)
        r(2)=TI(2,k)-TI(2,nat)
        r(3)=TI(3,k)-TI(3,nat)
        if(ask_equiv(r,BCELL,3,tiny_equiv)) then
           write(*,'(a17,i2,a27,i2)')'ERROR! Your atom ',  &
                nat,' is equivalent to atom ',k
           good=.false.
        end if
     end do
  end if
  return
end subroutine check_atoms

subroutine ask(iX,iY,iZ,Point,Radius,DIRC,iask)
!....................................................................
! Replies with:
!       iask='y'  -  if the grid point (iX,iY,iZ) is inside a sphere
!                    of the Radius drawn around the Point;
!       iask='n'  -  if not.
!....................................................................
!  DIRC - direct lattice vectors
!....................................................................
use param
implicit none
character iask*1
integer i,j,k,iX,iY,iZ
real*8 Point(3),r(3),x(3),DIRC(3,3),Radius,dist
!
! _______ get the point (iX,iY,iZ) in Cartesian coordinates as folows:
  r(1) = iX*DIRC(1,1)/NGX + iY*DIRC(2,1)/NGY + iZ*DIRC(3,1)/NGZ
  r(2) = iX*DIRC(1,2)/NGX + iY*DIRC(2,2)/NGY + iZ*DIRC(3,2)/NGZ
  r(3) = iX*DIRC(1,3)/NGX + iY*DIRC(2,3)/NGY + iZ*DIRC(3,3)/NGZ
!
!_______ a distance, dist, between r() and 27 equivalent spheres 
! with the centers put in x(), is calculated. The latter are built 
! from Point (the center of the input sphere) by adding the smallest 
! 27 direct lattice vectors:
!
  do i=-1,1
     do j=-1,1
        do k=-1,1
           x(1) = i*DIRC(1,1)+j*DIRC(2,1)+k*DIRC(3,1) + Point(1)
           x(2) = i*DIRC(1,2)+j*DIRC(2,2)+k*DIRC(3,2) + Point(2)
           x(3) = i*DIRC(1,3)+j*DIRC(2,3)+k*DIRC(3,3) + Point(3)
           dist=sqrt((r(1)-x(1))**2+(r(2)-x(2))**2+(r(3)-x(3))**2)
           if(dist.lt.Radius) then
              iask='y'
              return
           end if
        end do
     end do
  end do
  iask='n'
  return
end subroutine ask

subroutine ask_box(iX,iY,iZ,Center,SIDES,RECIP,DIRC,iask)
!.........................................................................
! Replies with:
!  iask='y'  -  if the grid point (iX,iY,iZ) is inside a parallelepiped
!               (the box) drawn around the Center with sides along the
!               directions of the vectors SIDES
!  iask='n'  -  if not.
!.........................................................................
!  DIRC(#,xyz) - direct lattice vectors
!  SIDES(#,xyz) - directions of the box sides
!  RECIP(#,xyz) - "reciprocal" vectors associated with SIDES, i.e.
!                    RECIP(#1)*SIDES(#2)=delta(#1,#2) (dot product)
!.........................................................................
use param
implicit none
character iask*1
integer iX,iY,iZ
real*8 Center(3),DIRC(3,3),r(3),SIDES(3,3),RECIP(3,3)
!
!_______ get the point (iX,iY,iZ) in Cartesian coordinates as folows:
      r(1) = iX*DIRC(1,1)/NGX + iY*DIRC(2,1)/NGY + iZ*DIRC(3,1)/NGZ
      r(2) = iX*DIRC(1,2)/NGX + iY*DIRC(2,2)/NGY + iZ*DIRC(3,2)/NGZ
      r(3) = iX*DIRC(1,3)/NGX + iY*DIRC(2,3)/NGY + iZ*DIRC(3,3)/NGZ
!
!_______ ask if the point r() is inside the box
      call ask_box1(r,Center,SIDES,RECIP,DIRC,iask)
end subroutine ask_box

subroutine ask_box1(r,Center,SIDES,RECIP,DIRC,iask)
!.........................................................................
! Replies with:
!  iask='y'  -  if the Cartesian grid point r() is inside a parallelepiped
!               (the box) drawn around the Center with sides along the
!               directions of the vectors SIDES
!  iask='n'  -  if not.
!.........................................................................
!  DIRC(#,xyz) - direct lattice vectors
!  SIDES(#,xyz) - directions of the box sides
!  RECIP(#,xyz) - "reciprocal" vectors associated with SIDES, i.e.
!                    RECIP(#1)*SIDES(#2)=delta(#1,#2) (dot product)
!.........................................................................
use param
implicit none
character iask*1
real*8,parameter :: tiny = 0.0001
real*8 Center(3),DIRC(3,3),r(3),x(3),corner(3),SIDES(3,3)
real*8 rins(3),RECIP(3,3),d0,d1,c1,c2,c3
integer i,j,k

!
!_______ the vector to the 1st corner of the box
!
      corner(1)=Center(1)-0.5*(SIDES(1,1)+SIDES(2,1)+SIDES(3,1))
      corner(2)=Center(2)-0.5*(SIDES(1,2)+SIDES(2,2)+SIDES(3,2))
      corner(3)=Center(3)-0.5*(SIDES(1,3)+SIDES(2,3)+SIDES(3,3))
!
!_______ part of r() inside the box
!
      rins(1)=r(1)-corner(1)
      rins(2)=r(2)-corner(2)
      rins(3)=r(3)-corner(3)
!
!_______ The vector rins() is inside the box fixed by vectors SIDES if
!        in the expansion of rins() wrt SIDES all 3 coefficients are
!        between 0.0  and 1.0
!_______ We run over all possible 27 images of the box to check if the
!        point rins() gets inside at least one of those.
!
      d0=-tiny
      d1=1.0+tiny
      do i=-1,1
         do j=-1,1
            do k=-1,1
               x(1) = i*DIRC(1,1)+j*DIRC(2,1)+k*DIRC(3,1) + rins(1)
               x(2) = i*DIRC(1,2)+j*DIRC(2,2)+k*DIRC(3,2) + rins(2)
               x(3) = i*DIRC(1,3)+j*DIRC(2,3)+k*DIRC(3,3) + rins(3)
               c1=x(1)*RECIP(1,1)+x(2)*RECIP(1,2)+x(3)*RECIP(1,3)
               if(c1.ge.d0 .and. c1.le.d1) then
                  c2=x(1)*RECIP(2,1)+x(2)*RECIP(2,2)+x(3)*RECIP(2,3)
                  if(c2.ge.d0 .and. c2.le.d1) then
                     c3=x(1)*RECIP(3,1)+x(2)*RECIP(3,2)+x(3)*RECIP(3,3)
                     if(c3.ge.d0 .and. c3.le.d1) then
                        iask='y'
                        return
                     end if
                  end if
               end if
            end do
         end do
      end do
      iask='n'
end subroutine ask_box1

subroutine ask_box2(r,Center,SIDES,RECIP,iask)
!.........................................................................
! Replies with:
!  iask='y'  -  if the grid point r() is inside a parallelepiped
!               (the box) drawn around the Center with sides along the
!               directions of the vectors SIDES
!  iask='n'  -  if not.
!.........................................................................
!  SIDES(#,xyz) - directions of the box sides
!  RECIP(#,xyz) - "reciprocal" vectors associated with SIDES, i.e.
!                    RECIP(#1)*SIDES(#2)=delta(#1,#2) (dot product)
!.........................................................................
! the only difference with ask_box1() is that images are not used.
!.........................................................................
use param
implicit none
character iask*1
real*8,parameter :: tiny=0.0001
real*8 Center(3),r(3),x(3),corner(3),SIDES(3,3),RECIP(3,3),d0,d1,c1,c2,c3

!
!_______ the vector to the 1st corner of the box
!
      corner(1)=Center(1)-0.5*(SIDES(1,1)+SIDES(2,1)+SIDES(3,1))
      corner(2)=Center(2)-0.5*(SIDES(1,2)+SIDES(2,2)+SIDES(3,2))
      corner(3)=Center(3)-0.5*(SIDES(1,3)+SIDES(2,3)+SIDES(3,3))
!
!_______ part of r() inside the box
!
      x(1)=r(1)-corner(1)
      x(2)=r(2)-corner(2)
      x(3)=r(3)-corner(3)
!
!_______ The vector x() is inside the box fixed by vectors SIDES if
!        in the expansion of x() wrt SIDES all 3 coefficients are
!        between 0.0  and 1.0

      d0=-tiny
      d1=1.0+tiny
      c1=x(1)*RECIP(1,1)+x(2)*RECIP(1,2)+x(3)*RECIP(1,3)
      if(c1.ge.d0 .and. c1.le.d1) then
         c2=x(1)*RECIP(2,1)+x(2)*RECIP(2,2)+x(3)*RECIP(2,3)
         if(c2.ge.d0 .and. c2.le.d1) then
            c3=x(1)*RECIP(3,1)+x(2)*RECIP(3,2)+x(3)*RECIP(3,3)
            if(c3.ge.d0 .and. c3.le.d1) then
               iask='y'
               return
            end if
         end if
      end if
      iask='n'
end subroutine ask_box2

subroutine show_atoms()
use param
use atoms
implicit none
integer k,i,j,j1,i1
real*8 x(3)
      k=0
      do i=1,NSPEC
         write(*,'(3a,i5)')'.... Species <',Species(i),'> ......> ',i
         write(*,'(a)')' Tot   # |----------  Fractional ---------|' &
              //'-------------  Cartesian ---------|'
         do j=1,NspN(i)
            k=k+1
            do i1=1,3
               x(i1)=0.0
               do j1=1,3
                  x(i1)=x(i1)+TI(j1,k)*BCELL(i1,j1)
               end do
            end do
            write(*,'(2(x,i3),x,3(x,f10.5),2x,3(x,f10.5),3x,a,4x,a)') &
                 k,j,(x(i1),i1=1,3),(TI(j1,k),j1=1,3)
         end do
      end do
      write(*,*)'Hit ENTER when done ...'
      read(*,*)
end subroutine show_atoms

subroutine show_atoms_tagged(listat,ngroup)
use param
use atoms
implicit none
logical listat(NIONS,ngroup)
character cha
integer k,i,j,j1,ngroup
real*8 a1,a2,a3
     
   k=0
   do i=1,NSPEC
      write(*,'(3a,i5)')'.... Species <',Species(i),'> ......> ',i
      write(*,'(a)')' Tot   # |----------  Fractional ---------|' &
              //'-------------  Cartesian ---------| Tag'
      do j=1,NspN(i)
         k=k+1
         a1=TI(1,k)*BCELL(1,1)+TI(2,k)*BCELL(1,2)+TI(3,k)*BCELL(1,3)
         a2=TI(1,k)*BCELL(2,1)+TI(2,k)*BCELL(2,2)+TI(3,k)*BCELL(2,3)
         a3=TI(1,k)*BCELL(3,1)+TI(2,k)*BCELL(3,2)+TI(3,k)*BCELL(3,3)
         cha=' '
         do j1=1,ngroup
            if(listat(k,j1)) cha='Y'
         end do
         write(*,'(2(x,i3),x,3(x,f10.5),2x,3(x,f10.5),3x,a)') &
                    k,j,a1,a2,a3,(TI(j1,k),j1=1,3),cha
      end do
   end do
   write(*,*)'Hit ENTER when done ...'
   read(*,*)
end subroutine show_atoms_tagged

subroutine show_atoms_in_plane()
!...............................................................................
!     Atomic positions with respect to the plane specified by the normal
! (vers1x,vers1y,vers1z) and the cnetral point (aCENTX,aCENTY,aCENTZ)
! are given for all atoms using the local coordinate system with:
!          X axis = (vers2x,vers2y,vers2z)
!          Y axis = (vers3x,vers3y,vers3z)
!...............................................................................
use param
use atoms
use menu
implicit none
real*8,parameter :: tiny=0.0001
real*8 x(3),cnt(3)
integer k,i,j,i1,j1
character cha*10
!
      cnt(1)=aCENTX
      cnt(2)=aCENTY
      cnt(3)=aCENTZ
!
      k=0
      do i=1,NSPEC
         write(*,'(3a,i5)')'.... Species <',Species(i),'> ......> ',i
         write(*,'(a)')' Tot   # |----------  Cartesian ----------|' &
              //'     X-plane    Y-plane    out-of-plane'
         do j=1,NspN(i)
            k=k+1

!______________ position x() of atom k in the coordinate system fixed to the
!          plane central point (aCENTX,aCENTY,aCENTZ) and axes X,Y,and Z=normal
!          (x(3) shows out-of-plane distance)
!
            x(1)=(TI(1,k)-aCENTX)*vers2x+(TI(2,k)-aCENTY)*vers2y+ &
                                          (TI(3,k)-aCENTZ)*vers2z
            x(2)=(TI(1,k)-aCENTX)*vers3x+(TI(2,k)-aCENTY)*vers3y+ &
                                          (TI(3,k)-aCENTZ)*vers3z
            x(3)=(TI(1,k)-aCENTX)*vers1x+(TI(2,k)-aCENTY)*vers1y+ &
                                          (TI(3,k)-aCENTZ)*vers1z

            if(abs(x(3)).lt.tiny) then
               cha='    NO    '
            else
               write(cha,'(f10.5)') x(3)
            end if
            write(*,'(2(x,i3),x,3(x,f10.5),2x,2(x,f10.5),x,a)') &
                 k,j,(TI(j1,k),j1=1,3),(x(i1),i1=1,2),cha
         end do
      end do
      write(*,*)'Hit ENTER when done ...'
      read(*,*)
end subroutine show_atoms_in_plane

subroutine at_in_box(Cent,Sides,Rec,Pos,iType,N_at,nt,Err)
!.................................................................
! Positions Pos() and species types iType() of all atoms fitted
! into the box (with Cent, Sides, Rec) are calculated.
!.................................................................
! nt - number of atoms found in the box
!.................................................................
use param
use atoms
implicit none
real*8 Pos(3,N_at),Cent(3),Sides(3,3),Rec(3,3),x(3)
integer iType(N_at),i,j,k,i1,j1,nat,nt,iErr,N_at
character iask
logical Err

      Err=.false.
      nat=0
      nt=0
      DO i1=1,NSPEC
         DO j1=1,NspN(i1)
            nat=nat+1

            do i=-1,1
               do j=-1,1
                  do k=-1,1
                     x(1)=i*DIRC(1,1)+j*DIRC(2,1)+k*DIRC(3,1)+TI(1,nat)
                     x(2)=i*DIRC(1,2)+j*DIRC(2,2)+k*DIRC(3,2)+TI(2,nat)
                     x(3)=i*DIRC(1,3)+j*DIRC(2,3)+k*DIRC(3,3)+TI(3,nat)
                     call ask_box2(x,Cent,Sides,Rec,iask)

                     if(iask.eq.'y') then
                        nt=nt+1
                        if(nt.eq.N_at) then
                           write(*,*) &
                                'ERROR! # of atoms in the box > ',N_at
                           Err=.true.
                           return
                        end if
                        Pos(1,nt)=x(1)
                        Pos(2,nt)=x(2)
                        Pos(3,nt)=x(3)
                        iType(nt)=i1
                     end if

                  end do
               end do
            end do

         END DO
      END DO
      write(*,*)'Number of atoms in the box found = ',nt
end subroutine at_in_box

subroutine reducn(x,DIRC,BCELL)
!..................................................................
!    It returns the point x() to the 0-th unit cell by extracting
! integer parts from the coefficients in the decomposition of x()
! with respect to DIRC and by making them positive
!..................................................................
implicit none
real*8 x(3),DIRC(3,3),BCELL(3,3),r(3)
real*8,parameter :: tiny=0.000001
integer i
!
!________ get coefficients r() in the decomposition of x() with
!         respect to DIRC;
!   These coefficients are iX/NGX, iY/NGY, iZ/NGZ, where (iX,iY,iZ)
! might be larger than NGX, NGY,NGZ and noninteger; this depends on
! the point x() which has nothing to do with the grid.
!
      r(1)=BCELL(1,1)*x(1)+BCELL(1,2)*x(2)+BCELL(1,3)*x(3)
      r(2)=BCELL(2,1)*x(1)+BCELL(2,2)*x(2)+BCELL(2,3)*x(3)
      r(3)=BCELL(3,1)*x(1)+BCELL(3,2)*x(2)+BCELL(3,3)*x(3)
!
!________ extract the integer part from the coefficients and make
!         the coefficients positive if they are not
!
      do i=1,3
         if(abs(r(i)-nint(r(i))).lt.tiny) r(i)=nint(r(i))
         r(i)=r(i)-int(r(i))+1
         r(i)=r(i)-int(r(i))
      end do
!
!________ return the point r() to the Cartesian representation, x()
!
      x(1)=r(1)*DIRC(1,1)+r(2)*DIRC(2,1)+r(3)*DIRC(3,1)
      x(2)=r(1)*DIRC(1,2)+r(2)*DIRC(2,2)+r(3)*DIRC(3,2)
      x(3)=r(1)*DIRC(1,3)+r(2)*DIRC(2,3)+r(3)*DIRC(3,3)
end subroutine reducn

subroutine normalize(x,y,z)
implicit none
real*8 znorm,x,y,z
   znorm=sqrt(x*x+y*y+z*z)
   x=x/znorm
   y=y/znorm
   z=z/znorm
end subroutine normalize

subroutine interpolate(R,BCELL,denval,grid)
!.......................................................
!	This subroutine does the internal interpolation
! of the point R within the eight cube vertices of the grid.
!.......................................................
! Every point in the grid R(i,j,k) is defined as:
!
!   R(i,j,k)=DIRC(1)*i/NGX+DIRC(2)*j/NGY+DIRC(3)*k/NGZ
!
! where i = 0,...,NGX-1; j = 0,...,NGY-1;
!       k = 0,...,NGZ-1  -  are used to numerate
! the grid.
!.......................................................
use param
implicit none
real*8 R(3), x(3), BCELL(3,3)
integer intx(3), cube(2,2,2,3),i1,i2,i3
real*8 grid(NGX,NGY,NGZ), frac(3),t5,t6,denval,t1,t2,t3,t4
integer,dimension(3) :: nmesh
      nmesh(1)=NGX ; nmesh(2)=NGY ; nmesh(3)=NGZ
!
!________ get coefficients x() in the decomposition of r() with
!         respect to DIRC;
!   These coefficients are iX/NGX, iY/NGY, iZ/NGZ, where (iX,iY,iZ)
! are exactly between 0 and NGX-1, 0 and NGY-1, etc. because of
! reducn() called previously
!
      x(1)=BCELL(1,1)*r(1)+BCELL(1,2)*r(2)+BCELL(1,3)*r(3)
      x(2)=BCELL(2,1)*r(1)+BCELL(2,2)*r(2)+BCELL(2,3)*r(3)
      x(3)=BCELL(3,1)*r(1)+BCELL(3,2)*r(2)+BCELL(3,3)*r(3)
!
!....... define position of x() in the net of grid points
!
      do i1=1,3
         x(i1)=nmesh(i1)*x(i1)
         frac(i1)=x(i1)-int(x(i1))
         intx(i1)=int(x(i1))
      end do
!
!___________ put x() inside a cube with vertices on the net and get
! its numbers in cube(i1,i2,i3,component), where (i1,i2,i3) give all
! possible 8 vertices of the cube (because of grid(iX+1,iY+1,iZ+1)
! defined for indices 1,...NGX, etc., the array cube gets +1):
!
      do i1=1,2
         do i2=1,2
            do i3=1,2
               cube(i1,i2,i3,1)=intx(1)+i1
               cube(i1,i2,i3,2)=intx(2)+i2
               cube(i1,i2,i3,3)=intx(3)+i3
               if(cube(i1,i2,i3,1).eq.(NGX+1)) cube(i1,i2,i3,1)=1
               if(cube(i1,i2,i3,2).eq.(NGY+1)) cube(i1,i2,i3,2)=1
               if(cube(i1,i2,i3,3).eq.(NGZ+1)) cube(i1,i2,i3,3)=1
            end do
         end do
      end do
!
!____________ interpolate the density in x() from that on the cube
!             vectices
!
      call lineint(cube,grid,1,1,1,2,1,1,frac(1),t1)
      call lineint(cube,grid,1,2,1,2,2,1,frac(1),t2)
      call lineint(cube,grid,1,1,2,2,1,2,frac(1),t3)
      call lineint(cube,grid,1,2,2,2,2,2,frac(1),t4)
      t5=t1+frac(2)*(t2-t1)
      t6=t3+frac(2)*(t4-t3)
      denval=t5+frac(3)*(t6-t5)
end subroutine interpolate

subroutine lineint(cube,grid,p1,p2,p3,q1,q2,q3,frac,val)
!.......................................................
!	A small subroutine to do linear interpolation between
!	two cube vertices
!.......................................................
use param
implicit none
integer cube(2,2,2,3),p1,p2,p3,q1,q2,q3
real*8 grid(NGX,NGY,NGZ),val,t1,t2,frac
      t1=grid(cube(p1,p2,p3,1),cube(p1,p2,p3,2),cube(p1,p2,p3,3))
      t2=grid(cube(q1,q2,q3,1),cube(q1,q2,q3,2),cube(q1,q2,q3,3))
      val=t1+frac*(t2-t1)
end subroutine lineint

subroutine vector3(DIRC)
!.....................................................................
! The plane case:
!   Choses a normal vector (vers1x,vers1y,vers1z) to the plane as well
! as two orthonormal vectors (vers2x,...) and (vers3x,...) in the
! plane.
!.....................................................................
! 'icase' gives the method used to specify the plane:
!  icase = 1 - by a normal vector;
!          2 - by 3 points lying in the plane.
!  angstr - gives a method of specifying coordinates
!.....................................................................
use menu
implicit none
real*8 Rbb(3),Rcc(3),DIRC(3,3),a
character answer
real*8,parameter :: tiny=0.00001
!
!____________ plane: a vector along the normal to the plane
!
!_________________ if the plane is chosen by the normal vector
      if(icase.eq.1) then
 7        write(*,*)'Give a vector (x,y,z) orthogonal to your plane:'
          read (*,*,err=7)  vers1x, vers1y, vers1z
!_________________ if the plane is chosen by three points
      else if(icase.eq.2) then
 21       write(*,*)'Give the 1st (A) point:'
          call givepoint(Ra(1),Ra(2),Ra(3),angstr)
          write(*,*)'Give the 2nd (B) point:'
          call givepoint(Rb(1),Rb(2),Rb(3),angstr)
          Rbb(1)=Rb(1)-Ra(1)
          Rbb(2)=Rb(2)-Ra(2)
          Rbb(3)=Rb(3)-Ra(3)
          write(*,*)'Give the 3rd (C) point:'
          call givepoint(Rc(1),Rc(2),Rc(3),angstr)
          Rcc(1)=Rc(1)-Ra(1)
          Rcc(2)=Rc(2)-Ra(2)
          Rcc(3)=Rc(3)-Ra(3)
          vers1x=Rbb(2)*Rcc(3)-Rbb(3)*Rcc(2)
          vers1y=Rbb(3)*Rcc(1)-Rbb(1)*Rcc(3)
          vers1z=Rbb(1)*Rcc(2)-Rcc(1)*Rbb(2)
          a=vers1x**2+vers1y**2+vers1z**2
          if(a.lt.tiny) then
             write(*,*)'Error! Your points lie on a line!'
             write(*,*)'Return/Try again (r,R/other char)?'
             read(*,'(a)') answer
             if(answer.ne.'r' .and. answer.ne.'R') go to 21
             return
          end if
       end if
!
!__________ normalize the normal and choose (rather arbitrarily)
!           two perpendicular vectors in the plane
!
      call normalize(vers1x,vers1y,vers1z)
      if (vers1x.eq.0.and.vers1y.eq.0) then
           vers2x=0.0
           vers2y=-1.0
           vers2z=0.0
      else
           vers2x=vers1y
           vers2y=-vers1x
           vers2z=0.0
      end if
      call normalize(vers2x,vers2y,vers2z)
      if (vers1x.eq.0.and.vers1y.eq.0) then
           vers3x=1.0
           vers3y=0.0
           vers3z=0.0
      else
           vers3x=-vers1x*vers1z
           vers3y=-vers1y*vers1z
           vers3z=vers1x**2+vers1y**2
      end if
      call normalize(vers3x,vers3y,vers3z)
end subroutine vector3

subroutine centralP(DIRC)
!.....................................................................
!     Asks for the central point (acentx,acenty,acentz) in the plane
!.....................................................................
!  icase = 1 - by a normal vector;
!          2 - by 3 points lying in the plane.
!  angstr - gives a method of specifying coordinates
!.....................................................................
use menu
implicit none
real*8 DIRC(3,3),alpha
!
!............. give the central point of the plain: ................
!
!______ if icase=2 (the plane was given by 3 points) and icase1=2,
!     the center of the triangle is chosen as the center of the plane
!
!______ in other cases the center is chosen in a general way and
!     (if icase=2 and icase1=1) the projection of it is made
!     to the plane with respect to the first point A (Ra):
!
      write(*,*)'Give CENTRAL POINT on the plane:'
      call givepoint(aCENTX,aCENTY,aCENTZ,angstr)
!
!________ removing a part normal to the plane (i.e. the projection)
      if(icase.eq.2) then
        write(*,*)'WARNING: The projection of this point on the'
        write(*,*)'      plane with respect to your A-point, '
        write(*,'(a3,2(f10.5,a1),f10.5,a2)') &
                '  (',Ra(1),',',Ra(2),',',Ra(3),'),'
        write(*,*)'      will be taken as the central point.'
        alpha=(acentx-Ra(1))*vers1x + (acenty-Ra(2))*vers1y + &
              (acentz-Ra(3))*vers1z
        acentx=acentx-alpha*vers1x
        acenty=acenty-alpha*vers1y
        acentz=acentz-alpha*vers1z
        write(*,*)'Your central point (in Angstroms) is given by:'
        write(*,'(a3,2(f10.5,a1),f10.5,a2)') &
                        '  (',acentx,',',acenty,',',acentz,').'
      end if
end subroutine centralP

subroutine givepoint2(aCENTX,aCENTY,angstr)
!...............................................................
!   Choose a point for the plot on the plane
!...............................................................
use param
use atoms
implicit none
real*8 x(2),aCENTX,aCENTY
character angstr*12
integer number
!
!__________ give coordinates in Angstroms
      if(angstr.eq.'<Angstroms> ') then
 5       write(*,*)'Give Cartesian coordinates as X,Y in Angstroms:'
         READ(*,*,err=5) aCENTX,aCENTY
!
!__________ give position via atomic number
      else if(angstr.eq.'<AtomNumber>') then
 30      write(*,*)'Give atomic number:'
         read(*,*,err=30) number
         if(number.lt.1.or.number.gt.NIONS) go to 31
         aCENTX=TI(1,number)
         aCENTY=TI(2,number)
         return
 31      write(*,*)'Error! Wrong atomic number! Try again!'
         go to 30
!
!__________ give through fractional coordinates
      else
 6       write(*,*)'Give fractional coordinates along A1,A2:'
         READ(*,*,err=11) x(1),x(2)
         go to 20
 11      write(*,*)'Error! Try again!'
         go to 6
!_________ transform to Cartesian coordinates
 20      aCENTX=x(1)*DIRC(1,1)+x(2)*DIRC(2,1)
         aCENTY=x(1)*DIRC(1,2)+x(2)*DIRC(2,2)
         return
      end if
end subroutine givepoint2

subroutine givepoint(aCENTX,aCENTY,aCENTZ,angstr)
!...............................................................
!   Choose a point for the plot
!...............................................................
use param
use atoms
implicit none
real*8 x(3),aCENTX,aCENTY,aCENTZ
character angstr*12
integer number
!
!__________ give coordinates in Angstroms
      if(angstr.eq.'<Angstroms> ') then
 5       write(*,*)'Give Cartesian coordinates as X,Y,Z in Angstroms:'
         READ(*,*,err=5) aCENTX,aCENTY,aCENTZ
!
!__________ give position via atomic number
      else if(angstr.eq.'<AtomNumber>') then
 30      write(*,*)'Give atomic number:'
         read(*,*,err=30) number
         if(number.lt.1.or.number.gt.NIONS) go to 31
         aCENTX=TI(1,number)
         aCENTY=TI(2,number)
         aCENTZ=TI(3,number)
         return
 31      write(*,*)'Error! Wrong atomic number! Try again!'
         go to 30
!
!__________ give through fractional coordinates
      else
 6       write(*,*)'Give fractional coordinates along A1,A2,A3:'
         READ(*,*,err=11) x(1),x(2),x(3)
         go to 20
 11      write(*,*)'Error! Try again!'
         go to 6
!_________ transform to Cartesian coordinates
 20      aCENTX=x(1)*DIRC(1,1)+x(2)*DIRC(2,1)+x(3)*DIRC(3,1)
         aCENTY=x(1)*DIRC(1,2)+x(2)*DIRC(2,2)+x(3)*DIRC(3,2)
         aCENTZ=x(1)*DIRC(1,3)+x(2)*DIRC(2,3)+x(3)*DIRC(3,3)
      end if
end subroutine givepoint
