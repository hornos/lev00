!===============================================================
! (1) INVERS of A | <==== this is what is here: TWO routines!!!
! (2) A*X=B       |
!===============================================================

subroutine linear(A,B,INDX,VV,N,NP,Err)
!......................................................................
!   Linear system of equations A*X=B is solved by the LU-decomposition.
! INPUT:
!            A(NP,NP)  - real matrix on the left
!            B(N) - vector on the right
!            NP - the physical dimension of A and B
!            N  - their actual dimension.
! OUTPUT:
!            B(N) - vector of solution
!......................................................................
! The initial matrix A and vector B will be destroyed at the end.
!......................................................................
! Two working arrays are used for the calculation:
!    INDX - an integer vector needed for the pivoting;
!    VV   - a real vector used to store the scaling.
!......................................................................
! Was taken from the book:
!    "Numerical recipies. The art of scientific computing (Fortran
!     version)." by W.H.Press et al (Cambridge Univ.Press, Cambridge,
!     1989), page 38.
!......................................................................

implicit none
real*8 A(NP,NP),VV(N),B(N),D
integer INDX(N),iErr,N,NP
logical Err
      Err=.false.
!________ decompose the matrix A just once
      call LUDCMPr(A,N,NP,INDX,D,VV,iErr)
      if(iErr.eq.1) then
         Err=.true.
         return
      end if
!________ solve the equation A*X=B.
      call LUBKSBr(A,N,NP,INDX,B)
end subroutine linear

subroutine invers(A,Y,INDX,VV,N,NP,Err)
!......................................................................
!  Invers Y(NP,NP) of the real matrix A(NP,NP)
!......................................................................
! NP - the physical dimension of the square matrices A and Y,
! N  - their actual dimension.
!......................................................................
!..... The initial matrix A will be destroyed at the end...............
! Two working arrays are used for the calculation:
!    INDX - an integer vector needed for the pivoting;
!    VV   - a real vector used to store the scaling.
!......................................................................
! Was taken from the book:
!    "Numerical recipies. The art of scientific computing (Fortran
!     version)." by W.H.Press et al (Cambridge Univ.Press, Cambridge,
!     1989), page 38.
!......................................................................
implicit none
integer N,NP,INDX(N),i,j,iErr
real*8 Y(NP,NP),A(NP,NP),VV(N),D
logical Err
      Err=.false.
!________ set the identity matrix
      do i=1,N
         do j=1,N
            Y(i,j)=0.0
         end do
         Y(i,i)=1.0
      end do
!________ decompose the matrix A just once
      call LUDCMPr(A,N,NP,INDX,D,VV,iErr)
      if(iErr.eq.1) then
         Err=.true.
         return
      end if
!________ find inverse by columns j of the Y matrix
      do j=1,N
         call LUBKSBr(A,N,NP,INDX,Y(1,j))
      end do
end subroutine invers

subroutine ludcmpR(a,n,np,indx,d,vv,iErr)
implicit none
integer n,np
real*8, parameter :: tiny=1.0e-20
real*8 a(np,np),vv(n),aamax,sum,dum,d
integer indx(n),iErr,i,j,k,imax
      iErr=0
      d=1.
!
      do i=1,n
         aamax=0.
         do j=1,n
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
         end do
         if (aamax.eq.0.) then
            write(*,*)' FATAL! Singular matrix!'
            iErr=1
            return
         end if
         vv(i)=1./aamax
      end do
!
!.......... the principal loop over columns j=1,...,n
      do j=1,n
!
         if (j.gt.1) then
            do i=1,j-1
               if (i.gt.1)then
                  sum=a(i,j)
                  do k=1,i-1
                     sum=sum-a(i,k)*a(k,j)
                  end do
                  a(i,j)=sum
               endif
            end do
         end if
!
         aamax=0.
!
         do i=j,n
            sum=a(i,j)
            if (j.gt.1)then
               do k=1,j-1
                  sum=sum-a(i,k)*a(k,j)
               end do
               a(i,j)=sum
            endif
            dum=vv(i)*abs(sum)
            if (dum.ge.aamax) then
               imax=i
               aamax=dum
            endif
         end do
!______________ interchanging of rows
         if (j.ne.imax)then
            do k=1,n
               dum=a(imax,k)
               a(imax,k)=a(j,k)
               a(j,k)=dum
            end do
            d=-d
            vv(imax)=vv(j)
         endif
         indx(j)=imax
!
         if(j.ne.n)then
            if(a(j,j).eq.0.) a(j,j)=tiny
            dum=1./a(j,j)
            do i=j+1,n
               a(i,j)=a(i,j)*dum
            end do
         endif

      end do
      if(a(n,n).eq.0.) a(n,n)=tiny
end subroutine ludcmpR

subroutine lubksbR(a,n,np,indx,b)
implicit none
integer n,np,ii,i,ll,j
real*8 a(np,np),b(n),sum
integer indx(n)
      ii=0
!....... forward substitution
      do i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         if (ii.ne.0)then
            do j=ii,i-1
               sum=sum-a(i,j)*b(j)
            end do
         else if (sum.ne.0.) then
            ii=i
         end if
         b(i)=sum
      end do
!...... backsubstitution
      do i=n,1,-1
         sum=b(i)
         if(i.lt.n)then
            do j=i+1,n
               sum=sum-a(i,j)*b(j)
            end do
         end if
         b(i)=sum/a(i,i)
      end do
end subroutine lubksbR

