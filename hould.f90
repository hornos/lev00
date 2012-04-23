subroutine diag(F,A,E,N,Prnt)
implicit none
!....... diagonalization of a square matrix F
real*8 F(N,N),A(N),E(N),epsH
logical Prnt
integer N,j,i
!.......................................................................
!................... I N P U T   F O R   H O U L D:.....................
!  F  - initial square matrix (N*N) which must be diagonalized;
!       only left thriangular part of it must actually be filled.
!  E  - 1-dim. vector of eigenvalues
!  One temporary 1-dim. vector AU of the dimension N
!  N - dimension of the square Fock matrix
!  epsH - precision
!................... O U T P U T  F R O M  H O U L D:..................
!  F -  eigenvectors will be allocated
!.......................................................................

      epsH=1.0E-7
      CALL HOULD (E,F,N,epsH,A)

!........... Write eigenvalues and eigenvectors
      if(Prnt) then
         write (*,*) 'Eigenvalues => ',(E(i),I=1,N)
         do j=1,N
            write(*,*) (F(j,i),i=1,N)
         end do
      end if
end subroutine diag

SUBROUTINE HOULD (D,Z,N,EPS,E)
implicit none
real*8 D(N),E(N),Z(N,N),EPS,f,g,h,p,di,b,em,r,c,s,ei
integer N,i,j,l,in,k,m,im1,im,k1

      do IN=2,N
         I=N+2-IN
         L=I-2
         F=Z(I,I-1)
         G=0.
         if(l.lt.1) go to 20
         do K=1,L
            G=G+Z(I,K)**2
         end do
20       H=F*F+G
         IF(G.GE.1.E-14) GO TO 30
         E(I)=F
         H=0.
         GO TO 100
30       L=L+1
         G= SQRT(H)
         IF(F.LT.0) GO TO 40
         G=-G
40       E(I)=G
         H=H-F*G
         Z(I,I-1)=F-G
         F=0.
         do J=1,L
            Z(J,I)=Z(I,J)/H
            G=0.
            do K=1,J
               G=G+Z(J,K)*Z(I,K)
            end do
            K1=J+1
            IF(K1.GT.L) GO TO 65
            do K=K1,L
               G=G+Z(K,J)*Z(I,K)
            end do
65          E(J)=G/H
            F=F+G*Z(J,I)
         end do
         P=F/(H+H)
         do J=1,L
            F=Z(I,J)
            G=E(J)-P*F
            E(J)=G
            do K=1,J
               Z(J,K)=Z(J,K)-F*E(K)-G*Z(I,K)
            end do
         end do
100      D(I)=H
      end do
!     TRIDIAGONALIZACIJAS BEIGAS
      D(1)=0.
      E(1)=0.
      DO I=1,N
         L=I-1
         DI=D(I)
         IF(DI.EQ.0.) GO TO 145
         DO J=1,L
            G=0.
            DO K=1,L
               G=G+Z(I,K)*Z(K,J)
            end do
            DO K=1,L
               Z(K,J)=Z(K,J)-G*Z(K,I)
            end do
         end do
145      D(I)=Z(I,I)
         Z(I,I)=1.
         IF(L.LT.1) GO TO 155
         DO J=1,L
            Z(I,J)=0.
            Z(J,I)=0.
         end do
155   end do
!     SAFORMETA P MATRICA
      DO I=2,N
         E(I-1)=E(I)
      end do
      E(N)=0.
      B=0.
      F=0.
      DO L=1,N
         J=0
         H=( ABS(D(L))+ ABS(E(L)))*EPS
         IF(B.LT.H) B=H
         DO M=L,N
            EM= ABS(E(M))
            IF(EM.LE.B) GO TO 190
         end do
190      IF(M.EQ.L) GO TO 320
200      IF(J.NE.30) GO TO 220
         GO TO 450
220      J=J+1
         EM=E(L)
         P=(D(L+1)-D(L))/(EM+EM)
         R= SQRT(P*P+1.)
         IF(P) 230,240,240
230      EM=P-R
         GO TO 250
240      EM=P+R
250      H=D(L)-E(L)/EM
         DO  I=L,N
            D(I)=D(I)-H
         end do
         F=F+H
         P=D(M)
         C=1.
         S=0.
         IM1=M-1
         IF(IM1.LT.L) GO TO 315
         DO IM=L,IM1
            I=IM1+L-IM
            EI=E(I)
            G=C*EI
            H=C*P
            IF( ABS(P)- ABS(EI)) 280,270,270
270         C=EI/P
            R= SQRT(C*C+1.)
            E(I+1)=S*P*R
            S=C/R
            C=1./R
            GO TO 290
280         C=P/EI
            R= SQRT(C*C+1.)
            E(I+1)=S*EI*R
            S=1./R
            C=C/R
290         DI=D(I)
            P=C*DI-S*G
            DI=C*G+S*DI
            D(I+1)=H+S*DI
            DO K=1,N
               H=Z(K,I+1)
               EI=Z(K,I)
               Z(K,I+1)=S*EI+C*H
               Z(K,I)=C*EI-S*H
            end do
         end do
315      E(L)=S*P
         D(L)=C*P
         EI= ABS(E(L))
         IF(EI.GT.B) GO TO 200
320      D(L)=D(L)+F
      end do
!     D-IPASVERTIBAS,Z-IPASVEKTORI
      I=0
  340 I=I+1
      IF(I.GE.N) GO TO 360
      H=D(I)
      F=D(I+1)
      IF(H.LE.F) GO TO 340
      D(I)=F
      D(I+1)=H
      DO K=1,N
         H=Z(K,I)
         Z(K,I)=Z(K,I+1)
         Z(K,I+1)=H
      end do
      IF(I.NE.1) I=I-2
      GO TO 340
360   CONTINUE
450   CONTINUE
end subroutine hould

