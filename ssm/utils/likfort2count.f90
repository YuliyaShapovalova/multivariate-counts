	   subroutine likfort2count(n,T,L,x,beta,phi,seta,likfin,hpsmT,hpT,hptraj,hpsmtraj,&
	   hpmean)
	   use random, only: random_order
	   	 
	   implicit none 
       integer, intent(in):: n, T, L
       double precision, intent(in):: seta(n,n), beta(n,1), phi(n,n)
       integer, intent(in):: x(n,T)
       double precision,intent(out):: likfin,hpsmT(n,L), hpT(n,L), hptraj(n,T),&
       hpsmtraj(n,T), hpmean(n,T) 
       integer:: j, k(L),xh(n,L),i
       integer:: Factorial
       double precision:: LogFactorial
       double precision:: hp(n,L),hpsm(n,L,T),w(L,T+1), Lik(T+1), rn(n,L), wmax, wmin
       double precision:: muhp(n,1), betam(n,L), lambda(n,L), wh(n,L),hpall(n,L,T)
       double precision:: wsm(L,T+1), wnorm(L,T+1)
       
       w(:,:)=dble(1)/dble(L)
       muhp=reshape([0.0,0.0],[n,1])
       hp(:,:)=0d0
       Lik(1)=0d0
       likfin=0d0
       
       ! Main SMC loop
       do j=2,(T+1)
          call resample(w(:,j-1),L,L,k)
          call rmultnorm(n, L, muhp, seta, rn)
      
          hp(:,:)=hp(:,k)
          hp=matmul(phi,hp)+rn
       
          call comb(L,beta(1,1),beta(2,1),betam)
          call combint(L,x(1,j-1),x(2,j-1),xh)
       
          lambda=betam(:,:)*exp(hp(:,:))
          call dpoisson(xh,L,n,lambda,wh) 
          w(:,j)=sum(wh,dim=1)
       wmax=maxval(w(:,j))
       w(:,j)=exp(w(:,j)-wmax)
       wnorm(:,j)=w(:,j)/(sum(w(:,j),dim=1))
       Lik(j)=Lik(j-1)+wmax+log(sum(w(:,j)))-log(dble(L))
       hpall(:,:,j-1)=hp(:,:)
       end do
       hpT(:,:)=hpall(:,:,T)    
       wmin=sum(wh(2,:))
       likfin=Lik(T+1)
       
       !smoothing
       wsm(:,T+1)=w(:,T+1) 
       do j=T,1,-1
       ! smoothes weights
       wsm(:,j)=wnorm(:,j)*sum(wsm(:,j+1)/(sum(w(:,j))))
       !
              ! smoothed distribution
         do i=1,n
           hpsm(i,:,j)=sum(hpall(i,:,j)*wsm(:,j))
         enddo   
       enddo 
       hpsmT(:,:)=hpsm(:,:,T)
          do i=1,n
            do j=1,T
             hpsmtraj(i,j)=sum(hpsm(i,:,j))/L
             hpmean(i,j)=sum(hpall(i,:,j))/L
           enddo
           enddo  
                do i=1,n
            do j=1,T
             hptraj(i,j)=sum(hpall(i,:,j-1)*wnorm(:,j))
           enddo
        enddo      
         
	   end subroutine likfort2count
        
        INTEGER FUNCTION Factorial(n)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        INTEGER:: i, ANS
        ANS=1
        DO i=1,n
          ANS=ANS*i
        END DO
        Factorial=ANS
        END FUNCTION Factorial
        
        DOUBLE PRECISION FUNCTION LogFactorial(n)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        INTEGER:: i
        DOUBLE PRECISION:: ANS
        ANS=0
        DO i=1,n
          ANS=ANS + log(dble(i))
        END DO
        LogFactorial=ANS
        END FUNCTION LogFactorial
        
        subroutine dpoisson(x,L,n,par,d)
        double precision, intent(in):: par(n,L)
        integer, intent(in):: n
        integer, intent(in):: L
        integer, intent(in):: x(n,L) 
        double precision:: a(n,L), b(n,L), c(n,L)
        integer:: i, j, m
        integer:: Factorial
        double precision:: LogFactorial
        double precision, intent(out):: d(n,L)
        a(:,:)=1
        do j=1,n
            do i=1,L 
        a(j,i)=-LogFactorial(x(j,i))    
        b(j,i)=dble(x(j,i))*dble(log(par(j,i)))    
        c(j,i)=-par(j,i)
        d(j,i)=a(j,i)+b(j,i)+c(j,i) 
             end do    
        end do
        
        end subroutine dpoisson
		
		REAL FUNCTION FindDet(matrix, n)
		IMPLICIT NONE
		double precision, DIMENSION(n,n) :: matrix
		INTEGER, INTENT(IN) :: n
		double precision :: m, temp
		INTEGER :: i, j, k, l
		LOGICAL :: DetExists = .TRUE.
		l = 1
		!Convert to upper triangular form
		DO k = 1, n-1
			IF (matrix(k,k) == 0) THEN
				DetExists = .FALSE.
				DO i = k+1, n
					IF (matrix(i,k) /= 0) THEN
						DO j = 1, n
							temp = matrix(i,j)
							matrix(i,j)= matrix(k,j)
							matrix(k,j) = temp
						END DO
						DetExists = .TRUE.
						l=-l
						EXIT
					ENDIF
				END DO
				IF (DetExists .EQV. .FALSE.) THEN
					FindDet = 0
					return
				END IF
			ENDIF
			DO j = k+1, n
				m = matrix(j,k)/matrix(k,k)
				DO i = k+1, n
					matrix(j,i) = matrix(j,i) - m*matrix(k,i)
				END DO
			END DO
		END DO
	
		!Calculate determinant by finding product of diagonal elements
		FindDet = l
		DO i = 1, n
			FindDet = FindDet * matrix(i,i)
		END DO
	
	END FUNCTION FindDet
        
        subroutine comb(length,a,b,d)
            !!create a matrix 5*n with elements a,b,c in a column
            integer:: i
            integer:: length
            double precision::a,b,c,e,f
            double precision::d(2,length)

            do i=1,length
            d(1,i)=a
            d(2,i)=b
            enddo
        end subroutine comb
        
        subroutine combint(length,a,b,d)
            !!create a matrix 5*n with elements a,b,c in a column
            integer:: i
            integer:: length
            integer::a,b,c,e,f
            integer::d(2,length)

            do i=1,length
            d(1,i)=a
            d(2,i)=b
            enddo
        end subroutine combint
        
      subroutine rmultnorm(k, N, mu, S, Z)
      use random
      
      integer, intent(in):: k, N
      double precision, intent(in):: mu(k,1), S(k,k)
      double precision, intent(out):: Z(k,N)
      integer:: i, M
      integer:: IDIST, ISEED(4)
      double precision:: y(k,N), x(k*N), cholS(k,k), S1(k,k)
      
      !call init_random_seed()
      IDIST=3
      ISEED=591
      call DLARNV( IDIST, ISEED, k*N, X ) 
      Y=reshape(x,[k,N])
      Z=matmul(S,Y) 
      end subroutine rmultnorm
			
             
        SUBROUTINE SETAL(P,N,A,Q,W)
		 INTEGER, intent(in):: N
		 INTEGER, intent(out):: A(N)
		 double precision, intent(out):: Q(N)
		 double precision, intent(in):: P(N)
		 INTEGER:: I,J,NN,NP,S,W(N)
		 NN=0
		 NP=N+1
		 DO I=1,N
		   Q(I)=N*P(I)
		   IF (Q(I) .LT. 1.0) THEN
			 NN=NN+1
			 W(NN)=I
		   ELSE
			 NP=NP-1
			 W(NP)=I
		   ENDIF
      END DO
	   DO S=1,N-1
			I=W(S)
			J=W(NP)
			A(I)=J
			Q(J)=Q(J)+Q(I)-1.0
			IF(Q(J) .LT. 1.0) NP=NP+1
      END DO
			A(W(N))=W(N)
	   DO I=1,N
        Q(I)=Q(I)+I-1            	
	   END DO
	   END SUBROUTINE SETAL
	   
	   SUBROUTINE ALRV(A,Q,N,R,ANSW)
	   implicit none
	   !use random
	   
	   INTEGER, intent(in):: A(N), N, R
	   integer,parameter :: seed = 86456
	   double precision, intent(in):: Q(N)
	   INTEGER:: I, J
	   double precision::  U
	   INTEGER, intent(out):: ANSW(R)
	   
	   call srand(seed)
	   do J=1,R
	     !call RANDOM_NUMBER(uu)
	     U=dble(N*rand())
	     I=1+INT(U)
	     IF (U .LE. Q(I)) THEN
	       ANSW(J)=I
	     ELSE
	       ANSW(J)=A(I) 
	     ENDIF
	     end do
	   END SUBROUTINE ALRV
	   
		!subroutine resample(w,N,R,indx)
        ! 
		!use random
		!implicit none
		! integer, intent(in):: N,R
		! double precision, intent(in):: w(N)
		! double precision:: P(N), Q(N)
		! INTEGER:: A(N), windx(N), order(R)
		! integer, intent(out):: indx(R)
	
		! P=w/sum(w)
		! print*,'P',P
		! !P=real(w/size(w))
		!!create alias table
		! call SETAL(P,N,A,Q,windx)
		! !print*,'A2',A
		! !print*,'Q2',Q
		! !create index vector
		 
		! CALL ALRV(A,Q,N,R,indx)
	   	! call random_order(order, R)
        ! indx=indx(order)
		
		! end subroutine resample
		 
	  !Systematic resampling	 
	  subroutine resample(w,N,K,indx)
	  use random
	  use nrutil

	  implicit none
	  
	  !N - is, K - want
	  integer, intent(in):: N, K
	  double precision, intent(in):: w(N)
	  integer:: i, j, cumsum
	  double precision::  T(K), u, unod
	  double precision:: n1, arth, Q(N), w1(N)
	  double precision:: a1, an, d
	  double precision :: a(K)
	  integer:: order(K)
	  integer, intent(out):: indx(K)

	  w1 = w / sum(w)
      call cumsm(N,w1,Q)
	  !Q = cumsum(w1)


      d=(1.0-1.0/dble(K))/(dble(K)-1.0)
      a(1)=0.0
      do i=2,K
      a(i)=a(i-1)+d
      end do

	  call RANDOM_NUMBER(u)
	  !unod=unod/K
	  T = a + dble(u)/dble(K)

	  i = 1
	  j = 1
	 		  !was do while(i<=N .and. j<=M)
		 
	   do while(i .LE. K .and. j .LE. N )
	   do while (Q(j) .LT. T(i)) 
	     j=j+1
        end do
          indx(i)=j	 
	      i=i+1  
       end do
      !call random_order(order, K)
      !indx=indx(order)
	  end subroutine resample
	  
	  subroutine cumsm(m,a,b)
        integer, intent(in):: m
        double precision, intent(in):: a(m) 
        double precision, intent(out):: b(m)
        integer:: i
        
        b(1)=a(1)
                
        do i=2,m
            b(i)=a(i)+b(i-1)
        end do
		
       end subroutine cumsm
             
             subroutine inverse(a,c,n)
			!============================================================
			! Inverse matrix
			! Method: Based on Doolittle LU factorization for Ax=b
			! Alex G. December 2009
			!-----------------------------------------------------------
			! input ...
			! a(n,n) - array of coefficients for matrix A
			! n      - dimension
			! output ...
			! c(n,n) - inverse matrix of A
			! comments ...
			! the original matrix a(n,n) will be destroyed 
			! during the calculation
			!===========================================================
			implicit none 
			integer n
			double precision a(n,n), c(n,n)
			double precision L(n,n), U(n,n), b(n), d(n), x(n)
			double precision coeff
			integer i, j, k

			! step 0: initialization for matrices L and U and b
			! Fortran 90/95 aloows such operations on matrices
			L=0.0
			U=0.0
			b=0.0

			! step 1: forward elimination
			do k=1, n-1
			   do i=k+1,n
				  coeff=a(i,k)/a(k,k)
				  L(i,k) = coeff
				  do j=k+1,n
					 a(i,j) = a(i,j)-coeff*a(k,j)
				  end do
			   end do
			end do

			! Step 2: prepare L and U matrices 
			! L matrix is a matrix of the elimination coefficient
			! + the diagonal elements are 1.0
			do i=1,n
			  L(i,i) = 1.0
			end do
			! U matrix is the upper triangular part of A
			do j=1,n
			  do i=1,j
				U(i,j) = a(i,j)
			  end do
			end do

			! Step 3: compute columns of the inverse matrix C
			do k=1,n
			  b(k)=1.0
			  d(1) = b(1)
			! Step 3a: Solve Ld=b using the forward substitution
			  do i=2,n
				d(i)=b(i)
				do j=1,i-1
				  d(i) = d(i) - L(i,j)*d(j)
				end do
			  end do
			! Step 3b: Solve Ux=d using the back substitution
			  x(n)=d(n)/U(n,n)
			  do i = n-1,1,-1
				x(i) = d(i)
				do j=n,i+1,-1
				  x(i)=x(i)-U(i,j)*x(j)
				end do
				x(i) = x(i)/u(i,i)
			  end do
			! Step 3c: fill the solutions x(n) into column k of C
			  do i=1,n
				c(i,k) = x(i)
			  end do
			  b(k)=0.0
			end do
			end subroutine inverse
			
			
			subroutine init_random_seed()
            use iso_fortran_env, only: int64
            implicit none
            integer, allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8), pid
            integer(int64) :: t
          
            call random_seed(size = n)
            allocate(seed(n))
            ! First try if the OS provides a random number generator
            open(newunit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
            if (istat == 0) then
               read(un) seed
               close(un)
            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(t)
               if (t == 0) then
                  call date_and_time(values=dt)
                  t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24_int64 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
               end if
               pid = getpid()
               t = ieor(t, int(pid, kind(t)))
               do i = 1, n
                  seed(i) = lcg(t)
               end do
            end if
            call random_seed(put=seed)
          contains
            ! This simple PRNG might not be good enough for real work, but is
            ! sufficient for seeding a better PRNG.
            function lcg(s)
              integer :: lcg
              integer(int64) :: s
              if (s == 0) then
                 s = 104729
              else
                 s = mod(s, 4294967296_int64)
              end if
              s = mod(s * 279470273_int64, 4294967291_int64)
              lcg = int(mod(s, int(huge(0), int64)), kind(0))
            end function lcg
          end subroutine init_random_seed