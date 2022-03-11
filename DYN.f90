PROGRAM GREENBERG                     
  USE, INTRINSIC                :: ISO_Fortran_env
  IMPLICIT NONE
  INTEGER, PARAMETER            :: NEQN=319            ! NUMBER OF ROIS (REGIONS OF INTEREST)
  INTEGER, PARAMETER            :: NTIME=2000, NT=0    ! NT=NUMBER OF T'S
  REAL*8,  DIMENSION(31)        :: T
  INTEGER, PARAMETER            :: p1=1000, p2=1000
  INTEGER, PARAMETER            :: r1=6, r2=361
  REAL*8,  DIMENSION(NEQN)      :: SC_P, SC_F          ! STATE VARIABLE S={0,0,1} FOR TWO TIME STEPS, F=FUTURE, P=PRESENT 
  INTEGER, DIMENSION(NEQN)      :: SD_P, SD_F          ! DISCRET VARIABLE S={-1,0,1} FOR TWO TIME STEPS
  REAL*8,  DIMENSION(NEQN,NEQN) :: M                   ! WEIGHTED HUMAM CONNECTOME 
  INTEGER                       :: I, J, TIME_STEP, SEMENTE, T_STEP, IPRINT
  REAL*8                        :: A, B, C
  INTEGER                       :: I1, J1, I2, K, L, MARCADOR
  INTEGER, DIMENSION(NEQN)      :: ARRAY, ARRAY_BINARY
  CHARACTER*20                  :: OUTFILE1, OUTFILE2, OUTFILE3, filename
  CHARACTER*5                   :: X1, X2
  CHARACTER(len=8)              :: fmt1, fmt2           !format descriptor
  INTEGER                       :: N_FMRI, I_T, I_THR, N_THR
  
  !--------T PARAMETER-------------------------------------------
	T(1)=0.0000d0
	T(2)=0.0132d0
	T(3)=0.0264d0
	T(4)=0.0396d0
	T(5)=0.0528d0
	T(6)=0.0660d0
	T(7)=0.0792d0
	T(8)=0.0924d0
	T(9)=0.1056d0
	T(10)=0.1105d0
	T(11)=0.1155d0
	T(12)=0.1188d0
	T(13)=0.1221d0
	T(14)=0.1287d0
	T(15)=0.1320d0
	T(16)=0.1353d0
	T(17)=0.1386d0
	T(18)=0.1419d0
	T(19)=0.1452d0
	T(20)=0.1485d0
	T(21)=0.1518d0
	T(22)=0.1551d0
	T(23)=0.1584d0
	T(24)=0.1617d0
	T(25)=0.1650d0
	T(26)=0.1683d0
	T(27)=0.1716d0
	T(28)=0.1782d0
	T(29)=0.1848d0
	T(30)=0.1914d0
	T(31)=0.1980d0

  N_THR = 31      ! number of T's  
  N_FMRI= 10      ! number of trials
  !--------------------------------------------------------------
  fmt1 = '(I3.3)' 
  fmt2 = '(F5.3)' 
  
  OPEN(UNIT=150, FILE='con_mat.txt', STATUS='OLD',    &
       ACCESS='SEQUENTIAL', FORM='FORMATTED', ACTION='READ' )
  
  DO J=1,NEQN
     READ (150,*)  (M(J,I),I=1,NEQN)
  ENDDO
 
  DO I_THR=1,N_THR
     
     write (X2,fmt2) T(I_THR)
     
     DO I_T=1,N_FMRI
        
        write (X1,fmt1) I_T                 ! converting integer to string using a 'internal file'
        filename = 'DYN_T='//trim(X2)//'_N='//trim(X1)//'.txt'
        
        OPEN(UNIT=I_T,FILE=filename)
        
        !---------starting the dynamics-------------------------------
        SC_P=0.0D0
        !-------------------------------------------------------------
        CALL SYSTEM_CLOCK(SEMENTE)
        CALL SRAND(SEMENTE)
        
        DO T_STEP=0,NT
           
           B=0.0D0
           !SD_P=-1
           DO I=1,NEQN
              SD_P(I)=-1 + FLOOR(RAND()*2)
           ENDDO
           
           DO TIME_STEP=1,NTIME
              
              DO I=1,NEQN
                 !  Q=0, R=-1, E=1
                 !  Q --> E
                 IF(SD_P(I) .EQ. 0) THEN
                    IF(1+floor(rand()*p1) .LT. r1) THEN     ! probability = r1/p1   
                       SC_F(I)=1.0D0
                       SD_F(I)=1
                    ELSE
                       A=0.0D0
                       DO J=1,NEQN
                          IF(M(I,J) .NE. 0.0D0 .AND. SD_P(J) .EQ. 1) A=A+M(I,J) 
                       ENDDO
                       IF(A .GE. T(I_THR)) THEN
                          SC_F(I)=1.0D0
                          SD_F(I)=1
                       ELSE
                          SC_F(I)=0.0D0
                          SD_F(I)=0
                       ENDIF
                    ENDIF
                 ENDIF
                 
                 ! E --> R
                 IF(SD_P(I) .EQ. 1) THEN
                    SC_F(I)=0.0D0       !WITH PROBABILITY 1
                    SD_F(I)=-1          !WITH PROBABILITY 1
                 ENDIF
                 
                 ! R --> Q
                 IF(SD_P(I) .EQ. -1) THEN
                    IF(1+floor(rand()*p2) .LT. r2) THEN    ! probability = r1/p1
                       SC_F(I)=0.0D0
                       SD_F(I)=0
                    ELSE
                       SC_F(I)=0.0D0
                       SD_F(I)=-1
                    ENDIF
                 ENDIF
                 
              ENDDO
              
              WRITE(I_T,*) (SC_P(IPRINT),IPRINT=1,NEQN)
              !-------------------------------------------------------------------------------------
              !-------------------------------------------------------------------------------------
              SD_P=SD_F
              SC_P=SC_F
           ENDDO
        END DO
        CLOSE(I_T)
     END DO
  END DO
END PROGRAM GREENBERG

