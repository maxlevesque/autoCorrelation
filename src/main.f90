! This program computes the autocorrelation of a list of sites during a trajectory in time.
! It is written by Maximilien Levesque while in postdoc in the group of Mathieu Salanne
! at UPMC Univ Paris 06, CNRS, ESPCI, UMR 7195, PECSA, F-75005 Paris, France.

PROGRAM autoCorrelation

    IMPLICIT NONE
    
    CHARACTER(LEN("acf.out")) :: outputFile = "acf.out"
    INTEGER :: Nat,i,nbTimeStepsInTraj,iostat,dt,t
    INTEGER, PARAMETER :: x=1, y=2, z=3
    DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: r ! position of site i at timestep t
    DOUBLE PRECISION :: time1, time0
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: acf
    CHARACTER(LEN=300) :: arg,trajectoryFileName,algo
  
    CALL CPU_TIME(time0)

    ! read all arguments that MUST be given at execution
    CALL readArguments(Nat,trajectoryFileName,algo)
    i =NbOfLinesInTraj()     ! deduce the number of timesteps in the trajectory from the number of lines in the trajectory file
    CALL testConsistencyOfAtomNumber(i,Nat)
    nbTimeStepsInTraj = i/Nat

    PRINT*,'I found',nbTimeStepsInTraj,' timesteps in your trajectory called ',trim(adjustl(trajectoryFileName))
    PRINT*,'Please be patient... Everything looks fine...'

    ! read vector of all sites i at all timesteps t
    ALLOCATE( r(Nat,nbTimeStepsInTraj,x:z), SOURCE=0.d0 )
    CALL opentraj
    DO t=1,nbTimeStepsInTraj
        IF( mod(t,1000)==1 .AND. t/=1) PRINT*,"I've read ",t-1," timesteps among ",nbTimeStepsInTraj
        DO i=1,Nat
            READ(10,*) r(i,t,x), r(i,t,y), r(i,t,z)
        END DO
    END DO
    CALL closetraj
    PRINT*,"I've read all the trajectories :). Be patient..."

    ! compute autocorrelation function acf(dt)= <v_i(t).v_i(t+dt)>_{i,t}   where . is the scalar product
    ALLOCATE( acf(0:nbTimeStepsInTraj-1), SOURCE=0.d0 )

    IF (TRIM(ADJUSTL(algo))=="bruteforce") THEN
        CALL bruteforce(r,acf)
    ELSE IF (TRIM(ADJUSTL(algo))=="fourierspace") THEN
!~         CALL fourierspace(r,acf)
    ELSE
        STOP "I do not understand the algorithm you ask for. I only recognize bruteforce and fourierspace. STOP"
    END IF

    ! acf(t) will be written in file unit 11
    OPEN(11,FILE=outputfile)
        DO dt = 0, nbTimeStepsInTraj-1
            WRITE(11,*) dt, acf(dt)
        END DO
    CLOSE(11)

    CALL CPU_TIME(time1)
    IF ((time1-time0)>120.d0) THEN
        PRINT*,"-- Finished in ",NINT((time1-time0)/60.d0)," min. GGHF ;) --"
    ELSE
        PRINT*,"-- Finished in ",NINT(time1-time0)," s. GGHF ;) --"
    END IF
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CONTAINS
    
    
        SUBROUTINE testConsistencyOfAtomNumber(i,Nat)
            INTEGER, INTENT(IN) :: i,Nat
            IF( MODULO(i,Nat)/=0 ) THEN
                PRINT*,"There is an inconsistency between the number of lines in file ",TRIM(ADJUSTL(trajectoryFileName))
                PRINT*,"and the number of atoms, ",Nat
                PRINT*,"since modulo(nlines,Nat)/=0"
                STOP "I stop -- :("
            END IF
        END SUBROUTINE testConsistencyOfAtomNumber
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        SUBROUTINE opentraj
            CALL inquireFileExistence(trajectoryFileName)
            ! read positions
            OPEN(10, FILE=trajectoryFileName,STATUS='old',IOSTAT=iostat)
            IF (iostat /= 0) THEN
                WRITE(*,*) 'File open failed for',trajectoryFileName
                WRITE(*,*) 'Error code is ', iostat
                STOP
            END IF
        END SUBROUTINE opentraj
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
        
        SUBROUTINE closetraj
            CLOSE(10)
        END SUBROUTINE closetraj
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        SUBROUTINE inquireFileExistence(fileName)
            CHARACTER(LEN=*), INTENT(IN) :: fileName
            INTEGER, PARAMETER :: stderr = 0
            LOGICAL :: exist
            INQUIRE(FILE=fileName, EXIST=exist)
            IF( .NOT. exist) THEN
                WRITE(stderr,*) "YOUR ERROR (not mine ;): The file ",fileName," does not exist. It should."
                STOP
            END IF
        END SUBROUTINE inquireFileExistence
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        FUNCTION NbOfLinesInTraj()
            INTEGER :: NbOfLinesInTraj
            CALL opentraj
            ! computes the number of lines in traj.in and deduces the number of timesteps
            NbOfLinesInTraj = -1
            DO WHILE (iostat == 0)
                READ(10,*,IOSTAT=iostat)
                NbOfLinesInTraj = NbOfLinesInTraj + 1
            END DO
            CALL closetraj
        END FUNCTION NbOfLinesInTraj
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
        
    !~     function NbOfLinesInTraj(filename)
    !~         character(len=*), intent(in) :: filename
    !~         integer :: NbOfLinesInTraj
    !~         character(len=180) :: cmd, msg
    !~         character(len=*), parameter :: tmpfilename = "000098767612398712309.TMP"
    !~         cmd="cat "//trim(adjustl(filename))//" | wc -l > "//tmpfilename
    !~         call execute_command_line(trim(adjustl(cmd)), wait=.true.)
    !~         call inquireFileExistence(tmpfilename)
    !~         open(86,file=tmpfilename)
    !~         read(86,*)NbOfLinesInTraj
    !~         close(86)
    !~         call execute_command_line("rm "//tmpfilename)
    !~     end function
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        SUBROUTINE readArguments(Nat,trajectoryFileName,algo)
            IMPLICIT NONE
            INTEGER, INTENT(OUT) :: Nat
            CHARACTER(LEN=*),INTENT(OUT) :: trajectoryFileName,algo
            CALL GET_COMMAND_ARGUMENT(1,arg,STATUS=i)
                IF (i<0) THEN
                    STOP "STOP. The length of the argument is too big for me :( "
                ELSE IF (i>0) THEN
                    STOP "Argument retrieval failed. You should execute the program with the number of atoms as 1st argument"
                END IF
                READ(arg,*) Nat
                PRINT*,"Number of atoms: ",Nat
            CALL GET_COMMAND_ARGUMENT(2,arg,STATUS=i)
                IF (i<0) THEN
                    STOP "STOP. The length of the argument is too big for me :( "
                ELSE IF (i>0) THEN
                    STOP "Argument retrieval failed. You should execute the program with the filename as 2nd argument"
                END IF
                READ(arg,*) i
                SELECT CASE (i)
                CASE (1)
                    algo='bruteforce'
                CASE (2)
                    algo="fourierspace"
                CASE DEFAULT
                    STOP "I did not understand the argument for the algo. Should be 1 (bruteforce) or 2 (fourierspace)"
                END SELECT
                PRINT*,"I'll use algorithm ",TRIM(ADJUSTL(algo))
            CALL GET_COMMAND_ARGUMENT(3,arg,STATUS=i)
                IF (i<0) THEN
                    STOP "STOP. The length of the argument is too big for me :( "
                ELSE IF (i>0) THEN
                    STOP "Argument retrieval failed. You should execute the program with the algorithm as 3nd argument"
                END IF
                trajectoryFileName = TRIM(ADJUSTL(arg))
                PRINT*,"I'll read the trajectory from ",TRIM(ADJUSTL(trajectoryFileName))
        END SUBROUTINE readArguments
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        SUBROUTINE bruteforce(r,acf)
            IMPLICIT NONE
            DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: r
            DOUBLE PRECISION, DIMENSION(0:), INTENT(OUT) :: acf
            DOUBLE PRECISION, DIMENSION(SIZE(r,2),SIZE(r,3)) :: ri ! position of site i at timestep t
            INTEGER :: i,dt,nt,Nat,nbTimeStepsInTraj
            INTEGER, PARAMETER :: x=1, y=2, z=3
            DOUBLE PRECISION :: time0, time1, remainingTimeInSec
            CALL CPU_TIME(time0)
            Nat = SIZE(r,1)   !r(Nat,nbTimeStepsInTraj,x:z)
            nbTimeStepsInTraj = SIZE(r,2)
            DO i= 1, Nat
                ri = r(i,:,:)
                DO dt = 0, nbTimeStepsInTraj-1
                    nt = nbTimeStepsInTraj-dt
                    acf(dt)= acf(dt)+&
                                (SUM(    ri(1:nt,x)*ri(1+dt:nt+dt,x)  &
                                        +ri(1:nt,y)*ri(1+dt:nt+dt,y)  &
                                        +ri(1:nt,z)*ri(1+dt:nt+dt,z))  )/DBLE(nt)
                END DO
    
                CALL CPU_TIME(time1)
                IF(MODULO(i,MAX(INT(Nat*0.1),1))==0) THEN
                    remainingTimeInSec = DBLE(Nat-i)*(time1-time0)/DBLE(i)
                    IF (remainingTimeInSec>60.d0) THEN
                        PRINT*,'Estimated remaining time ≈ ',NINT(remainingTimeInSec/60.d0),' min'
                    ELSE
                        PRINT*,'Estimated remaining time ≈ ',NINT(remainingTimeInSec),' s'
                    END IF
                END IF
            END DO
            acf = acf/DBLE(Nat)
        END SUBROUTINE bruteforce

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END PROGRAM autoCorrelation
