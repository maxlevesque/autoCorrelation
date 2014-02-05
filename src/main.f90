! This program computes the autocorrelation of a list of sites during a trajectory in time.
! It is written by Maximilien Levesque while in postdoc in the group of Mathieu Salanne
! at UPMC Univ Paris 06, CNRS, ESPCI, UMR 7195, PECSA, F-75005 Paris, France.

PROGRAM autoCorrelation

    IMPLICIT NONE
    
    CHARACTER(LEN("acf.out")) :: outputFile = "acf.out"
    INTEGER :: Nat,i,nbTimeStepsInTraj,iostat,dt,t,nt
    INTEGER, PARAMETER :: x=1, y=2, z=3
    DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: r ! position of site i at timestep t
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ri
    DOUBLE PRECISION :: time1, time0, remainingTimeInSec
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: acf
    CHARACTER(LEN=300) :: arg, trajectoryFileName

    ! read all arguments that MUST be given at execution
    call readArguments(Nat,trajectoryFileName)

    ! deduce the number of timesteps in the trajectory from the number of lines in the trajectory file
    nbTimeStepsInTraj = NbOfLinesInTraj()/Nat

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

!~     ! compute autocorrelation function acf(dt)= <v_i(t).v_i(t+dt)>_{i,t}   where . is the scalar product
    ALLOCATE( ri(nbTimeStepsInTraj,x:z), SOURCE=0.d0 )
    ALLOCATE( acf(0:nbTimeStepsInTraj-1), SOURCE=0.d0 )
    CALL CPU_TIME(time0)

    PRINT*,"I'm now brutally computing the velocity autocorrelation function, which is its auto-cross-correlation."
    DO i= 1, Nat
        ri = r(i,:,:)
        DO dt = 0, nbTimeStepsInTraj-1
            nt = nbTimeStepsInTraj-dt
            acf(dt)=acf(dt)+&
                        (SUM(    ri(1:nt,1)*ri(1+dt:nt+dt,1)  &
                                +ri(1:nt,2)*ri(1+dt:nt+dt,2)  &
                                +ri(1:nt,3)*ri(1+dt:nt+dt,3))  )/DBLE(nt)
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

    DEALLOCATE(r,ri)

    
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


    CONTAINS

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
    
    SUBROUTINE closetraj
        CLOSE(10)
    END SUBROUTINE closetraj

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

    
    
    SUBROUTINE readArguments(Nat,trajectoryFileName)
        INTEGER, INTENT(OUT) :: Nat
        CHARACTER(LEN=*),INTENT(OUT) :: trajectoryFileName
        CALL GET_COMMAND_ARGUMENT(1,arg,STATUS=i)
        IF (i<0) THEN
            STOP "STOP. The length of the argument is too big for me :( "
        ELSE IF (i>0) THEN
            STOP "Argument retrieval failed. You should execute the program with the number of atoms as argument, e.g. in ./acf 10 "
        END IF
        READ(arg,*) Nat
        CALL GET_COMMAND_ARGUMENT(2,arg,STATUS=i)
        IF (i<0) THEN
            STOP "STOP. The length of the argument is too big for me :( "
        ELSE IF (i>0) THEN
            STOP "Argument retrieval failed. You should execute the program with the number of atoms as argument, e.g. in ./acf 10 "
        END IF
        trajectoryFileName = TRIM(ADJUSTL(arg))
    END SUBROUTINE readArguments
  
END PROGRAM
