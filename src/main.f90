! This program computes the autocorrelation of a list of sites during a trajectory in time.
! It is written by Maximilien Levesque while in postdoc in the group of Mathieu Salanne
! at UPMC Univ Paris 06, CNRS, ESPCI, UMR 7195, PECSA, F-75005 Paris, France.

program autoCorrelation

    IMPLICIT NONE
    
    character(len("acf.out")) :: outputFile = "acf.out"
    integer :: Nat,i,nbTimeStepsInTraj,iostat,dt,t,nt,d
    integer, parameter :: x=1, y=2, z=3
    double precision, dimension(:,:,:), allocatable :: r ! position of site i at timestep t
    double precision, dimension(:,:), allocatable :: ri
    double precision :: lx, ly, lz, diffx, diffy, diffz, rc, dx2, dy2, dz2, r0, r1, time1, time0
    double precision, dimension(:), allocatable :: acf
    character(len=300) :: arg, trajectoryFileName
    double precision, dimension(x:z) :: l
    logical :: doagain
    double precision :: acf_dt_i, acf_dt_i_t

    ! read all arguments that MUST be given at execution
    call readArguments(Nat,trajectoryFileName)

    ! deduce the number of timesteps in the trajectory from the number of lines in the trajectory file
    nbTimeStepsInTraj = NbOfLinesInTraj(trajectoryFileName)/Nat

    print*,'I found',nbTimeStepsInTraj,' timesteps in your trajectory called ',trim(adjustl(trajectoryFileName))
    print*,'Please be patient... Everything looks fine...'

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
    do i= 1, Nat
        ri = r(i,:,:)
        do dt = 0, nbTimeStepsInTraj-1
            nt = nbTimeStepsInTraj-dt
            acf(dt)=acf(dt)+(sum(ri(1:nt,1)*ri(1+dt:nt+dt,1)+ri(1:nt,2)*ri(1+dt:nt+dt,2)+ri(1:nt,3)*ri(1+dt:nt+dt,3)))/dble(nt)
        end do
        call cpu_time(time1)
        if(mod(i,100)==1) print*,'Remaining time ≈ ',nint(dble(Nat-i)*(time1-time0)/dble(i)/60.d0),' min'
    end do
    acf = acf/dble(Nat)

    deallocate(r,ri)

    
    ! acf(t) will be written in file unit 11
    open(11,file=outputfile)
    do dt = 0, nbTimeStepsInTraj-1
        write(11,*) dt, acf(dt)
    end do
    close(11)

    PRINT*,"-- Finished. GGHF ;) --"


    CONTAINS

    subroutine opentraj
        call inquireFileExistence(trajectoryFileName)
        ! read positions
        open(10, file=trajectoryFileName,status='old',iostat=iostat)
        if (iostat /= 0) then
            write (*,*) 'File open failed for',trajectoryFileName
            write (*,*) 'Error code is ', iostat
            stop
        end if
    end subroutine
    
    subroutine closetraj
        close(10)
    end subroutine

    subroutine inquireFileExistence(fileName)
        character(len=*), intent(in) :: fileName
        integer, parameter :: stderr = 0
        logical :: exist
        inquire(file=fileName, exist=exist)
        if( .not. exist) then
            write(stderr,*) "YOUR ERROR (not mine ;): The file ",fileName," does not exist. It should."
            stop
        end if
    end subroutine

    function NbOfLinesInTraj(filename)
        character(len=*), intent(in) :: filename
        integer :: NbOfLinesInTraj
        call opentraj
        ! computes the number of lines in traj.in and deduces the number of timesteps
        NbOfLinesInTraj = -1
        do while (iostat == 0)
            read(10,*,iostat=iostat)
            NbOfLinesInTraj = NbOfLinesInTraj + 1
        end do
        call closetraj
    end function
    
    
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
