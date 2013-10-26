! This program computes the autocorrelation of a list of sites during a trajectory in time.
! It is written by Maximilien Levesque while in postdoc in the group of Mathieu Salanne
! at UPMC Univ Paris 06, CNRS, ESPCI, UMR 7195, PECSA, F-75005 Paris, France.

program autoCorrelation

    implicit none
    character(len("acf.out")) :: outputFile = "acf.out"
    integer :: Nat
    integer :: i, nbTimeStepsInTraj, iostat, dt, t, nt, d
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

    print*,'You have' ,nbTimeStepsInTraj,' time steps in your trajectory file ',trim(adjustl(trajectoryFileName))
    print*,'Please be patient. Everything seems fine... Multiorigin effect! ;)'

    ! read vector of all sites i at all timesteps t
    allocate( r(Nat,nbTimeStepsInTraj,x:z) )
    call opentraj
    do t = 1, nbTimeStepsInTraj
        if( mod(t,1000)==1 ) print*,"READING timestep ",t," over ",nbTimeStepsInTraj
        do i = 1, Nat
            read(10,*) r(i,t,x), r(i,t,y), r(i,t,z)
        end do
    end do
    call closetraj

!~     ! compute autocorrelation function acf(dt)= <v_i(t).v_i(t+dt)>_{i,t}   where . is the scalar product
    allocate( ri(nbTimeStepsInTraj,x:z) )
    ri = 0.d0
    allocate( acf(0:nbTimeStepsInTraj-1) )
    acf = 0.d0
    call cpu_time(time0)

    do i= 1, Nat
        ri = r(i,:,:)
        do dt = 0, nbTimeStepsInTraj-1
            nt = nbTimeStepsInTraj-dt
            acf(dt)=acf(dt)+(sum(ri(1:nt,1)*ri(1+dt:nt+dt,1)+ri(1:nt,2)*ri(1+dt:nt+dt,2)+ri(1:nt,3)*ri(1+dt:nt+dt,3)))/dble(nt)
        end do
        call cpu_time(time1)
        if(mod(i,100)==1) print*,'Estimated remaining time = ',nint(dble(Nat-i)*(time1-time0)/dble(i)/60.d0),' min'
    end do
    acf = acf/dble(Nat)

    deallocate(r,ri)

    ! acf(t) will be written in file unit 11
    open(11,file=outputfile)
    do dt = 0, nbTimeStepsInTraj-1
        write(11,*) dt, acf(dt)
    end do
    close(11)

    print*,"-- Everything seems OK -- ;)"

    contains

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

    
    
    subroutine readArguments(Nat,trajectoryFileName)
        integer, intent(out) :: Nat
        character (len=*), intent(out) :: trajectoryFileName
        call get_command_argument(1,arg,status=i)
        if( i < 0 ) then
            stop "STOP. The length of the argument is too big for me :( "
        else if ( i > 0 ) then
            stop "Argument retrieval failed. You should execute the program with the number of atoms as argument, e.g. in ./acf 10 "
        end if
        read(arg,*) Nat
    
        call get_command_argument(2,arg,status=i)
        if( i < 0 ) then
            stop "STOP. The length of the argument is too big for me :( "
        else if ( i > 0 ) then
            stop "Argument retrieval failed. You should execute the program with the number of atoms as argument, e.g. in ./acf 10 "
        end if
        trajectoryFileName = trim(adjustl(arg))
    
    end subroutine
  
end program
