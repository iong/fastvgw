program clustergs
    use vgw
    use vgwspfm
    use vgwfm
    use utils, only: M_PI
    implicit none
!    include 'mkl_service.fi'

    character(256) :: arg, coords, waste

    integer :: i, Natom, status

    double precision :: endtime(3), begtime(3), mass, taustop
    double precision ::  deBoer=0.1, kT=0.01d0, RC=100d0, BL=100d0, sigma0 = 2.749, epsilon0 = 35.6d0

    double precision, allocatable :: r0(:,:)
    double precision :: U(4), Uinf(4), rt(3)

    call get_command_argument(1, coords)

    open(33, file=trim(coords))
    read (33, *) Natom


    allocate(r0(3,Natom))
    U = 0
    read(33,IOSTAT=status,FMT=*) deBoer, U(1)
    if (status > 0) then
        U(1) = 0d0
    else if (status<0) then
        write (*,*) 'read(33,*) failed with status', status
        stop
    end if

    read(33, *) (waste, r0(:,i), i=1,Natom)
    close(33)

    ! deBoer can be supplied as a second argument on the command line
    if (command_argument_count() > 1) then
        call get_command_argument(2, arg)
        read(arg, *) deBoer
    end if


    r0 = r0*1.05d0

    !r0 = r0 / sigma0

    taustop = 1d0/kT
    mass = 1/deBoer**2
    !call vgwinit(natom, 'LJ', RCUTOFF=100d0, massx=mass)
    !call vgw0(r0, BL, taustop, U(2, i), rt(1, i))
    !call vgwcleanup()

    !call vgwfminit(natom, 'LJ', massx=mass, RCUTOFF=100d0)
    !call vgw0fm(r0, BL, taustop, U(3), rt(2))
    !call vgwfmcleanup()

    call vgwspfminit(natom, 'LJ', RCUTOFF=100d0, massx=mass)
    call vgw0spfm(r0, BL, taustop, U(4))
    call vgwspfmcleanup()

   ! U(2:4) = U(2:4) * epsilon0

	print *, Natom, U(4)/real(Natom)

end program clustergs
