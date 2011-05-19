program dimer
    use vgw
    !use vgwfm
    implicit none
    integer :: i, Npts=50, natom=2
    double precision :: r(3,2), Uc(3), rmin=2.8, rmax=10, kT=20.0!, RC=8.0, mass = 2.0
    double precision, allocatable :: Y(:), FX(:,:), InvMeff(:,:,:), SqrtMeff(:,:,:)
    double precision :: LNP, ENRG
    
    allocate(Y(1+21*Natom), FX(3,Natom), InvMeff(3,3,Natom), SqrtMeff(3,3,Natom))
    
    MASS = 2.0 * 0.020614788876D0
    NGAUSS = 4
    LJA(1) = 1.038252215127D0
    LJA(2) = 0.5974039109464D0
    LJA(3) = 0.196476572277834D0
    LJA(4) = 0.06668611771781D0
    LJC(1) = 96609.488289873d0
    LJC(2) = 14584.62075507514d0
    LJC(3) = -365.460614956589d0
    LJC(4) = -19.5534697800036d0
    !the rest
    rc=8
    rtol=1d-4
    atol=1d-4
    taumin=1d-4

    call vgwinit(natom, 100d0)
    !call vgwfminit(natom, 'pH2-4g')
    
    open(33,file='VGW_SP_FM.dat')
    
    r=0; Uc=0
    do i = 1,Npts
        r(1,2) = rmin + real(i-1)*(rmax - rmin)/real(Npts-1)
        !call vgw0(r, 100d0, 1.0/kT, 0d0, Uc(1))
        call vgw0(r, Uc(1), 1.0/kT, 0d0, Y)
        !call vgw0fm(r, 100d0, 1.0/kT, 0d0, Uc(2))
        call vgwquenchspb(natom,MASS/0.020614788876D0,r,FX,LNP,Uc(3),&
                      ENRG,1d0/kT,1d-4,100d0,1d-4,RC,Y, InvMeff,SqrtMeff)
        write (33, '(4F12.6)') r(1,2), Uc
    end do
    close(33)
end program dimer
