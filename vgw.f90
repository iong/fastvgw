module vgw
    use utils
    implicit none
    private
    public :: vgwinit, vgw0, vgw0v, get_fx, init_fmcorr
    
    integer :: Natom, Nmax, maxthreads, Ncorr=200
    real*8 :: BL, rcorrmin=2.7, rcorrmax=8.0, drcorr, Ucorrtot
    real*8, dimension(10) :: LJA, LJC
    integer :: NGAUSS
    integer, allocatable :: NBIDX(:,:), NNB(:)
    
    real*8 :: T, gama, gamap
    real*8, allocatable :: Q(:,:), G(:,:,:), Qnk(:,:,:), gamak(:,:), &
                            QP(:,:), GP(:,:,:), fmcorr(:), Ucorr(:)
                            
    real*8, allocatable :: UPV(:,:,:), UPM(:,:,:,:), TRUXXG(:), U(:)
    
    real*8 :: invmass, RC, TAUMIN, mass
    logical :: finished
    integer :: tid=0, nthr=1, thread_start, thread_stop, nnbmax
!$OMP THREADPRIVATE(tid, nthr, thread_start, thread_stop, nnbmax)

contains

subroutine vgwinit(Nmax_, species, M, rcutoff)
    use omp_lib
    implicit none
    integer, intent(in) :: Nmax_
    character(*), intent(in) :: species
    real*8, intent(in), optional :: M, rcutoff

    Nmax = Nmax_
    maxthreads = omp_get_max_threads()
    allocate(NNB(Nmax), NBIDX(Nmax,Nmax), upv(3,Nmax,0:maxthreads-1), &
        upm(3,3,Nmax,0:maxthreads-1), TRUXXG(0:maxthreads-1), &
        U(0:maxthreads-1),Ucorr(0:maxthreads-1), &
        Q(3,Nmax), G(3,3,Nmax), Qnk(3,3,Nmax), gamak(3,Nmax), QP(3,Nmax), &
        GP(3,3,Nmax), fmcorr(Ncorr))
    
    
    if (species=='pH2-3g') then
        NGAUSS=3
        LJA(1:3) = (/ 0.669311, 0.199426, 0.092713/)
        LJC(1:3) = (/ 29380.898517, -303.054026, -40.574585 /)
        mass = 2.0
        rc = 8.0
        TAUMIN=1d-4
    else if (species=='pH2-4g') then
            NGAUSS=4
        LJA(1:4) = (/ 1.038252215127D0, 0.5974039109464D0, 0.196476572277834D0, &
                    0.06668611771781D0 /)
        LJC(1:4) = (/ 96609.488289873d0, 14584.62075507514d0, -365.460614956589d0, &
                    -19.5534697800036d0 /)
        mass = 2.0
        rc = 8.0
        TAUMIN=1d-4
    end if
    
    if (present(M)) then
        mass = M
    end if
    if (present(rcutoff)) then
        rc = rcutoff
    end if
        
    mass = mass*0.020614788876D0
    invmass = 1.0/mass
    Ucorrtot = 0.d0
end subroutine

subroutine init_gaussians(q0, tau, mm)
    REAL*8, intent(in) :: Q0(:,:), tau
    logical, intent(in) :: mm
    real*8, parameter :: E3(3,3) = reshape( (/1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0/), (/3, 3/) )
    
    if (mm) then
        call Upot_UX_tau0(Q0)
    else
        call Upot_tau0(Q0)
    end if

!$OMP SINGLE
    gama = -tau*sum(U(0:nthr-1))
    Q(:,1:Natom) = Q0

    G(:,:,1:Natom) = spread(tau*invmass*E3, 3, Natom)
    if (mm) then
        Qnk(:,:,1:Natom) = spread(E3, 3, Natom)
        gamak(:,1:Natom) = - tau * UPV(:,1:Natom,0)
    end if
    
    Ucorrtot = sum(Ucorr(0:nthr-1))
!$OMP END SINGLE
end subroutine

function get_fx() result(fx)
    real*8 :: fx(3, Natom)
    fx = - gamak(3,Natom)/T
end function

subroutine init_fmcorr(taumax)
    use vgwm
    implicit none
    double precision, intent(in) :: taumax
    double precision :: rl(3,2), Ul(2)
    integer :: i
    
    call INITIALIZE_VGWM(2,1d-4,1d-5,0d0,0d0)
    
    rl = 0.0
    drcorr = (rcorrmax - rcorrmin) / (Ncorr - 1)
    
    do i = 1,Ncorr
        rl(1,2) = rcorrmin + real(i-1)*(rcorrmax - rcorrmin)/real(Ncorr-1)
        call vgw0(rl, 100d0, taumax, 0d0, Ul(1))
        call VGWMQUENCH(rl, taumax, Ul(2))
        fmcorr(i) = Ul(2) - Ul(1)
    end do
    
    fmcorr = fmcorr - fmcorr(Ncorr)
    
    call CLEANUP_VGWM()
end subroutine init_fmcorr

include 'interaction_lists.f90'
include 'potential_energy.f90'
include 'vgw0.f90'
!include 'vgw1.f90'
include 'rhss0.f90'
!include 'rhss1.f90'

end module vgw

