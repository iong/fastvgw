program gibbs3h
    implicit none
    real*8, parameter :: M_PI = 3.141592653d0
    real*8 :: rhcore
    real*8 :: Vtot, V(2), bl(2), kT, beta, U0(4), xstep(2), &
            logVstep, Vstep
    real*8 :: Z, Um(2), Vm(2), Nm(2), pm(2), rho0(2), rhom(2)
    real*8, allocatable :: rs(:,:, :)
    real*8, allocatable, target :: Umat(:,:,:), Umatnew(:,:,:)
    real*8, allocatable :: Ucacheline(:,:), y6cache(:)
    integer :: Ntot, N(2), Nswap=0, Nvol
    integer :: npickatt(2), npickattmax, nswapacc = 0
    real*8 :: swapdest(3,2)
    integer :: nxtrials(2)=0, nxacc(2)=0, nvoltrials=0, nvolacc=0
    integer :: imc, NMC=10000000, jmc, iliq,igas, icycle, Ncalibrate=500000
    character(LEN=256) :: arg


    real*8 :: rn

    call get_command_argument(1, arg)
    read (arg, *) N(1)
    call get_command_argument(2, arg)
    read (arg, *) N(2)
    call get_command_argument(3, arg)
    read (arg, *) rho0(1)
    call get_command_argument(4, arg)
    read (arg, *) rho0(2)
    call get_command_argument(5, arg)
    read (arg, *) kT
    
    write (*,*) 'kT =', kT

    if (command_argument_count() >= 6) then
        call get_command_argument(6, arg)
        read(arg, *) Nswap
    end if


    Ntot = sum(N)
    allocate(Umat(Ntot, Ntot, 4), Umatnew(Ntot, Ntot, 4), rs(3, Ntot, 2), &
        Ucacheline(Ntot,2), y6cache(Ntot))

    V = N/rho0
    Vtot = sum(V)
    Vstep=0.01*minval(V)

    bl = V**(1d0/3d0)
    xstep = minval(0.3/bl)

    if (minval(bl) < 5.0) then
        write (*,'("Box #",I1," is too small. L = ",F8.4, &
            " < 5*sigma = ",F8.4)') minloc(bl), minval(bl), 5.0
        write (*,*) 'Consider increasing the number of particles or decreasing the density'
    end if

    Nvol = 1!minval(N)/20
    write (*,*) 'Nvol =', Nvol

    beta = 1.0/kT

    rhcore = 1.0/(0.5 + 0.5*sqrt(1 - kT*log10(1d-15)))**(1d0/6d0)

    call populate_cube2(rs(:,1:N(1), 1))
    call populate_cube2(rs(:,1:N(2), 2))
    call init_interaction_matrix(bl, U0)

    do imc=1,2000000
        call random_number(rn)
        jmc=int(rn*Ntot) + 1
        call mc_move(jmc)
        call cumulate()
        if (mod(imc,100000) == 0) then
            write (*,'("NMC =",I10,", N =",2I5,", V =",2F9.2,", rho =",2F9.6,", U =",2F15.6,", p =",2F8.4,", Nswap = ",I0)') &
                imc, N, V, N/V, U0(1:2), pm/Z, Nswap
        end if
    end do

    do imc=1,100000
        call random_number(rn)
        jmc=1+int(rn*(Ntot + Nvol))
        if (jmc<= Ntot) then
            call mc_move(jmc)
        else
            call mc_vol()
        endif
        if (mod(imc,100000) == 0) then
            write (*,*) "xstep =", xstep, ", Vstep=", Vstep
        end if
    end do

    Z = 0; Um = 0; Nm = 0; pm = 0; Vm = 0; rhom = 0
    do imc=1,NMC
        call random_number(rn)
        jmc=1+int(rn*(Ntot + Nswap + Nvol))
        if (jmc<= Ntot) then
            call mc_move(jmc)
        elseif (jmc <= Ntot + Nswap) then
            call mc_swap()
        elseif (jmc <= Ntot + Nswap + Nvol) then
            call mc_vol()
        endif
        call cumulate()
        if (mod(imc,100000) == 0) then
            write (*,'("NMC =",I10,", N =",2I5,", V =",2F9.2,", rho =",2F9.6,", U =",2F15.6,", p =",2F8.4,", Nswap = ",I0)') &
                imc, N, V, N/V, U0(1:2), pm/Z, Nswap
        end if
        if (mod(imc,1000000) == 0) then
            call dump_data()
        end if
    enddo
contains

subroutine populate_cube2(rs)
    real*8, intent(out) :: rs(:,:)
    logical, allocatable :: occupied(:)
    integer :: np, nu, i, j, lattpos
    np = size(rs, 2)
    nu = ceiling(real(np)**(1.0/3.0))

    allocate(occupied(nu*nu*nu))
    occupied = .FALSE.

    do i=1,np
        do
            call random_number(rn)
            lattpos = int(rn*nu**3)
            if (.not.occupied(lattpos+1)) exit
        end do
        occupied(lattpos+1) = .TRUE.

        do j=1,3
            rs(j, i) = real(mod(lattpos, nu))
            lattpos = lattpos / nu
        end do
    end do

    rs = (rs+0.5)/real(nu)

    deallocate(occupied)
end subroutine

pure function min_image(r)
    implicit none
    real*8, intent(in) :: r(3)
    real*8 :: min_image(3)
    integer :: i

    do i=1,3
        if (r(i) > 0.5) then
            min_image(i) = r(i) - 1.0
        elseif (r(i) < -0.5) then
            min_image(i) = r(i) + 1.0
        else
            min_image(i) = r(i)
        end if
    end do
end function

pure function too_close(rsj, rs, blj)
    real*8, intent(in) :: rsj(3), rs(:,:), blj
    logical :: too_close
    real*8 :: drsq, rshc_sq, drs(3)
    integer :: i

    too_close = .FALSE.
    rshc_sq = (rhcore / blj)**2
    do i=1,size(rs, 2)
        drs = rsj - rs(:,i)
        drsq = sum(min_image(drs)**2)
        if (drsq < rshc_sq) then
            too_close = .TRUE.
            return
        end if
    enddo
end function


subroutine mc_move(jin)
    implicit none
    integer, intent(in) :: jin
    integer, parameter :: nstepcheck=200
    real*8 :: rsj(3), dU, dr(3), p, dvir
    integer :: ibox, j, i

    ibox = 1
    j = jin
    if (j > N(1)) then
        ibox = 2
        j = j - N(1)
    end if

    nxtrials(ibox) = nxtrials(ibox) + 1

    call random_number(dr)
    rsj = rs(:, j, ibox) + xstep(ibox) * (dr-0.5)
    rsj = rsj - floor(rsj)

    call dU_particle_move(j, ibox, rsj, dU)
    p = exp(-beta * dU)
    call random_number(rn)
    if (p>rn) then
        nxacc(ibox) = nxacc(ibox)  + 1
        rs(:, j, ibox) = rsj
        call dU_accept_move(j, ibox, dvir)
        U0(ibox) = U0(ibox) + dU
        U0(ibox+2) = U0(ibox+2) + dvir
    end if

    if (mod(nxtrials(ibox), nstepcheck) == 0) then
        if (nxacc(ibox) > (0.4 * nstepcheck)) then
            xstep(ibox) = xstep(ibox)*1.08351d0
        elseif (nxacc(ibox) < (0.2 * nstepcheck)) then
            xstep(ibox) = xstep(ibox)/1.04792d0
        end if
        xstep(ibox) = min(xstep(ibox), 0.25)
        nxacc(ibox) = 0
    end if
end subroutine

subroutine mc_swap()
    implicit none
    integer :: isrc, idest, jkill
    real*8 :: pacc, pacc0, rsn(3), rso(3), U0new(4)

    call random_number(rn)
    isrc = 1 + int(2.0*rn)
    idest = 3 - isrc

    call random_number(rsn)
    if (too_close(rsn, rs(:,1:N(idest), idest), bl(idest)) ) return

    call random_number(rn)
    jkill = 1 + int(rn*N(isrc))

    rso = rs(:,jkill, isrc)
    rs(:,jkill,isrc) = rs(:,N(isrc),isrc)

    rs(:,N(idest) + 1,idest) = rsn

    pacc0 = ( V(idest) * real(N(isrc)) ) / ( V(isrc) * real(N(idest) + 1) )

    N(isrc) = N(isrc) - 1
    N(idest) = N(idest) + 1

    call init_interaction_matrix(bl, U0new)

    pacc = pacc0 * exp( - beta * sum(U0new(1:2) - U0(1:2)) )

    call random_number(rn)
    if (pacc > rn) then
        nswapacc = nswapacc + 1
    else
        N(isrc) = N(isrc) + 1
        N(idest) = N(idest) - 1
        rs(:,N(isrc),isrc) = rso
        call init_interaction_matrix(bl, U0)
    end if
end subroutine


subroutine mc_vol()
    implicit none
    integer, parameter :: nstepcheck=50
    real*8 :: U0new(4), Vn(2), bln(2), logp, Uctrl(4)

    nvoltrials = nvoltrials + 1

    call random_number(rn)
    Vn(1) = V(1) + Vstep*(rn-0.5)
    Vn(2) = Vtot - Vn(1)
    bln = Vn**(1.0/3.0)

    ! ensure bl > 5\sigma
    if (minval(bln) < 5.0) return

    Uctrl = total_energy(bln)
    call init_interaction_matrix(bln, U0new)
    if (maxval(abs(U0new - Uctrl)) > 1d-10) then
        write (*,*) 'Large error!'
        write (*,'(4ES16.8)') U0new - Uctrl
    end if
    logp = sum ( -beta * (U0new(1:2) - U0(1:2)) + real(N)*(log(Vn/V)) )
    call random_number(rn)
    if (logp>log(rn)) then
        nvolacc = nvolacc + 1
        V = Vn
        bl = bln
        U0 = U0new
    else
        call init_interaction_matrix(bl, U0)
    end if

    if (mod(nvoltrials, nstepcheck) == 0) then
        if (nvolacc > 0.55 * nstepcheck) then
            Vstep = Vstep*1.18351d0
        elseif (nvolacc < 0.45 * nstepcheck) then
            Vstep = Vstep/1.14792d0
        end if
        !Vstep = max(Vstep, 0.001*Vtot)
        nvolacc = 0
    end if
end subroutine

function total_energy(bln) result(utot)
    real*8, intent(in) :: bln(2)
    real*8 :: utot(4)
    real*8:: lrcorr(2), blsq, d, dsq
    integer :: i, j, ibox, k

    utot = 0.0
    do ibox=1,2
     do i=1,N(ibox)-1
         blsq = bln(ibox)**2
         do j=i+1, N(ibox)
             dsq = 0.0
             do k=1,3
                 d = abs(rs(k, i, ibox) - rs(k, j, ibox))
                 if (d > 0.5) d = 1.0 - d
                 dsq = dsq + d*d
             end do
             if (dsq > 0.25) cycle

             dsq = dsq * blsq

             utot(ibox) = utot(ibox) + 1.0/dsq**6 - 1.0/dsq**3
             utot(ibox+2) = utot(ibox+2) + 12.0/dsq**6 - 6.0/dsq**3
         end do
     end do
    end do
    utot = utot * 4.0
    lrcorr = real(N**2)/bln**3*8.0/9.0*M_PI*(2.0**9/bln**9 - 24.0/bln**3)
    utot(1:2) = utot(1:2) + lrcorr
end function

subroutine ULJp(rsj, rs, blj, ULJ, VLJ)
    implicit none
    real*8, intent(in) :: rsj(3), rs(:,:), blj
    real*8, intent(out) :: ULJ(:)
    real*8, intent(out), optional :: VLJ(:)
    real*8 :: drs(3, size(rs, 2)), drs0(3)
    integer :: Nl, i

    Nl = size(rs, 2)

    do i=1,Nl
        drs0 = rsj - rs(:,i)
        drs(:,i) = min_image(drs0)
    enddo
    y6cache(1:Nl) = 1.0 / sum(drs**2, 1)**3
    forall (i=1:Nl, y6cache(i) < 64.0)
        y6cache(i) = 0.0
    end forall
    y6cache(1:Nl) = 1.0/(blj**6) * y6cache(1:Nl)
    ULJ = 4d0*(y6cache(1:Nl)**2 - y6cache(1:Nl))
    if (present(VLJ)) then
        VLJ(1:Nl) = 12.0*ULJ(1:Nl) + 24.0*y6cache(1:Nl)
    end if
end subroutine

subroutine update_virial_cache(Nbox)
    integer, intent(in) :: Nbox
    Ucacheline(1:Nbox, 2) = 12.0*Ucacheline(1:Nbox, 1) &
                                + 24.0*y6cache(1:Nbox)
end subroutine

subroutine init_interaction_matrix(bln, Utot)
    implicit none
    real*8, intent(in) :: bln(2)
    real*8, intent(out) :: Utot(4)
    integer :: i, ibox, j
    real*8 :: lrcorr(2)

    Umat = 0d0
    do ibox=1,2
        do i=1, N(ibox)-1
            call ULJp(rs(:, i, ibox), rs(:,i+1:N(ibox), ibox), bln(ibox), &
                    Umat(i+1:N(ibox), i, ibox), &
                    Umat(i+1:N(ibox), i, ibox + 2))
        end do
    end do

    Utot = 0d0
    do i=1,4
        ibox = 1 + mod(i-1,2)
        do j=1,N(ibox)-1
            Utot(i) = Utot(i) + sum(Umat(j+1:N(ibox), j, i))
        end do
    end do

    do i=1,4
        ibox = 1 + mod(i-1,2)
        do j=2,N(ibox)
            Umat(1:j-1,j, i) = Umat(j, 1:j-1, i)
        end do
    end do

    lrcorr = real(N**2)/bln**3*8.0/9.0*M_PI*(2.0**9/bln**9 - 24.0/bln**3)
    Utot(1:2) = Utot(1:2) + lrcorr
end subroutine

subroutine dU_particle_move(j, ibox, rsj, du)
    implicit none
    integer, intent(in) :: j, ibox
    real*8, intent(in) :: rsj(3)
    real*8, intent(out) :: dU

    Ucacheline = 0d0
    call ULJp(rsj, rs(:,1:N(ibox), ibox), bl(ibox), &
            Ucacheline(1:N(ibox), 1))
    Ucacheline(j,1) = 0d0

    dU =  sum(Ucacheline( 1:N(ibox), 1)) - sum(Umat( 1:N(ibox) , j, ibox))
end subroutine 

subroutine dU_accept_move(p, pbox, dvir)
    implicit none
    integer, intent(in) :: p, pbox
    real*8, intent(out) :: dvir
    integer :: i, ibox

    call update_virial_cache(N(pbox))
    Ucacheline(p,2) = 0.0
    dvir =  sum( Ucacheline( 1:N(pbox), 2) ) &
                                    - sum(Umat( 1:N(pbox) , p, pbox + 2))

    do i=1,2
        ibox = pbox + (i-1)*2
        Umat(1:N(pbox),p, ibox) = Ucacheline(1:N(pbox), i)
        Umat(p,1:N(pbox), ibox) = Ucacheline(1:N(pbox), i)
    end do
end subroutine

subroutine cumulate()
    integer :: i
    Z = Z + 1.0
    do i=1,2
        if (N(i) ==0) cycle
        Um(i) = Um(i) + U0(i)/N(i)
        Nm(i) = Nm(i) + N(i)
        Vm(i) = Vm(i) + V(i)
        rhom(i) = rhom(i) + N(i)/V(i)
        pm(i) = pm(i) + (N(i)*kT + U0(i+2)/3.0)/V(i)
    end do
end subroutine

subroutine dump_data()
    character(256) :: fname

    write(fname, '("Um_",F4.2,".dat")') kT
    open(30,file=trim(fname))
    write(30,'(11F16.5)') kT, rhom/Z, pm/Z, Um/Z, Nm/Z, Vm/Z
    close(30)
end subroutine
end program gibbs3h
