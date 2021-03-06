SUBROUTINE RHSSFM(NEQ, T, Y, YP)
!    use omp_lib
    IMPLICIT NONE
    integer, intent(in) :: NEQ!, IPAR(:)
    double precision, intent(in) :: T, Y(NEQ)!, RPAR(:)
    double precision, intent(out) :: YP(NEQ)
    INTEGER :: J,I1,I2,IG, I1_3, I2_3
    REAL*8 AG(3,3), &
            DETA,DETAG,QZQ,U12, &
            G12(3,3),A(3,3), &
            Zq(3), Z(3,3),Q12(3)
    real*8 :: UXY0(3,3), UX0(3)

    if (y(3*Natom+1) == 0d0) then
        call rhss_zero_time(NEQ, y, yp)
        return
    end if

    call unpack_y(y, Q, G, gama)

!    do I1=1,Natom-1
	call unpack_y(y, Q, G, gama) 

    U = 0; UX = 0; UXY = 0;
    do I1=1,Natom-1
        I1_3 = 3*(I1-1) + 1
        DO I2=I1+1,Natom
            I2_3 = 3*(I2-1) + 1
            Q12 = Q(I1_3:I1_3+2) - Q(I2_3:I2_3+2)
            Q12 = min_image(Q12, bl)
            G12=G(I1_3 : I1_3 + 2, I1_3 : I1_3 + 2) &
                + G(I2_3 : I2_3 + 2, I2_3 : I2_3 + 2) &
                - G(I2_3 : I2_3 + 2, I1_3 : I1_3 + 2) &
                - G(I1_3 : I1_3 + 2, I2_3 : I2_3 + 2)

            call detminvm(G12, DETA, A)
            DETA = 1.0d0/DETA

            UX0 = 0d0; UXY0 = 0d0
            DO IG=1,NGAUSS      ! BEGIN SUMMATION OVER GAUSSIANS
                AG = A
                do J=1,3
                    AG(J,J)=LJA(IG)+AG(J,J)
                end do

                call detminvm(AG, DETAG, Z)
                Z = - LJA(IG)**2 * Z

                do J=1,3
                    Z(J,J) = Z(J,J) + LJA(IG)
                end do

                Zq = matmul(Z, Q12) ! R = -2.0*Zq
                qZq = dot_product(Q12, Zq) 

                U12 = SQRT(DETA/DETAG)*EXP(-qZq)*LJC(IG)
                U = U + U12

                UX0 = UX0 - 2d0*U12*Zq
                do J=1,3
                    UXY0(:,J) = UXY0(:,J) + 2d0*U12*(2d0*Zq*Zq(J) - Z(:,J))
                end do
            end do ! IG
            
! Avoid load and store as much as possbile. Store now, process as a stream later. Much faster.
            UX(I1_3 : I1_3 + 2) = UX(I1_3 : I1_3 + 2) + UX0
            UX(I2_3 : I2_3 + 2) = UX(I2_3 : I2_3 + 2) - UX0
            UXY(I1_3 : I1_3 + 2, I1_3 : I1_3 + 2) = UXY(I1_3 : I1_3 + 2, I1_3 : I1_3 + 2) + UXY0
            UXY(I2_3 : I2_3 + 2, I2_3 : I2_3 + 2) = UXY(I2_3 : I2_3 + 2, I2_3 : I2_3 + 2) + UXY0
            
            !if (sum(q12**2) <= rfullmatsq) then
                UXY(I1_3 : I1_3 + 2, I2_3 : I2_3 + 2) = -UXY0
                UXY(I2_3 : I2_3 + 2, I1_3 : I1_3 + 2) = -UXY0
            !end if
        end do ! I2
    end do ! I1

    
    !QP = -matmul(G, UX)
    call dsymm('L', 'L', 3*Natom, 1, -1d0, G, 3*Natom, UX, 3*Natom, &
        0d0, QP, 3*Natom)

    !GU = matmul(G, UXY)
    call dsymm('L', 'L', 3*Natom, 3*Natom, 1d0, G, 3*Natom, UXY, 3*Natom, &
        0d0, GU, 3*Natom)

    !GP = -matmul(GU, G)
    call dsymm('R', 'L', 3*Natom, 3*Natom, -1d0, G, 3*Natom, GU, 3*Natom, &
        0d0, GP, 3*Natom)

    do J=1,3*Natom
        GP(J,J) = GP(J,J) + invmass
    end do

!    do I1=1,Natom-1
!        I1_3 = 3*(I1-1) + 1
!        DO I2=I1+1,Natom
!            I2_3 = 3*(I2-1) + 1
!            Q12 = Q(I1_3:I1_3+2) - Q(I2_3:I2_3+2)
!            Q12 = min_image(Q12, bl)
!            if (sum(q12**2) > rfullmatsq) then
!                GP(I1_3 : I1_3 + 2, I2_3 : I2_3 + 2) = 0
!                GP(I2_3 : I2_3 + 2, I1_3 : I1_3 + 2) = 0
!           end if
!        end do
!    end do
    
    gamap = -(0.25d0 * sum(UXY*G) + U)/real(Natom)

	call pack_y(QP, GP, gamap, yp)
END SUBROUTINE RHSSFM

subroutine rhss_zero_time(NEQ, y, yp)
    integer, intent(in) :: NEQ
    double precision, intent(in) :: y(:)
    double precision, intent(out) :: yp(:)

    double precision :: qij(3), qi(3), qj(3), rsq
    integer :: i, j

    yp = 0d0

    j = 3*Natom + 1
    do i=1,3*Natom
        yp(j) = invmass
        j = j + 3*Natom - i + 1
    end do

    if (NEQ > 3*Natom + nlg + 1) then
        j = 3*Natom + nlg
        do i=1,Natom
            yp(j + 9*(i - 1) + 1) = 1d0
            yp(j + 9*(i - 1) + 5) = 1d0
            yp(j + 9*i          ) = 1d0
        end do
    end if

    U=0d0

    DO I=1,Natom-1
        qi = y(3*I-2:3*I)
        DO J=I+1,Natom
                qj = y(3*J-2 : 3*J)
                qij = qi - qj
                rsq = sum(min_image(qij, BL)**2)
                U = U + sum(LJC(1:NGAUSS)*EXP(-LJA(1:NGAUSS)*rsq))
        ENDDO
    ENDDO

    yp(NEQ) = -U/real(Natom)
end subroutine
