module fea

    !! This module contains procedures that are common to FEM-analysis

    implicit none
    save
    private
    public :: displ, initial, buildload, buildstiff, enforce, recover, eigen

contains

!
!--------------------------------------------------------------------------------------------------
!
    subroutine initial

        !! This subroutine is mainly used to allocate vectors and matrices

        use fedata
        use link1
        use plane42rect

        integer, dimension(mdim) :: edof
        integer :: e, i, maxdist, dist, nen

        maxdist = 0
        dist = 0

        ! Hint for continuum elements
        ! This subroutine computes the number of global equation,
        ! half bandwidth, etc and allocates global arrays.

        ! Calculate number of equations
        neqn = 3*nn

        if (.not. banded) then
            allocate (kmat(neqn, neqn))
        else
            do e = 1, ne
                do i = 1, element(e)%numnode
                    edof(3*i-2) = 3 * element(e)%ix(i) - 2
                    edof(3*i-1) = 3 * element(e)%ix(i) - 1
                    edof(3*i  ) = 3 * element(e)%ix(i)
                end do
                dist = maxval(edof) - minval(edof)

                if (dist > maxdist) then
                    maxdist = dist
                end if

            end do

            bw = maxdist + 1
            allocate (kmat (bw, neqn))

        end if
        allocate (p(neqn), d(neqn))
        allocate (strain(ne, 3), stress(ne, 3))
        allocate (principal_stresses(ne, 3), sigmavm(ne))

        ! Initial stress and strain
        strain = 0
        stress = 0
        principal_stresses = 0
        sigmavm = 0


    end subroutine initial
!
!--------------------------------------------------------------------------------------------------
!
    subroutine displ

        !! This subroutine calculates displacements

        use fedata
        use numeth
        use processor

        integer :: e, i,l
        real(wp), dimension(:), allocatable :: plotval
        real(wp) :: comp

        theta = pi


        ! Build load-vector
        call buildload

        ! Build stiffness matrix
        call buildstiff

        ! Remove rigid body modes
        call enforce

        if (.not. banded) then
            ! Factor stiffness matrix
            call factor(kmat)
            ! Solve for displacement vector
            call solve(kmat, p)
        else
            ! Factor banded stiffness matrix
            call bfactor(kmat)
            ! Solve for banded displacement vector
            call bsolve(kmat, p)
        end if

        ! Transfer results
        d(1:neqn) = p(1:neqn)
        !print*, 'd'
        !print*, d

        ! Recover stress
        call recover

        ! compute global compliance
        comp = dot_product(d,p)
        print*, 'compliance', comp

        ! Output results
        call output


        allocate (plotval(ne))
        do e = 1, ne

            ! original statement
            !if (element(e)%id == 1) then
            !    plotval(e) = stress(e,1)
            !else if (element(e)%id == 2) then
            !    plotval(e) = stress(e,1)
            !end if

            if (element(e)%id == 1) then
                plotval(e) = sigmavm(e)
            else if (element(e)%id == 2) then
                plotval(e) = sigmavm(e)
            end if
        end do

        ! Plot deformed structure
        !call plotmatlabdef('Deformed')

        call plotmatlabdefplate('Deformed')

        ! Plot element values
        !call plotmatlabeval('Stresses',plotval)

        ! Plot principal stresses
        !call plotmatlabevec('Principal Stresses',principal_stresses(:,1), principal_stresses(:,2), principal_stresses(:,3))

    end subroutine displ
!
!--------------------------------------------------------------------------------------------------
!
    subroutine eigen(neig)

        !! This subroutine calculates eigenmodes

        use fedata
        use numeth
        use processor

        integer :: i, j, kmax, k
        real(wp), dimension(:), allocatable :: plotval
        real(wp) :: comp, estop, rk, cj
        integer, intent(in) :: neig

        real(wp), dimension(neqn) :: xeigen, yeigen, xeigennew, yeigenprev, eigenvec, zeigenvec
        real(wp), dimension(neqn, neig) :: deigen, zeigen
        real(wp), dimension(neqn, 1) :: matrix1
        real(wp), dimension(neig) :: lambda, freq
        call stopwatch('star')

        ! Build load-vector
        call buildload

        ! Build stiffness matrix
        call buildstiff

        ! Remove rigid body modes
        call enforce

        deigen = 0

        kmax = 100
        estop = 1e-12

        ! factorize K
        if (.not. banded) then
            ! Factor stiffness matrix
            call factor(kmat)
        else
            ! Factor banded stiffness matrix
            call bfactor(kmat)
        end if


        do i = 1, neig

            ! initial guess
            xeigen = 1
            yeigen = 0
            zeigen = 0

            ! compute Y
            call mmul(xeigen,yeigen)

            print*, 'i'
            print*, i

            ! compute z
            do j = 1, i-1

                !print*, 'j'
                !print*, j

                call mmul(reshape(deigen(:,j),shape(xeigen)),zeigenvec)
                zeigen(:,j) = zeigenvec

                !print*, 'zeigen'
                !print*, zeigen

            end do


            do k = 1, kmax

                yeigenprev = yeigen

                ! compute xeigen
                if (.not. banded) then
                    call solve(kmat, yeigen)
                else
                    call bsolve(kmat, yeigen)
                end if
                xeigennew(1:neqn) = yeigen(1:neqn)

                ! compute cj and orthogonalize xeigen
                do j = 1, i-1
                    cj = dot_product(xeigennew, reshape(zeigen(:,j),shape(xeigen)))

                    xeigennew = xeigennew - cj*reshape(deigen(:,j),shape(xeigen))

                    !print*, 'cj'
                    !print*, cj
                    !print*, '(deigen(:,j)'
                    !print*, deigen(:,j)

                end do

                call mmul(xeigennew,yeigen)

                rk = sqrt(dot_product(xeigennew, yeigen))

                yeigen = yeigen/rk

                if (norm2(xeigennew - xeigen)/ norm2(xeigennew) < estop) exit

                xeigen = xeigennew

            end do

            deigen(:,i) = xeigennew/rk

            lambda(i) = dot_product(xeigennew, yeigenprev)/rk**2
            freq(i)=sqrt(lambda(i))

        end do

        print*, neqn
        print*, 'lambda'
        print*, lambda
        print*, 'freq'
        print*, freq

        ! neqn = 410
        !print*,'kmat'
        !do l = 1,size(kmat,1)
        !    print"(24(f6.2,tr1))",kmat(l,1:neqn)
        !end do



        call stopwatch('stop')

        call plotmatlabeig('Eigenmode', 10.0d0, 5.0d0*deigen(:,1),[1.0d0,0.010d0])

    end subroutine eigen
!
!--------------------------------------------------------------------------------------------------
!
    subroutine buildload

        !! This subroutine builds the global load vector

        use fedata
        use plane42rect

        integer :: i,j, eface, l, m
        ! Hint for continuum elements:
        real(wp), dimension(mdim) :: xe
        real(wp), dimension(mdim) :: re
        integer, dimension(mdim) :: edof
        real(wp) :: fe, thk

        integer :: dof, en, nen

        !print*, 'loads'
        !do i = 1,neqn
        !    print*, loads(i, 1:4)
        !end do

        ! Build load vector
        p(1:neqn) = 0

        do i = 1, np
            select case(int(loads(i, 1)))
            case( 1 )
            	! Build nodal load contribution
            	dof = loads(i, 2)*3-3+loads(i,3)
                p(dof) = loads(i,4)

            case( 2 )

            	! Build uniformly distributed surface (pressure) load contribution
            	en = loads(i,2)
                nen = element(en)%numnode
            	eface = loads(i,3)
            	fe = loads(i,4)
            	thk = mprop(element(en)%mat)%thk

            	! Find coordinates and degrees of freedom
                do j = 1, nen
                    xe(3*j-2) = x(element(en)%ix(j),1)
                    xe(3*j-1) = x(element(en)%ix(j),2)
                    xe(3*j  ) = x(element(en)%ix(j),3)

                    !edof(2*j-1) = 2 * element(en)%ix(j) - 1
                    !edof(2*j)   = 2 * element(en)%ix(j)

                    edof(3*j-2) = 3 * element(en)%ix(j) - 2
                    edof(3*j-1) = 3 * element(en)%ix(j) - 1
                    edof(3*j  ) = 3 * element(en)%ix(j)
                end do

                call shell41_re(xe,eface, fe, thk, re)

                do l = 1, size(re)
                    p(edof(l)) = p(edof(l)) + re(l)
                end do

            case default
                print *, 'ERROR in fea/buildload'
                print *, 'Load type not known'
                stop
            end select
        end do

        !print*, 'p'
        !do m = 1,size(p)
        !    print*, p(m)
        !end do

    end subroutine buildload
!
!--------------------------------------------------------------------------------------------------
!
    subroutine buildstiff

        !! This subroutine builds the global stiffness matrix from
        !! the local element stiffness matrices

        use fedata
        use link1
        use plane42rect

        integer :: e, i, j, k, m, l, c, r, krow, kcol
        integer :: nen

        ! Hint for system matrix in band form:
        ! integer :: irow, icol

        integer, dimension(mdim) :: edof
        real(wp), dimension(mdim) :: xe
        real(wp), dimension(mdim, mdim) :: ke
        real(wp), dimension(mdim, mdim) :: me

        ! Hint for modal analysis:
        ! real(wp), dimension(mdim, mdim) :: me
        real(wp) :: young, youngy, shear, area
        ! Hint for modal analysis and continuum elements:
        real(wp) :: nu, dens, thk

        print*,'theta',theta

        ! Reset stiffness matrix
        if (.not. banded) then
            kmat = 0
        else
            kmat = 0

        end if

        do e = 1, ne

            ! Find coordinates and degrees of freedom
            nen = element(e)%numnode

            do i = 1, nen
                 !xe(2*i-1) = x(element(e)%ix(i),1)
                 !xe(2*i  ) = x(element(e)%ix(i),2)

                 xe(3*i-2) = x(element(e)%ix(i),1)
                 xe(3*i-1) = x(element(e)%ix(i),2)
                 xe(3*i  ) = x(element(e)%ix(i),3)
                 edof(3*i-2) = 3 * element(e)%ix(i) - 2
                 edof(3*i-1) = 3 * element(e)%ix(i) - 1
                 edof(3*i  ) = 3 * element(e)%ix(i)
            end do

            ! Gather material properties and find element stiffness matrix
            select case( element(e)%id )
            case( 1 )
                 young = mprop(element(e)%mat)%young
                 area  = mprop(element(e)%mat)%area
                 call link1_ke(xe, young, area, ke)
            case( 2 )
                 young = mprop(element(e)%mat)%young
                 youngy = mprop(element(e)%mat)%youngy
                 shear = mprop(element(e)%mat)%shear

                 nu  = mprop(element(e)%mat)%nu
                 thk  = mprop(element(e)%mat)%thk

                 call shell41_ke(xe, young, youngy, shear, theta, nu, thk, dens, ke, me)
            end select

            ! Assemble into global matrix
            if (.not. banded) then
                do i = 1, 3*nen
                    do j = 1, 3*nen
                        kmat(edof(i), edof(j)) = kmat(edof(i), edof(j)) + ke(i, j)
                    end do
                end do

            else
                do i = 1,3*nen
                    do j = 1,3*nen
                        if (edof(j) <= edof(i)) then
                            krow = - (edof(j) - edof(i)) +1
                            kmat(krow,edof(j)) = kmat(krow,edof(j)) + ke(i,j)
                        end if
                    end do
                end do

            end if
        end do


        !print*,'kmat'
        !do l = 1,size(kmat,1)
        !    print"(24(f12.4,tr1))",kmat(l,1:neqn)
        !end do

    end subroutine buildstiff

!
!--------------------------------------------------------------------------------------------------
!

    subroutine mmul(xeigen,yeigen)

        !! This subroutine builds the global stiffness matrix from
        !! the local element stiffness matrices

        use fedata
        use link1
        use plane42rect

        real(wp), dimension(:), intent(in) :: xeigen
        real(wp), dimension(:), intent(out) :: yeigen

        integer, dimension(mdim) :: edof
        real(wp), dimension(mdim) :: xe
        real(wp), dimension(mdim, mdim) :: ke
        real(wp), dimension(mdim, mdim) :: me
        real(wp) :: young, youngy, shear, area
        real(wp) :: nu, dens, thk
        integer :: e, i, nen, j, k, h
        integer :: idof

        real(wp), dimension(mdim) :: xeigene, yeigene


        xeigene = 0
        yeigen = 0
        yeigene = 0

        do e = 1, ne

            ! Find coordinates and degrees of freedom
            nen = element(e)%numnode

            do i = 1, nen
                 xe(3*i-2) = x(element(e)%ix(i),1)
                 xe(3*i-1) = x(element(e)%ix(i),2)
                 xe(3*i  ) = x(element(e)%ix(i),3)
                 edof(3*i-2) = 3 * element(e)%ix(i) - 2
                 edof(3*i-1) = 3 * element(e)%ix(i) - 1
                 edof(3*i  ) = 3 * element(e)%ix(i)
            end do

            ! Gather material properties and find element stiffness matrix
            select case( element(e)%id )
            case( 1 )
                 young = mprop(element(e)%mat)%young
                 area  = mprop(element(e)%mat)%area
                 call link1_ke(xe, young, area, ke)
            case( 2 )
                 young = mprop(element(e)%mat)%young
                 youngy = mprop(element(e)%mat)%youngy
                 shear = mprop(element(e)%mat)%shear
                 nu  = mprop(element(e)%mat)%nu
                 thk  = mprop(element(e)%mat)%thk
                 dens  = mprop(element(e)%mat)%dens

                 call shell41_ke(xe, young, youngy, shear, theta, nu, thk, dens, ke, me)
            end select

            !yeigen (edof) = yeigen(edof) + matmul(me, xeigen(edof))

            do j = 1, mdim
                xeigene(j) =  xeigen(edof(j))
            end do


            yeigene = matmul(me, xeigene)

            do k = 1, mdim
                yeigen(edof(k)) = yeigen(edof(k)) + yeigene(k)
            end do


            do h = 1, nb
                idof = int(2*(bound(h,1)-1) + bound(h,2))
                yeigen(idof) = 0
            end do

        end do


    end subroutine mmul

!
!--------------------------------------------------------------------------------------------------
!

     subroutine enforce

        !! This subroutine enforces the support boundary conditions

        use fedata

        integer :: i, idof, m, l, r, c, n, a, limit, posit, j
        real(wp) :: penal
        integer :: ito, e

        real(wp), dimension(neqn) :: kevector
        real(wp), dimension(bw,neqn) :: kmatUnbound



        ! Correct for supports
        if (.not. banded) then
            if (.not. penalty) then
                do i = 1, nb
                    idof = int(3*(bound(i,1)-1) + bound(i,2))
                    p(1:neqn) = p(1:neqn) - kmat(1:neqn, idof) * bound(i, 3)
                    p(idof) = bound(i, 3)
                    kmat(1:neqn, idof) = 0
                    kmat(idof, 1:neqn) = 0
                    kmat(idof, idof) = 1
                end do
            else
                penal = penalty_fac*maxval(kmat)
                do i = 1, nb
                    idof = int(3*(bound(i,1)-1) + bound(i,2))
                    kmat(idof, idof) = kmat(idof, idof) + penal
                    p(idof) = penal * bound(i, 3)
                end do
            end if

        else
            kmatUnbound = 0
            kmatUnbound = kmat
            kevector=0

            do i = 1, nb
                    idof = int(3*(bound(i,1)-1) + bound(i,2))

                if (idof > bw) then
                    ito = bw-1
                else
                    ito = idof-1
                end if

                do e = 1, ito
                    kevector(idof-e)=kmatUnbound((1+e), (idof - e))
                    kmat((1+e), (idof - e)) = 0
                end do

                if (idof > neqn-bw) then
                    kevector(idof:neqn)=kmat(1:(neqn-idof+1), idof)
                else
                    kevector(idof:(idof+bw-1))=kmat(1:bw, idof)
                end if

                p(1:neqn) = p(1:neqn) - kevector(1:neqn) * bound(i, 3)
                p(idof) = bound(i, 3)
                kmat(1:bw, idof) = 0
                kmat(1, idof) = 1

                end do
        end if


        !print*,'kmat after enforce'
        !do l = 1,size(kmat,1)
        !    print"(24(f10.4,tr1))",kmat(l,1:neqn)
        !end do

        !print*, 'p after enforce'
        !do m = 1,size(p)
        !    print*, m, p(m)
        !end do

    end subroutine enforce
!
!--------------------------------------------------------------------------------------------------
!
    subroutine recover

        !! This subroutine recovers the element stress, element strain,
        !! and nodal reaction forces

        use fedata
        use link1
        use plane42rect

        integer :: e, i, nen, l, h
        integer :: edof(mdim)
        real(wp), dimension(mdim) :: xe, de
        real(wp), dimension(mdim, mdim) :: ke, me
        real(wp) :: young, youngy, shear, area

        ! Hint for continuum elements:
        real(wp):: nu, dens, thk
        real(wp), dimension(3) :: estrain, estress, eprincipal_stresses
        real(wp) :: esigmavm

        ! Reset force vector
        p = 0

        do e = 1, ne

            ! Find coordinates etc...
            nen = element(e)%numnode
            do i = 1,nen
                xe(3*i-2) = x(element(e)%ix(i), 1)
                xe(3*i-1) = x(element(e)%ix(i), 2)
                xe(3*i  ) = x(element(e)%ix(i), 3)
                edof(3*i-2) = 3 * element(e)%ix(i) - 2
                edof(3*i-1) = 3 * element(e)%ix(i) - 1
                edof(3*i  ) = 3 * element(e)%ix(i)
                de(3*i-2) = d(edof(3*i-2))
                de(3*i-1) = d(edof(3*i-1))
                de(3*i  ) = d(edof(3*i  ))
            end do

            ! Find stress and strain
            select case( element(e)%id )
            case( 1 )
                ! truss problem
                young = mprop(element(e)%mat)%young
                area  = mprop(element(e)%mat)%area

                call link1_ke(xe, young, area, ke)
                p(edof(1:2*nen)) = p(edof(1:2*nen)) + matmul(ke(1:2*nen,1:2*nen), de(1:2*nen))

                call link1_ss(xe, de, young, estress, estrain)
                stress(e, 1:3) = estress
                strain(e, 1:3) = estrain

            case( 2 )
                 ! continuum problem
                 young = mprop(element(e)%mat)%young
                 youngy = mprop(element(e)%mat)%youngy
                 shear = mprop(element(e)%mat)%shear
                 nu = mprop(element(e)%mat)%nu
                 thk = mprop(element(e)%mat)%thk


                 call shell41_ke(xe, young, youngy, shear, theta, nu, thk, dens, ke, me)
                 p(edof(1:3*nen)) = p(edof(1:3*nen)) + matmul(ke(1:3*nen,1:3*nen), de(1:3*nen))

                 call shell41_ss(xe, de, young, youngy, shear, theta, nu, thk, estress, estrain, eprincipal_stresses, esigmavm)
                 stress(e, 1:3) = estress
                 strain(e, 1:3) = estrain
                 principal_stresses(e, 1:3) = eprincipal_stresses
                 sigmavm(e) = esigmavm

            end select
        end do

    end subroutine recover

end module fea
