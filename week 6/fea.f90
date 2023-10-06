module fea

    !! This module contains procedures that are common to FEM-analysis

    implicit none
    save
    private
    public :: displ, initial, buildload, buildstiff, enforce, recover

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
        neqn = 2*nn

        if (.not. banded) then
            allocate (kmat(neqn, neqn))
        else
            do e = 1, ne
                do i = 1, element(e)%numnode
                    edof(2*i-1) = 2 * element(e)%ix(i) - 1
                    edof(2*i)   = 2 * element(e)%ix(i)
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

        ! Initial stress and strain
        strain = 0
        stress = 0
    end subroutine initial
!
!--------------------------------------------------------------------------------------------------
!
    subroutine displ

        !! This subroutine calculates displacements

        use fedata
        use numeth
        use processor

        integer :: e, i
        real(wp), dimension(:), allocatable :: plotval

        ! Build load-vector
        call buildload

        ! Build stiffness matrix
        call buildstiff

        ! Remove rigid body modes
        call enforce


        !print*, 'ix'
        !print*, size(ix)
        !do i = 1, sixe(IX)
        !    print'(f5.2)', IX(i,:)
        !end do

        if (.not. banded) then
            ! Factor stiffness matrix
            call factor(kmat)
            ! Solve for displacement vector
            call solve(kmat, p)
        else
            print *, 'ERROR in fea/displ'
            print *, 'Band form not implemented -- you need to add your own code here'
            stop
        end if

        ! Transfer results
        d(1:neqn) = p(1:neqn)

        ! Recover stress
        call recover

        ! Output results
        call output


        allocate (plotval(ne))
        do e = 1, ne
            if (element(e)%id == 1) then
                plotval(e) = stress(e,1)
            else if (element(e)%id == 2) then
                plotval(e) = stress(e,1)
                !print*, e, plotval(e)
            end if
        end do

        ! Plot deformed structure
        call plotmatlabdef('Deformed')

        ! Plot element values
        call plotmatlabeval('Stresses',plotval)

    end subroutine displ
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

        ! Build load vector
        p(1:neqn) = 0

        do i = 1, np
            select case(int(loads(i, 1)))
            case( 1 )
            	! Build nodal load contribution
            	dof = loads(i, 2)*2-2+loads(i,3)
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
                    xe(2*j-1) = x(element(en)%ix(j),1)
                    xe(2*j  ) = x(element(en)%ix(j),2)
                    edof(2*j-1) = 2 * element(en)%ix(j) - 1
                    edof(2*j)   = 2 * element(en)%ix(j)
                end do


                call plane42rect_re(xe,eface, fe, thk, re)


                do l = 1, size(re)
                    p(edof(l)) = p(edof(l)) + re(l)
                end do

            case default
                print *, 'ERROR in fea/buildload'
                print *, 'Load type not known'
                stop
            end select
        end do

        !do m = 1,size(p)
        !    print*, m, p(m)
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

        ! Hint for modal analysis:
        ! real(wp), dimension(mdim, mdim) :: me
        real(wp) :: young, area
        ! Hint for modal analysis and continuum elements:
        real(wp) :: nu, dens, thk

        ! Reset stiffness matrix
        if (.not. banded) then
            kmat = 0
        else
            kmat = 0
            ! something else?

        end if

        do e = 1, ne

            ! Find coordinates and degrees of freedom
            nen = element(e)%numnode

            do i = 1, nen
                 xe(2*i-1) = x(element(e)%ix(i),1)
                 xe(2*i  ) = x(element(e)%ix(i),2)
                 edof(2*i-1) = 2 * element(e)%ix(i) - 1
                 edof(2*i)   = 2 * element(e)%ix(i)
            end do

            ! Gather material properties and find element stiffness matrix
            select case( element(e)%id )
            case( 1 )
                 young = mprop(element(e)%mat)%young
                 area  = mprop(element(e)%mat)%area
                 call link1_ke(xe, young, area, ke)
            case( 2 )
                 young = mprop(element(e)%mat)%young
                 nu  = mprop(element(e)%mat)%nu
                 thk  = mprop(element(e)%mat)%thk

                 call plane42rect_ke(xe, young, nu, thk, ke)

            end select

            ! Assemble into global matrix
            if (.not. banded) then
                do i = 1, 2*nen
                    do j = 1, 2*nen
                        kmat(edof(i), edof(j)) = kmat(edof(i), edof(j)) + ke(i, j)
                    end do
                end do

                !do m = 1, size(kmat, dim = 1)
                !    print *, kmat(m,:)
                !end do

                ! Hint: Can you eliminate the loops above by using a different Fortran array syntax?

            else
                do c = 1,2*nen
                    do r = 1,2*nen
                        krow = edof(r) - edof(c) +1
                        kcol = edof(c)

                        if (krow > 0) then
                            kmat(krow,kcol) = kmat(krow,kcol) + ke(r,c)

                        end if
                    end do
                end do

            end if
        end do


        !print*,'kmat'
        !do l = 1,size(kmat,1)
        !    print"(24(f6.2,tr1))",kmat(l,1:neqn)
        !end do

    end subroutine buildstiff
!
!--------------------------------------------------------------------------------------------------
!
    subroutine enforce

        !! This subroutine enforces the support boundary conditions

        use fedata

        integer :: i, idof, m, l, r, c, n, a, limit, posit
        real(wp) :: penal

        real(wp), dimension(neqn) :: kevector


        ! Correct for supports
        if (.not. banded) then
            if (.not. penalty) then
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    p(1:neqn) = p(1:neqn) - kmat(1:neqn, idof) * bound(i, 3)
                    p(idof) = bound(i, 3)
                    kmat(1:neqn, idof) = 0
                    kmat(idof, 1:neqn) = 0
                    kmat(idof, idof) = 1
                end do
            else
                penal = penalty_fac*maxval(kmat)
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    kmat(idof, idof) = kmat(idof, idof) + penal
                    p(idof) = penal * bound(i, 3)
                end do
            end if
        else
            do i = 1, nb
                idof = int(2*(bound(i,1)-1) + bound(i,2))
                p(idof) = bound(i, 3)


                print*,'kmat'
                do l = 1,size(kmat,1)
                    print"(24(f6.2,tr1))",kmat(l,1:neqn)
                end do

                kevector = 0

                if (idof <= bw) then
                    posit = 1
                    do n = 1, idof
                        c = n
                        r = idof + 1 - n
                        kevector(n) = kmat(r,c)
                        kmat(r,c) = 0
                    end do

                else
                    posit = idof - bw
                    do n = 1, bw
                        c = idof - bw + n
                        r = bw + 1 - n
                        kevector(n+posit) = kmat(r,c)
                        kmat(r,c) = 0
                    end do
                end if

                print*,'kevector'
                print'(f6.2)',kevector

                ! find and store column elements
                if (idof < neqn-bw+2) then
                    do a = 1,bw-1
                        kevector(n+a-2+posit) = kmat(a+1,idof)
                    end do
                else
                    do a = 1,(neqn-idof+1)
                        print*,a
                        kevector(n+a-2+posit) = kmat(a+1,idof)
                    end do

                end if

                p(1:neqn) = p(1:neqn) - kevector * bound(i, 3)
                kmat(:,idof) = 0
                kmat(1,idof) = 1

            end do
        end if



    print*,'kmat'
    do l = 1,size(kmat,1)
        print"(24(f6.2,tr1))",kmat(l,1:neqn)
    end do

    print*,'kevector'
    print'(f6.2)',kevector


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

        integer :: e, i, nen
        integer :: edof(mdim)
        real(wp), dimension(mdim) :: xe, de
        real(wp), dimension(mdim, mdim) :: ke
        real(wp) :: young, area
! Hint for continuum elements:
        real(wp):: nu, dens, thk
        real(wp), dimension(3) :: estrain, estress

        ! Reset force vector
        p = 0

        do e = 1, ne

            ! Find coordinates etc...
            nen = element(e)%numnode
            do i = 1,nen
                xe(2*i-1) = x(element(e)%ix(i), 1)
                xe(2*i)   = x(element(e)%ix(i), 2)
                edof(2*i-1) = 2 * element(e)%ix(i) - 1
                edof(2*i)   = 2 * element(e)%ix(i)
                de(2*i-1) = d(edof(2*i-1))
                de(2*i)   = d(edof(2*i))
            end do

            ! Find stress and strain
            select case( element(e)%id )
            case( 1 )
                young = mprop(element(e)%mat)%young
                area  = mprop(element(e)%mat)%area

                call link1_ke(xe, young, area, ke)
                p(edof(1:2*nen)) = p(edof(1:2*nen)) + matmul(ke(1:2*nen,1:2*nen), de(1:2*nen))

                call link1_ss(xe, de, young, estress, estrain)
                stress(e, 1:3) = estress
                strain(e, 1:3) = estrain

            case( 2 )
                 young = mprop(element(e)%mat)%young
                 nu  = mprop(element(e)%mat)%nu
                 thk  = mprop(element(e)%mat)%thk

                 call plane42rect_ke(xe, young, nu, thk, ke)
                 p(edof(1:2*nen)) = p(edof(1:2*nen)) + matmul(ke(1:2*nen,1:2*nen), de(1:2*nen))

                 call plane42rect_ss(xe, de, young, nu, estress, estrain)
                 stress(e, 1:3) = estress
                 strain(e, 1:3) = estrain

            end select
        end do
    end subroutine recover

end module fea
