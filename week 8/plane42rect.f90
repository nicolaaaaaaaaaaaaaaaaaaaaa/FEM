module plane42rect

    !! This module contains subroutines specific to the plane42 element
    !!
    !! The plane42 element has 4 nodes. Each node has 2 degrees-of-freedom,
    !! namely, displacement along the \(x\)- and \(y\)-coordinate directions.
    !!
    !! The nodes are numbered counter-clockwise as indicated below. The
    !! figure also shows how the element edges are labelled. For example,
    !! edge 1 is the edge between element node 1 and 2.
    !!
    !!       N4    E3    N3
    !!          o------o
    !!          |      |
    !!       E4 |      | E2
    !!          |      |
    !!          o------o
    !!       N1    E1    N2
    !!
    !!
    !! `N1` = element node 1, `N2` = element node 2, etc
    !! `E1` = element face 1, `E2` = element face 2, etc

    use types
    implicit none
    save

    private
    public :: plane42rect_ke, plane42rect_re, plane42rect_ss

contains

    subroutine plane42rect_ke(xe, young, nu, thk, ke)

        !! This subroutine constructs the stiffness matrix for
        !! a rectangular 4-noded quad element.

        real(wp), intent(in) :: young
            !! Young's Modulus for this element
        real(wp), intent(in) :: nu
            !! Poisson's Ratio for this element
        real(wp), intent(in) :: thk
            !! Thickness of this element
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed configuration
            !!
            !! * `xe(1:2)` = \((x,y)\)-coordinates of element node 1
            !! * `xe(3:4)` = \((x,y)\)-coordinates of element node 1
            !! * `xe(5:6)` = \((x,y)\)-coordinates of element node 2
            !! * `xe(7:8)` = \((x,y)\)-coordinates of element node 2
            !!
            !! See also [[plane42rect]]
        real(wp), dimension(:,:), intent(out) :: ke
            !! Stiffness matrix

        real(wp) :: cmat(3,3), fact, aa, bb
        real(wp) :: t1, t2, t3, t4, t5, t6, t9, t11, t13, t15, &
                        t16, t17, t18, t19, t22, t25, t29, t33, &
                        t34, t38, t42, t43, t47, t51, t55, t56, &
                        t57, t58, t59, t62, t66, t70, t74, t78, &
                        t82, t86, t90, t94, t98, t102, t106, t110
        real(wp) :: d11, d12, d13, d22, d23, d33

        real(wp) :: coord, w1, w2, eta, xi, volume, ng
        real(wp), dimension(1,2) :: g, weight

        real(wp), dimension(2,8) :: n
            ! shape function
        real(wp), dimension(3,8) :: bmat
        real(wp), dimension(2,2) :: jac
        real(wp) :: detjac


        real(wp) :: i,j            ! variables for integration

        ke = 0
        volume = 0


        ! Gaussian points
        coord = 1/sqrt(3.0)
        g(1,1) = +coord
        g(1,2) = -coord


        ! weights
        weight(1,1) = 1
        weight(1,2) = 1


        ! constitutive matrix (plane stress)
        cmat = 0
        fact = young/(1-nu**2)
        cmat(1,1) = fact
        cmat(1,2) = fact*nu
        cmat(2,1) = fact*nu
        cmat(2,2) = fact
        cmat(3,3) = fact*(1-nu)/2


        ng = 2

        do i = 1, ng
            do j = 1, ng
                eta = g(1,i)
                xi = g(1,j)
                w1 = weight(1,i)
                w2 = weight(1,j)

                call shape_42(xe,xi,eta,n,bmat,jac,detjac)

                ke = ke + w1*w2*thk*matmul(matmul(transpose(bmat),cmat),bmat)*detjac
                volume = volume + thk*detjac*w1*w2

            end do
        end do

    end subroutine plane42rect_ke
!
!--------------------------------------------------------------------------------------------------
!
    subroutine shape_42(xe,xi,eta,n,bmat,jac,detjac)

        real(wp), dimension(:), intent(in) :: xe

        real(wp), intent(in) :: xi
        real(wp), intent(in) :: eta
            ! coordinates of Gauss point

        real(wp), dimension(2,8), intent(out) :: n
            ! shape function
        real(wp), dimension(3,8), intent(out) :: bmat
        real(wp), dimension(2,2), intent(out) :: jac
        real(wp), intent(out) :: detjac

        real(wp), dimension(2,4) :: nder
        real(wp), dimension(3,4) :: l
        real(wp), dimension(4,2) :: xemat
        real(wp), dimension(2,2) :: gam
        real(wp), dimension(4,4) :: gamtilda
        real(wp), dimension(4,8) :: ntilda
        real(wp):: n1, n2, n3, n4


        n = 0
        bmat = 0
        jac = 0
        detjac = 0

        ! calculate L
        l(1,1) = 1
        l(1,2) = 0
        l(1,3) = 0
        l(1,4) = 0
        l(2,1) = 0
        l(2,2) = 0
        l(2,3) = 0
        l(2,4) = 1
        l(3,1) = 0
        l(3,2) = 1
        l(3,3) = 1
        l(3,4) = 0

        ! calculate jac
        nder(1,1) = - (1-eta)
        nder(1,2) = + (1-eta)
        nder(1,3) = + (1+eta)
        nder(1,4) = - (1+eta)
        nder(2,1) = - (1-xi)
        nder(2,2) = - (1+xi)
        nder(2,3) = + (1+xi)
        nder(2,4) = + (1-xi)
        nder = 0.25*nder

        xemat(1,1) = xe(1)
        xemat(1,2) = xe(2)
        xemat(2,1) = xe(3)
        xemat(2,2) = xe(4)
        xemat(3,1) = xe(5)
        xemat(3,2) = xe(6)
        xemat(4,1) = xe(7)
        xemat(4,2) = xe(8)


        jac = matmul(nder, xemat)

        ! calculate norm(jac)
        detjac = jac(1,1)*jac(2,2) - jac(2,1)*jac(1,2)


        ! calculate Gamma
        gam(1,1) = jac(2,2)
        gam(1,2) = -jac(1,2)
        gam(2,1) = -jac(2,1)
        gam(2,2) = jac(1,1)
        gam = 1/detjac*gam

        gamtilda(1,1) = gam(1,1)
        gamtilda(1,2) = gam(1,2)
        gamtilda(1,3) = 0
        gamtilda(1,4) = 0
        gamtilda(2,1) = gam(2,1)
        gamtilda(2,2) = gam(2,2)
        gamtilda(2,3) = 0
        gamtilda(2,4) = 0
        gamtilda(3,1) = 0
        gamtilda(3,2) = 0
        gamtilda(3,3) = gam(1,1)
        gamtilda(3,4) = gam(1,2)
        gamtilda(4,1) = 0
        gamtilda(4,2) = 0
        gamtilda(4,3) = gam(2,1)
        gamtilda(4,4) = gam(2,2)


        ntilda(1,1) = nder(1,1)
        ntilda(1,2) = 0
        ntilda(1,3) = nder(1,2)
        ntilda(1,4) = 0
        ntilda(1,5) = nder(1,3)
        ntilda(1,6) = 0
        ntilda(1,7) = nder(1,4)
        ntilda(1,8) = 0
        ntilda(2,1) = nder(2,1)
        ntilda(2,2) = 0
        ntilda(2,3) = nder(2,2)
        ntilda(2,4) = 0
        ntilda(2,5) = nder(2,3)
        ntilda(2,6) = 0
        ntilda(2,7) = nder(2,4)
        ntilda(2,8) = 0
        ntilda(3,1) = 0
        ntilda(3,2) = nder(1,1)
        ntilda(3,3) = 0
        ntilda(3,4) = nder(1,2)
        ntilda(3,5) = 0
        ntilda(3,6) = nder(1,3)
        ntilda(3,7) = 0
        ntilda(3,8) = nder(1,4)
        ntilda(4,1) = 0
        ntilda(4,2) = nder(2,1)
        ntilda(4,3) = 0
        ntilda(4,4) = nder(2,2)
        ntilda(4,5) = 0
        ntilda(4,6) = nder(2,3)
        ntilda(4,7) = 0
        ntilda(4,8) = nder(2,4)

        bmat = matmul(matmul(l,gamtilda),ntilda)

        n1 = 0.25*(1-xi)*(1-eta)
        n2 = 0.25*(1+xi)*(1-eta)
        n3 = 0.25*(1+xi)*(1+eta)
        n4 = 0.25*(1-xi)*(1+eta)

        n(1,1) = n1
        n(1,2) = 0
        n(1,3) = n2
        n(1,4) = 0
        n(1,5) = n3
        n(1,6) = 0
        n(1,7) = n4
        n(1,8) = 0
        n(2,1) = 0
        n(2,2) = n1
        n(2,3) = 0
        n(2,4) = n2
        n(2,5) = 0
        n(2,6) = n3
        n(2,7) = 0
        n(2,8) = n4


    end subroutine



!
!--------------------------------------------------------------------------------------------------
!

    subroutine plane42rect_re(xe, eface, fe, thk, re)

        !! This subroutine computes the element load vector due
        !! to surface traction (traction is always perpendicular
        !! to element face).

        integer, intent(in) :: eface
            !! Element face where traction (pressure) is applied

        real(wp), intent(in) :: fe
            !! Value of surface traction (pressure)
        real(wp), intent(in) :: thk
            !! Thickness of this element
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed configuration (see also [[plane42rect_ke]])
        real(wp), intent(out) :: re(8)
            !! Element force vector
            !!
            !! * `re(1:2)` = \((f_x^1, f_y^1)\) force at element node 1 in \(x\)- and y-direction
            !! * `re(3:4)` = \((f_x^2, f_y^2)\) force at element node 1 in \(x\)- and y-direction
            !! * etc...

        real(wp), dimension(1,2) :: g, weight
        real(wp), dimension(2,1) :: jacvec
        real(wp) :: eta, xi, w1, w2, coord
        integer(wp) :: ng, i

        real(wp), dimension(2,8) :: n
        real(wp), dimension(3,8) :: bmat
        real(wp), dimension(2,2) :: jac
        real(wp) :: detjac
        real(wp), dimension(8,1) :: f

        ng = 2
        re = 0

        !coord = 1/sqrt(3.0)
        coord = 1
        g(1,1) = +coord
        g(1,2) = -coord

        weight(1,1) = 1.0
        weight(1,2) = 1.0

        if (eface == 1) then
            eta = -1
            do i = 1,ng
                xi = g(1,i)
                w1 = weight(1,i)

                call shape_42(xe,xi,eta,n,bmat,jac,detjac)

                jacvec(1,1) = - jac(1,2)
                jacvec(2,1) =  jac(1,1)

                re = re + reshape(w1*thk*fe*matmul(transpose(n),jacvec), shape(re))

            end do

        elseif (eface == 2) then
            xi = 1
            do i = 1,ng
                eta = g(1,i)
                w1 = weight(1,i)

                call shape_42(xe,xi,eta,n,bmat,jac,detjac)

                jacvec(1,1) = - jac(2,2)
                jacvec(2,1) =  jac(2,1)

                re = re + reshape(w1*thk*fe*matmul(transpose(n),jacvec), shape(re))

            end do



        elseif (eface == 3) then
            eta = 1
            do i = 1,ng
                xi = g(1,i)
                w1 = weight(1,i)

                call shape_42(xe,xi,eta,n,bmat,jac,detjac)

                jacvec(1,1) =  jac(1,2)
                jacvec(2,1) = - jac(1,1)

                re = re + reshape(w1*thk*fe*matmul(transpose(n),jacvec), shape(re))

            end do
        elseif (eface == 4) then
            xi = -1
            do i = 1,ng
                eta = g(1,i)
                w1 = weight(1,i)

                call shape_42(xe,xi,eta,n,bmat,jac,detjac)

                jacvec(1,1) =  jac(2,2)
                jacvec(2,1) = - jac(2,1)

                re = re + reshape(w1*thk*fe*matmul(transpose(n),jacvec), shape(re))

            end do

        end if

    end subroutine plane42rect_re
!
!--------------------------------------------------------------------------------------------------
!
    subroutine plane42rect_ss(xe, de, young, nu, estress, estrain, eprincipal_stresses, esigmavm)

        !! This subrotuine computes the element stress and strain (The location inside the element
        !! where stress and and strain is evaluated, is defined inside the subroutine).

        real(wp), intent(in) :: young
            !! Young's Modulus for this element
        real(wp), intent(in) :: nu
            !! Poisson's Ratio for this element
        real(wp), dimension(:), intent(in)  :: xe
            !! Nodal coordinates of this element in undeformed configuration (see also [[plane42rect_ke]])
        real(wp), dimension(:), intent(in)  :: de
            !! Displacement field
            !!
            !! * `de(1:2)` = displacement of element node 1 along \(x\)- and \(y\)-axis, respectively
            !! * `de(3:4)` = displacement of element node 2 along \(x\)- and \(y\)-axis, respectively
            !! * etc...
        real(wp), dimension(:), intent(out) :: estress
            !! Stress at a point inside the element
            !!
            !! * `estress(1)` = \(\sigma_{11}\)
            !! * `estress(2)` = \(\sigma_{22}\)
            !! * `estress(3)` = \(\sigma_{12}\)
        real(wp), dimension(:), intent(out) :: estrain
            !! Strain at a point inside the element
            !!
            !! * `estrain(1)` = \(\epsilon_{11}\)
            !! * `estrain(2)` = \(\epsilon_{22}\)
            !! * `estrain(3)` = \(\epsilon_{12}\)

        real(wp), dimension(:), intent(out) :: eprincipal_stresses
        real(wp), intent(out) :: esigmavm


        real(wp) :: xi, eta, detjac
        real(wp), dimension(2,8) :: n
        real(wp), dimension(3,8) :: bmat
        real(wp), dimension(2,2) :: jac

        real(wp) :: sigma1, sigma2, sinpsi, cospsi, psi

        real(wp) :: cmat(3, 3)
        real(wp) :: cmat1(3, 3)     ! for plane strain
        real(wp) :: cmat2(3, 3)     ! for plane stress

        real(wp) :: aa, bb, x, y, bmatmult, cmatmult1, cmatmult2


        ! Build strain-displacement matrix
        xi = 0
        eta = 0

        call shape_42(xe,xi,eta,n,bmat,jac,detjac)

        estrain = matmul(bmat, de)


        ! Build constitutive matrix 1 (plane stress)

        cmatmult1 = young/(1-nu**2)
        cmat1 = 0

        cmat1(1,1) = 1
        cmat1(1,2) = nu
        cmat1(1,3) = 0
        cmat1(2,1) = nu
        cmat1(2,2) = 1
        cmat1(2,3) = 0
        cmat1(3,1) = 0
        cmat1(3,2) = 0
        cmat1(3,3) = (1-nu)/2

        cmat1 = cmatmult1*cmat1


        ! Build constitutive matrix 2 (plane strain)

        !cmatmult2 = young/((1-2*nu)*(1+nu))
        cmat2 = 0

        cmat2(1,1) = 1-nu
        cmat2(1,2) = nu
        cmat2(1,3) = 0
        cmat2(2,1) = nu
        cmat2(2,2) = 1-nu
        cmat2(2,3) = 0
        cmat2(3,1) = 0
        cmat2(3,2) = 0
        cmat2(3,3) = (1-2*nu)/2

        cmat2 = cmatmult2*cmat2


        ! Compute element stress
        estress = matmul(cmat1, estrain)

        ! Compute principal stress and direction
        eprincipal_stresses = 0
        sigma1 = 0.5*(estress(1)+estress(2)) + sqrt(((0.5*(estress(1)-estress(2)))**2)+estress(3)**2)
        sigma2 = 0.5*(estress(1)+estress(2)) - sqrt(((0.5*(estress(1)-estress(2)))**2)+estress(3)**2)

        cospsi = (estress(1)-estress(2))/(sigma1-sigma2)
        sinpsi = -2*estress(3)/(sigma1-sigma2)
        psi = atan2(sinpsi, cospsi)/2

        ! output quantities
        eprincipal_stresses(1) = sigma1
        eprincipal_stresses(2) = sigma2
        eprincipal_stresses(3) = psi

        esigmavm = sqrt(sigma1**2 + sigma2**2 - sigma1*sigma2)

    end subroutine plane42rect_ss

end module plane42rect
