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
    public :: plane42rect_ke, plane42rect_re, plane42rect_ss, shell41_ke

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

        ke = 0

        aa = (xe(3)-xe(1))/2
        bb = (xe(8)-xe(2))/2

        ! build constitutive matrix (plane stress)
        cmat = 0
        fact = young/(1-nu**2)
        cmat(1,1) = fact
        cmat(1,2) = fact*nu
        cmat(2,1) = fact*nu
        cmat(2,2) = fact
        cmat(3,3) = fact*(1-nu)/2

        d11 = cmat(1,1)
        d12 = cmat(1,2)
        d13 = cmat(1,3)
        d22 = cmat(2,2)
        d23 = cmat(2,3)
        d33 = cmat(3,3)

        t1 = bb**2
        t2 = t1*d11
        t3 = 2*t2
        t4 = aa**2
        t5 = t4*d33
        t6 = 2*t5
        t9 = 3*d13*aa*bb
        t11 = 1/aa
        t13 = 1/bb
        t15 = (t3+t6+t9)*t11*t13/6
        t16 = d13*t1
        t17 = 4*t16
        t18 = t4*d23
        t19 = 4*t18
        t22 = 3*d12*aa*bb
        t25 = 3*aa*d33*bb
        t29 = (t17+t19+t22+t25)*t11*t13/12
        t33 = (-t3+t5)*t11*t13/6
        t34 = 2*t18
        t38 = (-t17+t34+t22-t25)*t11*t13/12
        t42 = (t2+t5+t9)*t11*t13/6
        t43 = 2*t16
        t47 = (t43+t34+t22+t25)*t11*t13/12
        t51 = (-t2+t6)*t11*t13/6
        t55 = (-t43+t19+t22-t25)*t11*t13/12
        t56 = d33*t1
        t57 = 2*t56
        t58 = t4*d22
        t59 = 2*t58
        t62 = 3*aa*d23*bb
        t66 = (t57+t59+t62)*t11*t13/6
        t70 = (-t17+t34-t22+t25)*t11*t13/12
        t74 = (-t57+t58)*t11*t13/6
        t78 = (t56+t58+t62)*t11*t13/6
        t82 = (-t43+t19-t22+t25)*t11*t13/12
        t86 = (-t56+t59)*t11*t13/6
        t90 = (t3+t6-t9)*t11*t13/6
        t94 = (t17+t19-t22-t25)*t11*t13/12
        t98 = (t2+t5-t9)*t11*t13/6
        t102 = (t43+t34-t22-t25)*t11*t13/12
        t106 = (t57+t59-t62)*t11*t13/6
        t110 = (t56+t58-t62)*t11*t13/6
        ke(1,1) = t15
        ke(1,2) = t29
        ke(1,3) = t33
        ke(1,4) = t38
        ke(1,5) = -t42
        ke(1,6) = -t47
        ke(1,7) = -t51
        ke(1,8) = -t55
        ke(2,1) = t29
        ke(2,2) = t66
        ke(2,3) = t70
        ke(2,4) = t74
        ke(2,5) = -t47
        ke(2,6) = -t78
        ke(2,7) = -t82
        ke(2,8) = -t86
        ke(3,1) = t33
        ke(3,2) = t70
        ke(3,3) = t90
        ke(3,4) = t94
        ke(3,5) = -t51
        ke(3,6) = -t82
        ke(3,7) = -t98
        ke(3,8) = -t102
        ke(4,1) = t38
        ke(4,2) = t74
        ke(4,3) = t94
        ke(4,4) = t106
        ke(4,5) = -t55
        ke(4,6) = -t86
        ke(4,7) = -t102
        ke(4,8) = -t110
        ke(5,1) = -t42
        ke(5,2) = -t47
        ke(5,3) = -t51
        ke(5,4) = -t55
        ke(5,5) = t15
        ke(5,6) = t29
        ke(5,7) = t33
        ke(5,8) = t38
        ke(6,1) = -t47
        ke(6,2) = -t78
        ke(6,3) = -t82
        ke(6,4) = -t86
        ke(6,5) = t29
        ke(6,6) = t66
        ke(6,7) = t70
        ke(6,8) = t74
        ke(7,1) = -t51
        ke(7,2) = -t82
        ke(7,3) = -t98
        ke(7,4) = -t102
        ke(7,5) = t33
        ke(7,6) = t70
        ke(7,7) = t90
        ke(7,8) = t94
        ke(8,1) = -t55
        ke(8,2) = -t86
        ke(8,3) = -t102
        ke(8,4) = -t110
        ke(8,5) = t38
        ke(8,6) = t74
        ke(8,7) = t94
        ke(8,8) = t106

        ke = ke*thk

    end subroutine plane42rect_ke
!
!--------------------------------------------------------------------------------------------------
!
subroutine shell41_ke(xe, young, nu, thk, ke)

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
        real(wp), dimension(:,:), intent(out) :: ke
            !! Stiffness matrix

        real(wp) :: cmat(3,3), bmat(3,12), aa, bb, EE, t, x, y


        real(wp) :: fact
        real(wp) :: t1, t2, t3, t4, t5, t6, t9, t11, t13, t15, &
                        t16, t17, t18, t19, t22, t25, t29, t33, &
                        t34, t38, t42, t43, t47, t51, t55, t56, &
                        t57, t58, t59, t62, t66, t70, t74, t78, &
                        t82, t86, t90, t94, t98, t102, t106, t110
        real(wp) :: d11, d12, d13, d22, d23, d33

        ke = 0

        EE = young
        t = thk

        aa = (xe(3)-xe(1))/2
        bb = (xe(8)-xe(2))/2


        Bmat(1,1) = 0.3D1 / 0.4D1 * x / aa ** 3 - 0.3D1 / 0.4D1 * x * y / aa ** 3 / bb
        Bmat(1,2) = -0.1D1 / aa / 0.4D1 + 0.3D1 / 0.4D1 * x / aa ** 2 + y / bb / aa / 0.4D1 - 0.3D1 / 0.4D1 * x * y / aa ** 2 / bb
        Bmat(1,4) = -0.3D1 / 0.4D1 * x / aa ** 3 + 0.3D1 / 0.4D1 * x * y / aa ** 3 / bb
        Bmat(1,5) = 0.1D1 / aa / 0.4D1 + 0.3D1 / 0.4D1 * x / aa ** 2 - y / bb / aa / 0.4D1 - 0.3D1 / 0.4D1 * x * y / aa ** 2 / bb
        Bmat(1,7) = -0.3D1 / 0.4D1 * x / aa ** 3 - 0.3D1 / 0.4D1 * x * y / aa ** 3 / bb
        Bmat(1,8) = 0.1D1 / aa / 0.4D1 + 0.3D1 / 0.4D1 * x / aa ** 2 + y /bb / aa / 0.4D1 + 0.3D1 / 0.4D1 * x * y / aa ** 2 / bb
        Bmat(1,10) = 0.3D1 / 0.4D1 * x / aa ** 3 + 0.3D1 / 0.4D1 * x * y / aa ** 3 / bb
        Bmat(1,11) = -0.1D1 / aa / 0.4D1 + 0.3D1 / 0.4D1 * x / aa ** 2 - y / bb / aa / 0.4D1 + 0.3D1 / 0.4D1 * x * y / aa ** 2 / bb
        Bmat(2,1) = 0.3D1 / 0.4D1 * y / bb ** 3 - 0.3D1 / 0.4D1 * x * y / bb ** 3 / aa
        Bmat(2,3) = -0.1D1 / bb / 0.4D1 + x / bb / aa / 0.4D1 + 0.3D1 / 0.4D1 * y / bb ** 2 - 0.3D1 / 0.4D1 * x * y / aa / bb ** 2
        Bmat(2,4) = 0.3D1 / 0.4D1 * y / bb ** 3 + 0.3D1 / 0.4D1 * x * y / bb ** 3 / aa
        Bmat(2,6) = -0.1D1 / bb / 0.4D1 - x / bb / aa / 0.4D1 + 0.3D1 / 0.4D1 * y / bb ** 2 + 0.3D1 / 0.4D1 * x * y / aa / bb ** 2
        Bmat(2,7) = -0.3D1 / 0.4D1 * y / bb ** 3 - 0.3D1 / 0.4D1 * x * y / bb ** 3 / aa
        Bmat(2,9) = 0.1D1 / bb / 0.4D1 + x / bb / aa / 0.4D1 + 0.3D1 / 0.4D1 * y / bb ** 2 + 0.3D1 / 0.4D1 * x * y / aa / bb ** 2
        Bmat(2,10) = -0.3D1 / 0.4D1 * y / bb ** 3 + 0.3D1 / 0.4D1 * x * y/ bb ** 3 / aa
        Bmat(2,12) = 0.1D1 / bb / 0.4D1 - x / bb / aa / 0.4D1 + 0.3D1 / 0.4D1 * y / bb ** 2 - 0.3D1 / 0.4D1 * x * y / aa / bb ** 2
        Bmat(3,1) = 0.1D1 / bb / aa - 0.3D1 / 0.4D1 * x ** 2 / aa ** 3 / bb - 0.3D1 / 0.4D1 * y ** 2 / bb ** 3 / aa
        Bmat(3,2) = 0.1D1 / bb / 0.4D1 + x / bb / aa / 0.2D1 - 0.3D1 / 0.4D1 * x ** 2 / aa ** 2 / bb
        Bmat(3,3) = 0.1D1 / aa / 0.4D1 + y / bb / aa / 0.2D1 - 0.3D1 / 0.4D1 * y ** 2 / aa / bb ** 2
        Bmat(3,4) = -0.1D1 / bb / aa + 0.3D1 / 0.4D1 * x ** 2 / aa ** 3 /bb + 0.3D1 / 0.4D1 * y ** 2 / bb ** 3 / aa
        Bmat(3,5) = 0.1D1 / bb / 0.4D1 - x / bb / aa / 0.2D1 - 0.3D1 / 0.4D1 * x ** 2 / aa ** 2 / bb
        Bmat(3,6) = -0.1D1 / aa / 0.4D1 - y / bb / aa / 0.2D1 + 0.3D1 / 0.4D1 * y ** 2 / aa / bb ** 2
        Bmat(3,7) = 0.1D1 / bb / aa - 0.3D1 / 0.4D1 * x ** 2 / aa ** 3 / bb - 0.3D1 / 0.4D1 * y ** 2 / bb ** 3 / aa
        Bmat(3,8) = -0.1D1 / bb / 0.4D1 + x / bb / aa / 0.2D1 + 0.3D1 / 0.4D1 * x ** 2 / aa ** 2 / bb
        Bmat(3,9) = -0.1D1 / aa / 0.4D1 + y / bb / aa / 0.2D1 + 0.3D1 / 0.4D1 * y ** 2 / aa / bb ** 2
        Bmat(3,10) = -0.1D1 / bb / aa + 0.3D1 / 0.4D1 * x ** 2 / aa ** 3 /bb + 0.3D1 / 0.4D1 * y ** 2 / bb ** 3 / aa
        Bmat(3,11) = -0.1D1 / bb / 0.4D1 - x / bb / aa / 0.2D1 + 0.3D1 / 0.4D1 * x ** 2 / aa ** 2 / bb
        Bmat(3,12) = 0.1D1 / aa / 0.4D1 - y / bb / aa / 0.2D1 - 0.3D1 / 0.4D1 * y ** 2 / aa / bb ** 2

        ke(1,1) = -0.10D2 * EE * (aa ** 4 - bb ** 2 * (nu - 0.7D1 / 0.2D1)* aa ** 2 / 0.5D1 + bb ** 4)&
            * t ** 3 / (0.120D3 * nu ** 2 - 0.120D3) / bb ** 3 / aa ** 3
        ke(1,2) = -t ** 3 * EE / aa ** 2 / bb * (0.4D1 * aa ** 2 * nu + aa** 2 + 0.10D2 * bb ** 2) / (nu ** 2 - 0.1D1) / 0.120D3
        ke(1,3) = t ** 3 * EE * (-(0.4D1 * nu + 0.1D1) * bb ** 2 - 0.10D2* aa ** 2) / (nu ** 2 - 0.1D1) / aa / bb ** 2 / 0.120D3
        ke(1,4) = -0.5D1 * (aa ** 4 + 0.2D1 / 0.5D1 * bb ** 2 * (nu - 0.7D1 / 0.2D1) * aa ** 2 - 0.2D1 * bb ** 4) * EE * t ** 3 / &
            (0.120D3 *nu ** 2 - 0.120D3) / bb ** 3 / aa ** 3
        ke(1,5) = EE * ((nu - 0.1D1) * aa ** 2 - 0.10D2 * bb ** 2) * t **3 / (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(1,6) = t ** 3 * EE * (-(-0.4D1 * nu - 0.1D1) * bb ** 2 - 0.5D1* aa ** 2) / (nu ** 2 - 0.1D1) / aa / bb ** 2 / 0.120D3
        ke(1,7) = 0.5D1 * EE * t ** 3 * (aa ** 4 + 0.2D1 / 0.5D1 * bb ** 2* (nu - 0.7D1 / 0.2D1) * aa ** 2 + bb ** 4) / &
            (0.120D3 * nu ** 2- 0.120D3) / bb ** 3 / aa ** 3
        ke(1,8) = -EE * t ** 3 * ((nu - 0.1D1) * aa ** 2 + 0.5D1 * bb ** 2) / (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(1,9) = t ** 3 * EE * (-(nu - 0.1D1) * bb ** 2 - 0.5D1 * aa ** 2) / (nu ** 2 - 0.1D1) / aa / bb ** 2 / 0.120D3
        ke(1,10) = t ** 3 * EE * (-0.2D1 * aa ** 2 * bb ** 2 * nu + 0.10D2* aa ** 4 + 0.7D1 * aa ** 2 * bb ** 2 - 0.5D1 * bb ** 4) &
            / (0.120D3 * nu ** 2 - 0.120D3) / bb ** 3 / aa ** 3
        ke(1,11) = t ** 3 * EE * (0.4D1 * aa ** 2 * nu + aa ** 2 - 0.5D1 *bb ** 2) / (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(1,12) = t ** 3 * EE * ((nu - 0.1D1) * bb ** 2 - 0.10D2 * aa **2) / (nu ** 2 - 0.1D1) / aa / bb ** 2 / 0.120D3
        ke(2,1) = -t ** 3 * EE / aa ** 2 / bb * (0.4D1 * aa ** 2 * nu + aa** 2 + 0.10D2 * bb ** 2) / (nu ** 2 - 0.1D1) / 0.120D3
        ke(2,2) = EE * ((nu - 0.1D1) * aa ** 2 - 0.5D1 * bb ** 2) * t ** 3/ (0.45D2 * nu ** 2 - 0.45D2) / aa / bb
        ke(2,3) = -t ** 3 * EE * nu / (0.12D2 * nu ** 2 - 0.12D2)
        ke(2,4) = -EE * ((nu - 0.1D1) * aa ** 2 - 0.10D2 * bb ** 2) * t **3 / (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(2,5) = -EE * ((nu - 0.1D1) * aa ** 2 + 0.10D2 * bb ** 2) * t **3 / (0.180D3 * nu ** 2 - 0.180D3) / aa / bb
        ke(2,7) = EE * t ** 3 * ((nu - 0.1D1) * aa ** 2 + 0.5D1 * bb ** 2)/ (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(2,8) = EE * ((nu - 0.1D1) * aa ** 2 - 0.5D1 * bb ** 2) * t ** 3/ (0.180D3 * nu ** 2 - 0.180D3) / aa / bb
        ke(2,10) = t ** 3 * EE * (0.4D1 * aa ** 2 * nu + aa ** 2 - 0.5D1 *bb ** 2) / (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(2,11) = -0.2D1 * EE * t ** 3 * ((nu - 0.1D1) * aa ** 2 + 0.5D1/ 0.2D1 * bb ** 2) / (0.90D2 * nu ** 2 - 0.90D2) / aa / bb
        ke(3,1) = t ** 3 * EE * (-(0.4D1 * nu + 0.1D1) * bb ** 2 - 0.10D2* aa ** 2) / (nu ** 2 - 0.1D1) / aa / bb ** 2 / 0.120D3
        ke(3,2) = -t ** 3 * EE * nu / (0.12D2 * nu ** 2 - 0.12D2)
        ke(3,3) = t ** 3 * EE * ((nu - 0.1D1) * bb ** 2 - 0.5D1 * aa ** 2)/ (nu ** 2 - 0.1D1) / aa / bb / 0.45D2
        ke(3,4) = t ** 3 * EE * (-(-0.4D1 * nu - 0.1D1) * bb ** 2 - 0.5D1* aa ** 2) / (nu ** 2 - 0.1D1) / aa / bb ** 2 / 0.120D3
        ke(3,6) = t ** 3 * EE * (-0.2D1 * (nu - 0.1D1) * bb ** 2 - 0.5D1 *aa ** 2) / (nu ** 2 - 0.1D1) / aa / bb / 0.90D2
        ke(3,7) = t ** 3 * EE * (bb ** 2 * nu + 0.5D1 * aa ** 2 - bb ** 2)/ (0.120D3 * nu ** 2 - 0.120D3) / aa / bb ** 2
        ke(3,9) = t ** 3 * EE * ((nu - 0.1D1) * bb ** 2 - 0.5D1 * aa ** 2)/ (nu ** 2 - 0.1D1) / aa / bb / 0.180D3
        ke(3,10) = t ** 3 * EE * (-bb ** 2 * nu + 0.10D2 * aa ** 2 + bb **2) / (0.120D3 * nu ** 2 - 0.120D3) / aa / bb ** 2
        ke(3,12) = t ** 3 * EE * (-(nu - 0.1D1) * bb ** 2 - 0.10D2 * aa **2) / (nu ** 2 - 0.1D1) / aa / bb / 0.180D3
        ke(4,1) = -0.5D1 * (aa ** 4 + 0.2D1 / 0.5D1 * bb ** 2 * (nu - 0.7D1 / 0.2D1) * aa ** 2 - 0.2D1 * bb ** 4) * EE * t ** 3 / &
            (0.120D3 *nu ** 2 - 0.120D3) / bb ** 3 / aa ** 3
        ke(4,2) = -EE * ((nu - 0.1D1) * aa ** 2 - 0.10D2 * bb ** 2) * t **3 / (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(4,3) = t ** 3 * EE * (-(-0.4D1 * nu - 0.1D1) * bb ** 2 - 0.5D1* aa ** 2) / (nu ** 2 - 0.1D1) / aa / bb ** 2 / 0.120D3
        ke(4,4) = -0.10D2 * EE * (aa ** 4 - bb ** 2 * (nu - 0.7D1 / 0.2D1)* aa ** 2 / 0.5D1 + bb ** 4) * t ** 3 / &
            (0.120D3 * nu ** 2 - 0.12D3) / bb ** 3 / aa ** 3
        ke(4,5) = t ** 3 * EE * (0.4D1 * aa ** 2 * nu + aa ** 2 + 0.10D2 *bb ** 2) / (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(4,6) = t ** 3 * EE * (-(0.4D1 * nu + 0.1D1) * bb ** 2 - 0.10D2* aa ** 2) / (nu ** 2 - 0.1D1) / aa / bb ** 2 / 0.120D3
        ke(4,7) = t ** 3 * EE * (-0.2D1 * aa ** 2 * bb ** 2 * nu + 0.10D2* aa ** 4 + 0.7D1 * aa ** 2 * bb ** 2 - 0.5D1 * bb ** 4) &
            / (0.120D3 * nu ** 2 - 0.120D3) / bb ** 3 / aa ** 3
        ke(4,8) = t ** 3 * EE * (-(0.4D1 * nu + 0.1D1) * aa ** 2 + 0.5D1 *bb ** 2) / (nu ** 2 - 0.1D1) / bb / aa ** 2 / 0.120D3
        ke(4,9) = t ** 3 * EE * ((nu - 0.1D1) * bb ** 2 - 0.10D2 * aa ** 2) / (nu ** 2 - 0.1D1) / aa / bb ** 2 / 0.120D3
        ke(4,10) = 0.5D1 * EE * t ** 3 * (aa ** 4 + 0.2D1 / 0.5D1 * bb **2 * (nu - 0.7D1 / 0.2D1) * aa ** 2 + bb ** 4) / &
            (0.120D3 * nu ** 2- 0.120D3) / bb ** 3 / aa ** 3
        ke(4,11) = EE * t ** 3 * ((nu - 0.1D1) * aa ** 2 + 0.5D1 * bb ** 2) / (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(4,12) = t ** 3 * EE * (-(nu - 0.1D1) * bb ** 2 - 0.5D1 * aa **2) / (nu ** 2 - 0.1D1) / aa / bb ** 2 / 0.120D3
        ke(5,1) = EE * ((nu - 0.1D1) * aa ** 2 - 0.10D2 * bb ** 2) * t **3 / (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(5,2) = -EE * ((nu - 0.1D1) * aa ** 2 + 0.10D2 * bb ** 2) * t **3 / (0.180D3 * nu ** 2 - 0.180D3) / aa / bb
        ke(5,4) = t ** 3 * EE * (0.4D1 * aa ** 2 * nu + aa ** 2 + 0.10D2 *bb ** 2) / (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(5,5) = EE * ((nu - 0.1D1) * aa ** 2 - 0.5D1 * bb ** 2) * t ** 3/ (0.45D2 * nu ** 2 - 0.45D2) / aa / bb
        ke(5,6) = t ** 3 * EE * nu / (0.12D2 * nu ** 2 - 0.12D2)
        ke(5,7) = t ** 3 * EE * (-(0.4D1 * nu + 0.1D1) * aa ** 2 + 0.5D1 *bb ** 2) / (nu ** 2 - 0.1D1) / bb / aa ** 2 / 0.120D3
        ke(5,8) = -0.2D1 * EE * t ** 3 * ((nu - 0.1D1) * aa ** 2 + 0.5D1 /0.2D1 * bb ** 2) / (0.90D2 * nu ** 2 - 0.90D2) / aa / bb
        ke(5,10) = -EE * t ** 3 * ((nu - 0.1D1) * aa ** 2 + 0.5D1 * bb **2) / (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(5,11) = EE * ((nu - 0.1D1) * aa ** 2 - 0.5D1 * bb ** 2) * t **3 / (0.180D3 * nu ** 2 - 0.180D3) / aa / bb
        ke(6,1) = t ** 3 * EE * (-(-0.4D1 * nu - 0.1D1) * bb ** 2 - 0.5D1* aa ** 2) / (nu ** 2 - 0.1D1) / aa / bb ** 2 / 0.120D3
        ke(6,3) = t ** 3 * EE * (-0.2D1 * (nu - 0.1D1) * bb ** 2 - 0.5D1 *aa ** 2) / (nu ** 2 - 0.1D1) / aa / bb / 0.90D2
        ke(6,4) = t ** 3 * EE * (-(0.4D1 * nu + 0.1D1) * bb ** 2 - 0.10D2* aa ** 2) / (nu ** 2 - 0.1D1) / aa / bb ** 2 / 0.120D3
        ke(6,5) = t ** 3 * EE * nu / (0.12D2 * nu ** 2 - 0.12D2)
        ke(6,6) = t ** 3 * EE * ((nu - 0.1D1) * bb ** 2 - 0.5D1 * aa ** 2)/ (nu ** 2 - 0.1D1) / aa / bb / 0.45D2
        ke(6,7) = t ** 3 * EE * (-bb ** 2 * nu + 0.10D2 * aa ** 2 + bb **2) / (0.120D3 * nu ** 2 - 0.120D3) / aa / bb ** 2
        ke(6,9) = t ** 3 * EE * (-(nu - 0.1D1) * bb ** 2 - 0.10D2 * aa **2) / (nu ** 2 - 0.1D1) / aa / bb / 0.180D3
        ke(6,10) = t ** 3 * EE * (bb ** 2 * nu + 0.5D1 * aa ** 2 - bb ** 2) / (0.120D3 * nu ** 2 - 0.120D3) / aa / bb ** 2
        ke(6,12) = t ** 3 * EE * ((nu - 0.1D1) * bb ** 2 - 0.5D1 * aa ** 2) / (nu ** 2 - 0.1D1) / aa / bb / 0.180D3
        ke(7,1) = 0.5D1 * EE * t ** 3 * (aa ** 4 + 0.2D1 / 0.5D1 * bb ** 2* (nu - 0.7D1 / 0.2D1) * aa ** 2 + bb ** 4) / &
            (0.120D3 * nu ** 2- 0.120D3) / bb ** 3 / aa ** 3
        ke(7,2) = EE * t ** 3 * ((nu - 0.1D1) * aa ** 2 + 0.5D1 * bb ** 2)/ (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(7,3) = t ** 3 * EE * (bb ** 2 * nu + 0.5D1 * aa ** 2 - bb ** 2)/ (0.120D3 * nu ** 2 - 0.120D3) / aa / bb ** 2
        ke(7,4) = t ** 3 * EE * (-0.2D1 * aa ** 2 * bb ** 2 * nu + 0.10D2* aa ** 4 + 0.7D1 * aa ** 2 * bb ** 2 - 0.5D1 * bb ** 4)&
            / (0.120D3 * nu ** 2 - 0.120D3) / bb ** 3 / aa ** 3
        ke(7,5) = t ** 3 * EE * (-(0.4D1 * nu + 0.1D1) * aa ** 2 + 0.5D1 *bb ** 2) / (nu ** 2 - 0.1D1) / bb / aa ** 2 / 0.120D3
        ke(7,6) = t ** 3 * EE * (-bb ** 2 * nu + 0.10D2 * aa ** 2 + bb **2) / (0.120D3 * nu ** 2 - 0.120D3) / aa / bb ** 2
        ke(7,7) = -0.10D2 * EE * (aa ** 4 - bb ** 2 * (nu - 0.7D1 / 0.2D1)* aa ** 2 / 0.5D1 + bb ** 4) * t ** 3 / &
            (0.120D3 * nu ** 2 - 0.120D3) / bb ** 3 / aa ** 3
        ke(7,8) = t ** 3 * EE * (0.4D1 * aa ** 2 * nu + aa ** 2 + 0.10D2 *bb ** 2) / (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(7,9) = t ** 3 * EE * (0.4D1 * bb ** 2 * nu + 0.10D2 * aa ** 2 +bb ** 2) / (0.120D3 * nu ** 2 - 0.120D3) / aa / bb ** 2
        ke(7,10) = -0.5D1 * (aa ** 4 + 0.2D1 / 0.5D1 * bb ** 2 * (nu - 0.7D1 / 0.2D1) * aa ** 2 - 0.2D1 * bb ** 4) * EE * t ** 3 /&
            (0.120D3* nu ** 2 - 0.120D3) / bb ** 3 / aa ** 3
        ke(7,11) = -EE * ((nu - 0.1D1) * aa ** 2 - 0.10D2 * bb ** 2) * t ** 3 / (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(7,12) = t ** 3 * EE * (-0.4D1 * bb ** 2 * nu + 0.5D1 * aa ** 2- bb ** 2) / (0.120D3 * nu ** 2 - 0.120D3) / aa / bb ** 2
        ke(8,1) = -EE * t ** 3 * ((nu - 0.1D1) * aa ** 2 + 0.5D1 * bb ** 2) / (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(8,2) = EE * ((nu - 0.1D1) * aa ** 2 - 0.5D1 * bb ** 2) * t ** 3/ (0.180D3 * nu ** 2 - 0.180D3) / aa / bb
        ke(8,4) = t ** 3 * EE * (-(0.4D1 * nu + 0.1D1) * aa ** 2 + 0.5D1 *bb ** 2) / (nu ** 2 - 0.1D1) / bb / aa ** 2 / 0.120D3
        ke(8,5) = -0.2D1 * EE * t ** 3 * ((nu - 0.1D1) * aa ** 2 + 0.5D1 /0.2D1 * bb ** 2) / (0.90D2 * nu ** 2 - 0.90D2) / aa / bb
        ke(8,7) = t ** 3 * EE * (0.4D1 * aa ** 2 * nu + aa ** 2 + 0.10D2 *bb ** 2) / (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(8,8) = EE * ((nu - 0.1D1) * aa ** 2 - 0.5D1 * bb ** 2) * t ** 3/ (0.45D2 * nu ** 2 - 0.45D2) / aa / bb
        ke(8,9) = -t ** 3 * EE * nu / (0.12D2 * nu ** 2 - 0.12D2)
        ke(8,10) = EE * ((nu - 0.1D1) * aa ** 2 - 0.10D2 * bb ** 2) * t **3 / (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(8,11) = -EE * ((nu - 0.1D1) * aa ** 2 + 0.10D2 * bb ** 2) * t ** 3 / (0.180D3 * nu ** 2 - 0.180D3) / aa / bb
        ke(9,1) = t ** 3 * EE * (-(nu - 0.1D1) * bb ** 2 - 0.5D1 * aa ** 2) / (nu ** 2 - 0.1D1) / aa / bb ** 2 / 0.120D3
        ke(9,3) = t ** 3 * EE * ((nu - 0.1D1) * bb ** 2 - 0.5D1 * aa ** 2)/ (nu ** 2 - 0.1D1) / aa / bb / 0.180D3
        ke(9,4) = t ** 3 * EE * ((nu - 0.1D1) * bb ** 2 - 0.10D2 * aa ** 2) / (nu ** 2 - 0.1D1) / aa / bb ** 2 / 0.120D3
        ke(9,6) = t ** 3 * EE * (-(nu - 0.1D1) * bb ** 2 - 0.10D2 * aa **2) / (nu ** 2 - 0.1D1) / aa / bb / 0.180D3
        ke(9,7) = t ** 3 * EE * (0.4D1 * bb ** 2 * nu + 0.10D2 * aa ** 2 +bb ** 2) / (0.120D3 * nu ** 2 - 0.120D3) / aa / bb ** 2
        ke(9,8) = -t ** 3 * EE * nu / (0.12D2 * nu ** 2 - 0.12D2)
        ke(9,9) = t ** 3 * EE * ((nu - 0.1D1) * bb ** 2 - 0.5D1 * aa ** 2)/ (nu ** 2 - 0.1D1) / aa / bb / 0.45D2
        ke(9,10) = t ** 3 * EE * (-0.4D1 * bb ** 2 * nu + 0.5D1 * aa ** 2- bb ** 2) / (0.120D3 * nu ** 2 - 0.120D3) / aa / bb ** 2
        ke(9,12) = t ** 3 * EE * (-0.2D1 * (nu - 0.1D1) * bb ** 2 - 0.5D1* aa ** 2) / (nu ** 2 - 0.1D1) / aa / bb / 0.90D2
        ke(10,1) = t ** 3 * EE * (-0.2D1 * aa ** 2 * bb ** 2 * nu + 0.10D2* aa ** 4 + 0.7D1 * aa ** 2 * bb ** 2 - 0.5D1 * bb ** 4)&
            / (0.120D3 * nu ** 2 - 0.120D3) / bb ** 3 / aa ** 3
        ke(10,2) = t ** 3 * EE * (0.4D1 * aa ** 2 * nu + aa ** 2 - 0.5D1 *bb ** 2) / (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(10,3) = t ** 3 * EE * (-bb ** 2 * nu + 0.10D2 * aa ** 2 + bb **2) / (0.120D3 * nu ** 2 - 0.120D3) / aa / bb ** 2
        ke(10,4) = 0.5D1 * EE * t ** 3 * (aa ** 4 + 0.2D1 / 0.5D1 * bb **2 * (nu - 0.7D1 / 0.2D1) * aa ** 2 + bb ** 4) / &
            (0.120D3 * nu ** 2- 0.120D3) / bb ** 3 / aa ** 3
        ke(10,5) = -EE * t ** 3 * ((nu - 0.1D1) * aa ** 2 + 0.5D1 * bb **2) / (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(10,6) = t ** 3 * EE * (bb ** 2 * nu + 0.5D1 * aa ** 2 - bb ** 2) / (0.120D3 * nu ** 2 - 0.120D3) / aa / bb ** 2
        ke(10,7) = -0.5D1 * (aa ** 4 + 0.2D1 / 0.5D1 * bb ** 2 * (nu - 0.7D1 / 0.2D1) * aa ** 2 - 0.2D1 * bb ** 4) * EE * t ** 3 / &
            (0.120D3* nu ** 2 - 0.120D3) / bb ** 3 / aa ** 3
        ke(10,8) = EE * ((nu - 0.1D1) * aa ** 2 - 0.10D2 * bb ** 2) * t **3 / (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(10,9) = t ** 3 * EE * (-0.4D1 * bb ** 2 * nu + 0.5D1 * aa ** 2- bb ** 2) / (0.120D3 * nu ** 2 - 0.120D3) / aa / bb ** 2
        ke(10,10) = -0.10D2 * EE * (aa ** 4 - bb ** 2 * (nu - 0.7D1 / 0.2D1) * aa ** 2 / 0.5D1 + bb ** 4) * t ** 3 / &
            (0.120D3 * nu ** 2 - 0.120D3) / bb ** 3 / aa ** 3
        ke(10,11) = -t ** 3 * EE / aa ** 2 / bb * (0.4D1 * aa ** 2 * nu +aa ** 2 + 0.10D2 * bb ** 2) / (nu ** 2 - 0.1D1) / 0.120D3
        ke(10,12) = t ** 3 * EE * (0.4D1 * bb ** 2 * nu + 0.10D2 * aa ** 2+ bb ** 2) / (0.120D3 * nu ** 2 - 0.120D3) / aa / bb ** 2
        ke(11,1) = t ** 3 * EE * (0.4D1 * aa ** 2 * nu + aa ** 2 - 0.5D1 *bb ** 2) / (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(11,2) = -0.2D1 * EE * t ** 3 * ((nu - 0.1D1) * aa ** 2 + 0.5D1/ 0.2D1 * bb ** 2) / (0.90D2 * nu ** 2 - 0.90D2) / aa / bb
        ke(11,4) = EE * t ** 3 * ((nu - 0.1D1) * aa ** 2 + 0.5D1 * bb ** 2) / (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(11,5) = EE * ((nu - 0.1D1) * aa ** 2 - 0.5D1 * bb ** 2) * t **3 / (0.180D3 * nu ** 2 - 0.180D3) / aa / bb
        ke(11,7) = -EE * ((nu - 0.1D1) * aa ** 2 - 0.10D2 * bb ** 2) * t ** 3 / (0.120D3 * nu ** 2 - 0.120D3) / bb / aa ** 2
        ke(11,8) = -EE * ((nu - 0.1D1) * aa ** 2 + 0.10D2 * bb ** 2) * t ** 3 / (0.180D3 * nu ** 2 - 0.180D3) / aa / bb
        ke(11,10) = -t ** 3 * EE / aa ** 2 / bb * (0.4D1 * aa ** 2 * nu +aa ** 2 + 0.10D2 * bb ** 2) / (nu ** 2 - 0.1D1) / 0.120D3
        ke(11,11) = EE * ((nu - 0.1D1) * aa ** 2 - 0.5D1 * bb ** 2) * t **3 / (0.45D2 * nu ** 2 - 0.45D2) / aa / bb
        ke(11,12) = t ** 3 * EE * nu / (0.12D2 * nu ** 2 - 0.12D2)
        ke(12,1) = t ** 3 * EE * ((nu - 0.1D1) * bb ** 2 - 0.10D2 * aa **2) / (nu ** 2 - 0.1D1) / aa / bb ** 2 / 0.120D3
        ke(12,3) = t ** 3 * EE * (-(nu - 0.1D1) * bb ** 2 - 0.10D2 * aa **2) / (nu ** 2 - 0.1D1) / aa / bb / 0.180D3
        ke(12,4) = t ** 3 * EE * (-(nu - 0.1D1) * bb ** 2 - 0.5D1 * aa **2) / (nu ** 2 - 0.1D1) / aa / bb ** 2 / 0.120D3
        ke(12,6) = t ** 3 * EE * ((nu - 0.1D1) * bb ** 2 - 0.5D1 * aa ** 2) / (nu ** 2 - 0.1D1) / aa / bb / 0.180D3
        ke(12,7) = t ** 3 * EE * (-0.4D1 * bb ** 2 * nu + 0.5D1 * aa ** 2- bb ** 2) / (0.120D3 * nu ** 2 - 0.120D3) / aa / bb ** 2
        ke(12,9) = t ** 3 * EE * (-0.2D1 * (nu - 0.1D1) * bb ** 2 - 0.5D1* aa ** 2) / (nu ** 2 - 0.1D1) / aa / bb / 0.90D2
        ke(12,10) = t ** 3 * EE * (0.4D1 * bb ** 2 * nu + 0.10D2 * aa ** 2+ bb ** 2) / (0.120D3 * nu ** 2 - 0.120D3) / aa / bb ** 2
        ke(12,11) = t ** 3 * EE * nu / (0.12D2 * nu ** 2 - 0.12D2)
        ke(12,12) = t ** 3 * EE * ((nu - 0.1D1) * bb ** 2 - 0.5D1 * aa **2) / (nu ** 2 - 0.1D1) / aa / bb / 0.45D2






        ! build constitutive matrix (plane stress)
        cmat = 0
        fact = young/(1-nu**2)
        cmat(1,1) = fact
        cmat(1,2) = fact*nu
        cmat(2,1) = fact*nu
        cmat(2,2) = fact
        cmat(3,3) = fact*(1-nu)/2

        d11 = cmat(1,1)
        d12 = cmat(1,2)
        d13 = cmat(1,3)
        d22 = cmat(2,2)
        d23 = cmat(2,3)
        d33 = cmat(3,3)

        t1 = bb**2
        t2 = t1*d11
        t3 = 2*t2
        t4 = aa**2
        t5 = t4*d33
        t6 = 2*t5
        t9 = 3*d13*aa*bb
        t11 = 1/aa
        t13 = 1/bb
        t15 = (t3+t6+t9)*t11*t13/6
        t16 = d13*t1
        t17 = 4*t16
        t18 = t4*d23
        t19 = 4*t18
        t22 = 3*d12*aa*bb
        t25 = 3*aa*d33*bb
        t29 = (t17+t19+t22+t25)*t11*t13/12
        t33 = (-t3+t5)*t11*t13/6
        t34 = 2*t18
        t38 = (-t17+t34+t22-t25)*t11*t13/12
        t42 = (t2+t5+t9)*t11*t13/6
        t43 = 2*t16
        t47 = (t43+t34+t22+t25)*t11*t13/12
        t51 = (-t2+t6)*t11*t13/6
        t55 = (-t43+t19+t22-t25)*t11*t13/12
        t56 = d33*t1
        t57 = 2*t56
        t58 = t4*d22
        t59 = 2*t58
        t62 = 3*aa*d23*bb
        t66 = (t57+t59+t62)*t11*t13/6
        t70 = (-t17+t34-t22+t25)*t11*t13/12
        t74 = (-t57+t58)*t11*t13/6
        t78 = (t56+t58+t62)*t11*t13/6
        t82 = (-t43+t19-t22+t25)*t11*t13/12
        t86 = (-t56+t59)*t11*t13/6
        t90 = (t3+t6-t9)*t11*t13/6
        t94 = (t17+t19-t22-t25)*t11*t13/12
        t98 = (t2+t5-t9)*t11*t13/6
        t102 = (t43+t34-t22-t25)*t11*t13/12
        t106 = (t57+t59-t62)*t11*t13/6
        t110 = (t56+t58-t62)*t11*t13/6
        ke(1,1) = t15
        ke(1,2) = t29
        ke(1,3) = t33
        ke(1,4) = t38
        ke(1,5) = -t42
        ke(1,6) = -t47
        ke(1,7) = -t51
        ke(1,8) = -t55
        ke(2,1) = t29
        ke(2,2) = t66
        ke(2,3) = t70
        ke(2,4) = t74
        ke(2,5) = -t47
        ke(2,6) = -t78
        ke(2,7) = -t82
        ke(2,8) = -t86
        ke(3,1) = t33
        ke(3,2) = t70
        ke(3,3) = t90
        ke(3,4) = t94
        ke(3,5) = -t51
        ke(3,6) = -t82
        ke(3,7) = -t98
        ke(3,8) = -t102
        ke(4,1) = t38
        ke(4,2) = t74
        ke(4,3) = t94
        ke(4,4) = t106
        ke(4,5) = -t55
        ke(4,6) = -t86
        ke(4,7) = -t102
        ke(4,8) = -t110
        ke(5,1) = -t42
        ke(5,2) = -t47
        ke(5,3) = -t51
        ke(5,4) = -t55
        ke(5,5) = t15
        ke(5,6) = t29
        ke(5,7) = t33
        ke(5,8) = t38
        ke(6,1) = -t47
        ke(6,2) = -t78
        ke(6,3) = -t82
        ke(6,4) = -t86
        ke(6,5) = t29
        ke(6,6) = t66
        ke(6,7) = t70
        ke(6,8) = t74
        ke(7,1) = -t51
        ke(7,2) = -t82
        ke(7,3) = -t98
        ke(7,4) = -t102
        ke(7,5) = t33
        ke(7,6) = t70
        ke(7,7) = t90
        ke(7,8) = t94
        ke(8,1) = -t55
        ke(8,2) = -t86
        ke(8,3) = -t102
        ke(8,4) = -t110
        ke(8,5) = t38
        ke(8,6) = t74
        ke(8,7) = t94
        ke(8,8) = t106

        ke = ke*thk

    end subroutine shell41_ke
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
        real(wp) :: aa, bb, nface(2,8), f(2)

        aa = (xe(3)-xe(1))/2
        bb = (xe(8)-xe(2))/2

        nface = 0
        f = 0
        if (eface == 1) then
            nface(1,1) = aa
            nface(1,3) = aa
            nface(2,2) = aa
            nface(2,4) = aa
            f(2) = -fe
        elseif (eface == 2) then
            nface(1,3) = bb
            nface(1,5) = bb
            nface(2,4) = bb
            nface(2,6) = bb
            f(1) = fe
        elseif (eface == 3) then
            nface(1,5) = aa
            nface(1,7) = aa
            nface(2,6) = aa
            nface(2,8) = aa
            f(2) = fe
        elseif (eface == 4) then
            nface(1,1) = bb
            nface(1,7) = bb
            nface(2,2) = bb
            nface(2,8) = bb
            f(1) = -fe
        endif
        re = matmul(transpose(nface), f) * thk
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

        real(wp) :: sigma1, sigma2, sinpsi, cospsi, psi

        real(wp) :: bmat(3, 8), cmat(3, 3)
        real(wp) :: cmat1(3, 3)     ! for plane strain
        real(wp) :: cmat2(3, 3)     ! for plane stress

        real(wp) :: aa, bb, x, y, bmatmult, cmatmult1, cmatmult2

        ! Build strain-displacement matrix

        aa = (xe(3)-xe(1))/2
        bb = (xe(8)-xe(2))/2

        y = bb
        x = -aa

        bmatmult = 1/(4*aa*bb)

        bmat = 0
        bmat(1,1) = -(bb-y)
        bmat(1,2) = 0
        bmat(1,3) = (bb-y)
        bmat(1,4) = 0
        bmat(1,5) = (bb+y)
        bmat(1,6) = 0
        bmat(1,7) = -(bb+y)
        bmat(1,8) = 0
        bmat(2,1) = 0
        bmat(2,2) = -(aa-x)
        bmat(2,3) = 0
        bmat(2,4) = -(aa+x)
        bmat(2,5) = 0
        bmat(2,6) = (aa+x)
        bmat(2,7) = 0
        bmat(2,8) = (aa-x)
        bmat(3,1) = -(aa-x)
        bmat(3,2) = -(bb-y)
        bmat(3,3) = -(aa+x)
        bmat(3,4) =  (bb-y)
        bmat(3,5) = (aa+x)
        bmat(3,6) = (bb+y)
        bmat(3,7) = (aa-x)
        bmat(3,8) = -(bb+y)

        bmat = bmatmult*bmat


        ! Compute element strain
        estrain = matmul(bmat, de)


        ! Build constitutive matrix 1 (plane strain)

        cmatmult1 = young/(1-nu**2)
        cmat1 = 0

        cmat1(1,1) = 1
        cmat1(1,2) = nu
        cmat1(1,3) = 0
        cmat1(2,1) = nu
        cmat1(2,2) = nu
        cmat1(2,3) = 0
        cmat1(3,1) = 0
        cmat1(3,2) = 0
        cmat1(3,3) = (1-nu)/2

        cmat1 = cmatmult1*cmat1


        ! Build constitutive matrix 2 (plane stress)

        cmatmult2 = young/((1-2*nu)*(1+nu))
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
        sigma1 = 0.5*(estress(1)+estress(2)) + sqrt((0.5*(estress(1)-estress(2))**2)+estress(3)**2)
        sigma2 = 0.5*(estress(1)+estress(2)) - sqrt((0.5*(estress(1)-estress(2))**2)+estress(3)**2)

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
