!=============================================================================================
! Calculation of RMSD by a the quaternion-based characteristic polynomial (QCP) algorithm of Theobald [1].
! 
! [1] Theobald DL. Rapid calculation of RMSDs using a quaternion-based characteristic polynomial. 
!     Acta Cryst., A61:478, 2005.  doi:10.1107/50108767305015266
!
! Written by John D. Chodera <jchodera@gmail.com>, Dill lab, UCSF, 2006.
!=============================================================================================
module theobald_rmsd

  use numeric_kinds ! for defined precisions

  implicit none

  private
  public :: ls_rmsd
       
  !=============================================================================================
  ! Constants.
  !=============================================================================================
  
  logical, parameter :: debug = .FALSE. ! Flag to enable printing of debugging information.

  !=============================================================================================
  ! Data types.
  !=============================================================================================

  !=============================================================================================
  ! Private module data.
  !=============================================================================================

contains

  !=============================================================================================
  ! Compute the rigid-body aligned RMSD between two sets of vectors.
  !=============================================================================================
  function ls_rmsd(N, x_original_nk, y_original_nk)

    ! Parameters.
    integer, intent(in)                       :: N
      ! The number of vectors in the set to be aligned.
    real(sp), dimension(N,3), intent(in) :: x_original_nk, y_original_nk
      ! Coordinate sets for which the LS-RMSD is to be computed.
      ! x_nk(n,k) is the kth dimension (k=1..3) of atom n (n=1..N)

    ! Return value.
    real(sp) :: ls_rmsd 
      ! The computed rigid-body aligned RMSD between x and y.

    ! Local variables.
    real(sp), dimension(3) :: x_centroid_k, y_centroid_k
      ! x_centroid_k(k) is the kth dimension of the centroid of x_original_nk.
      ! y_centroid_k(k) is the kth dimension of the centroid of y_original_nk.
    real(sp), dimension(N,3) :: x_nk, y_nk
      ! x_original_nk and y_original_nk translated to their centroids.
    integer :: nIndex
      ! Loop indices.
    real(sp) :: G_x, G_y
      ! Inner products of structures x and y.
    real(sp), dimension(3,3) :: M
      ! Inner product matrix M (Eq. 4)
    real(sp), dimension(4,4) :: K
      ! 4x4 symmetric key matrix
    real(sp) :: C_4, C_3, C_2, C_1, C_0
    real(sp) :: lambda, lambda_old
    integer :: iteration
    integer, parameter :: maxits = 50
    real(sp), parameter :: tolerance = 1.0e-6
    real(sp) :: lambda2, a, b
    real(sp) :: rmsd2
    
    ! Compute centroids of x_original_nk and y_original_nk.
    x_centroid_k = sum(x_original_nk,1) / real(N,sp)
    y_centroid_k = sum(y_original_nk,1) / real(N,sp)

    ! Translate x and y to their centroids.
    forall(nIndex = 1:N)
       x_nk(nIndex,:) = x_original_nk(nIndex,:) - x_centroid_k(:)
       y_nk(nIndex,:) = y_original_nk(nIndex,:) - y_centroid_k(:)
    end forall

    ! Comptue inner products of individual structures.
    G_x = sum(x_nk**2)
    G_y = sum(y_nk**2)
    
    ! Compute the inner product matrix M (Eq. 4 from Ref. [1]).
    M = matmul(transpose(y_nk), x_nk)

    ! Form the 4x4 symmetric key matrix K.
    K(1,1) = M(1,1) + M(2,2) + M(3,3)
    K(1,2) = M(2,3) - M(3,2)
    K(1,3) = M(3,1) - M(1,3)
    K(1,4) = M(1,2) - M(2,1)
    K(2,1) = K(1,2)
    K(2,2) = M(1,1) - M(2,2) - M(3,3)
    K(2,3) = M(1,2) + M(2,1)
    K(2,4) = M(3,1) + M(1,3)
    K(3,1) = K(1,3)
    K(3,2) = K(2,3)
    K(3,3) = -M(1,1) + M(2,2) - M(3,3)
    K(3,4) = M(2,3) + M(3,2)
    K(4,1) = K(1,4)
    K(4,2) = K(2,4)
    K(4,3) = K(3,4)
    K(4,4) = -M(1,1) - M(2,2) + M(3,3)

    ! Compute coefficients of the characteristic polynomial.
    ! TODO: Can replace these with formulas that save a few FLOPS (see Eq. 7).
    C_4 = 1.0
    C_3 = 0.0
    C_2 = - 2.0 * sum(M**2)
    C_1 = - 8.0 * det(M)
    C_0 = det(K)
    
    ! Construct inital guess at lambda using upper bound.
    lambda = (G_x + G_y) / 2.0

    ! Iterate Newton-Raphson scheme.
    do iteration = 1,maxits
       ! Save old iterate.
       lambda_old = lambda

       ! Update lambda using Newton iteration.
       ! This is an optimized version of
       ! lambda = lambda_old - P(lambda) / (dP/dlambda)
       lambda2 = lambda_old*lambda_old
       b = (lambda2 + C_2)*lambda_old
       a = b + C_1
       lambda = lambda_old - (a*lambda_old + C_0) / (2.0*lambda2*lambda_old + b + a)

       ! Check for convergence.
       if (abs(lambda - lambda_old) < abs(tolerance*lambda)) exit
    end do
    
    ! Return RMSD.
    ! If rmsd2 is slightly negative, they are likely the same point -- rmsd should be zero.
    rmsd2 = (G_x + G_y - 2.0 * lambda) / real(N,dp);
    ls_rmsd = 0.0
    if(rmsd2 > 0) ls_rmsd = sqrt(rmsd2)
    ! write(*,*) 'rmsd = ', ls_rmsd, ' after ', iteration, ' iterations'

!    if(isnan(ls_rmsd)) then       
!       write(*,*) 'rmsd = ', ls_rmsd, ' after ', iteration, ' iterations'
!       write(*,*) 'G_x + G_y = ', G_x + G_y, ', - 2 * lambda = ', -2*lambda, ', sum = ', G_x + G_y - 2*lambda
!       !stop
!       ls_rmsd = -1
!    end if

  end function ls_rmsd

  !=============================================================================================
  ! Compute the determinant of matrix A recursively.
  !
  ! TODO: Replace 'submatrix' stuff with a specialized case for 4x4?
  !=============================================================================================
  recursive real(sp) function det(A) result(res)

    ! Parameters.
    real(sp), dimension(:,:), intent(in) :: A

    ! Local variables.
    integer :: i, pow
    real(sp) ::  d;
    real(sp), pointer :: B(:,:)
    character*1 :: c
    integer :: n

    n = size(A,1)
    
    d=0.0
    if (n==2) then
       res = A(1,1)*A(2,2)-A(1,2)*A(2,1)
       return
    else if (n==3) then
       d = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
       d = d-A(1,2)*(A(2,1)*A(3,3)-A(3,1)*A(2,3))
       d = d+A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
       res = d
       return
    else if (n==4) then
       res = A(1,1)*det(A((/2,3,4/),(/2,3,4/))) &
            -A(1,2)*det(A((/2,3,4/),(/1,3,4/))) &
            +A(1,3)*det(A((/2,3,4/),(/1,2,4/))) &
            -A(1,4)*det(A((/2,3,4/),(/1,2,3/))) 
       return
    else
       write(*,*) 'det for matrices bigger than 4x4 not implemented'
       stop
    end if
    
  end function det

end module theobald_rmsd
