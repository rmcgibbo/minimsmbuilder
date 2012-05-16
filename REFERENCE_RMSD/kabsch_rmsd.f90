!=============================================================================================
! Calculation of RMSD by a fast variant of the method of Kabsch [1,2].
! 
! [1]	W.Kabsch. A solution for the best rotation to relate two sets of vectors. Acta Cryst., A32:922, 1976.
! [2]	W.Kabsch. A discussion of the solution for the best rotation to relate two sets of vectors. Acta Cryst., A34:827, 1978.
!
! Written by John D. Chodera <jchodera@gmail.com>, Dill lab, UCSF, 2006.
!=============================================================================================
! TODO:
! - Copy LAPACK wrapper code eig() here, marking storage space as SAVE for efficiency?
!=============================================================================================
module kabsch_rmsd

  use numeric_kinds ! for defined precisions
  !use eispack, only : rs ! for symmetric eigenproblem
  use eigs_lapack, only : eig ! for symmetric eigenproblem

  implicit none

  private
  public :: ls_rmsd, ls_align
       
  !=============================================================================================
  ! Constants.
  !=============================================================================================
  
  logical, parameter :: debug = .false. ! Flag to enable printing of debugging information.

  !=============================================================================================
  ! Data types.
  !=============================================================================================

  !=============================================================================================
  ! Private module data.
  !=============================================================================================

contains

  !=============================================================================================
  ! Compute the rigid-body aligned RMSD between two sets of vectors.
  !
  ! TODO: Could optimize this by replacing call to eigenvalue decomposition rs() by cubic root
  ! method.
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
    real(dp), dimension(3) :: x_centroid_k, y_centroid_k
      ! x_centroid_k(k) is the kth dimension of the centroid of x_original_nk.
      ! y_centroid_k(k) is the kth dimension of the centroid of y_original_nk.
    real(dp), dimension(N,3) :: x_nk, y_nk
      ! x_original_nk and y_original_nk translated to their centroids.
    integer :: nIndex, kIndex
      ! Loop indices.
    real(dp) :: E0
      ! The fixed component of the residual, defined in step (a) at the end of Ref. [2].
    real(dp), dimension(3,3) :: R
      ! Matrix defined in Eqs. (3) and (4) of Ref. [2].
    real(dp), dimension(3) :: mu_k
      ! Eigenvalues of R'R, in descending order.
    real(dp), dimension(3,3) :: a
      ! Eigenvectors of R'R, stored as columns.
    real(dp), dimension(3,3) :: b
      ! Vectors deduced from a.
    integer :: ierr
      ! error code trapping for eigenvalue decomposition
    real(dp) :: sigma_3
      ! Factor needed in Step (d) of Ref. [2].
    integer :: i
      ! Loop indices
    real(dp) :: E
      ! Error function (eq. 1 of ref. [1])
    
    ! Step (a) from Ref. [2].

    ! Compute centroids of x_original_nk and y_original_nk.
    x_centroid_k = sum(x_original_nk,1) / real(N,dp)
    y_centroid_k = sum(y_original_nk,1) / real(N,dp)

    ! Translate x and y to their centroids.
    forall(nIndex = 1:N)
       x_nk(nIndex,:) = x_original_nk(nIndex,:) - x_centroid_k(:)
       y_nk(nIndex,:) = y_original_nk(nIndex,:) - y_centroid_k(:)
    end forall

    ! Compute E_0, defined in step (a) at the end of Ref. [2].
    ! E_0 = (1/2) \sum_n (|x_n|^2 + |y_n|^2)
    E0 = 0.5 * sum(sum(x_nk**2, 2) + sum(y_nk**2, 2))

    ! Compute the matrix R (Eq. 3 of Ref. [2]).
    ! NOTE: These variables have been promoted to double precision because it was found that
    ! additional precision was needed to obtain an accurate R that does not have negative
    ! eigenvalues.
    R = matmul(transpose(real(y_nk,dp)), real(x_nk,dp))

    ! Step (b) from Ref. [2].
    ! Compute the eigenvalues mu_k (in descending order) and eigenvectors a_k of the symmetric normal matrix (R'R).
    call eig(3, matmul(transpose(R), R), mu_k, a)    

    ! If any of the mu_k are slightly negative, set it to zero and emit a warning.
    ! This may indicate a precision problem.
    if( any(mu_k < 0.0) ) then
       if(debug) then
          write(*,*) 'WARNING: An element of mu_k is negative.'
          write(*,*) 'mu_k = ', mu_k
          write(*,*) 'This may indicate a precision problem.'
       end if       
       
       where(mu_k < 0.0) mu_k = 0.0
    end if
    
    ! Set a_3 = a_1 x a_2 to ensure a right-handed coordinate system.
    a(:,3) = cross(a(:,1), a(:,2))

    ! Step (c) from Ref. [2].
    ! Determine R a_k, k = 1, 2, 3, and normalize the first two to obtain b_1 and b_2.
    b(:,1:2) = matmul(R, a(:,1:2))
    b(:,1) = b(:,1) / norm( b(:,1) )
    b(:,2) = b(:,2) / norm( b(:,2) )
    ! Set b_3 = b_1 x b_2.
    b(:,3) = cross( b(:,1), b(:,2) )
    
    ! Step (d) from Ref. [2].
    ! Set sigma_3 = -1 if b3' R a_3 < 0, otherwise sigma_3 = +1.
    sigma_3 = + 1.0
    if(dot_product(b(:,3), matmul(R, a(:,3))) < 0) sigma_3 = -1.0
    
    ! Compute residual error.
    ! E = E_0 - sqrt(mu_1) - sqrt(mu_2) - sigma_3 sqrt(mu_3)
    E = E0 - sqrt(mu_k(1)) - sqrt(mu_k(2)) - sigma_3 * sqrt(mu_k(3))   

    ! Compute RMSD from residual error.
    if(E > 0) then
       ls_rmsd = sqrt(2.0 / real(N,dp) * E)
    else
       ! Occassionally, the error is small and negative if the two structures are identical.
       ls_rmsd = 0.0
    end if
    
    ! Sanity check on NaN.
    if(isnan(ls_rmsd)) then
!    gfortran chokes on this -- the following might work with it...    
!    if(.not. (ls_rmsd <= 0 .or. ls_rmsd > 0)) then
       write(*,*) 'Error: ls_rmsd is NaN.'
       write(*,*) 'x_centroid = ', x_centroid_k
       write(*,*) 'y_centroid = ', y_centroid_k
       write(*,*) 'centered atom coordinates:'
       do i = 1, N
          write(*,*) i, x_nk(i,:), '    ', y_nk(i,:)
       end do
       write(*,*) 'E0 = ', E0
       write(*,*) 'R = '
       do i = 1, 3
          write(*,*) R(i,:)
       end do
       write(*,*) 'ierr = ', ierr
       write(*,*) 'mu_k = ', mu_k
       write(*,*) 'a = '
       do i = 1, 3
          write(*,*) a(i,:)
       end do
       write(*,*) 'b = '
       do i = 1, 3
          write(*,*) b(i,:)
       end do

       stop
    end if

  end function ls_rmsd

  !=============================================================================================
  ! Compute the rotation matrix for mapping one set of vectors onto another.
  ! Both sets of vectors should already have been translated to their barycenters.
  !=============================================================================================
  subroutine compute_rotation_matrix(N, x_nk, y_nk, U)

    ! Parameters.
    integer, intent(in)                       :: N
      ! The number of vectors in the set to be aligned.
    real(sp), dimension(N,3), intent(in) :: x_nk
      ! Coordinate set to be rotated and translated.
      ! x_nk(n,k) is the kth dimension (k=1..3) of atom n (n=1..N)    
      ! Must already be translated so that origin is at barycenter.
    real(sp), dimension(N,3), intent(in) :: y_nk
      ! Coordinate set for which x_nk is to be aligned to.
      ! Must already be translated so that origin is at barycenter.
    real(sp), dimension(3,3), intent(out) :: U
      ! Rotation matrix, to be applied to x_nk to optimally align to y_nk.
   
    ! Local variables.
    real(dp), dimension(3,3) :: R
      ! Matrix defined in Eqs. (3) and (4) of Ref. [2].
    real(dp), dimension(3) :: mu_k
      ! Eigenvalues of R'R, in descending order.
    real(dp), dimension(3,3) :: a
      ! Eigenvectors of R'R, stored as columns.
    real(dp), dimension(3,3) :: b
      ! Vectors deduced from a.
    integer :: ierr
      ! error code trapping for eigenvalue decomposition

    ! Step (a) from Ref. [2].

    ! Compute the matrix R  (Eq. 3 of Ref. [2]).
    ! NOTE: These variables have been promoted to double precision because it was found that
    ! additional precision was needed to obtain an accurate R that does not have negative
    ! eigenvalues.
    R = matmul(transpose(real(y_nk,dp)), real(x_nk,dp))

    ! Step (b) from Ref. [2].
    ! Compute the eigenvalues mu_k (in descending order) and eigenvectors a_k of the symmetric normal matrix (R'R).
    call eig(3, matmul(transpose(R), R), mu_k, a)    
    
    ! If any of the mu_k are slightly negative, set it to zero and emit a warning.
    ! This may indicate a precision problem.
    if( any(mu_k < 0.0) ) then
       if(debug) then
          write(*,*) 'WARNING: An element of mu_k is negative.'
          write(*,*) 'mu_k = ', mu_k
          write(*,*) 'This may indicate a precision problem.'
       end if       
       
       where(mu_k < 0.0) mu_k = 0.0
    end if
        
    ! Set a_3 = a_1 x a_2 to ensure a right-handed coordinate system.
    a(:,3) = cross(a(:,1), a(:,2))

    ! Step (c) from Ref. [2].
    ! Determine R a_k, k = 1, 2, 3, and normalize the first two to obtain b_1 and b_2.
    b(:,1:2) = matmul(R, a(:,1:2))
    b(:,1) = b(:,1) / norm( b(:,1) )
    b(:,2) = b(:,2) / norm( b(:,2) )
    ! Set b_3 = b_1 x b_2.
    b(:,3) = cross( b(:,1), b(:,2) )

    ! Step (d) from Ref. [2].
    ! Form U according to Eq. (7) to obtain the best rotation.
    U = matmul(b, transpose(a))    

  end subroutine compute_rotation_matrix

  !=============================================================================================
  ! Align one set of vectors to another by rigid body translation and rotation.
  !=============================================================================================
  subroutine ls_align(N, x_nk, y_target_nk, atom_indices)

    ! Parameters.
    integer, intent(in)                       :: N
      ! The number of vectors in the set to be aligned.
    real(sp), dimension(N,3), intent(inout) :: x_nk
      ! Coordinate set to be rotated and translated.
      ! x_nk(n,k) is the kth dimension (k=1..3) of atom n (n=1..N)    
    real(sp), dimension(N,3), intent(in) :: y_target_nk
      ! Coordinate set for which x_nk is to be aligned to.
    integer, dimension(:), intent(in), optional :: atom_indices
      ! List of atom indices to use for alignment.
   
    ! Local variables.
    real(sp), dimension(3) :: x_centroid_k, y_centroid_k
      ! x_centroid_k(k) is the kth dimension of the centroid of x_original_nk.
      ! y_centroid_k(k) is the kth dimension of the centroid of y_original_nk.
    real(sp), dimension(N,3) :: y_nk
      ! y_original_nk translated to its centroid.
    integer :: nindices
      ! Number of atom indices (if specified).
    integer :: nIndex
      ! Loop indices.
    real(sp), dimension(3,3) :: U
      ! Rotation matrix applied to centered x to obtain optimal alignment to y.
    real(sp), dimension(3) :: x
      ! temporary vector
    integer :: k
      ! index

    ! Step (a) from Ref. [2].

    ! Compute centroids of x_original_nk and y_original_nk.
    if(present(atom_indices)) then
       nindices = size(atom_indices,1)
       x_centroid_k = sum(x_nk(atom_indices,:),1) / real(nindices,sp)
       y_centroid_k = sum(y_target_nk(atom_indices,:),1) / real(nindices,sp)       
    else
       x_centroid_k = sum(x_nk,1) / real(N,sp)
       y_centroid_k = sum(y_target_nk,1) / real(N,sp)
    end if

    ! Translate x and y to their centroids.
    forall(nIndex = 1:N)
       x_nk(nIndex,:) = x_nk(nIndex,:) - x_centroid_k(:)
       y_nk(nIndex,:) = y_target_nk(nIndex,:) - y_centroid_k(:)
    end forall

    ! Compute rotation matrix for optimal LS-alignment using atom list.
    if(present(atom_indices)) then
       call compute_rotation_matrix(nindices, x_nk(atom_indices,:), y_nk(atom_indices,:), U)
    else 
       call compute_rotation_matrix(N, x_nk, y_nk, U)
    end if
    
    ! Apply U to x and transform to centroid of y.
    forall(nIndex = 1:N)
       x_nk(nIndex,:) = matmul(U, x_nk(nIndex,:)) + y_centroid_k(:)
    end forall

  end subroutine ls_align

  !=============================================================================================
  ! Compute the cross product of two vectors.
  !=============================================================================================
  function cross(a, b)
    
    ! Parameters.
    real(dp), dimension(3) :: a, b
      ! vectors a and b

    ! Return value.
    real(dp), dimension(3) :: cross
      ! cross = a x b

    ! Compute cross product.
    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
    
  end function cross

  !=============================================================================================
  ! Compute the norm of a vector.
  !=============================================================================================
  function norm(a)
    
    ! Parameters.
    real(dp), dimension(:) :: a
      ! vector whose norm is to be computed
    
    ! Return value.
    real(dp) :: norm
      ! the two-norm

    ! Compute two-norm.
    norm = sqrt(sum(a**2))
  end function norm
    
  !=============================================================================================
  ! Compute the determinant of matrix A recursively.
  !
  ! TODO: Replace 'submatrix' stuff with a specialized case for 4x4?
  !=============================================================================================
  recursive real(sp) function det(A) result(res)

    ! Parameters.
    real(dp), dimension(:,:), intent(in) :: A

    ! Local variables.
    integer :: i, pow
    real(dp) ::  d;
    real(dp), pointer :: B(:,:)
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

end module kabsch_rmsd
