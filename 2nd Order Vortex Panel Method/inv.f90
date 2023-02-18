subroutine inv(matrix,size,inver)
    implicit none
    integer*4 :: size,iter1,iter2
    real*8, dimension(size,size) :: matrix,Identity,inver
    real*8, dimension(size,2*size) :: Augmented
 
    ForAll(iter1 = 1:size, iter2 = 1:size) Identity(iter1,iter2) = (iter1/iter2)*(iter2/iter1)
 
    Augmented(:,:size) = matrix
    Augmented(:,size+1:2*size) = Identity
    
    ! Backward Elimination
    iter1 = 1
    do while (iter1 .le. 2*size-1)
        iter2 = iter1 + 1 
    do while (iter2 .le. size)
     Augmented(iter2,:) = Augmented(iter2,:) - (Augmented(iter2,iter1)/Augmented(iter1,iter1))*Augmented(iter1,:)
        iter2=iter2+1
    enddo
        iter1=iter1+1
    enddo
 
    ! Forward Elimination
    iter1 = 1
    do while (iter1 .le. size-1)
        iter2 = iter1 + 1
    do while (iter2 .le. size)
     Augmented(iter1,:) = Augmented(iter1,:) - (Augmented(iter1,iter2)/Augmented(iter2,iter2))*Augmented(iter2,:)
        iter2=iter2+1
    enddo
        iter1=iter1+1
    enddo

    ! Diagonal matrice to Identity matrice
    ForAll(iter1 = 1:size) Augmented(iter1,:) = Augmented(iter1,:)/Augmented(iter1,iter1)
 
    inver = Augmented(:,size+1:2*size)
return
end