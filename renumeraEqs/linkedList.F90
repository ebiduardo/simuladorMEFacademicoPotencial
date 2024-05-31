!http://www.mathcs.emory.edu/~cheung/Courses/561/Syllabus/6-Fortran/Progs/list01.f90

MODULE ListModule

    TYPE ListElem
     REAL                    :: value;
     TYPE(ListElem), POINTER :: next;
     contains
   procedure :: PrintList
    END TYPE ListElem

contains

   SUBROUTINE PrintList(head)
    IMPLICIT NONE
    class( ListElem ), target, intent(in) :: head

    type( ListElem ), pointer :: ptr
    ptr => head
    print *, "The list is: "
    DO WHILE ( associated(ptr) )
      write(*,'(E12.6,", ")', advance='NO') ptr%value
      ptr => ptr%next
    END DO
    print *
   END SUBROUTINE

! =======================================

   subroutine insertElement(head, v_)
    IMPLICIT NONE
    type( ListElem ), pointer :: head 
    real :: v_
    type( ListElem ), pointer :: newElem
    allocate(newElem)
    newElem%value =v_ 
    newElem%next => head 
    head => newElem
   END subroutine

! =======================================

   FUNCTION removeElement(head)
    IMPLICIT NONE
    type( ListElem ), pointer :: head
    type( ListElem ), pointer :: removeElement
    removeElement => head
    IF ( ASSOCIATED(head) ) THEN
       head => head%next
    END IF
   END FUNCTION

   END MODULE ListModule

! =======================================

   PROGRAM LinkedList
    USE ListModule
    type( ListElem ), pointer :: head
    type( ListElem ), pointer :: anElem
    real :: value
    integer, parameter :: N = 4

   nullify( head )                 ! Initialize list to point to no target.

! Add the N elements
   DO i = 1, N
      CALL random_number( value )
      print*, "inserting :",  value 
      call insertElement(head, value)
   END DO
   CALL PrintList(head)
   print *, "removing elements from list...."
   anElem => removeElement(head)
   print*, "removing :",  anElem%value 
   CALL PrintList(head)
   anElem => removeElement(head)
   print*, "removing :",  anElem%value 
   CALL PrintList(head)
   end program LinkedList

