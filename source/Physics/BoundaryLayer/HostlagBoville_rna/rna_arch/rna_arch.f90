module rna_arch
 implicit none
 private

  !variables of the rna_arch
    real, public, pointer            :: valores_normalizados(:,:)
    real, public, pointer            :: wh1(:,:)
    real, public, pointer            :: ws(:,:)
    real, public, pointer            :: bh1(:,:)
    real, public, pointer            :: bs(:,:)
    real, public, pointer            :: v1(:,:)

    double precision, public, pointer :: vh1(:,:)
    double precision, public, pointer :: yh1(:,:)
    double precision, public, pointer ::  vs(:,:)
    double precision, public, pointer ::  ys(:,:)
    double precision, public, pointer :: valores_denormalizados(:,:)
end module rna_arch

module rna_class
 use rna_arch

 implicit none
 private

   type, public :: rna_AI
      integer                     :: in_id
      character(50)               :: ch_id
      integer                     :: nlines
      integer                     :: ncolumns
      real, pointer               :: sp_weights(:,:) 
      double precision, pointer   :: dp_weights(:,:) 
      real, pointer               :: sp_entries(:,:) !entries for the rna
    contains
      procedure :: get => rna_get_entries
      procedure :: weights => rna_get_weights
   end type rna_AI
  contains
   subroutine rna_get_entries(this)
      class(rna_AI), intent(inout) :: this
      character(120)            :: path, finput
      logical                   :: lexist
      integer                   :: i
      integer                   :: j

      !path="/home/enver/work/Research/rna/" 
      path="/home/enver/work/Research/Implementing_rna/ArtificialNeuralNetworkForPBL/"
      finput=trim(path)//trim(this%ch_id)//'/entrada.dat'

   !open file if exist
     inquire(file=finput,exist=lexist)
     if(lexist) then
       this%in_id=226
       open(unit=this%in_id,file=finput,action='read',form='formatted')
     else
       print*, finput,'file do not exist'
       stop
     end if  

   !read file
      read(unit=this%in_id,fmt=*) this%nlines, this%ncolumns
    
   !allocate entries from the readed file
      allocate(this%sp_entries(max(1,this%nlines),max(1,this%ncolumns)))

   !read entries from file
      do i=1, this%nlines
        read(unit=this%in_id,fmt=*) (this%sp_entries(i,j), j=1,this%ncolumns)
      end do

   !close file
      close(unit=this%in_id)


    print*, trim(path)//trim(this%ch_id)//'/entrada.dat'
    print*, 'inside rna_get_entries subroutine'
    print*, ((this%sp_entries(i,j), i=1,this%nlines), j=1,this%ncolumns)
   end subroutine rna_get_entries

   subroutine rna_get_weights(this)
      class(rna_AI), intent(inout) :: this
      character(120)            :: path, finput
      logical                   :: lexist
      integer                   :: i
      integer                   :: j

      !path="/home/enver/work/Research/rna/" 
      path="/home/enver/work/Research/Implementing_rna/ArtificialNeuralNetworkForPBL/"
      finput=trim(path)//trim(this%ch_id)

   !open file if exist
     inquire(file=finput,exist=lexist)
     if(lexist) then
       this%in_id=226
       open(unit=this%in_id,file=finput,action='read',form='formatted')
     else
       print*, finput,'file do not exist'
       stop
     end if  

    !read file
      read(unit=this%in_id,fmt=*) this%nlines, this%ncolumns
    
    !allocate weights from the readed file
      allocate(this%sp_weights(max(1,this%nlines),max(1,this%ncolumns)))

    !read weights from file
      do i=1, this%nlines
        read(unit=this%in_id,fmt=*) (this%sp_weights(i,j), j=1,this%ncolumns)
      end do

    !close file
      close(unit=this%in_id)

      print*, 'inside rna_get_weights subroutine'
      print*, ((this%sp_weights(i,j), i=1,this%nlines), j=1,this%ncolumns)

   end subroutine rna_get_weights

end module rna_class
