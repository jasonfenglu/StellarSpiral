module ImageTool

type image
        double precision,allocatable    ::imdata(:,:)
        integer                         ::DataSize(2)
        double precision                ::Origin(2)
        double precision                ::Vertex(2)

endtype image

contains


subroutine InitImage(ImageObj,imdata,datasize,origin,vertex)
        type(image) ImageObj
        double precision,allocatable    ::imdata(:,:)
        integer                         ::DataSize(2)
        double precision                ::Origin(2)
        double precision                ::Vertex(2)
        
        allocate(ImageObj%imdata(DataSize(1),DataSize(2)))
        ImageObj%imdata   = imdata
        ImageObj%DataSize = DataSize
        ImageObj%Origin   = Origin
        ImageObj%Vertex   = Vertex


endsubroutine

function ToRealCoord(p,ImageObj)
        IMPLICIT NONE
        type (image)            ImageObj
        double precision        ::p(2)
        double precision        ::ToRealCoord(2)

        ToRealCoord =  &
      (ImageObj%Vertex - ImageObj%Origin)/(ImageObj%datasize - 1.d0) & 
        *(p - ImageObj%DataSize) + ImageObj%Origin
endfunction ToRealCoord     

function DataCoord(p,ImageObj)
        IMPLICIT NONE
        type (image)            ImageObj
        double precision        ::p(2)
        Integer                 ::DataCoord(2)

        DataCoord =  INT(&
       (REAL(ImageObj%DataSize) - 1.d0)/(ImageObj%Vertex - ImageObj%Origin)& 
        *(p - ImageObj%Origin) + REAL(ImageObj%DataSize) )
endfunction DataCoord     

function PixelValue()
        IMPLICIT NONE
        type (image)                    ImageObj
!       double precision                ::p(2)
        double precision                ::PixelValue
!       Integer                         ::pp(2)

!       write(*,*) 'p',p
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!'
!       write(*,*) DataCoord(p,ImageObj)
!       write(*,*)'pp',p,pp
        PixelValue = 1.d0

endfunction
        
subroutine DeallocateImage(ImageObj)
        type(image)     ImageObj

        deallocate(ImageObj%imdata)
endsubroutine

endmodule
