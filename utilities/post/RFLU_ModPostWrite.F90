! ******************************************************************************
!
! Purpose: Files of the post program (rflupost) for ParaView: unstructured-grid
!   files (.vtu), the multiblock file of a step (.vtm) and the series (.pvd).
!
! Description: A .vtu holds its arrays as raw appended binary data (VTK XML,
!   little endian, UInt64 block sizes): points and fields in double precision,
!   connectivity and offsets as Int64, cell types as UInt8. A field is either
!   on the cells (fluid, bins) or on the points (particles).
!
! ******************************************************************************

MODULE RFLU_ModPostWrite

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_postField, &
            AddPostField, &
            WriteUnstructuredVTU, &
            WriteStepVTM, &
            WriteSeriesPVD

! A named field: values(iComp,i) for the i-th cell or point
  TYPE t_postField
    CHARACTER(32) :: name
    INTEGER :: nComp
    REAL(8), ALLOCATABLE :: values(:,:)
  END TYPE t_postField

  CHARACTER(1), PARAMETER :: NEWLINE = ACHAR(10)

CONTAINS

! ******************************************************************************
!
! Purpose: Append a field to a list of fields.
!
! ******************************************************************************

  SUBROUTINE AddPostField(fields,nFields,name,values)

    TYPE(t_postField), ALLOCATABLE, INTENT(INOUT) :: fields(:)
    INTEGER, INTENT(INOUT) :: nFields
    CHARACTER(*), INTENT(IN) :: name
    REAL(8), INTENT(IN) :: values(:,:)

    TYPE(t_postField), ALLOCATABLE :: grown(:)
    INTEGER :: i

    IF ( .NOT. ALLOCATED(fields) ) ALLOCATE(fields(8))
    IF ( nFields == SIZE(fields) ) THEN
      ALLOCATE(grown(2*SIZE(fields)))
      DO i = 1,nFields
        grown(i)%name  = fields(i)%name
        grown(i)%nComp = fields(i)%nComp
        CALL MOVE_ALLOC(fields(i)%values,grown(i)%values)
      END DO ! i
      CALL MOVE_ALLOC(grown,fields)
    END IF ! nFields

    nFields = nFields + 1
    fields(nFields)%name  = name
    fields(nFields)%nComp = SIZE(values,1)
    ALLOCATE(fields(nFields)%values,SOURCE=values)

  END SUBROUTINE AddPostField

! ******************************************************************************
!
! Purpose: Write an unstructured-grid file.
!
! Input:
!   fileName      File to write
!   time, step    Written as field data TIME and STEP
!   points        Coordinates (3,nPoints)
!   connectivity  Point numbers of the cells, 0-based, all cells in a row
!   offsets       End of each cell in connectivity
!   cellTypes     VTK cell type of each cell
!   fields        Fields, fields(1:nFields)
!   onPoints      Whether the fields are on the points (else on the cells)
!
! ******************************************************************************

  SUBROUTINE WriteUnstructuredVTU(fileName,time,step,points,connectivity, &
                                  offsets,cellTypes,fields,nFields,onPoints)

    CHARACTER(*), INTENT(IN) :: fileName
    REAL(8), INTENT(IN) :: time
    INTEGER, INTENT(IN) :: step,nFields
    REAL(8), INTENT(IN) :: points(:,:)
    INTEGER(8), INTENT(IN) :: connectivity(:),offsets(:)
    INTEGER(1), INTENT(IN) :: cellTypes(:)
    TYPE(t_postField), INTENT(IN) :: fields(:)
    LOGICAL, INTENT(IN) :: onPoints

    CHARACTER(256) :: line
    CHARACTER(9) :: location
    INTEGER :: errorFlag,i,iFile
    INTEGER(8) :: offset

    OPEN(NEWUNIT=iFile,FILE=fileName,FORM='UNFORMATTED',ACCESS='STREAM', &
         STATUS='REPLACE',IOSTAT=errorFlag)
    IF ( errorFlag /= 0 ) THEN
      WRITE(*,'(A,A)') 'rflupost: cannot write ',TRIM(fileName)
      STOP 1
    END IF ! errorFlag

    CALL Put(iFile,'<?xml version="1.0"?>')
    CALL Put(iFile,'<VTKFile type="UnstructuredGrid" version="1.0" '// &
                   'byte_order="LittleEndian" header_type="UInt64">')
    CALL Put(iFile,'<UnstructuredGrid>')
    CALL Put(iFile,'<FieldData>')
    WRITE(line,'(A,ES24.16E3,A)') '<DataArray type="Float64" Name="TIME" '// &
          'NumberOfTuples="1" format="ascii">',time,'</DataArray>'
    CALL Put(iFile,TRIM(line))
    WRITE(line,'(A,I0,A)') '<DataArray type="Int64" Name="STEP" '// &
          'NumberOfTuples="1" format="ascii">',step,'</DataArray>'
    CALL Put(iFile,TRIM(line))
    CALL Put(iFile,'</FieldData>')
    WRITE(line,'(A,I0,A,I0,A)') '<Piece NumberOfPoints="',SIZE(points,2), &
          '" NumberOfCells="',SIZE(cellTypes),'">'
    CALL Put(iFile,TRIM(line))

! Offsets of the appended blocks: each is an 8-byte size and the data
    offset = 0
    CALL Put(iFile,'<Points>')
    CALL PutArray(iFile,'Float64','Points',3,offset)
    offset = offset + 8 + 8*SIZE(points,KIND=8)
    CALL Put(iFile,'</Points>')
    CALL Put(iFile,'<Cells>')
    CALL PutArray(iFile,'Int64','connectivity',1,offset)
    offset = offset + 8 + 8*SIZE(connectivity,KIND=8)
    CALL PutArray(iFile,'Int64','offsets',1,offset)
    offset = offset + 8 + 8*SIZE(offsets,KIND=8)
    CALL PutArray(iFile,'UInt8','types',1,offset)
    offset = offset + 8 + SIZE(cellTypes,KIND=8)
    CALL Put(iFile,'</Cells>')

    location = 'CellData'
    IF ( onPoints ) location = 'PointData'
    CALL Put(iFile,'<'//TRIM(location)//'>')
    DO i = 1,nFields
      CALL PutArray(iFile,'Float64',TRIM(fields(i)%name),fields(i)%nComp, &
                    offset)
      offset = offset + 8 + 8*SIZE(fields(i)%values,KIND=8)
    END DO ! i
    CALL Put(iFile,'</'//TRIM(location)//'>')

    CALL Put(iFile,'</Piece>')
    CALL Put(iFile,'</UnstructuredGrid>')

    WRITE(iFile) '<AppendedData encoding="raw">'//NEWLINE//'_'
    WRITE(iFile) 8*SIZE(points,KIND=8),points
    WRITE(iFile) 8*SIZE(connectivity,KIND=8),connectivity
    WRITE(iFile) 8*SIZE(offsets,KIND=8),offsets
    WRITE(iFile) SIZE(cellTypes,KIND=8),cellTypes
    DO i = 1,nFields
      WRITE(iFile) 8*SIZE(fields(i)%values,KIND=8),fields(i)%values
    END DO ! i
    WRITE(iFile) NEWLINE//'</AppendedData>'//NEWLINE
    CALL Put(iFile,'</VTKFile>')

    CLOSE(iFile)

  END SUBROUTINE WriteUnstructuredVTU

! ******************************************************************************
!
! Purpose: Write the multiblock file of a step: a block of the fluid regions,
!   and with particles a block of the particles and one of the bins.
!
! Input:
!   fileName       File to write (post/root/collection_<step>.vtm)
!   regionFiles    The fluid files, relative to the .vtm
!   particleFile   The particle file ('' for none)
!   binFile        The bin file ('' for none)
!
! ******************************************************************************

  SUBROUTINE WriteStepVTM(fileName,regionFiles,particleFile,binFile)

    CHARACTER(*), INTENT(IN) :: fileName,particleFile,binFile
    CHARACTER(*), INTENT(IN) :: regionFiles(:)

    INTEGER :: errorFlag,i,iFile

    OPEN(NEWUNIT=iFile,FILE=fileName,FORM='FORMATTED',STATUS='REPLACE', &
         IOSTAT=errorFlag)
    IF ( errorFlag /= 0 ) THEN
      WRITE(*,'(A,A)') 'rflupost: cannot write ',TRIM(fileName)
      STOP 1
    END IF ! errorFlag

    WRITE(iFile,'(A)') '<?xml version="1.0"?>'
    WRITE(iFile,'(A)') '<VTKFile type="vtkMultiBlockDataSet" version="1.0">'
    WRITE(iFile,'(A)') '<vtkMultiBlockDataSet>'
    WRITE(iFile,'(A)') '  <Block index="0" name="fluid">'
    DO i = 1,SIZE(regionFiles)
      WRITE(iFile,'(A,I0,A,I0,A,A,A)') '    <DataSet index="',i-1, &
            '" name="region ',i,'" file="',TRIM(regionFiles(i)),'"/>'
    END DO ! i
    WRITE(iFile,'(A)') '  </Block>'
    IF ( LEN_TRIM(particleFile) > 0 ) THEN
      WRITE(iFile,'(A)') '  <Block index="1" name="particles">'
      WRITE(iFile,'(A,A,A)') '    <DataSet index="0" file="', &
            TRIM(particleFile),'"/>'
      WRITE(iFile,'(A)') '  </Block>'
    END IF ! particleFile
    IF ( LEN_TRIM(binFile) > 0 ) THEN
      WRITE(iFile,'(A)') '  <Block index="2" name="bins">'
      WRITE(iFile,'(A,A,A)') '    <DataSet index="0" file="',TRIM(binFile),'"/>'
      WRITE(iFile,'(A)') '  </Block>'
    END IF ! binFile
    WRITE(iFile,'(A)') '</vtkMultiBlockDataSet>'
    WRITE(iFile,'(A)') '</VTKFile>'

    CLOSE(iFile)

  END SUBROUTINE WriteStepVTM

! ******************************************************************************
!
! Purpose: Write the series file: one entry per step, with its time.
!
! Input:
!   fileName    File to write (post/collection.pvd)
!   times       Time of each step
!   stepFiles   The .vtm of each step, relative to the .pvd
!
! ******************************************************************************

  SUBROUTINE WriteSeriesPVD(fileName,times,stepFiles)

    CHARACTER(*), INTENT(IN) :: fileName
    REAL(8), INTENT(IN) :: times(:)
    CHARACTER(*), INTENT(IN) :: stepFiles(:)

    INTEGER :: errorFlag,i,iFile

    OPEN(NEWUNIT=iFile,FILE=fileName,FORM='FORMATTED',STATUS='REPLACE', &
         IOSTAT=errorFlag)
    IF ( errorFlag /= 0 ) THEN
      WRITE(*,'(A,A)') 'rflupost: cannot write ',TRIM(fileName)
      STOP 1
    END IF ! errorFlag

    WRITE(iFile,'(A)') '<?xml version="1.0"?>'
    WRITE(iFile,'(A)') '<VTKFile type="Collection" version="1.0">'
    WRITE(iFile,'(A)') '<Collection>'
    DO i = 1,SIZE(times)
      WRITE(iFile,'(A,ES24.16E3,A,A,A)') '  <DataSet timestep="',times(i), &
            '" file="',TRIM(stepFiles(i)),'"/>'
    END DO ! i
    WRITE(iFile,'(A)') '</Collection>'
    WRITE(iFile,'(A)') '</VTKFile>'

    CLOSE(iFile)

  END SUBROUTINE WriteSeriesPVD

! ******************************************************************************
! Helpers: a line of XML, and the header of an appended array
! ******************************************************************************

  SUBROUTINE Put(iFile,text)

    INTEGER, INTENT(IN) :: iFile
    CHARACTER(*), INTENT(IN) :: text

    WRITE(iFile) text//NEWLINE

  END SUBROUTINE Put

  SUBROUTINE PutArray(iFile,dataType,name,nComp,offset)

    INTEGER, INTENT(IN) :: iFile,nComp
    CHARACTER(*), INTENT(IN) :: dataType,name
    INTEGER(8), INTENT(IN) :: offset

    CHARACTER(256) :: line

    WRITE(line,'(A,A,A,A,A,I0,A,I0,A)') '<DataArray type="',dataType, &
          '" Name="',name,'" NumberOfComponents="',nComp, &
          '" format="appended" offset="',offset,'"/>'
    CALL Put(iFile,TRIM(line))

  END SUBROUTINE PutArray

END MODULE RFLU_ModPostWrite
