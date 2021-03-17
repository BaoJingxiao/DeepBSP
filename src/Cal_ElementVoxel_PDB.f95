!基于受体和配体的PDB文件（去除非极性氢原子），计算ElementVoxel描述符
!使用AutoDock4中定义的原子半径，Na+和K+的原子半径均设为1.200
!受体8通道分别为：H、C、N、O、P、S、M(金属，包括：Na,NA,Mg,MG,K,Ca,CA,Mn,MN,Fe,FE,Zn,ZN)、A(所有元素)
!配体8通道分别为：H、C、N、O、P、S、X(卤素，包括：F,Cl,CL,Br,BR,I)、A(所有元素)
!将数值范围为[0,1]的描述符浮点数据，通过减去0.5后乘以2^15，转换到范围为[-2^14,2^14]的整型数据
!使用NetCDF4格式存储数据并启用压缩，压缩等级设为6

!读取计算参数文件
module ReadConfigure
    !文件名称与参数设置文件样式：(不包括行首的!)
    !#configures for calculate the descriptor of complex
    !#file path (should end with /):
    !/media/Data/Work/Test_CalDescriptor/2ewb_GB2_Diel1.0_noRst/
    !#receptor PDB file name:
    !REC.pdb
    !#ligand PDB file name:
    !LIG.pdb
    !#output data file name:
    !Descriptor.dat.nc
    !#element VDW radii file:
    !/home/paul/Fortran/Work/Cal_ElementVoxel_PDB/dat/ElementVDWRadii.dat
    !#grid size of X:
    !24.0
    !#grid size of Y:
    !24.0
    !#grid size of Z:
    !24.0
    !#grid offset:
    !1.0
    !#VDW distance cutoff:
    !10.0
    !#end
    !申明变量
    implicit none
    character*1024 :: RecPDBFile
    character*1024 :: LigPDBFile
    character*1024 :: OutFile
    character*1024 :: RadiiFile
    real*8 :: GridSizeX,GridSizeY,GridSizeZ,GridOffset
    real*8 :: VDWCutoff
    !子程序
    contains
        !读取计算参数文件
        subroutine ReadConfigureFile(ConfFile)
            !申明变量
            implicit none
            character*1024 :: ConfFile
            integer*4 :: ConfFileID
            character*1024 :: FilePath
            !读取参数
            ConfFileID=10
            open(ConfFileID,file=trim(ConfFile),status="old",action="read")
            read(ConfFileID,*)
            read(ConfFileID,*);read(ConfFileID,'(a)')FilePath
            read(ConfFileID,*);read(ConfFileID,'(a)')RecPDBFile
            read(ConfFileID,*);read(ConfFileID,'(a)')LigPDBFile
            read(ConfFileID,*);read(ConfFileID,'(a)')OutFile
            read(ConfFileID,*);read(ConfFileID,'(a)')RadiiFile
            read(ConfFileID,*);read(ConfFileID,'(e8.4)')GridSizeX
            read(ConfFileID,*);read(ConfFileID,'(e8.4)')GridSizeY
            read(ConfFileID,*);read(ConfFileID,'(e8.4)')GridSizeZ
            read(ConfFileID,*);read(ConfFileID,'(e8.4)')GridOffset
            read(ConfFileID,*);read(ConfFileID,'(e8.4)')VDWCutoff
            close(ConfFileID)
            RecPDBFile=trim(FilePath)//trim(RecPDBFile)
            LigPDBFile=trim(FilePath)//trim(LigPDBFile)
            OutFile=trim(FilePath)//trim(OutFile)
        end subroutine ReadConfigureFile
end module ReadConfigure

!读取元素的范德华参数文件
module ReadVDW
    !调用模块
    use ReadConfigure
    !申明变量
    implicit none
    integer*4 :: ElementNum
    character*2,allocatable,dimension(:) :: ElementName
    integer*4,allocatable,dimension(:) :: ElementIndex
    real*8,allocatable,dimension(:) :: ElementVDWRadii
    !子程序
    contains
        !读取元素的范德华参数文件
        subroutine ReadVDWRadiiFile()
            !申明变量
            implicit none
            integer*4 :: FileID
            character*20 :: Line
            integer*4 :: FlagEOF
            integer*4 :: Temp1
            !确定元素类型总数
            ElementNum=0
            FileID=11
            open(FileID,file=trim(RadiiFile),status="old",action="read")
            do
                read(FileID,'(a11)',iostat=FlagEOF)Line
                if (FlagEOF/=0) exit
                if (trim(Line)/="") then
                    ElementNum=ElementNum+1
                end if
            end do
            close(FileID)
            !初始化数组
            allocate(ElementName(ElementNum));ElementName=""
            allocate(ElementIndex(ElementNum));ElementIndex=0
            allocate(ElementVDWRadii(ElementNum));ElementVDWRadii=0.0D0
            !读取元素名称、序号和范德华半径
            Temp1=0
            open(FileID,file=trim(RadiiFile),status="old",action="read")
            do
                read(FileID,'(a11)',iostat=FlagEOF)Line
                if (FlagEOF/=0) exit
                if (trim(Line)/="") then
                    Temp1=Temp1+1
                    write(ElementName(Temp1),'(a2)')trim(Line(1:2))
                    read(Line(4:5),'(i2)')ElementIndex(Temp1)
                    read(Line(7:11),'(f5.3)')ElementVDWRadii(Temp1)
                end if
            end do
            close(FileID)
        end subroutine ReadVDWRadiiFile
end module ReadVDW

!读取受体配体PDB文件，确定受体配体的原子数量，原子元素名称、原子范德华半径的12次方和原子坐标(以配体几何中心为原点)
module RecLig
    !调用模块
    use ReadConfigure
    use ReadVDW
    !申明变量
    integer*4 :: AtomNum_Rec,AtomNum_Lig
    character*2,allocatable,dimension(:) :: AtomElement_Rec,AtomElement_Lig
    real*8,allocatable,dimension(:) :: AtomVDWRadii12_Rec,AtomVDWRadii12_Lig
    real*8,allocatable,dimension(:) :: Crd_Rec,Crd_Lig
    !子程序
    contains
        !读取受体配体PDB文件，确定受体配体的原子数量，原子元素名称、原子范德华半径的12次方和原子坐标(以配体几何中心为原点)
        subroutine ReadRecLigPDB()
            !申明变量
            implicit none
            integer*4 :: FileID
            integer*4 :: FlagEOF
            character*80 :: Line
            real*8 :: LigCrdCenter(3)
            integer*4 :: Temp1
            integer*4 :: Flag1
            integer*4 :: n
            !读取受体PDB文件并确定原子数量
            AtomNum_Rec=0
            FileID=11
            open(FileID,file=trim(RecPDBFile),status="old",action="read")
            do
                read(FileID,'(a80)',iostat=FlagEOF)Line
                if (FlagEOF/=0) exit
                if ((Line(1:4)=="ATOM").or.(Line(1:6)=="HETATM")) then
                    AtomNum_Rec=AtomNum_Rec+1
                end if
            end do
            close(FileID)
            !读取配体PDB文件并确定原子数量
            AtomNum_Lig=0
            FileID=11
            open(FileID,file=trim(LigPDBFile),status="old",action="read")
            do
                read(FileID,'(a80)',iostat=FlagEOF)Line
                if (FlagEOF/=0) exit
                if ((Line(1:4)=="ATOM").or.(Line(1:6)=="HETATM")) then
                    AtomNum_Lig=AtomNum_Lig+1
                end if
            end do
            close(FileID)
            !初始化数组
            allocate(AtomElement_Rec(AtomNum_Rec));AtomElement_Rec=""
            allocate(AtomElement_Lig(AtomNum_Lig));AtomElement_Lig=""
            allocate(AtomVDWRadii12_Rec(AtomNum_Rec));AtomVDWRadii12_Rec=0.0D0
            allocate(AtomVDWRadii12_Lig(AtomNum_Lig));AtomVDWRadii12_Lig=0.0D0
            allocate(Crd_Rec(AtomNum_Rec*3));Crd_Rec=0.0D0
            allocate(Crd_Lig(AtomNum_Lig*3));Crd_Lig=0.0D0
            !读取受体PDB文件并确定各原子的元素名称、坐标和范德华半径的12次方
            Temp1=0
            FileID=11
            open(FileID,file=trim(RecPDBFile),status="old",action="read")
            do
                read(FileID,'(a80)',iostat=FlagEOF)Line
                if (FlagEOF/=0) exit
                if ((Line(1:4)=="ATOM").or.(Line(1:6)=="HETATM")) then
                    Temp1=Temp1+1
                    write(AtomElement_Rec(Temp1),'(a2)')trim(adjustl(Line(77:78)))
                    read(Line(31:54),'(3(f8.3))')Crd_Rec((Temp1*3-2):(Temp1*3))
                    Flag1=0
                    do n=1,ElementNum,1
                        if (trim(adjustl(AtomElement_Rec(Temp1)))==trim(adjustl(ElementName(n)))) then
                            AtomVDWRadii12_Rec(Temp1)=ElementVDWRadii(n)**12
                            Flag1=1
                            exit
                        end if
                    end do
                    if (Flag1==0) then
                        write(*,'("Warning: Ignore receptor atom with element name: ",a2)')AtomElement_Rec(Temp1)
                    end if
                end if
            end do
            close(FileID)
            !读取配体PDB文件并确定各原子的元素名称、坐标和范德华半径的12次方
            Temp1=0
            FileID=11
            open(FileID,file=trim(LigPDBFile),status="old",action="read")
            do
                read(FileID,'(a80)',iostat=FlagEOF)Line
                if (FlagEOF/=0) exit
                if ((Line(1:4)=="ATOM").or.(Line(1:6)=="HETATM")) then
                    Temp1=Temp1+1
                    write(AtomElement_Lig(Temp1),'(a2)')trim(adjustl(Line(77:78)))
                    read(Line(31:54),'(3(f8.3))')Crd_Lig((Temp1*3-2):(Temp1*3))
                    Flag1=0
                    do n=1,ElementNum,1
                        if (trim(adjustl(AtomElement_Lig(Temp1)))==trim(adjustl(ElementName(n)))) then
                            AtomVDWRadii12_Lig(Temp1)=ElementVDWRadii(n)**12
                            Flag1=1
                            exit
                        end if
                    end do
                    if (Flag1==0) then
                        write(*,'("Warning: Ignore ligand atom with element name: ",a2)')AtomElement_Lig(Temp1)
                    end if
                end if
            end do
            close(FileID)
            !计算配体原子的几何中心
            LigCrdCenter(1)=sum(Crd_Lig(1::3))/AtomNum_Lig
            LigCrdCenter(2)=sum(Crd_Lig(2::3))/AtomNum_Lig
            LigCrdCenter(3)=sum(Crd_Lig(3::3))/AtomNum_Lig
            !以配体原子的几何中心为原点，计算受体和配体各原子的坐标
            Crd_Rec(1::3)=Crd_Rec(1::3)-LigCrdCenter(1)
            Crd_Rec(2::3)=Crd_Rec(2::3)-LigCrdCenter(2)
            Crd_Rec(3::3)=Crd_Rec(3::3)-LigCrdCenter(3)
            Crd_Lig(1::3)=Crd_Lig(1::3)-LigCrdCenter(1)
            Crd_Lig(2::3)=Crd_Lig(2::3)-LigCrdCenter(2)
            Crd_Lig(3::3)=Crd_Lig(3::3)-LigCrdCenter(3)
        end subroutine ReadRecLigPDB
end module RecLig

!计算格点
module CalGrid
    !调用模块
    use ReadConfigure
    use RecLig
    !申明变量
    integer*4 :: GridNumX,GridNumY,GridNumZ
    real*8,allocatable,dimension(:) :: GridCrdX,GridCrdY,GridCrdZ
    real*8,allocatable,dimension(:,:,:) :: Voxel_H_Rec,Voxel_H_Lig            !受体配体中的H
    real*8,allocatable,dimension(:,:,:) :: Voxel_C_Rec,Voxel_C_Lig            !受体配体中的C
    real*8,allocatable,dimension(:,:,:) :: Voxel_N_Rec,Voxel_N_Lig            !受体配体中的N
    real*8,allocatable,dimension(:,:,:) :: Voxel_O_Rec,Voxel_O_Lig            !受体配体中的O
    real*8,allocatable,dimension(:,:,:) :: Voxel_P_Rec,Voxel_P_Lig            !受体配体中的P
    real*8,allocatable,dimension(:,:,:) :: Voxel_S_Rec,Voxel_S_Lig            !受体配体中的S
    real*8,allocatable,dimension(:,:,:) :: Voxel_M_Rec                        !受体中的金属:Na,Mg,K,Ca,Mn,Fe,Zn
    real*8,allocatable,dimension(:,:,:) :: Voxel_X_Lig                        !配体中的卤素:F,Cl,Br,I
    real*8,allocatable,dimension(:,:,:) :: Voxel_All_Rec,Voxel_All_Lig        !受体配体中的所有元素
    !子程序
    contains
        !计算格点的坐标
        subroutine GetGridCrd()
            !申明变量
            implicit none
            integer*4 :: n
            !确定XYZ各方向上格点的数量
            GridNumX=(floor((GridSizeX-GridOffset)/2.0D0/GridOffset)+1)*2
            GridNumY=(floor((GridSizeY-GridOffset)/2.0D0/GridOffset)+1)*2
            GridNumZ=(floor((GridSizeZ-GridOffset)/2.0D0/GridOffset)+1)*2
            !初始化数组
            allocate(GridCrdX(GridNumX));GridCrdX=0.0D0
            allocate(GridCrdY(GridNumY));GridCrdY=0.0D0
            allocate(GridCrdZ(GridNumZ));GridCrdZ=0.0D0
            !确定格点的坐标
            do n=1,GridNumX,1
                GridCrdX(n)=(n-GridNumX/2-0.5D0)*GridOffset
            end do
            do n=1,GridNumY,1
                GridCrdY(n)=(n-GridNumY/2-0.5D0)*GridOffset
            end do
            do n=1,GridNumZ,1
                GridCrdZ(n)=(n-GridNumZ/2-0.5D0)*GridOffset
            end do
        end subroutine GetGridCrd
        !计算受体和配体在各格点上的体素描述符
        subroutine CalGridVoxel()
            !申明变量
            implicit none
            real*8 :: MaxDist2
            real*8 :: TempDist2,TempDist12
            real*8 :: TempVoxel
            integer*4 :: n,x,y,z
            !初始化数组
            allocate(Voxel_H_Rec(GridNumZ,GridNumY,GridNumX));Voxel_H_Rec=0.0D0
            allocate(Voxel_H_Lig(GridNumZ,GridNumY,GridNumX));Voxel_H_Lig=0.0D0
            allocate(Voxel_C_Rec(GridNumZ,GridNumY,GridNumX));Voxel_C_Rec=0.0D0
            allocate(Voxel_C_Lig(GridNumZ,GridNumY,GridNumX));Voxel_C_Lig=0.0D0
            allocate(Voxel_N_Rec(GridNumZ,GridNumY,GridNumX));Voxel_N_Rec=0.0D0
            allocate(Voxel_N_Lig(GridNumZ,GridNumY,GridNumX));Voxel_N_Lig=0.0D0
            allocate(Voxel_O_Rec(GridNumZ,GridNumY,GridNumX));Voxel_O_Rec=0.0D0
            allocate(Voxel_O_Lig(GridNumZ,GridNumY,GridNumX));Voxel_O_Lig=0.0D0
            allocate(Voxel_P_Rec(GridNumZ,GridNumY,GridNumX));Voxel_P_Rec=0.0D0
            allocate(Voxel_P_Lig(GridNumZ,GridNumY,GridNumX));Voxel_P_Lig=0.0D0
            allocate(Voxel_S_Rec(GridNumZ,GridNumY,GridNumX));Voxel_S_Rec=0.0D0
            allocate(Voxel_S_Lig(GridNumZ,GridNumY,GridNumX));Voxel_S_Lig=0.0D0
            allocate(Voxel_M_Rec(GridNumZ,GridNumY,GridNumX));Voxel_M_Rec=0.0D0
            allocate(Voxel_X_Lig(GridNumZ,GridNumY,GridNumX));Voxel_X_Lig=0.0D0
            allocate(Voxel_All_Rec(GridNumZ,GridNumY,GridNumX));Voxel_All_Rec=0.0D0
            allocate(Voxel_All_Lig(GridNumZ,GridNumY,GridNumX));Voxel_All_Lig=0.0D0
            !当原子到格点的距离的平方大于MaxDist2时，忽略该原子格点的作用
            MaxDist2=VDWCutoff**2
            !计算受体在各格点上的体素描述符
            do n=1,AtomNum_Rec,1
                do x=1,GridNumX,1
                    do y=1,GridNumY,1
                        do z=1,GridNumZ,1
                            TempDist2=(Crd_Rec(n*3-2)-GridCrdX(x))**2+(Crd_Rec(n*3-1)-GridCrdY(y))**2+(Crd_Rec(n*3)-GridCrdZ(z))**2
                            if (TempDist2<1.0D-4) then
                                TempDist2=1.0D-4
                            end if
                            if (TempDist2<=MaxDist2) then
                                TempDist12=TempDist2**6
                                TempVoxel=1.0D0-exp(-AtomVDWRadii12_Rec(n)/TempDist12)
                                select case (trim(adjustl(AtomElement_Rec(n))))
                                    case ("H")
                                        if (Voxel_H_Rec(z,y,x)<TempVoxel) Voxel_H_Rec(z,y,x)=TempVoxel
                                    case ("C")
                                        if (Voxel_C_Rec(z,y,x)<TempVoxel) Voxel_C_Rec(z,y,x)=TempVoxel
                                    case ("N")
                                        if (Voxel_N_Rec(z,y,x)<TempVoxel) Voxel_N_Rec(z,y,x)=TempVoxel
                                    case ("O")
                                        if (Voxel_O_Rec(z,y,x)<TempVoxel) Voxel_O_Rec(z,y,x)=TempVoxel
                                    case ("P")
                                        if (Voxel_P_Rec(z,y,x)<TempVoxel) Voxel_P_Rec(z,y,x)=TempVoxel
                                    case ("S")
                                        if (Voxel_S_Rec(z,y,x)<TempVoxel) Voxel_S_Rec(z,y,x)=TempVoxel
                                    case ("Na","NA","Mg","MG","K","Ca","CA","Mn","MN","Fe","FE","Zn","ZN")
                                        if (Voxel_M_Rec(z,y,x)<TempVoxel) Voxel_M_Rec(z,y,x)=TempVoxel
                                end select
                                if (Voxel_All_Rec(z,y,x)<TempVoxel) Voxel_All_Rec(z,y,x)=TempVoxel
                            end if
                        end do
                    end do
                end do
            end do
            !计算配体在各格点上的体素描述符
            do n=1,AtomNum_Lig,1
                do x=1,GridNumX,1
                    do y=1,GridNumY,1
                        do z=1,GridNumZ,1
                            TempDist2=(Crd_Lig(n*3-2)-GridCrdX(x))**2+(Crd_Lig(n*3-1)-GridCrdY(y))**2+(Crd_Lig(n*3)-GridCrdZ(z))**2
                            if (TempDist2<1.0D-4) then
                                TempDist2=1.0D-4
                            end if
                            if (TempDist2<=MaxDist2) then
                                TempDist12=TempDist2**6
                                TempVoxel=1.0D0-exp(-AtomVDWRadii12_Lig(n)/TempDist12)
                                select case (trim(adjustl(AtomElement_Lig(n))))
                                    case ("H")
                                        if (Voxel_H_Lig(z,y,x)<TempVoxel) Voxel_H_Lig(z,y,x)=TempVoxel
                                    case ("C")
                                        if (Voxel_C_Lig(z,y,x)<TempVoxel) Voxel_C_Lig(z,y,x)=TempVoxel
                                    case ("N")
                                        if (Voxel_N_Lig(z,y,x)<TempVoxel) Voxel_N_Lig(z,y,x)=TempVoxel
                                    case ("O")
                                        if (Voxel_O_Lig(z,y,x)<TempVoxel) Voxel_O_Lig(z,y,x)=TempVoxel
                                    case ("P")
                                        if (Voxel_P_Lig(z,y,x)<TempVoxel) Voxel_P_Lig(z,y,x)=TempVoxel
                                    case ("S")
                                        if (Voxel_S_Lig(z,y,x)<TempVoxel) Voxel_S_Lig(z,y,x)=TempVoxel
                                    case ("F","Cl","CL","Br","BR","I")
                                        if (Voxel_X_Lig(z,y,x)<TempVoxel) Voxel_X_Lig(z,y,x)=TempVoxel
                                end select
                                if (Voxel_All_Lig(z,y,x)<TempVoxel) Voxel_All_Lig(z,y,x)=TempVoxel
                            end if
                        end do
                    end do
                end do
            end do
        end subroutine CalGridVoxel
end module CalGrid

!主程序
program main
    !调用模块
    use netcdf
    use ReadConfigure
    use ReadVDW
    use RecLig
    use CalGrid
    !申明变量
    implicit none
    character*1024 :: ConfFile
    !读取计算参数文件
    call getarg(1,ConfFile)
    !ConfFile="/media/Data/Dynamics/Test_CalDescriptor/Descriptor.config"
    call ReadConfigureFile(ConfFile)
    !读取元素的范德华半径参数文件
    call ReadVDWRadiiFile()
    !读取受体配体PDB文件，确定受体配体的原子数量，原子元素名称、原子范德华半径的12次方和原子坐标(以配体几何中心为原点)
    call ReadRecLigPDB()
    !计算格点的坐标
    call GetGridCrd()
    !计算受体和配体在各格点上的体素描述符
    call CalGridVoxel()
    !将计算的描述符输出至文件
    call GetOutFile()
end

!将计算的描述符输出至文件
subroutine GetOutFile()
    !调用模块
    use netcdf
    use ReadConfigure
    use CalGrid
    !申明变量
    implicit none
    integer*4 :: OutFileID
    integer*4 :: Flag_NF90
    integer*4 :: DimID_DS
    integer*4 :: DimID_GridX,DimID_GridY,DimID_GridZ
    integer*4 :: VarID_Data
    integer*4 :: Start_Data(4)
    integer*4 :: Count_Data(4)
    !定义输出文件
    Flag_NF90=nf90_create(trim(OutFile),NF90_CLOBBER+NF90_NETCDF4,OutFileID)
    Flag_NF90=nf90_def_dim(OutFileID,"DS",NF90_UNLIMITED,DimID_DS)
    Flag_NF90=nf90_def_dim(OutFileID,"GridX",GridNumX,DimID_GridX)
    Flag_NF90=nf90_def_dim(OutFileID,"GridY",GridNumY,DimID_GridY)
    Flag_NF90=nf90_def_dim(OutFileID,"GridZ",GridNumZ,DimID_GridZ)
    Flag_NF90=nf90_def_var(OutFileID,"Data",NF90_SHORT,(/DimID_GridZ,DimID_GridY,DimID_GridX,DimID_DS/),VarID_Data,deflate_level=6)
    Flag_NF90=nf90_put_att(OutFileID,NF90_GLOBAL,"Title","Element Voxel Data")
    Flag_NF90=nf90_put_att(OutFileID,NF90_GLOBAL,"GridNumX",GridNumX)
    Flag_NF90=nf90_put_att(OutFileID,NF90_GLOBAL,"GridNumY",GridNumY)
    Flag_NF90=nf90_put_att(OutFileID,NF90_GLOBAL,"GridNumZ",GridNumZ)
    Flag_NF90=nf90_put_att(OutFileID,NF90_GLOBAL,"GridOffset",GridOffset)
    Flag_NF90=nf90_enddef(OutFileID)
    !输出描述符至文件
    Start_Data=(/1,1,1,0/)
    Count_Data=(/GridNumZ,GridNumY,GridNumX,1/)
    !输出受体在各格点上的体素描述符
    Start_Data(4)=Start_Data(4)+1
    Flag_NF90=nf90_put_var(OutFileID,VarID_Data,int((Voxel_H_Rec-0.5D0)*(2**15)),Start_Data,Count_Data)
    Start_Data(4)=Start_Data(4)+1
    Flag_NF90=nf90_put_var(OutFileID,VarID_Data,int((Voxel_C_Rec-0.5D0)*(2**15)),Start_Data,Count_Data)
    Start_Data(4)=Start_Data(4)+1
    Flag_NF90=nf90_put_var(OutFileID,VarID_Data,int((Voxel_N_Rec-0.5D0)*(2**15)),Start_Data,Count_Data)
    Start_Data(4)=Start_Data(4)+1
    Flag_NF90=nf90_put_var(OutFileID,VarID_Data,int((Voxel_O_Rec-0.5D0)*(2**15)),Start_Data,Count_Data)
    Start_Data(4)=Start_Data(4)+1
    Flag_NF90=nf90_put_var(OutFileID,VarID_Data,int((Voxel_P_Rec-0.5D0)*(2**15)),Start_Data,Count_Data)
    Start_Data(4)=Start_Data(4)+1
    Flag_NF90=nf90_put_var(OutFileID,VarID_Data,int((Voxel_S_Rec-0.5D0)*(2**15)),Start_Data,Count_Data)
    Start_Data(4)=Start_Data(4)+1
    Flag_NF90=nf90_put_var(OutFileID,VarID_Data,int((Voxel_M_Rec-0.5D0)*(2**15)),Start_Data,Count_Data)
    Start_Data(4)=Start_Data(4)+1
    Flag_NF90=nf90_put_var(OutFileID,VarID_Data,int((Voxel_All_Rec-0.5D0)*(2**15)),Start_Data,Count_Data)
    !输出配体在各格点上的体素描述符
    Start_Data(4)=Start_Data(4)+1
    Flag_NF90=nf90_put_var(OutFileID,VarID_Data,int((Voxel_H_Lig-0.5D0)*(2**15)),Start_Data,Count_Data)
    Start_Data(4)=Start_Data(4)+1
    Flag_NF90=nf90_put_var(OutFileID,VarID_Data,int((Voxel_C_Lig-0.5D0)*(2**15)),Start_Data,Count_Data)
    Start_Data(4)=Start_Data(4)+1
    Flag_NF90=nf90_put_var(OutFileID,VarID_Data,int((Voxel_N_Lig-0.5D0)*(2**15)),Start_Data,Count_Data)
    Start_Data(4)=Start_Data(4)+1
    Flag_NF90=nf90_put_var(OutFileID,VarID_Data,int((Voxel_O_Lig-0.5D0)*(2**15)),Start_Data,Count_Data)
    Start_Data(4)=Start_Data(4)+1
    Flag_NF90=nf90_put_var(OutFileID,VarID_Data,int((Voxel_P_Lig-0.5D0)*(2**15)),Start_Data,Count_Data)
    Start_Data(4)=Start_Data(4)+1
    Flag_NF90=nf90_put_var(OutFileID,VarID_Data,int((Voxel_S_Lig-0.5D0)*(2**15)),Start_Data,Count_Data)
    Start_Data(4)=Start_Data(4)+1
    Flag_NF90=nf90_put_var(OutFileID,VarID_Data,int((Voxel_X_Lig-0.5D0)*(2**15)),Start_Data,Count_Data)
    Start_Data(4)=Start_Data(4)+1
    Flag_NF90=nf90_put_var(OutFileID,VarID_Data,int((Voxel_All_Lig-0.5D0)*(2**15)),Start_Data,Count_Data)
    !关闭输出文件
    Flag_NF90=nf90_close(OutFileID)
end subroutine GetOutFile
