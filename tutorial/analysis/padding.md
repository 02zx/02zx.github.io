# 用于补全RASPA中产生的pdb文件以供VMD可视化的程序

```fortran
program read_pdb
    implicit none
    character(80) :: line
    character(6) :: record
    integer :: atom_num, loop_n, max_n, frame,i
    character(4) :: atom_name
    character(3) :: resname
    real :: x, y, z
    
    integer :: iunit, ierr
    character(256) :: pdb_file
    
    ! Prompt the user for the PDB file name
    call getarg(1,pdb_file)
    call getarg(2, record)
    read(record, *, iostat=ierr) max_n
    
    ! Open the PDB file
    open(newunit=iunit, file=trim(pdb_file), status='old', action='read',iostat=ierr)
    if (ierr /= 0) then
        write(*,*) "Error opening the PDB file."
        stop
    end if
    frame=0
    ! Read each line of the PDB file
    do 
        read(iunit, '(A)',iostat=ierr) line
        if (ierr /= 0) exit
        if (index(line, "MODEL") /= 0) then
           write(*,*) "ITEM: TIMESTEP"
           write(*,'(I4)') frame
           frame=frame+1
           write(*,*) "ITEM: NUMBER OF ATOMS"
           write(*,'(I4)') max_n
           
           cycle  
        end if

        if (index(line, "CRYST1") /= 0) then  
           write(*,*) "ITEM: BOX BOUNDS ss ss ss"
           read(line, *) record, x,y,z,record
           write(*,'(2F8.3)') 0.000,x
           write(*,'(2F8.3)') 0.000,y
           write(*,'(2F8.3)') 0.000,z
           write(*,*) "ITEM: ATOMS id type xu yu zu vx"
           loop_n=1  
           cycle  
        end if

        if (index(line, "ENDMDL") /= 0) then  
           do i = loop_n, max_n  
              write(*,'(I3,I2,3F8.3,I3)')  i, 1, 0.000, 0.000, 0.000, 1
           end do  
           cycle 
        end if
        loop_n=loop_n+1
        read(line, *) record,atom_num,atom_name,resname, x,y,z,record,record,record
        write(*,'(I3,I2,3F8.3,I3)')  atom_num, 1, x, y, z, 0
    end do
    
    ! Close the PDB file
    close(iunit)
    
end program read_pdb

```

## 编译

```bash
gfortran padding.f90 -o padding.x
```

## 使用

```bash
#获得文件(traj.pdb)中的最大分子数
 N=`grep 'ATOM' traj.pdb |awk '{print $2}'|sort -g|tail -n 1`

#补全文件
./padding.x traj.pdb $N >out.lammpstrj
```

文件的输出格式为lammps, 其中vx 0为原始数据, vx 1为补全所添加的点
