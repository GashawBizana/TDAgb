import numpy as np
import subprocess

def generate_lammps_input(latP, xl, yl, zl, xh, yh, zh, gl, gu):
    input_file = f"""
# LAMMPS Input File for Grain Boundaries 
# Mark Tschopp, Dec2009 
# This file will generate a single Sigma5(310) STGB

# ---------- Initialize Simulation --------------------- 
clear 
units metal 
dimension 3 
boundary p p p 
atom_style atomic 

# ---------- Create Atomistic Structure --------------------- 

variable latP equal {latP}

lattice fcc ${{latP}}

variable        xl equal {xl}
variable        yl equal {yl}
variable        zl equal {zl}

variable        xh equal {xh}
variable        yh equal {yh}
variable        zh equal {zh}

variable        gl equal {gl}
variable        gu equal {gu}

region whole block 0.000000 ${{xh}} 0 ${{yh}} ${{zl}} ${{zh}} units box
create_box 2 whole 

variable        z_min equal zlo
variable        z_max equal zhi

region lower block INF INF INF  INF ${{z_min}} ${{gl}} units box
lattice fcc {latP} orient x  0  3 1 orient y  -1  0  0 orient z  0  -1  3
create_atoms 2 region lower 

region upper block INF INF INF INF ${{gu}} ${{z_max}} units box
lattice fcc {latP} orient x  0  3 1 orient y  -1  0  0 orient z  0  -1  3
create_atoms 2 region upper 

region middel block INF INF INF INF ${{gl}} ${{gu}} units box
lattice fcc {latP} orient x  0  9  13 orient y  -1  0  0 orient z  0  -13  9 
create_atoms 1 region middel 

group middel type 1 
group  outer type 2  

# ---------- Define Interatomic Potential --------------------- 
pair_style eam/alloy 
pair_coeff * * C:/Users/Zero/Documents/lammps_test/Al99.eam.alloy Al Al
neighbor 2.0 bin 
neigh_modify delay 10 check yes 
 
# ---------- Displace atoms and delete overlapping atoms --------------------- 
displace_atoms middel move 0 0 0 units lattice 
delete_atoms overlap 0.35 middel outer
 
# ---------- Define Settings --------------------- 
compute csym all centro/atom fcc
compute eng all pe/atom 
compute eatoms all reduce sum c_eng 

# ---------- Run Minimization --------------------- 
reset_timestep 0 
thermo 10 
thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms 
dump        1 all cfg 25 dump.sig5_minimization_*.cfg mass type xs ys zs c_csym c_eng fx fy fz
dump_modify     1 element Al Al
min_style cg 
minimize 1e-15 1e-15 5000 5000 
undump 1

# ---------- Run Minimization 2--------------------- 
# Now allow the box to expand/contract perpendicular to the grain boundary
reset_timestep 0 
thermo 10 
thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms 
fix 1 all box/relax y 0 vmax 0.001
min_style cg 
minimize 1e-15 1e-15 5000 5000 

# ---------- Calculate GB Energy --------------------- 
variable minimumenergy equal -3.360000
variable esum equal "v_minimumenergy * count(all)" 
variable xseng equal "c_eatoms - (v_minimumenergy * count(all))" 
variable gbarea equal "lx * lz * 2" 
variable gbe equal "(c_eatoms - (v_minimumenergy * count(all)))/v_gbarea" 
variable gbemJm2 equal ${{gbe}}*16021.7733 
variable gbernd equal round(${{gbemJm2}}) 
print "GB energy is ${{gbemJm2}} mJ/m^2" 
 
# ---------- Dump data into Data file ------------- 
reset_timestep 0 
dump        1 all cfg 10000 dump.al_sig5_310_*.cfg mass type xs ys zs c_csym c_eng fx fy fz
dump_modify     1 element Al Al
minimize 1e-15 1e-15 5000 5000
undump 1

write_data data_min.al_sig5_310_stgb

print "All done"
"""
    return input_file
def generate_and_wite_ori_files(LatP,VA,VB,gb_name,Temp):
    
    R_11=np.array([LatP/2,LatP/2,0])
    R_21=np.array([LatP/2,0,LatP/2])
    R_31=np.array([0,LatP/2,LatP/2])
    R_12=np.array([LatP/2,-LatP/2,0])
    R_22=np.array([LatP/2,0,-LatP/2])
    R_32=np.array([0,LatP/2,-LatP/2])
    R_13=np.array([-LatP/2,LatP/2,0])
    R_23=np.array([-LatP/2,0,LatP/2])
    R_33=np.array([0,-LatP/2,LatP/2])

    R_14=np.array([-LatP/2,-LatP/2,0])
    R_24=np.array([-LatP/2,0,-LatP/2])
    R_34=np.array([0,-LatP/2,-LatP/2])

    Rij=[R_11,R_21,R_31,R_12,R_22,R_32]
    
    XA,YA,ZA=np.array(VA[0:3]),np.array(VA[3:6]),np.array(VA[6:9])
    XB,YB,ZB=np.array(VB[0:3]),np.array(VB[3:6]),np.array(VB[6:9])
    oA=[]
    oB=[]

    for j in range(0,6):
        
        ujp=np.dot(Rij[j],XA.T)/np.linalg.norm(XA)
        vjp=np.dot(Rij[j],YA.T)/np.linalg.norm(YA)
        wjp=np.dot(Rij[j],ZA.T)/np.linalg.norm(ZA)
        
        oA.append([ujp,vjp,wjp])
        np.savetxt(gb_name+str(Temp)+"_RA.txt", oA, delimiter=' ',fmt='%.18f') 

    for j in range(0,6):
        
        ujp=np.dot(Rij[j],XB)/np.linalg.norm(XB)
        vjp=np.dot(Rij[j],YB)/np.linalg.norm(YB)
        wjp=np.dot(Rij[j],ZB)/np.linalg.norm(ZB)
        
        oB.append([ujp,vjp,wjp])
        np.savetxt(gb_name+str(Temp)+"_RB.txt", oB, delimiter=' ',fmt='%.18f') 

def gemerate_lammps_input_for_minimum(LatP,cellsize,min_z,max_z,cellorigin,VA,VB,gb_name,Temp):
    
    
    xl, yl, zl = 0, 0, -cellsize[2,2]
    xh, yh, zh = cellsize[0,0], cellsize[1,1], cellsize[2,2]
    gl,gu=cellorigin[2][0],cellsize[2,2]+cellorigin[2][0]
    lcut,ucut=min_z,max_z
    #lcutdel,ucutdel=15-np.abs(min_z)-0.1,0
    lcutdel,ucutdel=0.5,1
    overlapL=5#abs(abs(min_z))
    overlapU=5#abs(max_z)
    Folder='F:/CGBS_and_E/Computed_GBE_and_Structure_Al/Al_minEnergyGBstructs_2022/minEnergyGBstructs/'
    with open(gb_name+"_min.input", "w") as f:
        f.write('# Gashaw Bizana, NOV 2023. Script for creating GB\n')
        f.write('# ---------- Initialize Simulation --------------------- \n')
        f.write('clear\n')
        f.write('units metal\n')
        f.write('dimension 3 \n')
        f.write('boundary p p p\n') 
        f.write('atom_style atomic\n')
        f.write('# ---------- define variables for creating box ---------------------\n')
        f.write('variable latP equal 4.05\n')
        f.write('variable		xl equal '+str(xl)+'\n')
        f.write('variable		yl equal '+str(yl)+'\n')
        f.write('variable		zl equal '+str(zl)+'\n')
        f.write('variable		xh equal '+str(xh)+'\n')
        f.write('variable		yh equal '+str(yh)+'\n')
        f.write('variable		zh equal '+str(zh)+'\n')
        f.write('variable		gl equal '+str(gl)+'\n')
        f.write('variable		gu equal '+str(gu)+'\n')
        f.write('variable		lctemp equal '+str(lcut)+'\n')
        f.write('variable		uctemp equal '+str(ucut)+'\n')
        f.write('variable       overlapL equal '+str(overlapL)+'\n')
        f.write('variable       overlapU equal '+str(overlapU)+'\n')
        f.write('variable		lc equal ${lctemp}+${overlapL}\n')
        f.write('variable		uc equal ${uctemp}-${overlapU}\n')
        
        f.write('variable		lctempdel equal '+str(lcut)+'+'+str(lcutdel)+'\n')
        f.write('variable		uctempdel equal '+str(ucut)+'-'+str(ucutdel)+'\n')
        
        f.write('# ---------- creat sim box---------------------\n')
        f.write('region whole block 0 ${xh} 0 ${yh} ${gl} ${gu} units box\n')
        f.write('create_box 5 whole\n')
        f.write('mass 1 26.982\n')
        f.write('mass 2 26.982\n')
        f.write('mass 3 26.982\n')
        f.write('mass 4 26.982\n')
        f.write('mass 5 26.982\n')
        
        f.write('read_dump '+Folder+gb_name+'.out 0 x y z box yes add yes replace no \n')
        
        f.write('region upper block INF INF INF  INF ${uc} ${gu} units box\n')
        f.write('lattice fcc ${latP} orient x  '+" ".join(str(x) for x in VB[0:3])+' orient y '+ " ".join(str(x) for x in VB[3:6])+' orient z  '+" ".join(str(x) for x in VB[6:9])+'\n')
        f.write('create_atoms 4 region upper \n')
        
        f.write('region lower block INF INF INF  INF ${gl} ${lc} units box\n')
        f.write('lattice fcc ${latP} orient x  '+" ".join(str(x) for x in VA[0:3])+' orient y '+ " ".join(str(x) for x in VA[3:6])+' orient z  '+" ".join(str(x) for x in VA[6:9])+'\n')
        f.write('create_atoms 5 region lower \n')
        
        f.write('variable cutoff_lc equal ${lc}-15\n')
        f.write('variable cutoff_uc equal ${uc}+15\n')
        
        f.write('region lowertemp block INF INF INF  INF ${lctemp} ${lc} units box\n')
        
        f.write('region uppertemp block INF INF INF  INF ${uc} ${uctemp} units box\n')
        
        f.write('region lowertempOrg block INF INF INF  INF ${lctemp} ${lc} units box\n')
        
        f.write('region uppertempOrg block INF INF INF  INF ${uc} ${uctemp} units box\n')
        
        f.write('region org block INF INF INF  INF ${lc} ${uc} units box\n')
        
        f.write('region lowertempdel block INF INF INF  INF ${lctempdel} ${lc} units box\n')
        
        f.write('region uppertempdel block INF INF INF  INF ${uc} ${uctempdel} units box\n')
        
        f.write('#-------group ---------\n')
        f.write('group lowerA type 5\n')
        f.write('group  upperA type 4\n')
        f.write('group  top_orginal type 2\n')
        f.write('group  bottom_orginal type 1\n')
        f.write('group  middel region org\n')
        f.write('group  lt region lowertemp\n')
        f.write('group  ut region uppertemp\n')
        
        f.write('group  ltOrg region lowertempOrg\n')
        f.write('group  utOrg region uppertempOrg\n')
        
        f.write('group  lt1 intersect lt lowerA \n')
        f.write('group  ut1 intersect ut upperA\n')
        
        f.write('group  lt1Org intersect ltOrg bottom_orginal \n')
        f.write('group  ut1Org intersect utOrg top_orginal\n')
        
        f.write('group  ltdel region lowertempdel\n')
        f.write('group  utdel region uppertempdel\n')
        
        f.write('group  lt1del intersect ltdel lowerA \n')
        f.write('group  ut1del intersect utdel upperA\n')
        
        f.write('# ---------- Define Interatomic Potential --------------------- \n')
        f.write('pair_style eam/alloy \n')
        f.write('pair_coeff * * C:/Users/Zero/Documents/lammps_test/Al99.eam.alloy Al Al Al Al Al\n')
        f.write('neighbor 2.0 bin \n')
        f.write('neigh_modify delay 10 check yes\n')
        
        f.write(' # ----------shifting distances ---------- #\n')
        f.write('variable avTop_added_x equal xcm(ut1,x)\n')
        f.write('variable avTop_added_y equal xcm(ut1,y)\n')
        f.write('variable avTop_added_z equal xcm(ut1,z)\n')
        
        f.write('variable avTop_orginal_x equal xcm(ut1Org,x)\n')
        f.write('variable avTop_orginal_y equal xcm(ut1Org,y)\n')
        f.write('variable avTop_orginal_z equal xcm(ut1Org,z)\n')
        
        f.write('variable avBottom_added_x equal xcm(lt1,x)\n')
        f.write('variable avBottom_added_y equal xcm(lt1,y)\n')
        f.write('variable avBottom_added_z equal xcm(lt1,z)\n')
        
        f.write('variable avBottom_orginal_x equal xcm(lt1Org,x)\n')
        f.write('variable avBottom_orginal_y equal xcm(lt1Org,y)\n')
        f.write('variable avBottom_orginal_z equal xcm(lt1Org,z)\n')
        
        
        f.write('# ---------- Displace atoms and delete overlapping atoms ---------------------\n')
        
        f.write('variable      shift_x_top equal v_avTop_orginal_x-v_avTop_added_x\n')
        f.write('variable      shift_y_top equal v_avTop_orginal_y-v_avTop_added_y\n')
        f.write('variable      shift_z_top equal v_avTop_orginal_z-v_avTop_added_z\n')
        
        f.write('variable      shift_x_bottom equal v_avBottom_orginal_x-v_avBottom_added_x\n')
        f.write('variable      shift_y_bottom equal v_avBottom_orginal_y-v_avBottom_added_y\n')
        f.write('variable      shift_z_bottom equal v_avBottom_orginal_z-v_avBottom_added_z\n')
        
        
        
        f.write('displace_atoms upperA move ${shift_x_top} ${shift_y_top} 0  units box\n')
        f.write('displace_atoms lowerA move ${shift_x_bottom} ${shift_y_bottom} 0 units box\n')
        
        f.write('variable smallest_distance_u equal  (1.4142135623730951*${latP}/4)+2.04\n')
        f.write('variable smallest_distance_l equal  (1.4142135623730951*${latP}/4)+1\n')
        f.write('delete_atoms overlap ${smallest_distance_l} lowerA bottom_orginal\n')
        f.write('delete_atoms overlap ${smallest_distance_u} upperA top_orginal\n')
        
        #f.write('delete_atoms group lt1del \n')
        #f.write('delete_atoms group ut1del \n')
        
        
    
        
        f.write('# ---------- Define Settings ---------------------\n')
        f.write('compute csym all centro/atom fcc\n')
        f.write('compute eng all pe/atom \n')
        f.write('compute eatoms all reduce sum c_eng\n')
        f.write('# ---------- Run Minimization --------------------- \n')
        f.write('reset_timestep 0\n')
        f.write('thermo 1000\n')
        f.write('thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms\n')
        f.write('min_style cg ')
        f.write('minimize 1e-30 1e-30 1000000 1000000\n')
    
        
        f.write('# Now allow the box to expand/contract perpendicular to the grain boundary\n')
        f.write('reset_timestep 0 \n')
        f.write('thermo 1000\n')
        f.write('thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms\n')
        f.write('fix 1 lowerA box/relax z 4 x 4 y 4\n')
        f.write('unfix 1\n')
        f.write('fix 2 upperA box/relax z 4 x 4 y 4\n')
        f.write('unfix 2\n')
        #f.write('fix 3 all box/relax z ${smallest_distance} x ${smallest_distance} y ${smallest_distance}\n')
        #f.write('unfix 3\n')
      
        f.write('min_style cg\n ')
        f.write('minimize 1e-15 1e-15 50000 50000\n')
        
        f.write('write_data ''full_'+gb_name+'.out\n')
        
        
        f.write('print "All done"\n')
        
    f.close() 
        
    generate_and_wite_ori_files(LatP,VA,VB,gb_name,Temp)
    
    return gb_name+"_min.input"
    


def generate_lammps_input_for_scratch(LatP,cellsize,min_z,max_z,cellorigin,VA,VB,gb_name,Temp):
    
    
    xl, yl, zl = 0.0, 0.0, -cellsize[2,2]
    xh, yh, zh = cellsize[0,0], cellsize[1,1], cellsize[2,2]
    gl,gu=cellorigin[2][0]-0.001,cellsize[2,2]+cellorigin[2][0]+0.001
    lcut,ucut=min_z,max_z
    Folder='F:/CGBS_and_E/Computed_GBE_and_Structure_Al/Al_minEnergyGBstructs_2022/minEnergyGBstructs/'
    with open(gb_name+"_min.input", "w") as f:
        f.write('# Gashaw Bizana, NOV 2023. Script for creating GB\n')
        f.write('# ---------- Initialize Simulation --------------------- \n')
        f.write('clear\n')
        f.write('units metal\n')
        f.write('dimension 3 \n')
        f.write('boundary p p p\n') 
        f.write('atom_style atomic\n')
        f.write('# ---------- define variables for creating box ---------------------\n')
        f.write('variable latP equal 4.05\n')
        f.write('variable		xl equal '+str(xl)+'\n')
        f.write('variable		yl equal '+str(yl)+'\n')
        f.write('variable		zl equal '+str(zl)+'\n')
        #f.write('variable latP equal 4.05\n')
        f.write('variable		xh equal '+str(xh)+'\n')
        f.write('variable		yh equal '+str(yh)+'\n')
        f.write('variable		zh equal '+str(zh)+'\n')
        f.write('variable		gl equal '+str(gl)+'\n')
        f.write('variable		gu equal '+str(gu)+'\n')
        f.write('variable		lc equal '+str(lcut)+'\n')
        f.write('variable		uc equal '+str(ucut)+'\n')
        
        f.write('# ---------- creat sim box---------------------\n')
        f.write('region whole block 0 ${xh} 0 ${yh} ${gl} ${gu} units box\n')
        f.write('create_box 2 whole\n')
        f.write('mass 1 26.982\n')
        f.write('mass 2 26.982\n')
        
        
        #f.write('read_dump '+Folder+gb_name+'.out 0 x y z box yes add yes replace yes \n')
        
        f.write('region upper block INF INF INF  INF 0 ${gu} units box\n')
        f.write('lattice fcc 4.05 orient x  '+" ".join(str(x) for x in VB[0:3])+' orient y '+ " ".join(str(x) for x in VB[3:6])+' orient z  '+" ".join(str(x) for x in VB[6:9])+'\n')
        f.write('create_atoms 1 region upper \n')
        
        f.write('region lower block INF INF INF  INF ${gl} 0 units box\n')
        f.write('lattice fcc 4.05 orient x  '+" ".join(str(x) for x in VA[0:3])+' orient y '+ " ".join(str(x) for x in VA[3:6])+' orient z  '+" ".join(str(x) for x in VA[6:9])+'\n')
        f.write('create_atoms 2 region lower \n')
        

        
        
        f.write('region org block INF INF INF  INF ${lc} ${uc} units box\n')
        
        f.write('#-------group ---------\n')
        #f.write('group lower type 1\n')
        #f.write('group  upper type 2\n')
        f.write('group lowerA region lower\n')
        f.write('group  upperA region upper\n')
        f.write('# ---------- Define Interatomic Potential --------------------- \n')
        f.write('pair_style eam/alloy \n')
        f.write('pair_coeff * * C:/Users/Zero/Documents/lammps_test/Al99.eam.alloy Al Al\n')
        f.write('neighbor 2.0 bin \n')
        f.write('neigh_modify delay 10 check yes\n')
        f.write('# ---------- Displace atoms and delete overlapping atoms ---------------------\n')
        f.write('displace_atoms upperA move 0 0 0 units lattice \n')
        f.write('delete_atoms overlap 0.35 lowerA upperA \n')
        
        
        f.write('# ---------- Define Settings ---------------------\n')
        f.write('compute csym all centro/atom fcc\n')
        f.write('compute eng all pe/atom \n')
        f.write('compute eatoms all reduce sum c_eng\n')
    
        
        f.write('# ---------- Run Minimization --------------------- \n')
        f.write('reset_timestep 0\n')
        f.write('thermo 10\n')
        f.write('thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms\n')
        f.write('min_style cg\n ')
        f.write('minimize 1e-50 1e-50 50000 100000\n')
        
        f.write('# Now allow the box to expand/contract perpendicular to the grain boundary\n')
        f.write('fix 1 all box/relax z 0 vmax 0.001\n')
        
        # f.write('reset_timestep 0 \n')
        # f.write('thermo 10\n')
        # f.write('thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms\n')
        
        
        # f.write('min_style cg\n ')
        f.write('minimize 1e-50 1e-50 50000 100000\n')
        
        
        # f.write('minimize 1e-15 1e-15 5000 5000\n')
        
        f.write('write_data ''full_'+gb_name+'.out\n')
        
        
        f.write('print "All done"\n')
        
    f.close() 
        
    generate_and_wite_ori_files(LatP,VA,VB,gb_name,Temp)
    
    return gb_name+"_min.input"

Run_name=[]

    
for i in [81]:
    
    LatP=4.0932
    Temp=700
    cellsize=SimCell[i]
    min_z=MinZ[i]
    max_z=MaxZ[i]
    cellorigin=CellOrigin[i]
    VA=MA[i]
    VB=MB[i]
    gb_name=Name[i][:-4]
    print(gb_name)
    
    #lamp_run=generate_lammps_input_for_scratch(LatP,cellsize,min_z,max_z,cellorigin,VA,VB,gb_name,Temp)
    lamp_run=gemerate_lammps_input_for_minimum(LatP,cellsize,min_z,max_z,cellorigin,VA,VB,gb_name,Temp)
    Run_name.append(lamp_run)
    narg="F:/CGBS_and_E/Computed_GBE_and_Structure_Al/Gen_GB_min/"+lamp_run
    command = ["C:/Program Files/Microsoft MPI/Bin/mpiexec.exe", '-n', '4', 'lmp','-in',narg]
    subprocess.run(command, stdout=subprocess.PIPE)

  
    
    
    
    
        
        
        
        
        
        