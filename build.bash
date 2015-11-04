#!/bin/bash

F77=pgf77
F77FLAGS="-fastsse -mcmodel=medium -Bstatic  -Mextend"
INC="/home/vesso/molaris-build/source_official"

MOLARIS_OBJ=" \
process_line.o \
readlib.o \
parmset.o \
read_xyz.o \
killme.o \
fintime.o \
get_code.o \
let2code.o \
typ_to_num.o \
liberr.o \
check_exist_file.o \
check_blank.o \
read_brk.o \
expand_input.o \
trim_to_centroid.o \
trim_to_alpha.o \
check_ang.o \
number_in_lib.o \
find_h.o \
fndchr.o \
fix_resname.o \
fix_atomname.o \
brkerr.o \
lgoodrange.o \
word_count.o \
upper_case.o \
chop_string.o \
help.o \
topoin.o \
check_h_exist.o \
make_alib.o \
chk_nbparm.o \
topology_driver.o \
write_topology.o \
write_sequence.o \
center.o \
create_evb_dat.o \
get_charge_res.o \
getmcphilist.o \
hb_pair.o \
dist.o \
sybyl_lib.o \
connect.o \
addhyd.o \
make_hb.o \
his_h_swich.o \
exclude.o \
bondcode.o \
angcode.o \
torcode_gromos.o \
torcode.o \
itorcode.o \
write_topology_charge.o \
readevblib.o \
def_bond_g.o \
def_bond.o \
def_ang_g.o \
def_ang.o \
puthyd.o \
spray.o \
thrust.o \
get_sp_h.o \
get_wat_h.o \
write_evb.o \
evb_bond.o \
evb_angle.o \
evb_torsion.o \
evb_itorsion.o \
puthydb.o \
elen_for.o \
lcharged_residue.o \
check_env.o \
set_exclude.o \
dist_atoms.o \
find_center.o \
nblist.o \
bondp.o \
thetap.o \
phip.o \
phipsi.o \
nonbond.o \
ephps.o \
exphi.o \
get_xs.o \
enz_default_p.o"

MOLARIS_OBJ_PATH="/home/vesso/molaris-build/source_official/"

OBJ=" \
pretopq.o \
energy_e.o"

OBJ_PATH="/home/vesso/monte"

OBJ_LIST=""

for i in ${MOLARIS_OBJ}
do
   OBJ_LIST=${OBJ_LIST}" "${MOLARIS_OBJ_PATH}/$i
done

for i in ${OBJ}
do
   OBJ_LIST=${OBJ_LIST}" "${OBJ_PATH}/$i
done

${F77} ${F77FLAGS} -c pretopq.f -I${INC} ${OBJ_LIST}
${F77} ${F77FLAGS} -c energy_e.f -I${INC} ${OBJ_LIST}
${F77} ${F77FLAGS} -o assemble assemble.f -I${INC} ${OBJ_LIST}
