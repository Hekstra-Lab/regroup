run efvector.py

# Settings
bg_color white
set ray_opaque_background, true
set map_auto_expand_sym, 1
set mesh_radius, 0.01
set antialias, 1
set direct, 0.5
util.ray_shadows('none')
set cartoon_gap_cutoff, 0

# Load Files
load ../../../data/efx/5_ESF_refinement/8G50_efx_on.pdb, DHFR_A_overview
load ../../../data/efx/5_ESF_refinement/8G50_efx_on.pdb, DHFR_A
remove DHFR_A and not chain A
remove DHFR_A_overview and not chain A
load ../../../data/efx/5_ESF_refinement/8G50_efx_on.pdb, DHFR_B_overview
load ../../../data/efx/5_ESF_refinement/8G50_efx_on.pdb, DHFR_B
remove DHFR_B and not chain B
remove DHFR_B_overview and not chain B

hide everything

# General Reps
set valence, 0
show sticks, (resname NAP+FOL and not hydrogens)
set stick_radius, 0.35, resname NAP+FOL and not hydrogens
show cartoon

color skyblue, DHFR_A_overview and element C
color skyblue, DHFR_A and element C
color raspberry, DHFR_B_overview and element C
color raspberry, DHFR_B and element C

color palegreen, resname FOL and element C
color palecyan, resname NAP and element C

#----------------------------------------------------------------------#
# Hinge-bending overview

distance DHFR_A_overview///53/CA, DHFR_A_overview///23/CA
distance DHFR_B_overview///53/CA, DHFR_B_overview///23/CA
set dash_radius, 0.15

add_efvector("DHFR_A", 0, 1, 0, [  -5.943,  29.579,  16.975], invert_polarity=True, expansion=2.0, color=[218/255, 77/255, 155/255])

viewport 1200, 800
set_view (\
    -0.398880333,   -0.296888173,   -0.867608190,\
     0.678282142,    0.541189015,   -0.497034013,\
     0.617103219,   -0.786743760,   -0.014500832,\
     0.000000000,    0.000000000, -230.398727417,\
    15.268459320,   34.231487274,   24.942375183,\
   204.869674683,  255.927658081,  -20.000000000 )

disable DHFR_A
disable DHFR_B

ray
png pngs/hinge_distance.png, dpi=1200
scene 1_overview, store

disable dist*
disable efvector*
disable DHFR_A_overview
disable DHFR_B_overview
enable DHFR_A
enable DHFR_B

#----------------------------------------------------------------------#
# Folate carboxylate

# Add EF vector
add_efvector("DHFR_A", 0, 1, 0, [  10.596,  18.768,  21.739], name="efvector2")

# Add EF vector
add_efvector("DHFR_A", 0, 1, 0, [  10.596,  18.768,  21.739], invert_polarity=True, name="efvector3", red=False)


# Set up maps
load_mtz ../../../data/efx/5_ESF_refinement/8G50_efx_on.mtz, esfA, 2FOFCWT, PH2FOFCWT
load_mtz ../../../data/efx/5_ESF_refinement/8G50_efx_on.mtz, esfB, 2FOFCWT, PH2FOFCWT

align DHFR_B////CA, DHFR_A////CA
matrix_copy DHFR_B, esfB

# Set up maps
map_double esfA
map_double esfA
map_double esfB
map_double esfB

isomesh iso_esfA1, esfA, 1.5, DHFR_A and resname FOL and not hydrogens, carve=2
isomesh iso_esfB1, esfB, 1.5, DHFR_B and resname FOL and not hydrogens, carve=2

set transparency, 0.1, iso*

color marine, iso_esfA1
color firebrick, iso_esfB1
color skyblue, DHFR_A and element C
color raspberry, DHFR_B and element C

set cartoon_color, gray70, DHFR_A or DHFR_B

set_view (\
     0.838004112,    0.309828222,   -0.449152619,\
    -0.468235463,    0.830964327,   -0.300407797,\
     0.280158103,    0.462050080,    0.841426253,\
    -0.000218548,   -0.000512302,  -26.530130386,\
    14.557941437,   18.741893768,   22.230127335,\
     8.322886467,   44.656707764,  -20.000001907 )

# Write image
ray
png pngs/folate_esf.png, dpi=1200
scene 2_folate, store

disable efector*
disable iso_esf*

#----------------------------------------------------------------------#
# C-terminus

hide sticks, resname FOL+NAP
show sticks, polymer and i. 159+134+157 and not hydrogens

isomesh iso_esfA2, esfA, 1.5, DHFR_A and polymer and i. 159+134+157 and not hydrogens, carve=1.5
isomesh iso_esfB2, esfB, 1.5, DHFR_B and polymer and i. 159+134+157 and not hydrogens, carve=1.5
color marine, iso_esfA2
color firebrick, iso_esfB2

add_efvector("DHFR_A", 0, 1, 0, [  35.923,   9.529,  -4.958], name="efvector4")
add_efvector("DHFR_A", 0, 1, 0, [  35.923,   9.529,  -4.958], invert_polarity=True, name="efvector5", red=False)

set_view (\
    -0.040904526,   -0.487356156,    0.872231543,\
    -0.486216098,   -0.752923191,   -0.443493575,\
     0.872869074,   -0.442237496,   -0.206161961,\
     0.000113845,    0.000462661,  -38.168052673,\
    37.400810242,    9.428668976,    1.106899261,\
    16.316272736,   59.615058899,  -20.000001907 )

# Write image
ray
png pngs/cterminus_esf.png, dpi=1200
scene 3_cterminus, store

#----------------------------------------------------------------------#
# Tyr128

hide sticks
disable efvector*
disable iso_esf*
hide spheres
hide cartoon, i. 103-108

select sel_tyr128A, (DHFR_A and i. 125-128 and not hydrogens)
select sel_tyr128B, (DHFR_B and i. 125-128 and not hydrogens)
show sticks, sel_tyr128A or sel_tyr128B or (i. 100 and not hydrogens)

isomesh iso_esfA3, esfA, 1.5, sel_tyr128A, carve=1.5
isomesh iso_esfB3, esfB, 1.5, sel_tyr128B, carve=1.5
color marine, iso_esfA3
color firebrick, iso_esfB3

add_efvector("DHFR_A", 0, 1, 0, [  30.217,  20.977,  -6.401], name="efvector6")
add_efvector("DHFR_A", 0, 1, 0, [  30.217,  20.977,  -6.401], invert_polarity=True, name="efvector7", red=False)

set_view (\
    -0.451755583,    0.393076003,    0.800854623,\
     0.055028968,   -0.883692443,    0.464793921,\
     0.890406966,    0.254065871,    0.377581537,\
     0.000124430,   -0.000217713,  -41.287311554,\
    25.262060165,   18.502674103,    1.177760124,\
    21.331256866,   59.498729706,  -20.000001907 )

# Write iamge
ray
png pngs/tyr128_esf.png, dpi=1200
scene 4_tyr128, store


#----------------------------------------------------------------------#
# Met20

hide sticks
disable efvector*
disable iso_esf*
set sphere_scale, 0.2

show cartoon
show sticks, resname FOL+NAP and not hydrogens
select met20A, (DHFR_A and polymer and i. 20-22+27 and not hydrogens)
select met20B, (DHFR_B and polymer and i. 20-22+27 and not hydrogens)

show sticks, met20A or met20B or (polymer and i. 22 and not hydrogens)
show spheres, resname HOH and i. 315
set sphere_color, skyblue, DHFR_A
set sphere_color, raspberry, DHFR_B

isomesh iso_esfA4, esfA, 1.5, met20A or (DHFR_A and resname HOH and i. 315), carve=1.5
isomesh iso_esfB4, esfB, 1.5, met20B or (DHFR_B and resname HOH and i. 315), carve=1.5
color marine, iso_esfA4
color firebrick, iso_esfB4

add_efvector("DHFR_A", 0, 1, 0, [   20.951,  17.15,  14.450], name="efvector8")
add_efvector("DHFR_A", 0, 1, 0, [   20.951,  17.15,  14.450], invert_polarity=True, name="efvector9", red=False)
set_view (\
     0.467959702,    0.452750951,   -0.758946896,\
    -0.425864846,    0.868020415,    0.255238980,\
     0.774351358,    0.203771502,    0.599020064,\
    -0.000331290,   -0.000881859,  -49.149620056,\
    15.016478539,   19.625869751,   10.643830299,\
    37.717838287,   60.571060181,  -20.000001907 )

# Write
ray
png pngs/met20_esf.png, dpi=1200
scene 5_activesite, store

#----------------------------------------------------------------------#
# Closed I/II

disable efvector*
disable iso_esf*

select pro21A,  DHFR_A and polymer and i. 21-22 and not hydrogens
select pro21B,  DHFR_B and polymer and i. 21-22 and not hydrogens

hide sticks
hide cartoon
hide spheres
show sticks, pro21A or pro21B

isomesh iso_esfA5, esfA, 1.0, pro21A, carve=1.5
isomesh iso_esfB5, esfB, 1.0, pro21B, carve=1.5
color marine, iso_esfA5
color firebrick, iso_esfB5

set_view (\
    -0.254565001,   -0.793197811,   -0.553176165,\
     0.389857829,    0.439302146,   -0.809316158,\
     0.884967506,   -0.421689302,    0.197404698,\
    -0.000518752,    0.000080653,  -22.323699951,\
     8.719263077,   19.605463028,    8.958110809,\
    20.096973419,   24.616836548,  -20.000001907 )

# Write
ray
png pngs/trp22_closed_states.png, dpi=1200
scene 6_activesite_inset, store
