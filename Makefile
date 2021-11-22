
#         Makefile for SOSIE
#         ******************

# Version: trunk


# Create a symbolic link from macro/make.macro_* of your choice to ./make.macro !!!

include make.macro

NCDF_LIB=$(NETCDF_DIR)/lib
NCDF_INC=$(NETCDF_DIR)/include

#========================================================
#         You should not change anything below
#========================================================

# LIBRARY linking:
LD_NC  = -L$(NCDF_LIB) $(L_NCDF) ; # "L_NCDF" defined in make.macro ... (ex: "-lnetcdff")
LD_ALL = -L./lib -lsosie $(LD_NC)

# Disable implicit rules to speedup build
.SUFFIXES:
SUFFIXES :=
%.o:
%.mod:

LIB_SOSIE = lib/libsosie.a

OBJ = obj/io_ezcdf.o \
obj/mod_conf.o \
obj/mod_init.o \
obj/mod_manip.o \
obj/mod_grids.o \
obj/mod_bdrown.o \
obj/mod_drown.o \
obj/mod_akima_2d.o \
obj/mod_bilin_2d.o \
obj/mod_akima_1d.o \
obj/mod_scoord.o \
obj/mod_interp.o \
obj/mod_nemotools.o \
obj/mod_poly.o

#obj/mod_drown.o \

OBJ_I2GT = obj/io_ezcdf.o \
obj/mod_conf.o \
obj/mod_manip.o \
obj/mod_drown.o \
obj/mod_bdrown.o \
obj/mod_bilin_2d.o \
obj/mod_poly.o

#obj/mod_drown.o \

OBJ_I2HS = obj/io_ezcdf.o \
obj/mod_conf.o \
obj/mod_manip.o \
obj/mod_bdrown.o \
obj/mod_bilin_2d.o \
obj/mod_akima_1d.o \
obj/mod_poly.o

OBJ_CRS = obj/mod_nemo.o obj/crs.o obj/crsdom.o


# Modules to install in $INSTALL_DIR/include :
MOD_INST= mod/io_ezcdf.mod \
mod/mod_akima_2d.mod \
mod/mod_bdrown.mod

all: bin/sosie3.x bin/corr_vect.x bin/mask_drown_field.x bin/ij_from_lon_lat.x bin/create_angle_file.x

i2gt: bin/interp_to_ground_track.x
i2hs: bin/interp_to_hydro_section.x

crs: bin/nemo_coarsener.x

#bin/nemo_coarsener_2d.x bin/nemo_coarsener_3d.x


test: bin/test_stuffs.x

bin/sosie3.x: src/sosie.f90 $(LIB_SOSIE)
	@mkdir -p bin
	$(FC) $(FF_SAFE) src/sosie.f90 -o bin/sosie3.x $(LD_ALL)

bin/corr_vect.x: src/corr_vect.f90 $(LIB_SOSIE)
	$(FC) $(FF_SAFE) src/corr_vect.f90 -o bin/corr_vect.x $(LD_ALL)

bin/create_angle_file.x: src/create_angle_file.f90 $(LIB_SOSIE)
	$(FC) $(FF_SAFE) src/create_angle_file.f90 -o bin/create_angle_file.x $(LD_ALL)

bin/test_stuffs.x: src/test_stuffs.f90 $(LIB_SOSIE)
	$(FC) $(FF_SAFE) src/test_stuffs.f90 -o bin/test_stuffs.x $(LD_ALL)


bin/interp_to_ground_track.x: src/interp_to_ground_track.f90 $(OBJ_I2GT)
	@mkdir -p bin
	$(FC) $(FF_SAFE) $(OBJ_I2GT) src/interp_to_ground_track.f90 -o bin/interp_to_ground_track.x $(LD_NC)

bin/interp_to_hydro_section.x: src/interp_to_hydro_section.f90 $(OBJ_I2HS)
	@mkdir -p bin
	$(FC) $(FF_SAFE) $(OBJ_I2HS) src/interp_to_hydro_section.f90 -o bin/interp_to_hydro_section.x $(LD_NC)

OBJ_IJLL = obj/io_ezcdf.o obj/mod_conf.o obj/mod_manip.o
bin/ij_from_lon_lat.x: src/ij_from_lon_lat.f90 $(OBJ_IJLL)
	@mkdir -p bin
	$(FC) $(FF_SAFE) $(OBJ_IJLL) src/ij_from_lon_lat.f90 -o bin/ij_from_lon_lat.x $(LD_NC)


### CRS:
obj/mod_nemo.o: src/crs/mod_nemo.f90
	@mkdir -p obj
	@mkdir -p mod
	$(FC) $(FF) -I$(NCDF_INC) -c src/crs/mod_nemo.f90 -o obj/mod_nemo.o

obj/crs.o: src/crs/crs.f90 obj/mod_nemo.o
	$(FC) $(FF) -I$(NCDF_INC) -c src/crs/crs.f90 -o obj/crs.o

obj/crsdom.o: src/crs/crsdom.f90 obj/crs.o
	$(FC) $(FF) -I$(NCDF_INC) -c src/crs/crsdom.f90 -o obj/crsdom.o

bin/nemo_coarsener.x: src/crs/nemo_coarsener.f90 $(OBJ_CRS) obj/io_ezcdf.o obj/mod_manip.o
	@mkdir -p bin
	$(FC) $(FF) -I$(NCDF_INC) obj/io_ezcdf.o obj/mod_manip.o $(OBJ_CRS) src/crs/nemo_coarsener.f90 -o bin/nemo_coarsener.x $(LD_NC)






$(LIB_SOSIE): $(OBJ)
	@echo
	@mkdir -p lib bin
	ar -rv $(LIB_SOSIE) $(OBJ)
	ranlib $(LIB_SOSIE)
	@echo



# Objects :
# ---------

obj/io_ezcdf.o: src/io_ezcdf.f90
	@mkdir -p obj
	@mkdir -p mod
	$(FC) $(FF_SAFE) -I$(NCDF_INC) -c src/io_ezcdf.f90 -o obj/io_ezcdf.o

obj/mod_conf.o: src/mod_conf.f90 obj/io_ezcdf.o
	@mkdir -p obj
	@mkdir -p mod
	$(FC) $(FF) -c src/mod_conf.f90 -o obj/mod_conf.o

obj/mod_init.o: src/mod_init.f90 obj/mod_conf.o obj/mod_scoord.o
	$(FC) $(FF) -c src/mod_init.f90 -o obj/mod_init.o



obj/mod_grids.o: src/mod_grids.f90 obj/mod_conf.o obj/io_ezcdf.o obj/mod_manip.o
	$(FC) $(FF) -c src/mod_grids.f90 -o obj/mod_grids.o

obj/mod_interp.o: src/mod_interp.f90 obj/mod_conf.o obj/mod_manip.o obj/mod_grids.o obj/mod_bdrown.o obj/mod_drown.o obj/mod_akima_2d.o obj/mod_bilin_2d.o obj/mod_akima_1d.o obj/io_ezcdf.o obj/mod_nemotools.o
	$(FC) $(FF) -c src/mod_interp.f90 -o obj/mod_interp.o

obj/mod_manip.o: src/mod_manip.f90 obj/mod_conf.o obj/io_ezcdf.o
	$(FC) $(FF) -c src/mod_manip.f90 -o obj/mod_manip.o

obj/mod_drown.o: src/mod_drown.f90 obj/mod_conf.o
	$(FC) $(FF) -c src/mod_drown.f90 -o obj/mod_drown.o

obj/mod_bdrown.o: src/mod_bdrown.f90 obj/mod_conf.o
	$(FC) $(FF) -c src/mod_bdrown.f90 -o obj/mod_bdrown.o

obj/mod_akima_2d.o: src/mod_akima_2d.f90
	$(FC) $(FF) -c src/mod_akima_2d.f90 -o obj/mod_akima_2d.o

obj/mod_bilin_2d.o: src/mod_bilin_2d.f90 obj/io_ezcdf.o obj/mod_conf.o obj/mod_manip.o obj/mod_poly.o
	$(FC) $(FF) -I$(NCDF_INC) -c src/mod_bilin_2d.f90 -o obj/mod_bilin_2d.o

obj/mod_akima_1d.o: src/mod_akima_1d.f90
	$(FC) $(FF) -c src/mod_akima_1d.f90 -o obj/mod_akima_1d.o

obj/mod_scoord.o: src/mod_scoord.f90
	$(FC) $(FF) -c src/mod_scoord.f90 -o obj/mod_scoord.o

obj/mod_nemotools.o: src/mod_nemotools.f90
	$(FC) $(FF) -c src/mod_nemotools.f90 -o obj/mod_nemotools.o

obj/mod_poly.o: src/mod_poly.f90 obj/io_ezcdf.o
	$(FC) $(FF) -c src/mod_poly.f90 -o obj/mod_poly.o



OBJ_MSK_DRWN = obj/mod_conf.o obj/mod_manip.o obj/mod_bdrown.o obj/io_ezcdf.o
bin/mask_drown_field.x: src/mask_drown_field.f90 $(OBJ_MSK_DRWN)
	@mkdir -p ./bin
	$(FC) $(FF) $(OBJ_MSK_DRWN) -o bin/mask_drown_field.x src/mask_drown_field.f90 $(LD_NC)


install: all
	@mkdir -p $(INSTALL_DIR)/bin $(INSTALL_DIR)/lib $(INSTALL_DIR)/include
	rsync -avP bin/*.x      $(INSTALL_DIR)/bin/
	rsync -avP $(LIB_SOSIE) $(INSTALL_DIR)/lib/
	rsync -avP $(MOD_INST)  $(INSTALL_DIR)/include/

#	ifneq ("$(wildcard bin/nemo_coarsener.x)","")
#	rsync -avP bin/nemo_coarsener.x $(INSTALL_DIR)/bin/
#	endif

uninstall:
	rm -f $(INSTALL_DIR)/bin/*.x
	rm -f $(INSTALL_DIR)/lib/$(LIB_SOSIE) $(INSTALL_DIR)/bin/apply_zon_corr.x
	rm -f $(INSTALL_DIR)/include/io_ezcdf.mod

clean:
	rm -f bin/* $(LIB_SOSIE) *.out mod/*.mod *.x *~ *\# src/*~ src/*\# *.log *.tmp fort.*
	rm -rf bin mod lib obj

distclean: clean
	rm -f examples/*.nc4 examples/*.nc examples/ORCAX_to_ORCAY/*.nc examples/interp_to_ground_track/*.nc
	rm -f examples/data/*.drwn examples/ij_from_lon_lat/*.out
